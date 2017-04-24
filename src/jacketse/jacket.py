#!/usr/bin/env python
# encoding: utf-8
"""
jacket.py

Created by Rick Damiani and Andrew Ning on 2013-10-28.
Copyright (c) NREL. All rights reserved.

!!!!!TO CHECK Buckling LENGTH PASSED TO Tube.py: element length vs. joint2joint length!!!!
ALSO lENGTH OF LEGS may not be right.

All lines marked CJB+, CJB-, or CJBe have been added, "removed" (commented out), or edited respectively by Casey Broslawski. Summer 2016
"""
print("Modified JacketSE") #CJB+

import matplotlib.pyplot as plt #first in the list for peregrine's sake
import math
import copy
import numpy as np
import warnings
from openmdao.main.api import VariableTree,  Component, Assembly, set_as_top
from openmdao.main.datatypes.api import Int, Float, Array, Str, VarTree, Bool, Dict, Instance
from openmdao.lib.drivers.api import COBYLAdriver
#from openmdao.main.api import enable_console
#enable_console()

from scipy.optimize import fsolve
from commonse.SegIntersect import SegIntersect
from commonse.SegIntersect import CalcDist#_akima,
from commonse.Material import Material
from commonse.Tube import Tube
from commonse.RigidMember import RigidMember
#from commonse.environment import Wind, Wave, Soil
from commonse.SoilC import SoilC, SubgrReact,SoilPileStiffness
from commonse.EmbedLength import EmbedLength
from commonse.axisEqual3D import axisEqual3D
from commonse.Frustum import frustum
from commonse.utilities import sind, cosd

from frame3dd import Frame, NodeData, ReactionData, ElementData, Options, \
    StaticLoadCase

from printJacketres import main as printJacket

from VarTrees import JcktGeoOutputs,TwrGeoOutputs,RNAprops, Frame3DDaux
from loads import JcktLoad,LoadOutputs, WaterInputs, WindInputs
from Utilization import UtilAssembly, IEC_PSFS

from pySubDyn import pySubDynA, WriteInputFile, WriteDriver, RunSubDyn, ReadOutput, SubDynAOutputs #CJB+

import sys, os
pi = math.pi
tan=np.tan
pi=np.pi
atan=np.arctan

# -----------------
#  Helper Functions
# -----------------


#______________________________________________________________________________#
def Mat2Member(nmems,matins0):
    """This function takes a set of material properties, and a number of members
        and spits out an array of material objects, one per member.
        INPUTS: \n
         matins0:   VarTree(matInput(),iotype='in',desc='Material Data (name,E,nu,rho,fy,fyc') \n
         nmems=Int(iotype='in',desc='Number of Members')                                       \n
        OUTPUTS: \n
         matobs=Array(iotype='out',dtype=np.object,desc='Array of Material Objects')
    """
    matins=matins0.copy()
    maxno=1 #number of properties given in input, we need to adjust them to fit the nmems members
    matobs=np.empty(nmems,dtype=np.object)

    for (k,v) in matins.__dict__.items():
        if not k.startswith('_'):
            if v.size and (v.size < nmems):
                junk=nmems-v.size
                setattr(matins,k,np.hstack((v,v[-1].repeat(junk))))  #replicate the last value input, for all other members

    for i in range(0,nmems):

        if not matins.E.size:
            E=[]
        else:
            E=matins.E[i]
        if not matins.rho.size:
            rho=[]
        else:
            rho=matins.rho[i]
        if not matins.nu.size:
            nu=[]
        else:
            nu=matins.nu[i]
        if not matins.fy.size:
            fy=[]
        else:
            fy=matins.fy[i]
        if not matins.fyc.size:
            fyc=[]
        else:
            fyc=matins.fyc[i]
        if not matins.matname.size:
            matname=[]
        else:
            matname=matins.matname[i]

        matobs[i]=Material(matname=matname,E=E,nu=nu,\
                           rho=rho,fy=fy,fyc=fyc)

    return matobs

#________________________________________________________________________________#



# -----------------
#  Variable Trees
# -----------------

class MatInputs(VariableTree):
    matname =Array(np.array(['']), dtype=str, desc='Names for Materials')  #This makes an Array of variable length strings; NOTE IT MUST BE CALLED names not name as singular is reserved for openmdao!!
    E    =Array(np.array([]),    units='N/m',  dtype=np.float, desc='Young''s Moduli')
    nu   =Array(np.array([]),      units=None,   dtype=np.float, desc='Poisson''s coeffcients')
    rho  =Array(np.array([]),     units='kg/m**3',dtype=np.float, desc='Densities')
    fy   =Array(np.array([]),    units='N/m**2', dtype=np.float, desc='Tensile Yield Strengths')
    #fyc  =Array(fy.default_value, units='N/m**2', dtype=np.float, desc='Compressive Strengths')
    fyc=fy.clone(desc='Compressive Strengths')   #TO DO I would like this one to be filled in with fy stuff if fyc is not given

class JcktGeoInputs(VariableTree):
    """Basic Geometric Inputs need to build a Jacket Model"""

    # inputs
    nlegs    =    Int(4,     units=None, desc='Number of Legs [3 or 4]??')
    nbays    =    Int(4,     units=None, desc='Number of Bays')
    batter   =    Float(8.,  units=None, desc='2D batter [e.g., 10]')
    dck_botz =    Float(16., units='m',  desc='Deck Underside Elevation MSL')
    dck_width=    Float(     units='m',  desc='Deck Width')   #This needs to be checked. How can I take this value as a default coming from the Tower Db?  2*Db TO DO
    dck_widthfrac= Float(2.,  units=None,  desc='Deck Width Factor in fraction of Tower Db. Careful: If dck_width given, this is ignored.')
    weld2D   =    Float(0.2, units=None, desc='Weldment allowance in fraction of Leg OD')
    VPFlag   =    Bool(False,units=None, desc='Vertical Pile Flag [Y/N]: If True the Mudbrace is put at the expected joint with Piles, i.e. at the bottom of leg.')
    clamped  =    Bool(True, units=None, desc='Bottom of structure clamped or connected through springs: Note if AFflag=True in piles, then clamped will be reset to True.')
    AFflag   =    Bool(False,units=None, desc='Apparent Fixity Flag: if True AF is activated. Note if AFflag=True in piles, then clamped will be reset to True.')
    PreBuildTPLvl=Int(0,  units=None, iotype='in', desc='Level of Prebuild [1-5], 0 ignored; see JcktPreBuild component')

class TPlumpMass(VariableTree):
    """Basic Inertial Properties of Optional TP lumped Mass"""
    mass=Float(0.0, units='kg',      desc='TP lumped mass at intersection of diagonal braces, base of stem')    #TP mass
    I   =Array(np.zeros(6), dtype=np.float, units='kg*m**2',desc='TP [IXX,IYY,IZZ,IXY,IXZ,IYZ] @ base of stem')
    CMoff=Array(np.zeros(3), dtype=np.float,units='m', desc='TP lumped mass CM [x,y,z] offset from base of stem')       # TP lumped CMx,y,zoff [m]


# -----------------
#  Components
# -----------------
class LegGeoInputs(VariableTree):
    """Basic Geometric Inputs need to build Legs of Jacket"""
    Dleg_def = 1.2     #default value for Dleg
    tleg_def = 0.0254  #default value for tleg

    Dleg     = Array(np.array([Dleg_def]).repeat([JcktGeoInputs().nbays+1]), units='m', dtype=np.float, desc='Leg OD')
    tleg     = Array(np.array([tleg_def]).repeat([JcktGeoInputs().nbays+1]), units='m', dtype=np.float, desc='Leg Wall Thickness')
    Dleg0    = Float( units='m',  desc='Leg OD that will be replicated for all the leg stations. Optional in place of Dleg. If set it will trump and replace Dleg with a constant value.')
    tleg0    = Float( units='m',  desc='Leg Wall Thickness that will be replicated for all the leg stations. Optional in place of tleg. If set it will trump and replace tleg with a constant value.')

    legZbot  = Float(0.,        units='m', desc='Bottom Elevation of Leg from Sea-Floor')
    Kbuck    = Float(1,         units=None, desc='Leg Kbuckling Factor')   #Effective length factor
    legmatins= VarTree(MatInputs(),        desc="Leg Material Data")
    ndiv     = Int(1,        units=None,  desc='Number of FE elements per leg member')
#_____________________________________________________#
class LegGeoOutputs(VariableTree):
    """Basic Geometric Outputs needed to build Xbraces of Jacket"""
    LegObj   = Instance(Klass=Tube, desc='object of tube class')
    xyz      = Array(units='m', dtype=np.float,  desc='Leg Node Coordinates') #shape=(3,nlegs,(nbays+1)*ndiv+1),
    nNodesLeg= Int(  units=None, desc='Number of (1) Leg Nodes')
    njs=       Int(  units=None, desc='Number of Leg Joints')
    joints   = Array(units='m', dtype=np.float, desc='Leg Joints Coordinates') #shape=(3,nlegs,nbays+2),
#_____________________________________________________#

class PreLegBuild(Component):
    """This component adjusts the Dleg,tleg arrays if Dleg0, tleg0 are given."""
    #inputs
    nbays     =Int(iotype='in', units=None, desc='Number of Bays')
    leginputs =VarTree(LegGeoInputs(),iotype='in',  desc="Leg Input Data")

    #outputs
    prelegouts =VarTree(LegGeoInputs(),iotype='out',  desc="Updated Leg Input Data")

    def execute(self):
        #Simplify nomenclature
        nbays     =self.nbays

        self.prelegouts=self.leginputs #initialize

        #Next two useful in case Dleg0 is the only optimization variable in optimization problems
        if self.leginputs.Dleg0:
           self.prelegouts.Dleg=np.array([self.leginputs.Dleg0]*(nbays+1))
        if self.leginputs.tleg0:
           self.prelegouts.tleg=np.array([self.leginputs.tleg0]*(nbays+1))

        if any(self.prelegouts.Dleg == 0.) :
            self.prelegouts.Dleg=self.prelegouts.Dleg[0].repeat(self.prelegouts.Dleg.size)
        if any(self.prelegouts.tleg == 0.) :
            self.prelegouts.tleg=self.prelegouts.tleg[0].repeat(self.prelegouts.tleg.size)

#_____________________________________________________#

class legs(Component):
    """This component creates a leg object and returns legs'' xyz node coordinates.
        Assumes stump at the bottom of 1.5OD, stump at top defined by TP"""
    #inputs

    #some of these inputs can be taken from Jacket variable trees
    JcktPrms=VarTree(JcktGeoInputs(), iotype='in',  desc='Some Jacket Parameters used to build Legs')

    innr_ang=  Float(      iotype='in',   units='rad',desc='Angle between radial direction and base side, in rad')  #This comes from PreJckt
    bay_hs   = Array(      iotype='in',   units='m', dtype=np.float,  desc='Bay Heights') #This comes from PreJckt
    JcktH    = Float(      iotype='in',   units='m',  desc='Jacket Total Height from legbot to legtop')             #This comes from PreJckt
    dck_width= Float(      iotype='in',   units='m',  desc='Deck Width')   #This comes from PreJckt; note there is dck_width under JcktPrms, but we want an updated version via PreJcktBuild here

    leginputs =VarTree(LegGeoInputs(),iotype='in',  desc="Leg Input Data")
    legbot_stmph =Float( iotype='in',     units='m',desc='Bottom Leg Stump Vertical Height, must be >0 !!')  #Updated legbot_stmph : This should be set to 1.5 Dleg[0] if user did not set it, done in PreBuild

    wdepth   = Float(units='m',  iotype='in',desc='Water Depth') #CJB+
    wlevel   = Float(units='m',  iotype='in',desc='Water Level (total height of submerged+buried structure)') #CJB+
    z_floor   = Float(units='m',  iotype='in',desc='z_floor') #CJB water flag

    #outputs
    bay_hslegs   = Array(          iotype='out',   units='m', dtype=np.float,  desc='Bay Heights Augmented for stump')
    legouts  = VarTree(LegGeoOutputs(), iotype='out', desc='Leg Geometry Output')
    #LegObj   = Instance(Klass=Tube,iotype='out', desc='object of tube class')
    #xyz      = Array(          iotype='out',  units='m', dtype=np.float,  desc='Leg Node Coordinates') #shape=(3,nlegs,(nbays+1)*ndiv+1),
    #nNodesLeg= Int(iotype='out', units=None, desc='Number of (1) Leg Nodes')
    #njs=       Int(iotype='out', units=None, desc='Number of Leg Joints')
    #joints   = Array(iotype='out',units='m', dtype=np.float, desc='Leg Joints Coordinates') #shape=(3,nlegs,nbays+2),



    def execute(self):
        #Simplify nomenclature
        nbays     =self.JcktPrms.nbays
        nlegs     =self.JcktPrms.nlegs
        batter=self.JcktPrms.batter
        weld2D   =  self.JcktPrms.weld2D
        dck_width=self.dck_width
        Dleg=self.leginputs.Dleg
        tleg=self.leginputs.tleg

        ndiv=self.leginputs.ndiv
        legmatins =self.leginputs.legmatins
        stump_h=self.legbot_stmph    #Leg bottom stump vertical height

        wdepth=self.wdepth #CJB+
        wlevel=self.wlevel #CJB+
        z_floor=self.z_floor #CJB water flag

        # Nodes associated with joints starting from joint with pile, and including X-brace joints
        # Note: these are not necessarily the FEA model nodes if more than 1 element are
        #       used for each member, though for the time being they are, we will need TO CHANGE THIS
        # Nodes are counted CounterClockWise starting from x<0,y<0 quadrant, and bottom-up

        #Initialize output
        self.legouts.njs=nbays+2 #number of joints
        self.legouts.nNodesLeg=( (nbays+1)*ndiv+1 )  #For 1 leg
        self.legouts.xyz=np.zeros([3,nlegs,self.legouts.nNodesLeg])

        #First joint (connection to pilehead) coordinates
        self.legouts.xyz[0,0,0]=-dck_width/2. + (1+weld2D)*Dleg[-1]/2. - (self.JcktH)/batter
        #self.legouts.xyz[2,0,0]=self.leginputs.legZbot #CJB- Original code modified in the line below
        self.legouts.xyz[2,0,0]=self.leginputs.legZbot-wlevel #CJBe Successfully shifts the jacket (everything above the mudbraces) down
        #Add stumps to bay heights, at the bottom
        self.bay_hslegs=np.insert(self.bay_hs,0,stump_h) #np.insert(self.bay_hs,-1*wdepth,stump_h-wdepth)
        baydiv=self.bay_hslegs.repeat(ndiv)/ndiv

        #Other coordinates
        self.legouts.xyz[2,0,1:]=  self.legouts.xyz[2,0,0]+baydiv.cumsum() #z coordinates 1st leg
        #self.legouts.xyz[2,0,-1]+= Dleg[-1]/2.    #not sure why i had this, i remove it for now                 #z coordinates 1st leg
        #self.legouts.xyz[0,0,1:]=  self.legouts.xyz[0,0,0]+\
        #                        (self.legouts.xyz[2,0,1:]-self.leginputs.legZbot+wdepth)/batter    #x coordinates 1st leg #CJB- Original code modified in the line below
        self.legouts.xyz[0,0,1:]=  self.legouts.xyz[0,0,0]+\
                                (self.legouts.xyz[2,0,1:]-self.leginputs.legZbot+wlevel)/batter #x coordinates 1st leg #CJBe Add wlevel to the parenthesis to account for shift in self.legouts.xyz[2,0,1:]
        self.legouts.xyz[1,0,:]=self.legouts.xyz[0,0,0:]*np.tan(self.innr_ang) #y coordinates 1st leg
        #Now create nodes for other legs
        self.legouts.xyz[2,1:,:]=self.legouts.xyz[2,0,:]  #z,to replicate perhaps

        if nlegs==4:
            self.legouts.xyz[0,1:3,:]=-self.legouts.xyz[0,0,:] #x legs 2 and 3
            self.legouts.xyz[0,3,:]  = self.legouts.xyz[0,0,:] #x leg 4
            self.legouts.xyz[1,2:4,:]=-self.legouts.xyz[1,0,:] #y legs 3 and 4
            self.legouts.xyz[1,1,:]  = self.legouts.xyz[1,0,:] #y leg 2
        else:  #3 legs
            self.legouts.xyz[0,1,:]=-self.legouts.xyz[0,0,:] #x
            self.legouts.xyz[0,2,:]= 0. # x
            self.legouts.xyz[1,2,:]=-self.legouts.xyz[0,0,:]/np.cos(self.innr_ang) #y
            self.legouts.xyz[1,1,:]= self.legouts.xyz[1,0,:] #y

        #Now store the Joints separately as they will be needed for Xbrc connections
        self.legouts.joints=self.legouts.xyz[:,:,::ndiv]  #joints
        legjnts=self.legouts.joints.transpose([1,2,0])
        dists=(np.roll(legjnts,-1,axis=1)-legjnts)[:,:-1,:]
        #For the lengths of the various members; we should assign the joint-joint leg to it for buckling checks;
        Leg_lgths=np.sqrt( (dists**2).sum(axis=2))[0].repeat(ndiv)  #Note we just assign lengths for the members of leg 1 for consistency with legs, Xbrc

        #materials
        legmats=Mat2Member(self.legouts.njs-1,legmatins)  #First create joint-to-joint member material
        #now to create the Tube object one per member:
        #first: replicate Dleg, tleg appropriately toa ccoutn for members
        Dlegs=Dleg.repeat(ndiv)
        tlegs=tleg.repeat(ndiv)
        #then replicate equally the material objects
        legmats=legmats.repeat(ndiv)
        #if (self.legmats.size==1) :
        #    legmats=legmats.repeat(self.ndiv)
        #elif (self.legmats.size != self.joints.size):
        #    sys.exit('Leg material objects need to be as many as the joints or just 1')

        #Now need to get structural properties  HERE : LEG LENGTHS NEED TO BE CALCULATED!!!
        if Leg_lgths.size != Dlegs.size:
            sys.exit('Error: Check number of bays: Dlegs.size <> then leg number of members!!!')
        self.legouts.LegObj=Tube(Dlegs,tlegs,Leg_lgths,self.leginputs.Kbuck,legmats) #Initialize just like the parent class D,t, etc., then add a few more attributes

#_____________________________________________________#
class XBrcGeoInputs(VariableTree):
    """Basic Geometric Inputs need to build Xbraces of Jacket"""
    Dbrc_def = 1.2  #default value for Dbrc
    tbrc_def = 0.0254  #default value for tbrc
    Dbrc     = Array(np.array([Dbrc_def]).repeat([JcktGeoInputs().nbays]),    units='m', dtype=np.float,    desc='X-Brace OD')
    tbrc     = Array(np.array([tbrc_def]).repeat([JcktGeoInputs().nbays]),    units='m', dtype=np.float,      desc='X-Brace Wall Thickness')
    Dbrc0     = Float(  units='m', dtype=np.float,    desc='X-Brace OD. Optional in place of Dbrc. If set it will trump and replace Dbrc with a constant value.')
    tbrc0     = Float(  units='m', dtype=np.float,    desc='X-Brace Wall Thickness. Optional in place of tbrc. If set it will trump and replace tbrc with a constant value.')

    Kbuck    = Float(0.8,      units=None, desc='XBrc Kbuckling Factor')   #Effective length factor
    Xbrcmatins=VarTree(MatInputs(),        desc="Xbrc Material Input Data")
    ndiv     = Int(1,          units=None, desc='Number of FE elements per Xbrace member')
    precalc  =Bool(False,      units=None, desc='Flag to identify whether or not we should precalculate the Xbrc D and t based on rules of thumb')

class XBrcGeoOutputs(VariableTree):
    """Basic Geometric Outputs needed to build Xbraces of Jacket"""
    LLURObj   = Instance(Klass=Tube,                           desc='Object of Class Tube for LLUR braces')
    ULLRObj   = Instance(Klass=Tube,                           desc='Object of Class Tube for ULLR braces')
    joints   = Array(          units='m', dtype=np.float, desc='Xbrace Joint (intersection) Coordinates(3,nbays,nlegs=nfaces)')
    nodesLLUR1=Array(          units='m', dtype=np.float, desc='Xbraces'' intermediate Coordinates for LL-UR below Xjoints (3,(ndiv-1)*nbays,nlegs=nfaces)')
    nodesLLUR2=Array(          units='m', dtype=np.float, desc='Xbraces'' intermediate Coordinates for LL-UR above Xjoints (3,(ndiv-1)*nbays,nlegs=nfaces)')
    nodesULLR1=Array(          units='m', dtype=np.float, desc='Xbraces'' intermediate Coordinates for UL-LR above Xjoints (3,(ndiv-1)*nbays,nlegs=nfaces)')
    nodesULLR2=Array(          units='m', dtype=np.float, desc='Xbraces'' intermediate Coordinates for UL-LR below Xjoints (3,(ndiv-1)*nbays,nlegs=nfaces)')
    nNodesXbrc= Int(            units=None,                desc='Number of Xbrc Nodes (1 Face) not counting joints at legs')

class Xbraces(Component):
    """This component creates nodes for the X-braces """

    #inputs
    Xbrcinputs=VarTree(XBrcGeoInputs(),iotype='in',desc="Xbrace Input Data")
    nbays =Int(iotype='in', units=None, desc='Number of Bays')

    legjnts   =Array(dtype=np.float,units='m', iotype='in',desc='Leg Joint Coordinates')
    LegObj   = Instance(Klass=Tube,                iotype='in',desc='Object of tube class as coming from Component Legs')#optional, needed only if precalc=True
    wdepth   = Float(              units='m',  iotype='in',desc='Water Depth')#optional, needed only if precalc=True
    bay_hs   =Array(dtype=np.float,units='m',  iotype='in',desc='Bay Lenghts')#optional, needed only if precalc=True
    bay_bs   =Array(dtype=np.float,units='m',  iotype='in',desc='Bay Base Widths')#optional, needed only if precalc=True
    al_bat2D  =Float(              units='rad',iotype='in',desc='Batter Angle in 2D, angle between vertical and leg in projection (2D) [rad]')#optional, needed only if precalc=True
    innr_ang  =Float(              units='rad',iotype='in',desc='Angle between radial direction and base side, in rad')

    #outputs
    Xbrcouts=VarTree(XBrcGeoOutputs(), iotype='out', desc='Basic output data for Xbraces')
    HbrcD =  Float(units='m',iotype='out',desc='Recommended Hbrace OD')
    Hbrct =  Float(units='m',iotype='out',desc='Recommended Hbrace Wall Thickness')

    def execute(self):

        #Simplify nomenclature
        Dbrc=self.Xbrcinputs.Dbrc
        tbrc=self.Xbrcinputs.tbrc
        nbays=self.nbays
        ndiv=self.Xbrcinputs.ndiv
        nlegs=self.legjnts.shape[1]
        Xbrcmatins=self.Xbrcinputs.Xbrcmatins

        #Next two useful in case Dbrc0 is the only optimization variable in optimization problems
        if self.Xbrcinputs.Dbrc0:
           Dbrc=np.array([self.Xbrcinputs.Dbrc0]*nbays)
        if self.Xbrcinputs.tbrc0:
           tbrc=np.array([self.Xbrcinputs.tbrc0]*nbays)

        if any(Dbrc == 0.) :
            Dbrc=Dbrc[0].repeat(Dbrc.size)
        if any(tbrc == 0.) :
            tbrc=tbrc[0].repeat(tbrc.size)


        # Now let us take care of brace intersection nodes, X-braces
        junk=np.arange(1,nbays+1) #A1 points (LL), leg joints counters
        junk1=junk+1 #A2 points (UR), leg joints counters

        A1=self.legjnts[:,:,junk].reshape(3,-1)  #LLeft point
        A2=np.roll(self.legjnts[:,:,junk1],-1,axis=1).reshape(3,-1) #start from 2nd leg-->roll, UR point
        B1=self.legjnts[:,:,junk1].reshape(3,-1)  #ULeft point
        B2=np.roll(A1,-nbays,axis=1).reshape(3,-1)  #LRight point

        junk,junk1=SegIntersect(A1,A2,B1,B2)
        self.Xbrcouts.joints=junk.reshape(3,-1,nbays) # this works, (it was (nbays,-1) and it would not so to be in same format/shape as brc class definition)

        #2 variables that will be updated later
        A1=A1.reshape(3,-1,nbays)
        A2=A2.reshape(3,-1,nbays)

        LLUR1_lgths=(self.Xbrcouts.joints-A1)   #Length of semi braces LL2UR below intersection (also UL2LR above intersection)
        LLUR2_lgths=(self.Xbrcouts.joints-A2)   #Length of semi braces LL2UR above intersection (also UL2LR below intersection)

        #Initialize

        self.Xbrcouts.nodesLLUR1=np.array([])
        self.Xbrcouts.nodesLLUR2=np.array([])
        self.Xbrcouts.nodesULLR1=np.array([])
        self.Xbrcouts.nodesULLR2=np.array([])
            #This section could be skipped if ndiv=1, in which case only Xjoints are going to be used
        if ndiv != 1:

            B1=B1.reshape(3,-1,nbays)
            B2=B2.reshape(3,-1,nbays)

            delta0=np.arange(1,ndiv) #this will be used to calculate distance vectors from starting joint
            #LL to UR :
            #nodes below the intersection

            deltas=LLUR1_lgths/ndiv #vector distance between 2 consecutive nodes
            deltas =deltas*delta0[:,np.newaxis,np.newaxis,np.newaxis] #dimensions have increased by one to account for ndiv, so each plane is 1,2,..ndiv*deltas [x for bays, then face, then 1..ndiv]
            junk=(ndiv-1)*(ndiv >1)+1*(ndiv ==1)    #repeat factor for A1 (inner nodes)
            self.Xbrcouts.nodesLLUR1=(A1.repeat(junk,axis=0).reshape(3,junk,nlegs,nbays)+deltas.transpose(1,0,2,3))
            self.Xbrcouts.nodesLLUR1=self.Xbrcouts.nodesLLUR1.reshape(3,junk,-1).transpose(0,2,1).reshape(3,nlegs,-1)  #(xyz,rows=faces,cols=bay)
            #nodesLLUR1=(A1.repeat(ndiv-1,axis=0).reshape(3,2,4,4).transpose(1,0,2,3)+deltas).reshape(3,nbays,nlegs,-1)

            #nodes above the intersection

            deltas=LLUR2_lgths/ndiv #vector distance between 2 consecutive nodes
            deltas =deltas*delta0[:,np.newaxis,np.newaxis,np.newaxis] #dimensions have increased by one to account for ndiv, so each plane is 1,2,..ndiv*deltas [x for bays, then face, then 1..ndiv]
            deltas =np.roll(deltas,3,0)  #This is so that I keep the order going from intersect node to A2 joints in the upper part of the brace
            self.Xbrcouts.nodesLLUR2=(A2.repeat(junk,axis=0).reshape(3,junk,nlegs,nbays)+deltas.transpose(1,0,2,3))
            self.Xbrcouts.nodesLLUR2=self.Xbrcouts.nodesLLUR2.reshape(3,junk,-1).transpose(0,2,1).reshape(3,nlegs,-1)

            #nodesLLUR=np.hstack([nodesLLUR1,nodesLLUR2])

            #UL to LR :
            #nodes above the intersection
            deltas=(self.Xbrcouts.joints-B1)/ndiv #vector distance between 2 consecutive nodes
            deltas =deltas*delta0[:,np.newaxis,np.newaxis,np.newaxis] #dimensions have increased by one to account for ndiv, so each plane is 1,2,..ndiv*deltas [x for bays, then face, then 1..ndiv]
            self.Xbrcouts.nodesULLR1=(B1.repeat(junk,axis=0).reshape(3,junk,nlegs,nbays)+deltas.transpose(1,0,2,3))
            self.Xbrcouts.nodesULLR1=self.Xbrcouts.nodesULLR1.reshape(3,junk,-1).transpose(0,2,1).reshape(3,nlegs,-1)

            #nodes below the intersection
            deltas=(self.Xbrcouts.joints-B2)/ndiv #vector distance between 2 consecutive nodes
            deltas =deltas*delta0[:,np.newaxis,np.newaxis,np.newaxis] #dimensions have increased by one to account for ndiv, so each plane is 1,2,..ndiv*deltas [x for bays, then face, then 1..ndiv]
            deltas =np.roll(deltas,3,0)  #This is so that I keep the order going from intersect node to B2 joints in the lower part of the brace
            self.Xbrcouts.nodesULLR2=(B2.repeat(junk,axis=0).reshape(3,junk,nlegs,nbays)+deltas.transpose(1,0,2,3))
            self.Xbrcouts.nodesULLR2=self.Xbrcouts.nodesULLR2.reshape(3,junk,-1).transpose(0,2,1).reshape(3,nlegs,-1)


            #nodesULLR=np.hstack([nodesULLR1,nodesULLR2])

        LLUR1_lgths=np.sqrt((LLUR1_lgths**2).sum(axis=0))[0,:]
        LLUR2_lgths=np.sqrt((LLUR2_lgths**2).sum(axis=0))[0,:]
        LLUR_lgths=np.vstack((LLUR1_lgths,LLUR2_lgths)).transpose().flatten()
        ULLR_lgths=np.roll(LLUR_lgths.reshape(nbays,2),-1,axis=1).flatten()  #bottom up 1 face joint2joint (2 semibrace per brace)
        #Number of nodes for 1 face of the jacket, not including joints at legs
        self.Xbrcouts.nNodesXbrc= (ndiv-1)*4*nbays+nbays        #for 1 face of the jacket only: number of internal nodes per semi-brace * 4 semibraces *nbays + intersection joints
        nmems=ndiv*4*nbays   #number of X members (node to node) per face of jacket
        #materials
        Xbrcmats=Mat2Member(nbays,Xbrcmatins)  #First create joint-to-joint member material

        #now to create the Tube object one per member: DO IT ONLY FOR LLUR; ULLR is the same

        #We need to calculate Dbrcs, tbrcs if precalc was enabled

        if self.Xbrcinputs.precalc:
            Dbrc,tbrc,res=BrcFind(self.LegObj,self.bay_bs,self.bay_hs,self.innr_ang,self.al_bat2D,self.wdepth,mudbrc=False,Kbuck=self.Xbrcinputs.Kbuck,Klr_mx=70,minmass=True,PipeSch=True)

        #first: replicate Dbrcs, tbrcs appropriately to account for members
        Dbrcs=Dbrc.repeat(2*ndiv).flatten()  # LL2UR braces all bays, 1 face
        tbrcs=tbrc.repeat(2*ndiv).flatten()  # LL2UR braces all bays, 1 face
        #then replicate equally the material objects, same deal in terms of order, first LLUR all bays, then ULLR all bays
        Xbrcmats=Xbrcmats.repeat(2*ndiv).flatten() # LL2UR braces all bays, 1 face

        #Now I need the lengths of the various members; but we should assign the joint-joint leg to it;
        #Follow order of members as LL2UR all bays bottom to top first, then UL2LR all bays bottom to top

        LLUR_lgths=LLUR_lgths.repeat(ndiv) #LLUR members from bottom up, all bays 1 face
        ULLR_lgths=ULLR_lgths.repeat(ndiv) #ULLR members from bottom up, all bays 1 face

        #Now need to get structural properties
        self.Xbrcouts.LLURObj=Tube(Dbrcs,tbrcs,LLUR_lgths.flatten(),self.Xbrcinputs.Kbuck,Xbrcmats) #LL2UR: This is for members of 1 face only, with members of each bay first, bottom to top
        self.Xbrcouts.ULLRObj=Tube(Dbrcs,tbrcs,ULLR_lgths.flatten(),self.Xbrcinputs.Kbuck,Xbrcmats) #UL2LR: This is for members of 1 face only, with members of each bay first, bottom to top

        #Store recommended HbrcD and Hbrct
        self.HbrcD=self.Xbrcouts.LLURObj.D[-1]
        self.Hbrct=self.Xbrcouts.LLURObj.t[-1]

#_____________________________________________________#
class MudBrcGeoInputs(VariableTree):
    """Basic Geometric Inputs need to build Mudbraces of Jacket"""
    Dbrc_mud_def = 0.813  #default value for Dbrc_mud
    tbrc_mud_def = 0.03   #default value for tbrc_mud
    Dbrc_mud = Float(Dbrc_mud_def,  units='m', desc='Mud-Brace OD')
    tbrc_mud = Float(tbrc_mud_def,  units='m', desc='Mud-Brace Wall Thickness')
    Kbuck    = Float(0.8,      units=None, desc='Mud-Brc Kbuckling Factor')   #Effective length factor
    Mbrcmatins=VarTree(MatInputs(),        desc="Mud-brc Material Input Data")
    ndiv     = Int(1,          units=None, desc='Number of FE elements per Mud-Brace member')
    precalc  =Bool(False,      units=None, desc='Flag to identify whether or not we should precalculate the Mudbrc D and t based on rules of thumb')

class MudBrcGeoOutputs(VariableTree):
    """Basic Geometric Outputs needed to build Xbraces of Jacket"""
    brcObj   = Instance(Klass=Tube,                           desc='Object of Class Tube')
    joints   = Array(          units='m', dtype=np.float, desc='Mud brace Joint (with legs) Coordinates (3,nlegs=nfaces)')
    nodes    = Array(          units='m', dtype=np.float, desc='Mud braces'' intermediate Coordinates for nodes (3,ndiv-1,nlegs=nfaces)')
    nNodesMbrc= Int(            units=None,                desc='Number of Mudbrc Nodes (1 Face) not counting joints at legs')

class MudBraces(Component):
    """This component creates nodes for the MudBraces"""

    #inputs
    nlegs =Int(4,iotype='in', units=None, desc='Number of Legs')     #TO DO : Need to figure out to get this one from JcktINput

    Mbrcinputs=VarTree(MudBrcGeoInputs(),iotype='in',desc="Mud Brace Input Data")
    VPFlag   = Bool(units=None, iotype='in',desc='Vertical Pile Flag')
    PileFlag = Bool(units=None, iotype='in',desc='Piles [Y/N] Flag')
    legjnts   =Array(units='m', iotype='in', dtype=np.float,  desc='Leg Joint Coordinates')

    LegObj   = Instance(Klass=Tube,                iotype='in',desc='Object of tube class as coming from Component Legs')#optional, needed only if precalc=True
    wdepth   = Float(              units='m',  iotype='in',desc='Water Depth')#optional, needed only if precalc=True
    bay_bs   =Array(dtype=np.float,units='m',  iotype='in',desc='Bay Base Widths')#optional, needed only if precalc=True

    #outputs
    Mbrcouts=VarTree(MudBrcGeoOutputs(), iotype='out', desc='Basic output data for MudBrace')

    def execute(self):
        #Simplify nomenclature
        Dbrc=self.Mbrcinputs.Dbrc_mud
        tbrc=self.Mbrcinputs.tbrc_mud
        ndiv=self.Mbrcinputs.ndiv
        nlegs=self.nlegs
        Mbrcmatins=self.Mbrcinputs.Mbrcmatins
        legjnts=self.legjnts
        nNodesleg=legjnts.shape[2]
        # Now let us take care of intersection nodes, joints
        if self.VPFlag and self.PileFlag:
            self.Mbrcouts.joints=legjnts[:,:,0].transpose()  #Mud-Brace Joints, are the joints with piles if vertical piles,else keep it at Xbrace joint
            VPFlag=True
        else:
            self.Mbrcouts.joints=legjnts[:,:,1].transpose()  #Mud-Brace Joints, are the joints with piles if vertical piles,else keep it at Xbrace joint
            VPFlag=False
        #Find the intermediate nodes
        factors=np.arange(1,ndiv)

        dists=(np.roll(self.Mbrcouts.joints,-1,0)-self.Mbrcouts.joints)
        #For the lengths of the various members; we should assign the joint-joint leg to it for buckling checks;
        Mbrc_lgths=np.sqrt( (dists**2).sum(axis=1))[0].repeat(ndiv)  #Note we just assign lengths for the members of face 1 for consistency with legs, Xbrc

        nodes=self.Mbrcouts.joints + dists*factors[...,np.newaxis,np.newaxis]/ndiv #nodes excluding joints
        self.Mbrcouts.nodes=nodes.transpose(1,0,2).reshape(-1,3)  #This puts them in the right order CCW


        #Number of nodes for 1 face of the jacket, not including joints at legs
        self.Mbrcouts.nNodesMbrc= (ndiv-1)     #for 1 face of the jacket only: number of internal nodes per mud-brace
        nmems=ndiv*nlegs   #number of mud members (node to node) per face of jacket

        #materials
        Mbrcmats=Mat2Member(1,Mbrcmatins)  #First create joint-to-joint member material

        #now to create the Tube object one per member:
        #We need to calculate Dbrcs, tbrcs if precalc was enabled
        if self.Mbrcinputs.precalc:
            Dbrc,tbrc,res=BrcFind(self.LegObj,self.bay_bs,np.array([]),[],[],self.wdepth,mudbrc=True,VPFlag=VPFlag,Kbuck=self.Mbrcinputs.Kbuck,Klr_mx=70,minmass=True,PipeSch=True)

        #first: replicate Dbrcs, tbrcs appropriately to account for members
        Dbrcs=np.tile(Dbrc,ndiv) # 1 face
        tbrcs=np.tile(tbrc,ndiv) #1 face
        #then replicate equally the material objects, CCW
        Mbrcmats=Mbrcmats.repeat(ndiv).flatten()

        #Now need to get structural properties
        self.Mbrcouts.brcObj=Tube(Dbrcs,tbrcs,Mbrc_lgths.flatten(),self.Mbrcinputs.Kbuck,Mbrcmats) #This is for members of 1 face only, CCW

#_____________________________________________________#
class HBrcGeoInputs(VariableTree):
    """Basic Geometric Inputs need to build Hbraces of Jacket"""
    Dbrc_def = 0.813  #default value for Dbrc
    tbrc_def = 0.03   #default value for tbrc
    Dbrch = Float(Dbrc_def,  units='m', desc='Top H-Brace OD')
    tbrch = Float(tbrc_def,  units='m', desc='Top H-Brace Wall Thickness')
    Kbuck    = Float(0.8,      units=None, desc='Top H-Brace Kbuckling Factor')   #Effective length factor
    Hbrcmatins=VarTree(MatInputs(),        desc="Top H-Brace Material Input Data")
    ndiv     = Int(1,          units=None, desc='Number of FE elements per H-Brace member')
    precalc  =Bool(False,      units=None, desc='Flag to identify whether or not we should precalculate the Hbrc D and t set = Top Xbrace D and t')

class HBrcGeoOutputs(VariableTree):
    """Basic Geometric Outputs needed to build Xbraces of Jacket"""
    brcObj   = Instance(Klass=Tube,                           desc='Object of Class Tube')
    joints   = Array(          units='m', dtype=np.float, desc='Top H-brace Joint (with legs) Coordinates (3,nlegs=nfaces)')
    nodes    = Array(          units='m', dtype=np.float, desc='Top H-braces'' intermediate Coordinates for nodes (3,ndiv-1,nlegs=nfaces)')
    nNodesHbrc= Int(           units=None,                desc='Number of Hbrc Nodes (1 Face) not counting joints at legs')

class HBraces(Component):
    """This component creates nodes for the Horizontal Braces at the top, below the girder if wanted."""

    #inputs
    nlegs =Int(4,iotype='in', units=None, desc='Number of Legs')     #TO DO : Need to figure out to get this one from JcktINput

    Hbrcinputs=VarTree(HBrcGeoInputs(),iotype='in',desc="Top H-Brace Input Data")

    legjnts   =Array(units='m', iotype='in', dtype=np.float,  desc='Leg Joint Coordinates')

    HbrcD =  Float(units='m',iotype='in',desc='Recommended Hbrace OD from Xbraces')#optional, needed only if precalc=True
    Hbrct =  Float(units='m',iotype='in',desc='Recommended Hbrace Wall Thickness')#optional, needed only if precalc=True

    #outputs
    Hbrcouts=VarTree(HBrcGeoOutputs(), iotype='out', desc='Basic output data for Top H-brace')

    def execute(self):
        #Simplify nomenclature
               #Simplify nomenclature
        Dbrch=self.Hbrcinputs.Dbrch
        tbrch=self.Hbrcinputs.tbrch
        ndiv=self.Hbrcinputs.ndiv
        nlegs=self.nlegs
        Hbrcmatins=self.Hbrcinputs.Hbrcmatins
        legjnts=self.legjnts
        nNodesleg=legjnts.shape[2]
        # Now let us take care of joints

        self.Hbrcouts.joints=legjnts[:,:,-1].transpose()  #Top H-Brace Joints

        #Find the intermediate nodes

        factors=np.arange(1,ndiv)

        dists=(np.roll(self.Hbrcouts.joints,-1,0)-self.Hbrcouts.joints)
        #For the lengths of the various members; we should assign the joint-joint leg to it for buckling checks;
        Hbrc_lgths=np.sqrt( (dists**2).sum(axis=1))[0].repeat(ndiv)  #Note we just assign lengths for the members of face 1 for consistency with legs, Xbrc

        nodes=self.Hbrcouts.joints + dists*factors[...,np.newaxis,np.newaxis]/ndiv #nodes excluding joints
        self.Hbrcouts.nodes=nodes.transpose(1,0,2).reshape(-1,3)  #This puts them in the right order CCW


        #Number of nodes for 1 face of the jacket, not including joints at legs
        self.Hbrcouts.nNodesHbrc= (ndiv-1 +(ndiv==0))     #for 1 face of the jacket only: number of internal nodes per h-brace; this makes sure we account for ndiv=0 i.e. no hbrace
        nmems=ndiv*nlegs   #number of mud members (node to node) per face of jacket

        #materials
        Hbrcmats=Mat2Member(1,Hbrcmatins)  #First create joint-to-joint member material

        #now to create the Tube object one per member:
        if self.Hbrcinputs.precalc:
            Dbrch=self.HbrcD
            tbrch=self.Hbrct
        #first: replicate Dbrcs, tbrcs appropriately to account for members
        Dbrcs=np.tile(Dbrch,ndiv) # 1 face
        tbrcs=np.tile(tbrch,ndiv) #1 face
        #then replicate equally the material objects, CCW
        Hbrcmats=Hbrcmats.repeat(ndiv).flatten()

        #Now need to get structural properties
        self.Hbrcouts.brcObj=Tube(Dbrcs,tbrcs,Hbrc_lgths.flatten(),self.Hbrcinputs.Kbuck,Hbrcmats) #This is for members of 1 face only, CCW

#_____________________________________________________#
class TPGeoInputs(VariableTree):
    """Basic Geometric Inputs need to build Transition Piece of Jacket"""

    #Struts are inclined arms
    Dstrut = Float(       units='m', desc='TP Inclined Arm OD')
    tstrut = Float(       units='m', desc='TP Inclined Arm Wall Thickness')
    strutndiv= Int(1,     units=None, desc='Number of FE elements per strut member')
    TPstrtmatins=VarTree(MatInputs(),        desc="TP Strut Material Input Data")
    Kbuck_strut= Float(0.8,      units=None, desc='TP Strut Kbuckling Factor')   #Effective length factor

    #stumps connect to legs
    Dstump = Float(      units='m', desc='TP-leg Stumps OD')
    tstump = Float(      units='m', desc='TP-leg Stumps Wall Thickness')
    hstump = Float(      units='m', desc='TP-leg Stumps Length')
    stumpndiv= Int(1,    units=None, desc='Number of FE elements per stump member')
    TPstmpmatins=VarTree(MatInputs(),desc="TP Stump Material Input Data")
    Kbuck_stump=Float(1,units=None,desc='TP Stump Kbuckling Factor')   #Effective length factor

    #braces are the horizontal Xbraces for the TP (diagonals)
    Dbrc = Float(      units='m', desc='TP brace OD')
    tbrc = Float(      units='m', desc='TP brace Wall Thickness')
    brcndiv= Int(1,    units=None, desc='Number of FE elements per brace member')
    TPbrcmatins=VarTree(MatInputs(),desc="TP brace Material Input Data")
    Kbuck_brc= Float(0.8,      units=None, desc='TP Brace Kbuckling Factor')   #Effective length factor

    #girders
    Dgir = Float(      units='m',   desc='TP Girder OD')
    tgir = Float(      units='m',   desc='TP Girder Wall Thickness')
    girndiv= Int(1,    units=None,  desc='Number of FE elements per Girder member')
    TPgirdmatins=VarTree(MatInputs(),desc="TP Girder Material Input Data")
    Kbuck_gir= Float(0.8,units=None,desc='TP Girder Kbuckling Factor')   #Effective length factor


    #Main Stem is the representation of the main shell
    nstems= Int(1,    units=None, desc='Number of Stem SubMembers (Submember is of uniform geometry)')
    Dstem     = Array(    units='m', dtype=np.float,    desc='TP Stem OD')  #np.array([6.])
    tstem     = Array(  units='m', dtype=np.float,    desc='TP Stem Wall Thickness') #np.array([0.10])
    hstem     = Array(np.array([]),     units='m', dtype=np.float,    desc='TP Stem Submember Length')  #sum is TPlth
    stemndiv= Int(1,    units=None,  desc='Number of FE elements per Stem SubMember')
    TPstemmatins=VarTree(MatInputs(),desc="TP Stem Material Input Data")
    Kbuck_stem= Float(0.8,units=None,desc='TP Stem Kbuckling Factor')   #Effective length factor

    stemwallfactor=Float(1.5,units=None,iotype='in', desc='Stem wall factor over tower base wall thickness. Used only in case of PreBuildTPlvl>0.')

    gussets  =  Bool(True, iotype='in', units='-', desc='TP arrangement')  #For now let us just implement gusset type  #TO DO


#_____________________________________________________#
class TPGeoOutputs(VariableTree):
    """Basic Geometric Outputs needed to build Xbraces of Jacket"""
    girObj   = Instance(Klass=Tube,                    desc='Object of Class Tube')  #For girders
    brcObj   = Instance(Klass=Tube,                    desc='Object of Class Tube')  #For diagonals
    stmpObj  = Instance(Klass=Tube,                    desc='Object of Class Tube')  #For stumps
    strtObj = Instance(Klass=Tube,                    desc='Object of Class Tube')  #For struts
    stemObj  = Instance(Klass=Tube,                    desc='Object of Class Tube')  #For main cylinder

    brcjoints   = Array(units='m', dtype=np.float, desc='TP diagonal brace joints (with stumps) Coordinates (nlegs=nfaces,3)')
    stemjoints  = Array(units='m', dtype=np.float, desc='TP main cylinder  joints (with diagonal, with submember, and with struts) Coordinates (nstems+1,3)')

    brcnodes    = Array(units='m', dtype=np.float, desc='TP diagonal brace internal Coordinates excluding joints with stumps and stem (nlegs=nfaces,brcndiv-1,3) CCW order')
    girnodes    = Array(units='m', dtype=np.float, desc='TP girders internal Coordinates excluding joints with stumps (nlegs=nfaces,girndiv-1,3) CCW order')
    stmpnodes   = Array(units='m', dtype=np.float, desc='TP stumps internal Coordinates excluding joints with legs but including joints with braces (nlegs=nfaces,stmpndiv,3) CCW order')
    strtnodes   = Array(units='m', dtype=np.float, desc='TP struts internal Coordinates excluding joints with stumps and stem (nlegs=nfaces,strtndiv-1,3) CCW order')
    stemnodes   = Array(units='m', dtype=np.float, desc='TP main Cylinder Internal Coordinates including joints with stumps and stem (nlegs=nfaces,stemndiv+1,3) CCW order')

    nNodesBrc   = Int(0,units=None,                  desc='Number of Diagonal-Brace-nodes (1 Leg) not counting joints at stumps and stem')
    nNodesGir   = Int(0,units=None,                  desc='Number of Girder-nodes (1 Leg) not counting joints at stumps')
    nNodesStrt  = Int(0,units=None,                  desc='Number of Strut-nodes (1 Leg) not counting joints at stumps and stem')
    nNodesStmp  = Int(1,units=None,                  desc='Number of Stump-nodes (1 Leg) not counting joints at legs, but including joints with braces')
    nNodesStem  = Int(2,units=None,                  desc='Number of Stem-nodes including all joints')
    TPlumpedMass=Array(np.zeros([10]),dtype=np.float, desc='TP lumped mass, Ixx, Iyy, Izz,Ixy,Ixz,Iyz,Cmxoff,Cmyoff,Cmzoff from TPlumpMass properties in input')

    mass        =Float(units='kg', desc='TP structural mass, each subcomponent is assumed to have constant density along the various ndiv segments')

class PreBuildTP(Component):
    """This component takes care of precalculating some dimensions for the TP structure,\n
       especially in the case where we do not verify the TP members."""
    #inputs
    LegObj    = Instance(Klass=Tube,        iotype='in',desc='Object of tube class as coming from Component Legs')
    XBrcObj   = Instance(Klass=Tube,        iotype='in',desc='Object of tube class as coming from Component Xbraces')

    TwrDb         =Float(units='m',     iotype='in', desc='Tower Base OD')
    TwrDTRb       =Float(units='m',     iotype='in', desc='Tower Base DTR')

    TPinputs=VarTree(TPGeoInputs(),     iotype='in', desc='Basic input data for Transition Piece')
    PreBuildTPLvl=Int(    units=None,   iotype='in', desc="""Level of Prebuild [1-5]: 1=Stump & Strut from leg; 2[default]=1+(Gir=diagBrace=Xbrace); \n
                                                               3=(Gir=Brace=Xbrace); 4=diagBrace=girders. 5=1+4    in all cases: Stem from Twr; """)
    #outputs
    BuildTPouts  =VarTree(TPGeoInputs(), iotype='out', desc='Revised Input data for Transition Piece')

    def execute(self):
        self.BuildTPouts=self.TPinputs #here I transfer all the other inputs before expanding them with new info
        PreBuildLvl=self.PreBuildTPLvl #shorten name

         #Set central stem following tower base and bump up thikness
        self.BuildTPouts.Dstem=np.asarray([self.TwrDb]).repeat(self.TPinputs.hstem.size)
        self.BuildTPouts.tstem=np.asarray([self.TwrDb/self.TwrDTRb*self.TPinputs.stemwallfactor]).repeat(self.TPinputs.hstem.size)

        if PreBuildLvl==1 or PreBuildLvl==2 or PreBuildLvl==5:
            #Set struts equal to top of leg
            self.BuildTPouts.Dstrut=self.LegObj.D[-1]
            self.BuildTPouts.tstrut=self.LegObj.t[-1]
            self.BuildTPouts.Dstump=self.LegObj.D[-1]
            self.BuildTPouts.tstump=self.LegObj.t[-1]

        if PreBuildLvl==2 or PreBuildLvl==3:
            #set girders and diagonals equal to top braces
            self.BuildTPouts.Dgir=self.XBrcObj.D[-1]
            self.BuildTPouts.tgir=self.XBrcObj.t[-1]
            self.BuildTPouts.Dbrc=self.XBrcObj.D[-1]
            self.BuildTPouts.tbrc=self.XBrcObj.t[-1]
        if PreBuildLvl==4 or PreBuildLvl==5:
            self.BuildTPouts.Dbrc=self.TPinputs.Dgir
            self.BuildTPouts.tbrc=self.TPinputs.tgir

class TP(Component):

    #inputs
    TPinputs=VarTree(TPGeoInputs(),         iotype='in', desc='Basic input data for Transition Piece')
    TPlumpinputs=VarTree(TPlumpMass(),      iotype='in', desc='Basic Inertial Properties of TP Lumped Mass')  #Extra Inertia at the TP
    nlegs   =Int(4,units=None,              iotype='in', desc='Number of Legs')
    legjnts =Array(units='m',dtype=np.float,iotype='in', desc='Leg Joint Coordinates')

    #outputs
    TPouts  =VarTree(TPGeoOutputs(), iotype='out', desc='Basic output data for Transition Piece')
    TwrBsZ  = Float(units='m', iotype='out',        desc='Top Joint Coordinates, i.e. Tower Connection Z (w.r.t. mudline)')

    def execute(self):
        #Simplify Nomenclature
       legjnts=self.legjnts[:,:,-1].transpose()  #joints with legs
       hstump=self.TPinputs.hstump*(self.TPinputs.stumpndiv !=0)#Assign hstump=0 if stumpndiv  is 0 as well, in that case no need to have stumps
       gussets=self.TPinputs.gussets
       nlegs=self.nlegs
       nstems=self.TPinputs.nstems
       Dgir=self.TPinputs.Dgir
       tgir=self.TPinputs.tgir
       Dbrc=self.TPinputs.Dbrc
       tbrc=self.TPinputs.tbrc
       Dstrut=self.TPinputs.Dstrut
       tstrut=self.TPinputs.tstrut
       Dstem=self.TPinputs.Dstem
       tstem=self.TPinputs.tstem
       Dstmp=self.TPinputs.Dstump
       tstmp=self.TPinputs.tstump
       stumpndiv=self.TPinputs.stumpndiv*(hstump !=0.)  #Assign stumpndiv=0 if hstump is 0 as well, in that case no need to have stumps
       brcndiv=self.TPinputs.brcndiv
       girndiv=self.TPinputs.girndiv
       strutndiv=self.TPinputs.strutndiv
       stemndiv=self.TPinputs.stemndiv
       TPlth=self.TPinputs.hstem.sum()  #Overall length of TP
       girdmatins=self.TPinputs.TPgirdmatins    #Material data for braces (girders + diagonals)
       brcmatins=self.TPinputs.TPbrcmatins    #Material data for braces (girders + diagonals)
       stmpmatins=self.TPinputs.TPstmpmatins  #Material data for stumps
       stemmatins=self.TPinputs.TPstemmatins  #Material data for main cylinder
       strtmatins=self.TPinputs.TPstrtmatins  #Material data for struts


       if gussets:  #simulate rigid gussets with horizontal and slanted braces
             njs=3 #3 joints at the base, mid, and top of TP, others are the leg-hbrc joints
       else:
             njs=7   #creating leg stumps as well

        # CREATE TP NOW
       if gussets:

            #-------------Start with stumps-------------#
            brcjoints=legjnts#Initialize  in case hstump=0 or stumpndiv=0
            factors=hstump*np.arange(1,stumpndiv+1)/stumpndiv
            self.TPouts.stmpnodes=legjnts.repeat(stumpndiv,0)
            self.TPouts.stmpnodes[:,2]+=factors.reshape([1,-1]).repeat(nlegs,0).flatten() #z nodes
            if hstump:  #could not avoif this one
                brcjoints=self.TPouts.stmpnodes[stumpndiv-1::stumpndiv,:]  # also=Girder Joints
            dists=(brcjoints-legjnts) #top stump joints - legjoints
            stmp_lgths=np.sqrt( (dists**2).sum(axis=1))[0].repeat(stumpndiv)  #Note we just assign lengths for the members of "1 leg", for consistency with legs, Xbrc, Mbrc,Hbrc

            self.TPouts.nNodesStmp=stumpndiv  #No. of nodes excluding joint with leg, but incliding joint with brace; for 1 leg-stump only
            #now to create the Tube object one per member:
            stmpmats=Mat2Member(1,stmpmatins)  #First create joint-to-joint member material
            #first: replicate Dbrcs, tbrcs appropriately to account for members
            Ds=np.tile(Dstmp,stumpndiv) # 1 face
            ts=np.tile(tstmp,stumpndiv) #1 face
            #then replicate equally the material objects, CCW
            stmpmats=stmpmats.repeat(stumpndiv).flatten()
            #Now need to get structural properties
            self.TPouts.stmpObj=Tube(Ds,ts,stmp_lgths.flatten(),self.TPinputs.Kbuck_stump,stmpmats) #This is for members of 1 leg only, CCW

            #All Stumps mass
            stmp_mass=0. #Initialize
            if stumpndiv: #since we may not have stumps
                stmp_mass=(self.TPouts.stmpObj.Area*(stmp_lgths/stumpndiv).repeat(stumpndiv)*self.TPouts.stmpObj.mat[0].rho).sum()*nlegs
                       #----------------#

            #---------Now I need to go to braces: perimeter girders + diagonals---------#
            #Find the intermediate nodes of the girders since I have already found the corner joints
            factors=np.arange(1.,girndiv)/girndiv
            dists=(np.roll(brcjoints,-1,0)-brcjoints)
            #For the lengths of the various members; we should assign the joint-joint leg to it for buckling checks;
            gird_lgths=np.sqrt( (dists**2).sum(axis=1))[0].repeat(girndiv)  #Note we just assign lengths for the members of face 1 for consistency with legs, Xbrc, Mbrc,Hbrc
            girnodes=brcjoints + dists*factors[...,np.newaxis,np.newaxis] #nodes excluding joints
            self.TPouts.girnodes=girnodes.transpose(1,0,2).reshape(-1,3)  #This puts them in the right order CCW- Girder nodes

            self.TPouts.nNodesGir=girndiv-1  #No. of nodes excluding joint with stumps; for 1 leg/face only
            #now to create the Tube object one per member:
            girdmats=Mat2Member(1,girdmatins)  #First create joint-to-joint member material
            #first: replicate Dbrcs, tbrcs appropriately to account for members
            Ds=np.tile(Dgir,girndiv) # 1 face
            ts=np.tile(tgir,girndiv) #1 face
            #then replicate equally the material objects, CCW
            girdmats=girdmats.repeat(girndiv).flatten()
            #Now need to get structural properties
            self.TPouts.girObj=Tube(Ds,ts,gird_lgths.flatten(),self.TPinputs.Kbuck_gir,girdmats) #This is for members of 1 face only, CCW

            gir_mass=(self.TPouts.girObj.Area*(gird_lgths/girndiv).repeat(girndiv)*self.TPouts.girObj.mat[0].rho).sum()*nlegs #Total mass
                       #----------------#

            #----------Now I need to go to diagonal braces-----------#
            botjnt=np.array([0.,0.,brcjoints[0,2]])  #bottom of stem
            factors=np.arange(1.,brcndiv)/brcndiv
            dists=botjnt[np.newaxis,:]-brcjoints
            self.TPouts.brcnodes=(brcjoints+dists*factors[...,np.newaxis,np.newaxis]).reshape(-1,3) #nodes excluding joints
            brc_lgths=np.sqrt( (dists**2).sum(axis=1))[0].repeat(brcndiv)  #Note we just assign lengths for the members of face 1 for consistency with legs, Xbrc, Mbrc,Hbrc

            self.TPouts.nNodesBrc=brcndiv-1  #No. of nodes excluding joint with stumps and stem; for 1 leg/face only
            #now to create the Tube object one per member:
            brcmats=Mat2Member(1,brcmatins)  #First create joint-to-joint member material
            #first: replicate Dbrcs, tbrcs appropriately to account for members
            Ds=np.tile(Dbrc,brcndiv) # 1 face
            ts=np.tile(tbrc,brcndiv) #1 face
            #then replicate equally the material objects, CCW
            brcmats=brcmats.repeat(brcndiv).flatten()
            #Now need to get structural properties
            self.TPouts.brcObj=Tube(Ds,ts,brc_lgths.flatten(),self.TPinputs.Kbuck_brc,brcmats) #This is for members of 1 leg only, CCW

            brc_mass=(self.TPouts.brcObj.Area*(brc_lgths/brcndiv).repeat(brcndiv)*self.TPouts.brcObj.mat[0].rho).sum()*nlegs #Total mass
                    #----------------#

            #----------Now I need to go to Struts-----------#
            #First find top node
            topjnt=np.array([0.,0.,brcjoints[0,2]+TPlth])
            #Then find all intermediate nodes
            factors=np.arange(1.,strutndiv)/strutndiv
            dists=topjnt[np.newaxis,:]-brcjoints
            self.TPouts.strtnodes=(brcjoints+dists*factors[...,np.newaxis,np.newaxis]).reshape(-1,3) #nodes excluding joints
            strt_lgths=np.sqrt( (dists**2).sum(axis=1))[0].repeat(strutndiv)  #Note we just assign lengths for the members of face 1 for consistency with legs, Xbrc, Mbrc,Hbrc....

            self.TPouts.nNodesStrt=strutndiv-1  #No. of nodes excluding joint with stumps and stem; for 1 leg/face only
            #now to create the Tube object one per member:
            strtmats=Mat2Member(1,strtmatins)  #First create joint-to-joint member material
            #first: replicate Dbrcs, tbrcs appropriately to account for members
            Ds=np.tile(Dstrut,strutndiv) # 1 face
            ts=np.tile(tstrut,strutndiv) #1 face
            #then replicate equally the material objects, CCW
            strtmats=strtmats.repeat(strutndiv).flatten()
            #Now need to get structural properties

            #Let us adjust the length to account for stem!!!!!
            strt_lgths2=strt_lgths-np.sqrt(Dstem[0]**2/4.+((Dstem[0]/2.*TPlth)**2) / ((dists[0,0:2]**2).sum()))
            self.TPouts.strtObj=Tube(Ds,ts,strt_lgths2.flatten(),self.TPinputs.Kbuck_strut,strtmats) #This is for members of 1 leg only, CCW

            #total struts' mass
            strt_mass=(self.TPouts.strtObj.Area[0]*strt_lgths2[0]*self.TPouts.strtObj.mat[0].rho)*nlegs
            #strt_mass=(self.TPouts.strtObj.Area*(strt_lgths/strutndiv).repeat(strutndiv)*self.TPouts.strtObj.mat[0].rho).sum()*nlegs
                    #----------------#

            #----------Now Main Stem----------#

            stemjnts=np.zeros([nstems+1,3])
            stemjnts[0,:]=botjnt

            stemjnts[1:nstems+1,2]=botjnt[2]+self.TPinputs.hstem.cumsum()
            #Store base of tower elevation w.r.t. mudline for conenctions later to tower
            self.TwrBsZ=stemjnts[-1,2]

            factors=(self.TPinputs.hstem[np.newaxis,:]*np.arange(0.,stemndiv)[:,np.newaxis]/stemndiv).transpose()
            stemnodes=stemjnts[:-1,:].repeat(stemndiv,axis=0) #initialize internal + joint nodes
            stemnodes[:,2]+=factors.flatten()
            self.TPouts.stemnodes=np.vstack((stemnodes,stemjnts[-1,:]))  #add top joint and send it to output


            nmems_stem=nstems*stemndiv
            stem_lgths=np.tile(TPlth,nmems_stem)  #Note we just assign lengths (joint-to-joint) for all elements for consistency with legs, Xbrc, Mbrc,Hbrc....

            self.TPouts.nNodesStem=nmems_stem+1  #No. of nodes including joints, total for stem
            #now to create the Tube object one per member:
            stemmats=Mat2Member(nstems,stemmatins)  #First create joint-to-joint member material
            #first: replicate Dbrcs, tbrcs appropriately to account for members
            Ds=Dstem.repeat(stemndiv)
            ts=tstem.repeat(stemndiv)
            #then replicate equally the material objects, bottom to top
            stemmats=stemmats.repeat(stemndiv)
            #Now need to get structural properties
            self.TPouts.stemObj=Tube(Ds,ts,stem_lgths.flatten(),self.TPinputs.Kbuck_stem,stemmats) #This is for members of stem, bottom to top

            #total stem's mass
            stem_mass=(self.TPouts.stemObj.Area*(self.TPinputs.hstem/stemndiv).repeat(stemndiv)*self.TPouts.stemObj.mat[0].rho).sum()

            #Calculate the TP mass
            self.TPouts.mass=stmp_mass+stem_mass+gir_mass+brc_mass+strt_mass
                    #----------------#
       else:
            pass  #to complete for other case with fortress box

        #Now store the TP concentrated mass for later in the assembling
       TPlumpIxx=np.copy(self.TPlumpinputs.I[0])
       TPlumpIyy=np.copy(self.TPlumpinputs.I[1])
       TPlumpIzz=np.copy(self.TPlumpinputs.I[2])
       if not(self.TPlumpinputs.I[0]) and self.TPlumpinputs.mass:
            lonlegs=np.sqrt(((legjnts[0,:]-legjnts[1,:])**2).sum())/nlegs
            TPlumpIxx=self.TPlumpinputs.mass/16.*lonlegs**2  #This works for nlegs=4 not sure about nlegs=3 to check
            TPlumpIyy+=TPlumpIxx
            TPlumpIzz+=2.*TPlumpIxx
       self.TPouts.TPlumpedMass=np.array([self.TPlumpinputs.mass,TPlumpIxx,TPlumpIyy,TPlumpIzz,\
                                          self.TPlumpinputs.I[3],self.TPlumpinputs.I[4],self.TPlumpinputs.I[5], \
                                          self.TPlumpinputs.CMoff[0],self.TPlumpinputs.CMoff[1],self.TPlumpinputs.CMoff[2]])

#_____________________________________________________#

class TwrGeoInputs(VariableTree):
    """Basic Geometric Inputs needed to build Tower"""
    Db =        Float(  units='m',  desc='Tower Base Diameter')
    Dt =        Float(  units='m',  desc='Tower Top Diameter')
    Htwr =      Float(  units='m',  desc='Tower Length')
    Htwr2frac = Float(  units=None, desc='Uniform X-section Tower Length as a fraction of tower height.')

    DTRb     = Float(              units=None, desc='Diameter to thickness ration at the base and below Htwr2.')
    DTRt     = Float(              units=None, desc='Diameter to thickness ration at the top.')
    DTRsdiff = Bool(True,          units=None, desc='Flag to set DTRt=DTRb, mainly used to reduce optimization design variables. Note if set to False, it will trump user-set DTRt value.')
    TwrSecH  = Float(30.,          units='m',  desc='Length of Buckling Section- it will be assumed that every submember is associated with this length.')
    Kbuck    = Float(1,            units=None, desc='Tower Kbuckling Factor.')   #Effective length factor
    Twrmatins= VarTree(MatInputs(),            desc='Tower Material Data.')
    ndiv     = Array([10],         units=None, desc='Array[1 or 2]: Number of FE elements per Tower member: two members max (uniform + tapered) ; CMzOFF Rigid Member is not included here. Note if ndiv is smaller than [2,10],it will be adjusted to a [2,10] value.')
    DeltaZmax= Float(              units='m',  desc='Maximum DeltaZ desired for the tower FE lengths. If input then ndiv will be adjusted to verify deltaz<=deltaZmax.')

    #Optional: these are used in case the user inputs stations
    ztwr     = Array( dtype=np.float, units='m',    desc='Tower stations'' heights from base of tower')
    Dtwr     = Array( dtype=np.float, units='m',    desc='Tower stations'' diameters')
    ttwr     = Array( dtype=np.float, units='m',    desc='Tower stations'' wall thicknesses')
    TwrlumpedMass = Array(np.zeros([1,11]),dtype=np.float, desc='Concentrated masses along tower: first column z''s from base of tower; 2nd through 11th column: mass and Ixx,Iyy,Izz,Ixy,Ixz,Iyz,,CMxoff,CMyoff,CMzoff values')

# TwrGeoOutputs is in VarTrees.py
##class TwrGeoOutputs(VariableTree):
##    """Basic Geometric Outputs needed to build Tower of Jacket"""
##    TwrObj   =  Instance(Klass=Tube,                             desc='Object of Class Tube for Tower portion of Jacket')
##    Twr2RNAObj =Instance(Klass=RigidMember,                      desc='Object of Class RigidMember for Rigid Tower portion')
##    joints   = Array(np.array([]),units='m', dtype=np.float, desc='Pile Joint (with legs) Coordinates (3,nlegs=nfaces)')
##    nodes    = Array(np.array([]),units='m', dtype=np.float, desc='Tower ALL Nodes'' Coordinates (3,nNodes)')
##    nNodes   = Int(0,           units=None,                  desc='Number of Tower Nodes INCLUDING joint at TP')
##    mass =     Float(           units='kg',                  desc='Tower Mass')
##    TopMass  = Array(np.zeros([4]),         dtype=np.float, desc='Tower Top mass, Ixx, Iyy, Izz from RNA properties in input')
##    HH       = Float(            units='m', desc='Hub-Height')

class Tower(Component):
    """Tower model.  Assumes: \n
    EITHER  \n
        - Constant taper in both diameter and thickness \n
        - A uniform section at the base is allowed \n
    OR \n
        -z,D,t given by user at desired end nodes of elements. \n

    This component returns nodes and tube object. mass calculation MUST BE REVISED to allow for multiple materials. \n
    """
    # inputs
    Twrins  =VarTree(TwrGeoInputs(), iotype='in', desc='Basic Input data for Tower')
    HH       = Float(     units='m', iotype='in', desc='Hub-Height=Htwrbase+Htwr+CmZoff. Could be input in lieu of Htwr')
    RNAinputs=VarTree(RNAprops(), iotype='in', desc='Basic Inertial Properties of RNA')
    RigidTop=Bool(True,  iotype='in', desc='Connection To RNA to be considered Through a RigidMember(TRUE) or mathematically (FALSE)')
    nlegs   =Int(4,iotype='in', units=None, desc='Number of Legs')
    legjnts =Array(units='m', iotype='in', dtype=np.float,  desc='Leg Joint Coordinates')

    TwrBsZ =Float(0.,  iotype='in', units='m', desc='Elevation w.r.t. mudline of tower base joint (=TP top node z coordinate)')
    wdepth =Float(     iotype='in', units='m', desc='Water Depth')

    #outputs
    Twrouts  =VarTree(TwrGeoOutputs(), iotype='out', desc='Basic Output data for Tower')
    Dt=       Float(     iotype='out', units='m', desc='TowerTop OD used later by other components.')
    tt=       Float(     iotype='out', units='m', desc='TowerTop wall thickness used later by other components.')

    def execute(self):
        #Simplify nomenclature
        yawangle=self.RNAinputs.yawangle
        CMzoff=np.copy(self.RNAinputs.CMoff[2])
        Thzoff  =self.RNAinputs.Thoff[2]   #Float(units='m',    desc='Distance from hub centerline to Tower top along z. If left blank, it will be set = RNAInputs.CMzoff' )
        Db=self.Twrins.Db
        Dt=self.Twrins.Dt
        tb = Db/self.Twrins.DTRb
        if not(self.Twrins.DTRsdiff):
            self.Twrins.DTRt=self.Twrins.DTRb

        tt = Dt/self.Twrins.DTRt
        Htwr=self.Twrins.Htwr
        Htwr2frac=self.Twrins.Htwr2frac
        DeltaZmax=self.Twrins.DeltaZmax
        ztwr=self.Twrins.ztwr
        Dtwr=self.Twrins.Dtwr
        ttwr=self.Twrins.ttwr


        usrstations=False  #Initialize

        if Htwr2frac>1. or Htwr2frac<0.:
            sys.exit('Error: Htwr2frac >1 or Htwr2frac <0: Not allowed!!!')

        TwrBsZ=self.TwrBsZ
        TwrSecH=self.Twrins.TwrSecH
        ndiv=self.Twrins.ndiv
        Twrmatins =self.Twrins.Twrmatins
        Twrmats=Mat2Member(ndiv.size,Twrmatins)  #First create joint-to-joint member material


        RigidTop=(CMzoff != 0.) and self.RigidTop  #Add an extra joint and member only if rigid member requested and CMzoff<>0

        CMzoff=CMzoff*RigidTop  #From now on I do not need CMzoff if RigidTop so it would be set to 0

        #Calculate HH if not given as input, or Htwr in case
        self.Twrouts.HH=self.HH

        nNodes=ztwr.size  #Initialize also to check whether stations have been assigned already

        if not Thzoff:  #INitialize this one if absent
            Thzoff=self.RNAinputs.CMoff[2]

        if not(Htwr) and not(self.Twrouts.HH) and not(nNodes):  #Htwr not specified, HH is specified
            sys.exit('!!!You must specify either HH and Thzoff, or Htwr!!!')
        elif nNodes: #Stations are given, thus HTwr is given, need to check HH
            Htwr=ztwr[-1]-ztwr[0] #CJB I believe this line is OK as is
        elif not(Htwr):  #HH is specified
                Htwr=self.Twrouts.HH-Thzoff-TwrBsZ #+self.wdepth #self.RNAinputs.CMzoff
                #CJBe Comment out self.wdepth from the above line to consider z=0 at MSL.
        elif not(self.Twrouts.HH):            #Htwr is specified
            self.Twrouts.HH=TwrBsZ+Htwr+Thzoff #+self.wdepth #+self.RNAinputs.CMzoff
            #CJBe Comment out self.wdepth from the above line to consider z=0 at MSL.

        if abs(Htwr-(self.Twrouts.HH-Thzoff-TwrBsZ)) >0.00001: #Both are specified but they do not jive -self.RNAinputs.CMzoff #CJB- Originally: (self.Twrouts.HH-Thzoff-TwrBsZ+self.wdepth) Remove +self.wdepth
            sys.exit('HH and HTwr incompatible: HH=Thzoff+Htwr+TwrBsZ')

        Htwr2=Htwr2frac*Htwr
        Htwr2Flg=(Htwr2 != 0.)

        njs=2 +Htwr2Flg  +int(RigidTop)  #Number of joints for Tower assumed at base, top, HTwr2 and CMzoff


        if nNodes:
            usrstations=True
            if nNodes==Dtwr.size and nNodes==ttwr.size:  #This means user has input tower geometry via stations already

                print 'User-input stations will be used: Db,Dt,DTRb,DTRt are trumped!'

                Dt=Dtwr[-1]
                tt=ttwr[-1]


                self.Twrouts.nNodes = ((nNodes-1)*ndiv[-1] +1) +int(RigidTop)   #This includes TP top point joint
                self.Twrouts.nodes = np.zeros([3,self.Twrouts.nNodes])   #Initialize
                self.Twrouts.nodes[2,:] = TwrBsZ+ztwr

                self.Twrouts.nElems=self.Twrouts.nNodes - int(RigidTop) -1  #Number of elements in the flexible portion of tower

                self.Twrouts.nodes[2,-1]=self.Twrouts.nodes[2,-int(RigidTop)-1]+CMzoff   #In case put CMzoff at the top, if CMzoff=0 it does not change anything, good way to eliminate if

                #Tube properties, lower node properties
                Ds=Dtwr[:-1]
                ts=ttwr[:-1] #thicknesses non-rigid members Ni nodes

                #Now store the Joints separately, in this case they coincide with nodes
                self.Twrouts.joints=np.copy(self.Twrouts.nodes)
                njs=self.Twrouts.nNodes
                            #materials
                Twrmats2=Mat2Member(self.Twrouts.nElems,Twrmatins)  #First create joint-to-joint member material

                #Now store concentrated masses
                self.Twrouts.TwrlumpedMass =self.Twrins.TwrlumpedMass
                self.Twrouts.TwrlumpedMass[:,0] += TwrBsZ

            else:
                sys.exit('Inconsistencies found in ztwr,Dtwr, and ttwr: revise tower stations !!')

        else:                                      #______________This means user has input paramteric values_______________#

            #Node coordinates starting from deltazs
            if Htwr2Flg:
                Htwrtop=Htwr-Htwr2 #tapered length
                if DeltaZmax: #if this is input <>0 make sure you do not exceed it, first at the bottom
                    ndiv=np.max([ndiv,np.ceil([Htwr2/DeltaZmax,Htwrtop/DeltaZmax])],0)
                ndiv=np.max([ndiv,np.array([2,10])],0).astype(int)  #make sure we have at least 2 and 10 elements in the uniform and tapered sections
                #Now make sure there is not too much difference between deltaz in Htwr2 and htwr (1:3 max)
                dzbot=Htwr2/ndiv[0]
                dzbot_min=0.33*Htwrtop/ndiv[1]
                if  dzbot<dzbot_min:
                    ndiv[0]=np.ceil(Htwr2/dzbot_min)
                deltaZs=np.hstack((np.array([Htwr2/ndiv[0]]).repeat(ndiv[0]), np.array([Htwrtop/ndiv[1]]).repeat(ndiv[1]))).flatten()
            else:
                #This is to bypass the issue where the user might have input ndiv(2), although it should be ndiv(1), so ndiv[1] rules
                ndiv[0]=ndiv[-1]=np.max([ndiv[-1],np.ceil(Htwr/DeltaZmax)])
                ndiv=np.max([ndiv,np.array([10]).repeat(ndiv.size)],0).astype(int)  #make sure we have at least 10 elements in the tapered section (only section), we still assume it is a 2-element array
                deltaZs=np.array([Htwr/ndiv[-1]]).repeat(ndiv[-1]).flatten()

            deltaZcum=np.hstack((0.,deltaZs.cumsum()))

            self.Twrouts.nNodes= np.dot(ndiv,np.array([Htwr2Flg,1]))+1 +int(RigidTop)   #This includes TP top point joint

            self.Twrouts.nodes=np.zeros([3,self.Twrouts.nNodes])   #Initialize

            self.Twrouts.nodes[2,0:self.Twrouts.nNodes-int(RigidTop)]=  TwrBsZ + deltaZcum#z coordinates


            self.Twrouts.nElems=self.Twrouts.nNodes - int(RigidTop) -1  #Number of elements in the flexible portion of tower

            self.Twrouts.nodes[2,-1]=self.Twrouts.nodes[2,-int(RigidTop)-1]+CMzoff   #In case put CMzoff at the top, if CMzoff=0 it does not change anything, good way to eliminate if


            #Now create Tube objects: we are interpolating the real diameter between the node points, getting the average and live with it for each element

            Ds=np.zeros(self.Twrouts.nElems)*np.NaN #diameters at the start of the actual non-rigid members Ni nodes
            ts=Ds.copy() #thicknesses non-rigid members Ni nodes

            idxtwr2=np.nonzero(deltaZcum <= Htwr2)[0] #indices of constant D tower portion .[0] to make it into an array
            idx=range(idxtwr2[-1]+1,self.Twrouts.nElems) #indices above last section of tower with D=Db

            Ds[idx]=Db-(Db-Dt)*(self.Twrouts.nodes[2,idx]-self.Twrouts.nodes[2,0]-Htwr2)/ \
                        (Htwr-Htwr2)
            Ds[idxtwr2]=Db

            ts[idx]=tb-(tb-tt)*(self.Twrouts.nodes[2,idx]-self.Twrouts.nodes[2,0]-Htwr2)/ \
                        (Htwr-Htwr2)
            ts[idxtwr2]=Db/self.Twrins.DTRb

            print '{:d} nodes in the constant-OD segment of the tower'.format(idxtwr2.size)

            #Now store the Joints separately
            self.Twrouts.joints=np.zeros([3,njs]) #Initialize
            self.Twrouts.joints[2,1:njs]=TwrBsZ+ np.array([Htwr2,Htwr,RigidTop*(Htwr+CMzoff)]).T.nonzero()[0]  #joints

            #then replicate equally the material objects
            Twrmats2=np.empty(self.Twrouts.nElems,dtype=object)

            for ii in range(0,Htwr2Flg+1):
                idx1=ndiv.cumsum()[ii]
                idx0=idx1-ndiv[ii]
                Twrmats2[idx0:idx1]=np.tile(Twrmats[ii],ndiv[ii])    #note the usual repeat as used in legs for instance, does not work here, since I have different ndiv per jnt2jnt member

        #Now need to get structural properties
        self.Twrouts.TwrObj=Tube(Ds,ts,np.array([TwrSecH]).repeat(self.Twrouts.nElems),self.Twrins.Kbuck,Twrmats2) #Initialize just like the parent class D,t, etc., then add a few more attributes

        self.Twrouts.Twr2RNAObj=Tube(np.array([]),np.array([]),np.array([]),self.Twrins.Kbuck,np.array([])) #Initialize

        #Now store the RNA concentrated mass for later in the assembling
        self.Twrouts.TopMass=np.array([self.RNAinputs.mass,self.RNAinputs.I[0],self.RNAinputs.I[1],self.RNAinputs.I[2],\
                                       self.RNAinputs.I[3],self.RNAinputs.I[4],self.RNAinputs.I[5],\
                                       self.RNAinputs.CMoff[0],self.RNAinputs.CMoff[1],self.RNAinputs.CMoff[2]])

        #Note account for position of RNA for weight effects onto loading, not onto frequency
        InerTensor=np.zeros([3,3])
        InerTensor[0,:]=np.array([self.RNAinputs.I[0],self.RNAinputs.I[3],self.RNAinputs.I[4]])
        InerTensor[1,:]=np.array([self.RNAinputs.I[3],self.RNAinputs.I[1],self.RNAinputs.I[5]])
        InerTensor[2,:]=np.array([self.RNAinputs.I[4],self.RNAinputs.I[5],self.RNAinputs.I[2]])
        DIRCOS=np.array([ [cosd(yawangle) ,    -sind(yawangle), 0.],
                          [    sind(yawangle) , cosd(yawangle), 0.],
                          [    0.,                         0.,  1.]])
        InerTensor=np.dot(np.dot(DIRCOS,InerTensor),DIRCOS.T)

        self.Twrouts.TopMass_yaw=self.Twrouts.TopMass.copy()  #This accounts for yaw angle
        self.Twrouts.TopMass_yaw[1:7]=np.array([InerTensor[0,0],InerTensor[1,1],InerTensor[2,2],\
                                             InerTensor[0,1],InerTensor[0,2],InerTensor[1,2]])
        self.Twrouts.TopMass_yaw[7:]=np.dot(DIRCOS, np.array([self.RNAinputs.CMoff[0],self.RNAinputs.CMoff[1],self.RNAinputs.CMoff[2]]))

        self.Twrouts.Thoff_yaw=np.dot(DIRCOS, np.array([self.RNAinputs.Thoff[0],self.RNAinputs.Thoff[1],self.RNAinputs.Thoff[2]]))  #THIS SHOULD GO ELSEWHERE, but it is convenient since I ma rotating stuff
        self.Twrouts.rna_yawedcm=np.copy(self.Twrouts.TopMass_yaw[7:])
        #RIGID MEMBER
        if RigidTop:
            self.Twrouts.Twr2RNAObj=RigidMember(D=Dt,t=tt,L=self.RNAinputs.CMoff[2])
            self.Twrouts.TopMass[9]=0.  #Note remove CMzoff in case rigid member (should also be cmx,yoff but those are not used in teh rigid member yet)
            self.Twrouts.TopMass_yaw[9]=0.  #Note remove CMzoff in case rigid member (should also be cmx,yoff but those are not used in teh rigid member yet)

        #Store a couple of things that may be needed by other components: note top OD and t are not known elsewhere otherwise
        self.Dt=Dt
        self.tt=tt
        self.Twrouts.Htwr=Htwr
        self.Twrouts.Htwr2=Htwr2

        # Calculate mass
        # volume of constant section
        Vunif = Htwr2 * pi/4*(Db**2 - (Db-2*tb)**2)  #Uniform Section
        # volume of conical frustrum
        V0, CM0  = frustum(Db, Dt, Htwr-Htwr2)
        V01,CM01 = frustum(Db-2*tb, Dt-2*tt, Htwr-Htwr2)
        Vfrst = V0 - V01  #Frustum
        # mass
        self.Twrouts.mass = Twrmats[0].rho*Vunif + Twrmats[-1].rho*Vfrst

        if usrstations:
            deltaZs=(np.roll(ztwr,-1)-ztwr)[:-1]
            V0s, CM0s  = frustum(Dtwr[:-1], Dtwr[1:], deltaZs)
            V01s,CM01s = frustum( (Dtwr-2.*ttwr)[:-1], (Dtwr-2.*ttwr)[1:], deltaZs)
            Vfrst = V0s - V01s  #Frustum
            rhos=np.array([mat.rho for mat in self.Twrouts.TwrObj.mat])
            self.Twrouts.mass = (rhos*Vfrst).sum()

#_____________________________________________________#
class SoilGeoInputs(VariableTree):
    """Basic Soil Properties needed to assess Soil and Embedment Length"""
    zbots   = Array(units='m',     dtype=np.float, desc='Layers'' bottom z-coordinate from mudline')
    gammas  = Array(units='N/m**3',dtype=np.float, desc='Layers'' Unit Weight')
    cus     = Array(units='N/m**2',dtype=np.float, desc='Layers'' Undrained Shear Strength')
    phis    = Array(units='deg',   dtype=np.float, desc='Layers'' Friction Angles')
    delta   = Float(units='deg',                   desc='Pile-soil friction angle')
    sndflg  = Bool(units=None, desc='Flag indicating whether it is a cohesionless (sand,True) or cohesive (clay, False) soil.')
    SoilSF  = Float(  units=None,    iotype='in', desc='Safety factor for soil (affects pile capacity): From ABS BOWTI and API RP2A use 2.0')  #Changed from 1.25 that was given to us by Keystone, on 4/16/2015
    PenderSwtch= Bool(False, units=None, desc='Flag indicating whether Pender''s or MAtlock&Reese (Default) method is used.')
    plug    = Bool(False, units=None, desc='Flag indicating whether plugged or unplugged pile.')

class SoilGeoOutputs(VariableTree):
    """Basic Soil-Pile Stiffness Data"""
    Kx=       Float(units='N/m', desc='Lateral Stiffness Equivalent Spring Constant')
    Ktheta=   Float(units='N*m/rad', desc='Rotational Stiffness Equivalent Spring Constant')

class Soil(Component):
    """This component creates soil object."""

    #inputs
    SoilIns =VarTree( SoilGeoInputs(), iotype='in', desc='Basic Soil Parameters')

    #outputs
    Soilouts=VarTree(SoilGeoOutputs(), iotype='out', desc='Basic output data for Soil')
    SoilObj =Instance(Klass=SoilC        , iotype='out', desc='Object of Class Soil')

    def execute(self):  #Should add Andrew's stiffness calculation program
        self.SoilObj=SoilC(zbots=self.SoilIns.zbots,gammas=self.SoilIns.gammas,cus=self.SoilIns.cus,\
                           phis=self.SoilIns.phis,delta=self.SoilIns.delta, sndflg=self.SoilIns.sndflg,\
                           PenderSwtch=self.SoilIns.PenderSwtch,SoilSF=self.SoilIns.SoilSF,plug=self.SoilIns.plug)


#_____________________________________________________#
class PileGeoOutputs(VariableTree):
    """Basic Geometric Outputs needed to build Piles of Jacket"""
    PileObj     = Instance(Klass=Tube,                          desc='Object of Class Tube for Non-AF portion of Pile')
    AFObj       = Instance(Klass=Tube,                          desc='Object of Class Tube for AF portion of Pile')
    joints      = Array(np.array([]),units='m', dtype=np.float, desc='Pile Joint (with legs) Coordinates (3,nlegs=nfaces)')
    nodes       = Array(np.array([]),units='m', dtype=np.float, desc='Piles'' intermediate Coordinates for nodes (3,ndiv-1,nlegs=nfaces)')
    nNodesPile  = Int(0,             units=None,                desc='Number of Pile Nodes (1 Face) not counting joints at legs')
    #Dpileout    = Float(units='m', desc='Pile Outer Diameter revisted for battered piles')    #TO DO I would like this one to take Dleg[0] value if not given, or be adjusted for size if VPFlag=false
    #tpileout    = Float(units='m', desc='Pile Wall Thickness passed through or assumed = leg bottom t')     #TO DO I would like this one to take tleg[0] value if not given

class AFOutputs(VariableTree):
    """Apparent Fixity Output Data"""
    L_AF=      Float(0.,   iotype='out',units='m',    desc="Pile Apparent Fixity Length")
    Jxx_AF=    Float(      iotype='out',units='m**4', desc="Pile Apparent Fixity Jxx Area Moment of Inertia")
    t_AF=      Float(0.05, iotype='out',units='m',   desc='Apparent Fixity wall thickness')
    AFObj   = Instance(Klass=Tube,iotype='out',    desc='Object of Class Tube')

class AFixity(Component):
    """FROM BUILDMPTWR.PY: This 'function/component' calculates the apparent fixity parameters for the pile \n
    given the soil conditions and selected radius of gyration. It augments the \n
    object self (PileC) as well.  IT assumes we keep the pile 1st diameter\n
    and we play with length and wall  thickness.  Material density will need to \n
    changed as well to realize a weightless member. The cross section is a tubuluar one. \n
    INPUT    \n
    soil   - object of soil class

    """
    #inputs
    Soilinputs=VarTree(SoilGeoOutputs(),iotype='in',desc="Soil-pile-interaction Stiffness Data")
    MoverF=      Float(iotpe='in',units='1/m',desc='Ratio of Bending Moment to Force at Pile Head')
    r_gyr=       Float(iotpe='in',units='m',  desc='Embedment Pile Xsection Radius of Gyration')
    AFmatin=VarTree(MatInputs(),iotpe='in',   desc="Apparent Fixity Material Input Data")
    EmbedMass=   Float(iotpe='in', units='kg',desc='Embedment Pile Mass')
    D_AF=        Float(iotpe='in', units='m', desc='Apparent Fixity OD, usually equal to Embedment Pile Diameter')
    ndiv=        Int(1,iotpe='in', units=None,desc='Number of FE elements per AF member')

    #outputs
    AFouts=VarTree(AFOutputs(),iotype='out',units='m',    desc="Pile Apparent Fixity Data")

    def execute(self):
        Ktheta=self.Soilinputs.Ktheta
        Kx=self.Soilinputs.Kx
        ndiv=self.ndiv
        AFmats=Mat2Member(ndiv,self.AFmatin)

        L_AF=np.roots([1./3., 0.5*(self.MoverF-Ktheta/(Kx*self.MoverF)), -Ktheta/Kx])
        self.AFouts.L_AF=L_AF[L_AF>0][0]
        self.AFouts.Jxx_AF=self.AFouts.L_AF*Ktheta / self.MoverF * (self.MoverF+self.AFouts.L_AF/2.) / self.AFmatin.E

        Area_AF=self.AFouts.Jxx_AF / self.r_gyr**2

        self.AFouts.t_AF=np.roots([1., -self.AFouts.D_AF, 2./pi*Area_AF])
        self.AFouts.t_AF=self.t_AF[self.AFouts.t_AF>0][0]
        #Now adjust density to make sure we get the right apparent fixity mass
        AFmats.rho=self.EmbedMass/(Area_AF*self.AFouts.L_AF)
        self.AFouts.AF_Obj=Tube(self.AF_D,self.AFouts.AF_t,self.AFouts.L_AF,1.,AFmats)  #Kbuck=1. hardwired

#_____________________________________________________#
class PileGeoInputs(VariableTree):
    """Basic Geometric Inputs need to build Legs of Jacket"""
    Dpile    = Float(      units='m', desc='Pile Outer Diameter')    #TO DO I would like this one to take Dleg[0] value if not given, or be adjusted for size if VPFlag=false
    tpile    = Float(      units='m', desc='Pile Wall Thickness')     #TO DO I would like this one to take tleg[0] value if not given
    pileZtop = Float(      units='m', desc='Pile Head Elevation from Sea-Floor')
    Lp       = Float(      units='m', desc='Pile Embedment Length.')  #We are now putting Lp as an input and it will be verified in the optimizer for capacity.
    Kbuck    = Float(0.8,  units=None,desc='Pile Kbuckling Factor')   #Effective length factor
    Pilematins=VarTree(MatInputs(),   desc="Pile Material Input Data")
    ndiv     = Int(        units=None,desc='Number of FE elements per Pile')

class Piles(Component):
    """This component creates nodes for the Piles.
       For now just vertical piles have been verified.
       Battered piles would also be modifying properties inside leg bottom, which has not been done yet."""

    #inputs
    nlegs =Int(           iotype='in', units=None, desc='Number of Legs')     #TO DO : Need to figure out to get this one from JcktINput
    Pileinputs=VarTree(PileGeoInputs(),iotype='in',desc="Pile Input Data")
    LegbD    = Float(     iotype='in',units='m', desc='Leg Bottom Section Outer Diameter')
    Legbt    = Float(     iotype='in',units='m', desc='Leg Bottom Section Wall Thickness')
    VPFlag   = Bool(      iotype='in',units=None,desc='Vertical Pile Flag [Y/N]: If True the Mudbrace is put at the expected joint with Piles, i.e. at the bottom of leg.')
    AFflag   = Bool(      iotype='in',units=None,desc='Apparent Fixity Flag: if True AF is activated. Note if AFflag=True in piles, then clamped will be reset to True.')

    AFouts=VarTree(AFOutputs(), iotype='in',      units='m',  desc="Pile Apparent Fixity Data")
    legjnts   =Array( iotype='in', dtype=np.float,units='m',  desc='Leg Joint Coordinates')  #From Legs somehow

    #outputs
    Pileouts=VarTree(PileGeoOutputs(), iotype='out', desc='Basic output data for Piles')
    PileObjout = Instance(Klass=Tube,      iotype='out', desc='Object of Class Tube for Non-AF portion of Pile, valid also when ndiv=0 to calculate Lembed')
    PileFlag = Bool(      iotype='out',units=None,desc='Piles [Y/N] Flag, to be used by Mudbrc component to figure out where to put mudbrace. It is true if ndiv!=0.')

    def execute(self):
       ndiv=self.Pileinputs.ndiv
       self.PileFlag= (ndiv!=0)

       nlegs=self.nlegs
       VPFlag=self.VPFlag
       AFflag=self.AFflag

       Dpile=self.Pileinputs.Dpile
       tpile=self.Pileinputs.tpile
       if not Dpile: #empty data, assign leg data first
            Dpile=self.LegbD
       if not tpile: #empty data, assign leg data first
            tpile=self.Legbt

       if not(VPFlag): #battered piles
          #force Dp to follow Dleg in this case
         Dpile=np.min([Dpile, self.LegbD-2.*self.Legbt-3.5*0.0254])

       if Dpile !=self.Pileinputs.Dpile:
                    warnings.warn('Pile OD has been modified to account for Leg ID, Dp= {:f} [m]'.format(Dpile))


       #store something in case pile is not used in the model at all, so that we can still calculate Lembed
       self.PileObjout=Tube(Dpile,tpile,[],self.Pileinputs.Kbuck,Mat2Member(1,self.Pileinputs.Pilematins)) #This is for 1 member only, for non-AF portion of pile, to calculate Lembed

       if not(ndiv) : #so if user is requesting no piles at all, just return default empty stuff
           self.Pileouts.nNodesPile=0
           self.Pileouts.nodes=np.empty([nlegs,0,3])
           self.Pileouts.PileObj=Tube(np.array([]),np.array([]),np.array([]),self.Pileinputs.Kbuck,np.array([]))
           self.Pileouts.AFObj=self.Pileouts.PileObj

       else:
           #Simplify nomenclature
           legbotjnts=self.legjnts[:,:,0].transpose()  #joints with legs
           legZbot=legbotjnts[0,2]



           L_AF=self.AFouts.L_AF #initialize

           self.Pileouts.AFObj=Tube(np.array([]),np.array([]),np.array([]),self.Pileinputs.Kbuck,np.array([])) #Initialize
           if AFflag: #calculate Apparent Fixity
                self.Pileouts.AFObj=self.AFouts.AF_Obj  #This returns it only if component has been run

           #Let us start by finding nodes
           njoints=1+AFflag+(legZbot>0.) #no. of joints per pile:joints with legs, at sea bed and potentially below through AF
           pilejnts=np.zeros([nlegs,njoints,3])  #joints with legs, at sea bed and potentially below through AF
           pilejnts[:,-1,:]=legbotjnts

           if VPFlag:  #vertical piles; also enforce Dpile=Dleg at one point
                pilejnts=legbotjnts.repeat(njoints,axis=0)
                dists=-np.array([0.,legZbot,L_AF])[:njoints].cumsum()[::-1].reshape([1,-1])  #from bottom up
                pilejnts[:,2]+=dists.repeat(nlegs,0).flatten()  #from bottom up, then leg 1 through 4
                pilejnts=pilejnts.reshape([nlegs,njoints,3])
           else:

                #find direction for joints from leg1; also enforce Dpile=Dleg-2*tleg-3.5in
                dists=-(self.legjnts[:,:,1].transpose()-legbotjnts)  #vectors showing directions for pile joints, 1 per leg, pointing down
                delta=np.sqrt((dists**2).sum(axis=1))[0]  #distance between legbotjnts and joint above it on the leg
                deltaz=-dists[0,2]                         #distance between legbotjnts and joint above it on the leg along z; "-" so is positive
                dists2=np.array([0.,legZbot*delta/deltaz,L_AF])[:njoints].cumsum()[::-1].reshape([-1])  #from bottom up, distances along leg direction
                pilejnts=legbotjnts.repeat(njoints,axis=0).reshape([nlegs,njoints,3]) + dists[:,np.newaxis]*dists2[:,np.newaxis]/delta

           #Now find nodes
           botjnts=pilejnts[:,:-1,:]  #joints excluding leg joint
           factors=(np.arange(0.,ndiv)/ndiv).reshape([ndiv,1])  #reshape for later multiplication
           dists=(np.roll(pilejnts,-1,axis=1)-pilejnts)[:,0:-1,:] #vectors from joint to joint

           self.Pileouts.nodes=botjnts[:,:,np.newaxis]+dists[:,:,np.newaxis]*factors[np.newaxis,:]  #These include the AF and or bottom joints, intermediate joints, but not joints with legs; all legs

           pile_lgths=np.sqrt( (dists**2).sum(axis=2))[0].repeat(ndiv)  #Note we just assign lengths for the members of face 1 for consistency with legs, Xbrc, Mbrc,Hbrc....

           self.Pileouts.nNodesPile=(njoints-1)*ndiv  #No. of nodes excluding joint with legs, but including AF in case; for 1 leg/face only

                #now to create the Tube object one per member:
           pilemats=Mat2Member(njoints-1-AFflag,self.Pileinputs.Pilematins)  #First create joint-to-joint member material, for the non-AF portion only here:

                #first: replicate Dpile, tpile appropriately to account for elements
           Ds=np.tile(Dpile,ndiv) # 1 face
           ts=np.tile(tpile,ndiv) #1 face
                #then replicate equally the material objects, CCW
           pilemats=pilemats.repeat(ndiv).flatten()
                #Now need to get structural properties
           self.Pileouts.PileObj=Tube(Ds,ts,pile_lgths.flatten(),self.Pileinputs.Kbuck,pilemats) #This is for members of 1 leg only, for non-AF portion of pile

#________________________________________________#

class PreJcktBuild(Component):
    """Small Component Helper that takes some input from PrepJcket and Leg and Creates
    some other inputs for components: JcktBuild and Legs.

    Jacket model.  Assumes:
           BuildJacket() really does the building, here we set general and default main parmaeters.
           Additionally, the following assumptions apply: \n
        0. Round tube jacket construction only\n
        1. 3 or 4 legs\n
        2. X-braces only\n
        3. All joints are moment-coonections (weldments)\n
        4. Shear Deformation (Timoshenko) included\n
        5. The joint finite size is ignored for the time being, it will be modified in the future\n
        6. Properties are constant within a bay (D,t for legs and braces)\n
        7. The input file for frame3DD is in the matlab variant style (see Frame3DD documentation)\n
        8. Jacket leg has a stump about 1.5D long (vertical height) and at the top\n
           it connects to
        9. Legs, TP are created in BuildJacket to give a chance to read input from a file
        10. Constant brace angle with chord only
        11. If legZbot>0, a vertical pile stump to sea-bed is created with same D and t as leg foot as default
    """
    # inputs
    JcktGeoIn = VarTree(JcktGeoInputs(), iotype='in', desc='Geometric Inputs for BuildJckt')
    #RNAgeo    = VarTree(RNAprops(),     iotype='in', desc='Geometric Inputs for BuildJckt')
    #wdepth   = Float( units='m',iotype='in',   desc='Water Depth') #CJB- water flag   #TO DO, this should really come from environment
    #wlevel  = Float( units='m',iotype='in',   desc='Water Level') #CJB- water flag  #TO DO, this should really come from environment
    WaterNew  = VarTree(WaterInputs(),  iotype='in', desc='Water Data') #CJB water flag

    LegbD     =Float(    iotype='in',           desc="Leg Base OD")
    LegtD     =Float(    iotype='in',           desc="Leg Top OD")
    Legtt     =Float(    iotype='in',           desc="Leg Top t")

    legZbot   = Float(0.,iotype='in',units='m', desc='Bottom Elevation of Leg from Sea-Floor')
    legbot_stmphin =Float( iotype='in',units='m', desc='Bottom Leg Stump Length') #User's input

    TwrDb         =Float(units='m',     iotype='in', desc='Tower Base OD')
    TwrDTRb       =Float(units='m',     iotype='in', desc='Tower Base DTR')

    TPinputs=VarTree(TPGeoInputs(),     iotype='in', desc='Basic input data for Transition Piece')
    PreBuildTPLvl=Int(    units=None,   iotype='in', desc="""Level of Prebuild [1-5]: 1=Stump & Strut from leg; 2[default]=1+(Gir=diagBrace=Xbrace); \n
                                                              3=(Gir=Brace=Xbrace); 4=diagBrace=girders. 5=1+4    in all cases: Stem from Twr; """)

    Xbrcinputs=VarTree(XBrcGeoInputs(),iotype='in',desc="Xbrace Input Data. They will be slightly modified here.")

    # outputs: some parameters for geometry of the jacket
    legbot_stmph =Float(   iotype='out',  units='m', desc='Bottom Leg Stump Length')  #Updated legbot_stmph : This should be set to 1.5 Dleg[0] if user did not set it
    dck_width   =Float(    iotype='out',   units='m',  desc='Deck Width') # This is twice the Tower Based Db

    al_bat2D  =Float(    iotype='out', units='rad',desc='Batter Angle in 2D, angle between vertical and leg in projection (2D) [rad]')
    al_bat3D  =Float(    iotype='out', units='rad',desc='Batter Angle in 3D, used later for loading')
    beta3D    =Float(    iotype='out', units='rad',desc='3D angle between leg and xbrace to compare to NORSOK, in rad')
    beta2D    =Float(    iotype='out', units='rad',desc='2D angle (projection) between leg and xbrace, in rad')
    innr_ang  =Float(pi/4., iotype='out',units='rad',desc='Angle between radial direction and base side, in rad')
    bay_hs    =Array(dtype=np.float, iotype='out',units='m',  desc='Bay Lenghts')
    bay_bs    =Array(dtype=np.float, iotype='out',units='m',  desc='Bay Base Widths')
    JcktH     =Float(       iotype='out',units='m',  desc='Jacket height from leg bottom to leg top')
    wbas0     =Float(       iotype='out',units='m',  desc='Jacket width at bottom of first bay')
    wbase     =Float(       iotype='out',units='m',  desc='Jacket virtual width at mudline')

    #TP outputs
    BuildTPouts  =VarTree(TPGeoInputs(),  iotype='out', desc='Revised Input data for Transition Piece')
    Xbrcouts     =VarTree(XBrcGeoInputs(),iotype='out',desc="Revised Xbrace Input Data.")

    def execute(self):
        """This function returns main geometry of the Jacket. Node coordinates and connectivity
          This function builds up a model for Frame3DD based on Jacket object.
          It requires an object as input and spits out an augmented version of it.

        !!!!soil and apparent fixity still to fix!!!!"""
        if self.WaterNew.wlevel==0 and self.WaterNew.wdepth!=0: #CJB+
            self.WaterNew.wlevel=self.WaterNew.wdepth #CJB+
        if self.WaterNew.z_floor==0 and self.WaterNew.wdepth!=0: #CJB+
            self.WaterNew.z_floor=-self.WaterNew.wdepth #CJB+

        nbays=self.JcktGeoIn.nbays

        #X-brace precalculations: Next two useful in case Dbrc0 is the only optimization variable in optimization problems
        self.Xbrcouts=self.Xbrcinputs #initialize

        if self.Xbrcinputs.Dbrc0:
          self.Xbrcouts.Dbrc=np.array([self.Xbrcinputs.Dbrc0]*nbays)
        if self.Xbrcinputs.tbrc0:
           self.Xbrcouts.tbrc=np.array([self.Xbrcinputs.tbrc0]*nbays)

        if any(self.Xbrcouts.Dbrc == 0.) :
            self.Xbrcouts.Dbrc=self.Xbrcouts.Dbrc[0].repeat(self.Xbrcinputs.Dbrc.size)
        if any(self.Xbrcouts.tbrc == 0.) :
            self.Xbrcouts.tbrc=self.Xbrcouts.tbrc[0].repeat(self.Xbrcinputs.tbrc.size)

        #TP precalculations
        self.BuildTPouts=self.TPinputs #here I transfer all the other inputs before expanding them with new info
        PreBuildLvl=self.PreBuildTPLvl #shorten name

         #Set central stem following tower base and bump up thikness
        if PreBuildLvl:
            self.BuildTPouts.Dstem=np.asarray([self.TwrDb]).repeat(self.TPinputs.hstem.size)
            self.BuildTPouts.tstem=np.asarray([self.TwrDb/self.TwrDTRb*self.TPinputs.stemwallfactor]).repeat(self.TPinputs.hstem.size)

        if PreBuildLvl==1 or PreBuildLvl==2 or PreBuildLvl==5:
            #Set struts equal to top of leg
            self.BuildTPouts.Dstrut=self.LegtD
            self.BuildTPouts.tstrut=self.Legtt
            self.BuildTPouts.Dstump=self.LegtD
            self.BuildTPouts.tstump=self.Legtt

        if PreBuildLvl==2 or PreBuildLvl==3:
            #set girders and diagonals equal to top braces
            self.BuildTPouts.Dgir=self.Xbrcinputs.Dbrc[-1]
            self.BuildTPouts.tgir=self.Xbrcinputs.tbrc[-1]
            self.BuildTPouts.Dbrc=self.Xbrcinputs.Dbrc[-1]
            self.BuildTPouts.tbrc=self.Xbrcinputs.tbrc[-1]
        if PreBuildLvl==4 or PreBuildLvl==5:
            self.BuildTPouts.Dbrc=self.TPinputs.Dgir
            self.BuildTPouts.tbrc=self.TPinputs.tgir

        #Leg bottom stump height set here if we have leg data and the user did not put a different value
        self.legbot_stmph=self.legbot_stmphin  #default value
        if (self.LegbD) and not(self.legbot_stmphin):   #User can never select =0, since it is forbidden
            self.legbot_stmph=1.5 * self.LegbD

        #Deck width set here if we have tower data and the user did not put a different fixed value for the deck, but a value function of Db
        if (self.TwrDb) and not(self.JcktGeoIn.dck_width) and self.JcktGeoIn.dck_widthfrac:
            self.dck_width = self.JcktGeoIn.dck_widthfrac * self.TwrDb
        else:
            self.dck_width = self.JcktGeoIn.dck_width

        #shorten names of class Jacket's attributes and calculates a few more parameters
                          #nbays
        nlegs=self.JcktGeoIn.nlegs                    #legs
        batter=self.JcktGeoIn.batter                  #2D batter (=tg(angle w.r.t. vertical))

        self.al_bat2D=  np.arctan(1./batter)               #angle between vertical and leg in projection (2D) [rad]
        self.al_bat3D=np.arctan(np.sqrt(2.)/batter)        #angle between vertical and leg  (3D) [rad]

        dck_botz=self.JcktGeoIn.dck_botz              #deck-bottom height
        weld2D=self.JcktGeoIn.weld2D                  #fraction of chord OD used in the weldment

        wdpth=self.WaterNew.wdepth    #CJBe                         #water depth
        #pileZtop=self.JcktGeoIn.pileZtop              #top of pile z coordinate

        legZbot=self.legZbot                          #bottom of leg z coordinate

        #Calculate some basic properties THAT CONTAIN SOME ASSUMPTIONS ON GEOMETRY
        self.innr_ang=pi/6. * (nlegs==3) + pi/4. *(nlegs==4) #Angle between radial direction and base side

        self.JcktH= dck_botz -self.BuildTPouts.hstump + wdpth - legZbot  #Jckt Height, from legtop to legbot
        Hbays= self.JcktH-self.legbot_stmph#-self.legtop_stmp  #Jacket Height available for bays,
            #removing 1D for stump at bottom and 0.5D at top and bottom for welds to brace

        #other params
        self.wbase=self.dck_width-2.*self.BuildTPouts.Dstump/2.*(1.+weld2D) +2.*(wdpth+dck_botz)/batter #virtual width at sea bed
        self.wbas0=self.wbase-2*np.tan(self.al_bat2D)*(legZbot+self.legbot_stmph) #width at 1st horiz. brace  joint

         #Calculate bay width, height, brace angle, and angle between X-braces (larger of the two)
        self.bay_bs,self.bay_hs,self.beta2D=FindBrcAng(Hbays,nbays,self.wbas0,self.al_bat2D)
        self.beta3D, al_Xbrc=FindBeta3D(self.beta2D, self.al_bat2D, self.al_bat3D, self.innr_ang, self.wbas0) #[rad],[rad]

#______________________________________________________________________________#
#JcktGeoOutputs is in VarTrees.py
##class JcktGeoOutputs(VariableTree):
##    """Node Coordinates and Member Connectivity"""
##    nNodesJckt = Int(  units=None,                 desc='Total Number of nodes in the Jacket, no Tower')
##    nodes =      Array(units='m',  dtype=np.float, desc='Node''s coordinates in the Substructure')
##    radii =      Array(units='m',  dtype=np.float, desc='Node''s Radii')
##    Reacts   =   Array(units=None, dtype=int,      desc='Node IDs for reactions + Fixity values (1/0)')  #Fixed for the time being with all fixity (6 dofs per node fixed)
##    nmems  =     Int(  units=None,                 desc='Total Number of Members in the Jacket, no Tower')
##    mems   =     Array(units=None, dtype=int,      desc='Member Connectivity Node i Node j for every member (i.e.,element)')
##    props  =     Array(units='m',  dtype=np.float, desc='Jacket Member''s xsectional and material properties')
##    XnsfM  =     Array(            dtype=np.float, desc='Jacket Member''s Global to Local DIrection Cosine Matrices [3,3,nelems]')
##    TubeObjs=    Instance(Klass=Tube,                  desc='Object of Class Tube for Jacket Elements. Tube Objects, one per element')
##    jnt_masses=  Array(            dtype=np.float, desc='Jacket Concentrated Masses at Joints: (joint no, mass, IXX, IYY, IZZ)')

#______________________________________________________________________________#
class BuildGeometry(Component):
    """This function assembles the full structure made of individual components to come up with a general xyz node distribution.
    It also calculates the properties of the various members and sets up variables to be passed to Frame3DD caller.
    """

    #inputs
    Frame3DDprms=VarTree(Frame3DDaux(),     iotype='in',desc='Frame3DD auxiliary parameters')
    JcktGeoIn=   VarTree(JcktGeoInputs(),   iotype='in',desc='Jacket Geometry Basic Inputs')
    Pileouts=    VarTree(PileGeoOutputs(),  iotype='in',desc='Pile Geometry Data')   #From TP
    Legouts=     VarTree(LegGeoOutputs(),   iotype='in',desc='Leg Geometry Data')   #From Legs
    TPouts=      VarTree(TPGeoOutputs(),    iotype='in',desc='TP Geometry Data')   #From TP
    Xbrcouts=    VarTree(XBrcGeoOutputs(),  iotype='in',desc='Xbrace Node Data')   #From Xbraces
    Mbrcouts=    VarTree(MudBrcGeoOutputs(),iotype='in',desc='MudBrace Node data')   #From Mudbraces
    Hbrcouts=    VarTree(HBrcGeoOutputs(),  iotype='in',desc='Top H-Brace Node data')   #From Hbraces
    Twrouts=     VarTree(TwrGeoOutputs(),   iotype='in',desc='Tower Node data')   #From Tower

    #outputs
    JcktGeoOut = VarTree(JcktGeoOutputs(),       iotype='out', desc='Geometry of the Jacket -Node Coordinates and Member Connectivity')
    pillegZs    =Array(dtype=np.float, units='m',iotype='out',desc='Z coordinates of nodes of pile and leg #1''s nodes.')
    pillegDs    =Array(dtype=np.float, units='m',iotype='out',desc='ODs of pile and leg #1''s nodes.')
    twrDs       =Array(dtype=np.float, units='m',iotype='out',desc='ODs of Tower''s nodes.')
    Pilemems=    Array(dtype=int,     iotype='out',desc='Connectivity Array for all members of the pile portion only')
    Legmems =    Array(dtype=int,     iotype='out',desc='Connectivity Array for all members in the leg portion only')
    Twrmems     =Array(dtype=int,     iotype='out',desc='Connectivity Array for all members of the tower portion only (including rigid member at the top if requested)')
    TPmems      =Array(dtype=int,     iotype='out', desc='Connectivity Array for all members of the TP portion only')
    XjntIDs     =Array(dtype=int,     iotype='out',desc='X-Joint IDs, as for Frame3DD, used to do checks later')
    KjntIDs     =Array(dtype=int,     iotype='out',desc='K-Joint IDs, as for Frame3DD, used to do checks later')

    SDPropSet=Array(iotype='out', desc='Properties to be used in SubDyn') #CJB+ test

    def execute(self):
        #simplify nomenclature
        nlegs=self.JcktGeoIn.nlegs
        nbays=self.JcktGeoIn.nbays
        VPFlag=self.JcktGeoIn.VPFlag                            #whether piles are vertical(true) or battered (false)

        nNodesLeg=self.Legouts.nNodesLeg  #nodes per leg
        nNodesXbrc=self.Xbrcouts.nNodesXbrc
        nNodesMbrc=self.Mbrcouts.nNodesMbrc #nodes per face (internal only, no joints with legs)
        nNodesHbrc=self.Hbrcouts.nNodesHbrc #nodes per face (internal only, no joints with legs)
        nNodesPile=self.Pileouts.nNodesPile #nodes per leg  (internal only, no joints with legs)

        nNodesBrc=self.TPouts.nNodesBrc     #nodes per leg  (internal only, no joints with legs or stem)
        nNodesStrt=self.TPouts.nNodesStrt   #nodes per leg  (internal only, no joints with legs or stem)
        nNodesGir=self.TPouts.nNodesGir     #nodes per leg  (internal only, no joints with legs)
        nNodesStmp=self.TPouts.nNodesStmp   #nodes per leg  (internal +top joint with braces but no joint with leg)
        nNodesStem=self.TPouts.nNodesStem   #nodes for the stem (all of them including joints with struts and braces)

        nNodesTwr=self.Twrouts.nNodes-1    #number of nodes subtracting joint at TP

    # NODES FIRST

    #Total number of nodes
        self.JcktGeoOut.nNodesJckt = nlegs* (nNodesPile+nNodesLeg+nNodesXbrc+nNodesMbrc+nNodesHbrc+nNodesBrc+nNodesStrt+nNodesGir+nNodesStmp)+nNodesStem +nNodesTwr
        out1=np.arange(0,self.JcktGeoOut.nNodesJckt,dtype=int)+1 #node IDs  (+1 since python starts at 0, but frame3DD needs 1)
        self.JcktGeoOut.radii=np.zeros(self.JcktGeoOut.nNodesJckt)      #radii for frame3DD, this will be fixed at one point

    #Combine piles and legs such that the nodes are sequential along individual leg
        pilenodes=self.Pileouts.nodes            #this puts the nodes in order: 1. bottom up of pile 1, then 2,3,4
        if self.Pileouts.nodes.size:
            pilenodes=pilenodes.squeeze(axis=1)            #this puts the nodes in order: 1. bottom up of pile 1, then 2,3,4
        legnodes=self.Legouts.xyz.transpose([1,2,0])
        pilegnodes=np.concatenate((pilenodes,legnodes),axis=1).reshape([-1,3]) #this puts the nodes in order: 1. bottom2up of leg 1, then 2,3,4

    # Now start combining Legs and Xbraces: we need to find a way to map nodes to list of IDs

    #Add x-brace nodes: first LLUR ones
        LLUR1=self.Xbrcouts.nodesLLUR1
        LLUR2=self.Xbrcouts.nodesLLUR2
        ULLR1=self.Xbrcouts.nodesULLR1
        ULLR2=self.Xbrcouts.nodesULLR2
        Xndiv=0  # Initialize. This is how many nodes per semi-brace LLUR or ULLR, not counting Xjoints or Legjnts

        if LLUR1.size: #This means ndvi>1 for Xbraces  -COULD NOT AVOID THIS IF.....
            LLUR1=self.Xbrcouts.nodesLLUR1.transpose([1,2,0])   #LLUR1 nodes' coordinates
            LLUR2=self.Xbrcouts.nodesLLUR2.transpose([1,2,0])   #LLUR2 nodes' coordinates
            ULLR1=self.Xbrcouts.nodesULLR1.transpose([1,2,0])   #LLUR1 nodes' coordinates
            ULLR2=self.Xbrcouts.nodesULLR2.transpose([1,2,0])   #LLUR2 nodes' coordinates
            Xndiv=self.Xbrcouts.nodesLLUR1.shape[2]/nbays  #Update is how many nodes per semi-brace LLUR or ULLR, not counting Xjoints or Legjnts

        Xjnts=self.Xbrcouts.joints.transpose([1,2,0]).reshape([-1,3])   #Xbrc joints coordinates

    #Add mud and top H-brace
        Mbrcnds=self.Mbrcouts.nodes  #Mud brace intermediate nodes' coordinates, CCW starting from face 1
        Hbrcnds=self.Hbrcouts.nodes

    #Add TP stumps, stem, girders, diagonal braces, struts  in this order
        stmpnds=self.TPouts.stmpnodes
        stemnds=self.TPouts.stemnodes
        girnds= self.TPouts.girnodes
        brcnds=self.TPouts.brcnodes
        strtnds=self.TPouts.strtnodes

    #Add Tower nodes
        twrnds=self.Twrouts.nodes[:,1:].T  #Exclude TP joint at top of stem

        self.JcktGeoOut.nodes=np.vstack((pilegnodes,np.concatenate((LLUR1,LLUR2),2).reshape([-1,3]),\
                                     np.concatenate((ULLR1,ULLR2),2).reshape([-1,3]),Xjnts,\
                                     Mbrcnds,Hbrcnds,\
                                     stmpnds,stemnds,girnds,brcnds,strtnds, twrnds))   #xyz of all nodes

        #Now I need to define the logic for the member definition

        #Pile members
        nNodesPileg=nNodesLeg+nNodesPile  #Nodes in one pile/leg combo
        nMpPile=nNodesPile     #members per pile
        idx=out1[0:nNodesPile+1].reshape([1,-1]).repeat(nlegs,0)
        idx+= +np.arange(0,nlegs*nNodesPileg,nNodesPileg)[:,np.newaxis]
        self.Pilemems=idx.repeat(2,1)[:,1:-1].reshape([-1,2])

        #leg members(elements)
        nMpLeg=nNodesLeg-1     #members per leg (no pile)
        legndiv= (nNodesLeg-1)/(nbays+1)  #number of elements per joint2joint member
        totnNodeleg=nlegs*nNodesPileg
        #first find indices to skip, which are the pile nodes, build a filter
        idx2skip=idx[:,:-1]-1   #Indices in out1 of pile nodes;
        junk=np.ones(totnNodeleg,dtype=bool)     #Set all to 1 on the final
        junk[idx2skip]=False
        idx=out1[0:totnNodeleg][junk] #just leg IDs (no piles)
        #Additionally I need to remove from out1 the bottom and top node of each leg before repeating the array to form member connectivity
        idxsize=idx.size*2
        junk=np.ones(idxsize,dtype=bool)     #Set all to 1 on the final
        idx2skip=np.arange(0,idxsize,2*nNodesLeg)
        idx2skip=np.hstack((idx2skip-1,idx2skip))
        junk[idx2skip]=False
        self.Legmems=idx.repeat(2,0)[junk].reshape([-1,2])

        #XBrace Members-note account for presence of pile nodes with the offset
        LLUR1jntIDs=np.arange(0,totnNodeleg).reshape(nlegs,-1)[:,nNodesPile+legndiv:-1:legndiv] +1 #leg internal nodes IDs organized per leg; +1 for python convention
        ULLR1jntIDs=np.arange(0,totnNodeleg).reshape(nlegs,-1)[:,nNodesPile+2*legndiv::legndiv]+1  #leg internal nodes IDs organized per leg;+1 for python convention
        ULLR2jntIDs=np.roll(LLUR1jntIDs,-1,0)
        LLUR2jntIDs=np.roll(ULLR1jntIDs,-1,0)

        idxXj=totnNodeleg+Xndiv*4*nbays*nlegs  #Starting index for Xjnts
        XjntIDs=out1[idxXj:idxXj+nbays*nlegs].reshape([nlegs,-1])   #Indices of joints, organized per face and nbays

        nMpXbrc=2*(Xndiv+1)  #number of members per brace (leg joint to leg joint)
        totnNodeLLUR= 2*Xndiv*nbays*nlegs  #nodes belonging to LLUR in out1  (does not include any joints)
        LLURmems=np.empty([nMpXbrc*nbays*nlegs,2],dtype=int)
        ULLRmems=LLURmems.copy()

        #MudBrace Members
        Mndiv=self.Mbrcouts.nodes.shape[0]/nlegs  #this is how many nodes per Mud-brace , not counting joints with Legs
        nMpMBrc=(nNodesMbrc+1)  #number of members per Mudbrace
        Mudmems=np.empty([nlegs*nMpMBrc,2],dtype=int)  #Mud members
        idxMb=idxXj+nbays*nlegs  #Starting indices for mud brace nodes in out1
        mudjntIDs=out1[nNodesPile+(1-VPFlag*(nNodesPile!=0))*legndiv:totnNodeleg:nNodesPileg]  #Leg base joint IDs  : THIS LOGIC SHOULD BE SET BY MudBrcs

        #Top H-Brace Members
        Hndiv=self.Hbrcouts.nodes.shape[0]/nlegs  #this is how many nodes per Mud-brace , not counting joints with Legs
        nMpHBrc=nNodesHbrc+1 -(Hndiv==0)  #number of members per Hbrace, accounting for possiblity of no Hbraces too.
        Hmems=np.empty([nlegs*nMpHBrc,2],dtype=int)  #H members
        idxHb=idxMb+Mndiv*nlegs  #Starting indices for H-brace nodes in out1
        legTopjntIDs=out1[nNodesPileg-1:totnNodeleg:nNodesPileg] #Leg Top joint IDs

        #TP stumps members
        brcjnts=legTopjntIDs #Initialize in case hstump is 0, therefore no upper stumps
        totnNodeStmp=nlegs*nNodesStmp   #total numder of stump nodes, except for leg joints
        stmpndiv=nNodesStmp#=self.TPouts.stmpnodes.shape[0]/nlegs  #this is how many nodes per stump , not counting joints with Legs
        nMpStmp=nNodesStmp  #number of members per stump, no +1 on this one
        Stmpmems=np.empty([nlegs*nMpStmp,2],dtype=int)  #Stump members
        idxSp=idxHb+Hndiv*nlegs  #Starting index for stump nodes in out1
        idxSp2=idxSp+stmpndiv*nlegs  #Ending index for stump nodes in out1
        if nMpStmp:
            brcjnts=out1[idxSp+nNodesStmp-1:idxSp2:nNodesStmp] # Update these also=Girder Joints indices in out1

        #TP stem members follow stumps
        #nNodesStem  =this is how many nodes per stem ; all included
        nMpStem=nNodesStem-1  #number of members per stem
        Stemmems=np.empty([nMpStmp,2],dtype=int)  #Stem members
        idxSm=idxSp2  #Starting index for stem nodes in out1, bottom up
        idxSm2=idxSm+nNodesStem-1     #End index for stem nodes in out1, bottom up (-1 since both start and end joint included here)
        Stemmems=out1[idxSm:idxSm2+1].repeat(2)[1:-1].reshape(-1,2)

        #TP Girder members
        girdndiv=self.TPouts.girnodes.shape[0]/nlegs  #this is how many nodes per girder , not counting joints with stumps
        nMpGir=(nNodesGir+1)  #number of members per girder
        Girdmems=np.empty([nlegs*nMpGir,2],dtype=int)  #Girder members
        idxGr=idxSm2+1 #Starting index for girder nodes in out1
        #idxGr2=idxGr+nlegs*nNodesGir #End index for girder nodes in out1

        #TP Diagonal Brace members
        brcndiv=self.TPouts.brcnodes.shape[0]/nlegs  #this is how many nodes per brace , not counting joints with stumps and stem
        nMpBrc=(nNodesBrc+1)  #number of members per brace
        Brcmems=np.empty([nlegs*nMpBrc,2],dtype=int)  #Brace members
        idxBr=idxGr+girdndiv*nlegs    #Starting indices for brace nodes in out1

        #TP Strut members
        strutndiv=self.TPouts.strtnodes.shape[0]/nlegs  #this is how many nodes per struts , not counting joints with stumps and stem
        nMpStrt=(nNodesStrt+1)  #number of members per struts
        Strtmems=np.empty([nlegs*nMpStrt,2],dtype=int)  #Struts members
        idxSt=idxBr+brcndiv*nlegs    #Starting indices for strut nodes in out1

        for ii in range(0,nlegs):    #I should be able to vectorize this one, for now leave as is
            for jj in range(0,nbays):

                #one brace at a time:member index fo
                idx00=nMpXbrc*(jj+ii*nbays)
                idx11=idx00+nMpXbrc

                #index in out1 for LLUR
                idx0=totnNodeleg+Xndiv*2*(jj+ii*nbays)
                idx1=idx0+Xndiv
                idx2=idx1+Xndiv
                junk=np.hstack((LLUR1jntIDs[ii,jj],out1[idx0:idx1],  XjntIDs[ii,jj], out1[idx1:idx2], LLUR2jntIDs[ii,jj]))  #All nodes IDs for current brace
                LLURmems[idx00:idx11,:]=junk.repeat(2)[1:-1].reshape(-1,2)
                #index in out1 for ULLR
                idx0=idx0+totnNodeLLUR
                idx1=idx0+Xndiv
                idx2=idx1+Xndiv
                junk=np.hstack((ULLR1jntIDs[ii,jj],out1[idx0:idx1],  XjntIDs[ii,jj], out1[idx1:idx2], ULLR2jntIDs[ii,jj]))  #All nodes IDs for current brace
                ULLRmems[idx00:idx11,:]=junk.repeat(2)[1:-1].reshape(-1,2)
                #LLURnodes[ii,idx00:idx11,:]=np.vstack((LLUR1jnts[ii,jj,:],LLUR1[ii,idx0:idx1,:],Xjnts[ii,jj,:],LLUR2[ii,idx0:idx1,:],LLUR2jnts[ii,jj,:]))
                #ULLRnodes[ii,idx00:idx11,:]=np.vstack((ULLR1jnts[ii,jj,:],ULLR1[ii,idx0:idx1,:],Xjnts[ii,jj,:],ULLR2[ii,idx0:idx1,:],ULLR2jnts[ii,jj,:]))

            #use same ii loop for Mudbraces and H braces
            #member index
            idx00=nMpMBrc*ii
            idx11=idx00+nMpMBrc
            #indices into out1
            idx0=idxMb+ii*nNodesMbrc
            idx1=idx0+nNodesMbrc
            junk=np.hstack((mudjntIDs[ii],out1[idx0:idx1],mudjntIDs[np.mod(ii+1,nlegs)]))  #All nodes IDs for current brace
            Mudmems[idx00:idx11,:]=junk.repeat(2)[1:-1].reshape(-1,2)
            #member index
            idx00=nMpHBrc*ii
            idx11=idx00+nMpHBrc
            #indices into out1
            idx0=idxHb+ii*nNodesHbrc
            idx1=idx0+nNodesHbrc
            junk=np.hstack((legTopjntIDs[ii],out1[idx0:idx1],legTopjntIDs[np.mod(ii+1,nlegs)]))  #All nodes IDs for current brace
            Hmems[idx00:idx11,:]=junk.repeat(2)[1:-1].reshape(-1,2)

            #use same ii loop for TP stumps,girders,braces, and struts
            #Stump member index
            idx00=nMpStmp*ii
            idx11=idx00+nMpStmp
            #indices into out1
            idx0=idxSp+ii*nNodesStmp
            idx1=idx0+nNodesStmp
            junk=np.hstack((legTopjntIDs[ii],out1[idx0:idx1]))  #All nodes IDs for current stump
            Stmpmems[idx00:idx11,:]=junk.repeat(2)[1:-1].reshape(-1,2)

            #Girder member index
            idx00=nMpGir*ii
            idx11=idx00+nMpGir
            #indices into out1
            idx0=idxGr+ii*nNodesGir
            idx1=idx0+nNodesGir
            junk=np.hstack((brcjnts[ii],out1[idx0:idx1],brcjnts[np.mod(ii+1,nlegs)]))  #All nodes IDs for current brace
            Girdmems[idx00:idx11,:]=junk.repeat(2)[1:-1].reshape(-1,2)

            #Diagonal Brace member index
            idx00=nMpBrc*ii
            idx11=idx00+nMpBrc
            #indices into out1
            idx0=idxBr+ii*nNodesBrc
            idx1=idx0+nNodesBrc
            junk=np.hstack((brcjnts[ii],out1[idx0:idx1],out1[idxSm]))  #All nodes IDs for current brace
            Brcmems[idx00:idx11,:]=junk.repeat(2)[1:-1].reshape(-1,2)

            #Strut member index
            idx00=nMpStrt*ii
            idx11=idx00+nMpStrt
            #indices into out1
            idx0=idxSt+ii*nNodesStrt
            idx1=idx0+nNodesStrt
            junk=np.hstack((brcjnts[ii],out1[idx0:idx1],out1[idxSm2]))  #All nodes IDs for current strut
            Strtmems[idx00:idx11,:]=junk.repeat(2)[1:-1].reshape(-1,2)

        #Now Tower Members
        idxtwr=idxSt+nlegs*strutndiv   #starting index for tower nodes
        self.Twrmems=np.hstack((out1[idxSm2],out1[idxtwr:])).repeat(2,0)[1:-1].reshape(-1,2)
        nMpTwr=self.Twrmems.shape[0]  #Number of tower members, including rigid Twr2RNA


        #Store TP mems for later
        self.TPmems=np.vstack((Stmpmems,Stemmems,Girdmems,Brcmems,Strtmems))

        #Combine all members (elements)
        self.JcktGeoOut.mems=np.vstack((self.Pilemems,self.Legmems,LLURmems,ULLRmems,Mudmems,Hmems,self.TPmems,self.Twrmems))
        self.JcktGeoOut.nmems=self.JcktGeoOut.mems.shape[0]

        #_________________Now get the properties for all members__________________#

        #piles
        pileprops=np.empty([nMpPile*nlegs,10])
        junk=np.array([self.Pileouts.PileObj.Area,self.Pileouts.PileObj.Asx,self.Pileouts.PileObj.Asy,\
                self.Pileouts.PileObj.J0,self.Pileouts.PileObj.Jxx,self.Pileouts.PileObj.Jyy]).transpose() # [6]
        pileprops[:,0:6]=junk.reshape([-1+(nMpPile==0),junk.shape[0],6]).repeat(nlegs,axis=0).reshape([nlegs*nMpPile,6]) #repeated for all nlegs piles [nlegs,6]; first condition accounts for 0 pile requested
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.Pileouts.PileObj.mat])
        Gs=np.array([mat.G for mat in self.Pileouts.PileObj.mat])
        rhos=np.array([mat.rho for mat in self.Pileouts.PileObj.mat])
        rolls=np.zeros([nMpPile,1])
        pileprops[:,6:10]=np.hstack((Es.reshape([-1,1]).repeat(nlegs,axis=0),Gs.reshape([-1,1]).repeat(nlegs,axis=0),\
                                    rolls.repeat(nlegs,axis=0),rhos.reshape([-1,1]).repeat(nlegs,axis=0)))

        #legs
        legprops=np.empty([nMpLeg*nlegs,10])
        junk=np.array([self.Legouts.LegObj.Area,self.Legouts.LegObj.Asx,self.Legouts.LegObj.Asy,\
                self.Legouts.LegObj.J0,self.Legouts.LegObj.Jxx,self.Legouts.LegObj.Jyy]).transpose() # [6]
        legprops[:,0:6]=junk.reshape([-1,junk.shape[0],6]).repeat(nlegs,axis=0).reshape([nlegs*nMpLeg,6])  #repeated for all nlegs  [nlegs,6]
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.Legouts.LegObj.mat])
        Gs=np.array([mat.G for mat in self.Legouts.LegObj.mat])
        rhos=np.array([mat.rho for mat in self.Legouts.LegObj.mat])
        rolls=np.zeros([nMpLeg,1])
        legprops[:,6:10]=np.hstack((Es.reshape([-1,1]).repeat(nlegs,axis=0),Gs.reshape([-1,1]).repeat(nlegs,axis=0),\
                                   rolls.repeat(nlegs,axis=0),rhos.reshape([-1,1]).repeat(nlegs,axis=0)))

        #Xbraces
        #xbrc
        LLURprops=np.empty([nMpXbrc*nbays*nlegs,10])
        junk=np.array([self.Xbrcouts.LLURObj.Area,self.Xbrcouts.LLURObj.Asy,self.Xbrcouts.LLURObj.Asy,\
                self.Xbrcouts.LLURObj.J0,self.Xbrcouts.LLURObj.Jxx,self.Xbrcouts.LLURObj.Jyy]).transpose() # [nbays,6] LL-UR
        LLURprops[:,0:6]=junk.reshape([-1,junk.shape[0],6]).repeat(nlegs,axis=0).reshape([nlegs*nbays*nMpXbrc,6])  #repeated for all nlegs  [nlegs,6]

        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.Xbrcouts.LLURObj.mat])
        Gs=np.array([mat.G for mat in self.Xbrcouts.LLURObj.mat])
        rhos=np.array([mat.rho for mat in self.Xbrcouts.LLURObj.mat])
        rolls=np.zeros([nMpXbrc*nbays,1])
        LLURprops[:,6:10]=np.hstack((Es.reshape([-1,1]).repeat(nlegs,axis=0),Gs.reshape([-1,1]).repeat(nlegs,axis=0),rolls.repeat(nlegs,axis=0),rhos.reshape([-1,1]).repeat(nlegs,axis=0)))
        #For the ULLR braces, the only thing that needs to be adjusted is how the lengths are handled, but Xsectional properties are the same as LLUR
        ULLRprops=LLURprops.copy()

        #Mudbraces
        Mbrcprops=np.empty([nMpMBrc*nlegs,10])
        junk=np.array([self.Mbrcouts.brcObj.Area,self.Mbrcouts.brcObj.Asx,self.Mbrcouts.brcObj.Asy,\
                self.Mbrcouts.brcObj.J0,self.Mbrcouts.brcObj.Jxx,self.Mbrcouts.brcObj.Jyy]).transpose() # [6]
        Mbrcprops[:,0:6]=junk.reshape([-1,junk.shape[0],6]).repeat(nlegs,axis=0).reshape([nlegs*nMpMBrc,6])  #repeated for all nlegs  [nlegs,6]
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.Mbrcouts.brcObj.mat])
        Gs=np.array([mat.G for mat in self.Mbrcouts.brcObj.mat])
        rhos=np.array([mat.rho for mat in self.Mbrcouts.brcObj.mat])
        rolls=np.zeros([nMpMBrc,1])
        Mbrcprops[:,6:10]=np.hstack((Es.reshape([-1,1]).repeat(nlegs,axis=0),Gs.reshape([-1,1]).repeat(nlegs,axis=0),rolls.repeat(nlegs,axis=0),rhos.reshape([-1,1]).repeat(nlegs,axis=0)))

        #Top HBraces
        Hbrcprops=np.empty([nMpHBrc*nlegs,10])
        junk=np.array([self.Hbrcouts.brcObj.Area,self.Hbrcouts.brcObj.Asx,self.Hbrcouts.brcObj.Asy,\
                self.Hbrcouts.brcObj.J0,self.Hbrcouts.brcObj.Jxx,self.Hbrcouts.brcObj.Jyy]).transpose() # [6]
        Hbrcprops[:,0:6]=junk.reshape([-1+(nMpHBrc==0),junk.shape[0],6]).repeat(nlegs,axis=0).reshape([nlegs*nMpHBrc,6])  #repeated for all nlegs  [nlegs,6] the condition at the start accounts for cases with Hndiv=0, i.e. no hbraces requested
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.Hbrcouts.brcObj.mat])
        Gs=np.array([mat.G for mat in self.Hbrcouts.brcObj.mat])
        rhos=np.array([mat.rho for mat in self.Hbrcouts.brcObj.mat])
        rolls=np.zeros([nMpHBrc,1])
        Hbrcprops[:,6:10]=np.hstack((Es.reshape([-1,1]).repeat(nlegs,axis=0),Gs.reshape([-1,1]).repeat(nlegs,axis=0),rolls.repeat(nlegs,axis=0),rhos.reshape([-1,1]).repeat(nlegs,axis=0)))

        #TPstumps
        Stmpprops=np.empty([nMpStmp*nlegs,10])
        junk=np.array([self.TPouts.stmpObj.Area,self.TPouts.stmpObj.Asx,self.TPouts.stmpObj.Asy,\
                self.TPouts.stmpObj.J0,self.TPouts.stmpObj.Jxx,self.TPouts.stmpObj.Jyy]).transpose() # [6]
        Stmpprops[:,0:6]=junk.reshape([-1+(nMpStmp==0),junk.shape[0],6]).repeat(nlegs,axis=0).reshape([nlegs*nMpStmp,6])  #repeated for all nlegs  [nlegs,6]; condition at the start accounts for cases with no stumps stumpndiv=0=hstump, i.e. no stumps requested
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.TPouts.stmpObj.mat])
        Gs=np.array([mat.G for mat in self.TPouts.stmpObj.mat])
        rhos=np.array([mat.rho for mat in self.TPouts.stmpObj.mat])
        rolls=np.zeros([nMpStmp,1])
        Stmpprops[:,6:10]=np.hstack((Es.reshape([-1,1]).repeat(nlegs,axis=0),Gs.reshape([-1,1]).repeat(nlegs,axis=0),rolls.repeat(nlegs,axis=0),rhos.reshape([-1,1]).repeat(nlegs,axis=0)))

        #TPGirders
        Girprops=np.empty([nMpGir*nlegs,10])
        junk=np.array([self.TPouts.girObj.Area,self.TPouts.girObj.Asx,self.TPouts.girObj.Asy,\
                self.TPouts.girObj.J0,self.TPouts.girObj.Jxx,self.TPouts.girObj.Jyy]).transpose() # [6]
        Girprops[:,0:6]=junk.reshape([-1,junk.shape[0],6]).repeat(nlegs,axis=0).reshape([nlegs*nMpGir,6])  #repeated for all nlegs  [nlegs,6]
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.TPouts.girObj.mat])
        Gs=np.array([mat.G for mat in self.TPouts.girObj.mat])
        rhos=np.array([mat.rho for mat in self.TPouts.girObj.mat])
        rolls=np.zeros([nMpGir,1])
        Girprops[:,6:10]=np.hstack((Es.reshape([-1,1]).repeat(nlegs,axis=0),Gs.reshape([-1,1]).repeat(nlegs,axis=0),rolls.repeat(nlegs,axis=0),rhos.reshape([-1,1]).repeat(nlegs,axis=0)))

        #TPBraces
        Brcprops=np.empty([nMpBrc*nlegs,10])
        junk=np.array([self.TPouts.brcObj.Area,self.TPouts.brcObj.Asx,self.TPouts.brcObj.Asy,\
                self.TPouts.brcObj.J0,self.TPouts.brcObj.Jxx,self.TPouts.brcObj.Jyy]).transpose() # [6]
        Brcprops[:,0:6]=junk.reshape([-1,junk.shape[0],6]).repeat(nlegs,axis=0).reshape([nlegs*nMpBrc,6])  #repeated for all nlegs  [nlegs,6]
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.TPouts.brcObj.mat])
        Gs=np.array([mat.G for mat in self.TPouts.brcObj.mat])
        rhos=np.array([mat.rho for mat in self.TPouts.brcObj.mat])
        rolls=np.zeros([nMpBrc,1])
        Brcprops[:,6:10]=np.hstack((Es.reshape([-1,1]).repeat(nlegs,axis=0),Gs.reshape([-1,1]).repeat(nlegs,axis=0),rolls.repeat(nlegs,axis=0),rhos.reshape([-1,1]).repeat(nlegs,axis=0)))

        #TPstruts
        Strtprops=np.empty([nMpStrt*nlegs,10])
        junk=np.array([self.TPouts.strtObj.Area,self.TPouts.strtObj.Asx,self.TPouts.strtObj.Asy,\
                self.TPouts.strtObj.J0,self.TPouts.strtObj.Jxx,self.TPouts.strtObj.Jyy]).transpose() # [6]
        Strtprops[:,0:6]=junk.reshape([-1,junk.shape[0],6]).repeat(nlegs,axis=0).reshape([nlegs*nMpStrt,6])  #repeated for all nlegs  [nlegs,6]
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.TPouts.strtObj.mat])
        Gs=np.array([mat.G for mat in self.TPouts.strtObj.mat])
        rhos=np.array([mat.rho for mat in self.TPouts.strtObj.mat])
        rolls=np.zeros([nMpStrt,1])
        Strtprops[:,6:10]=np.hstack((Es.reshape([-1,1]).repeat(nlegs,axis=0),Gs.reshape([-1,1]).repeat(nlegs,axis=0),rolls.repeat(nlegs,axis=0),rhos.reshape([-1,1]).repeat(nlegs,axis=0)))

        #TPstem
        Stemprops=np.empty([nMpStem,10])
        junk=np.array([self.TPouts.stemObj.Area,self.TPouts.stemObj.Asx,self.TPouts.stemObj.Asy,\
                self.TPouts.stemObj.J0,self.TPouts.stemObj.Jxx,self.TPouts.stemObj.Jyy]).transpose() # [6]
        Stemprops[:,0:6]=junk
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.TPouts.stemObj.mat])
        Gs=np.array([mat.G for mat in self.TPouts.stemObj.mat])
        rhos=np.array([mat.rho for mat in self.TPouts.stemObj.mat])
        rolls=np.zeros([nMpStem,1])
        Stemprops[:,6:10]=np.hstack((Es.reshape([-1,1]),Gs.reshape([-1,1]),rolls,rhos.reshape([-1,1])))

        #Tower including rigid twr2Rna component at the top in case
        Twrprops=np.empty([nMpTwr,10])
        junk=np.array([self.Twrouts.TwrObj.Area,self.Twrouts.TwrObj.Asx,self.Twrouts.TwrObj.Asy,\
                self.Twrouts.TwrObj.J0,self.Twrouts.TwrObj.Jxx,self.Twrouts.TwrObj.Jyy]).transpose() # [6]
        #add E,G, rho : this should really be improved this is pretty bad
        Es=np.array([mat.E for mat in self.Twrouts.TwrObj.mat])
        Gs=np.array([mat.G for mat in self.Twrouts.TwrObj.mat])
        rhos=np.array([mat.rho for mat in self.Twrouts.TwrObj.mat])
        rolls=np.zeros([Es.size,1])
        junk1=np.hstack((Es.reshape([-1,1]),Gs.reshape([-1,1]),rolls,rhos.reshape([-1,1])))

        Twrprops=np.hstack([junk,junk1])#Flexible tower properties

        TwrRgdProps=np.empty([])  #Initialize

        if self.Twrouts.Twr2RNAObj.Area: #not empty, that means that the tower top is made of a rigid member
            junk=np.array([self.Twrouts.Twr2RNAObj.Area,self.Twrouts.Twr2RNAObj.Asx,self.Twrouts.Twr2RNAObj.Asy,\
                            self.Twrouts.Twr2RNAObj.J0,self.Twrouts.Twr2RNAObj.Jxx,self.Twrouts.Twr2RNAObj.Jyy]) # [6]
            junk1=np.hstack((self.Twrouts.Twr2RNAObj.mat.E,self.Twrouts.Twr2RNAObj.mat.G,0.,self.Twrouts.Twr2RNAObj.mat.rho))
            TwrRgdProps=np.hstack([junk,junk1])#Rigid tower properties
            Twrprops=np.vstack((Twrprops,TwrRgdProps)) #Augmented TwrProps

        #Assemble all properties
        self.JcktGeoOut.props=np.vstack((pileprops,legprops,LLURprops,ULLRprops,Mbrcprops,Hbrcprops,Stmpprops,Stemprops,Girprops,Brcprops,Strtprops,Twrprops))

        #Also store the Tube objects into one big Tube Object, which should have 1 per element                         (note: ULLRObj=LLURObj)
        #The following is horrible, but I could not come up with anything to create the massive tube object
        #ODs
        Ds=np.hstack(( np.tile(np.hstack((self.Pileouts.AFObj.D,self.Pileouts.PileObj.D)),nlegs), \
                       np.tile(self.Legouts.LegObj.D,nlegs), np.tile(self.Xbrcouts.LLURObj.D,nlegs),np.tile(self.Xbrcouts.LLURObj.D,nlegs),\
                       np.tile(self.Mbrcouts.brcObj.D,nlegs),np.tile(self.Hbrcouts.brcObj.D,nlegs), np.tile(self.TPouts.stmpObj.D,nlegs),\
                       self.TPouts.stemObj.D, np.tile(self.TPouts.girObj.D,nlegs),   np.tile(self.TPouts.brcObj.D,nlegs), \
                       np.tile(self.TPouts.strtObj.D,nlegs),\
                       self.Twrouts.TwrObj.D,self.Twrouts.Twr2RNAObj.D))

        ts=np.hstack(( np.tile(np.hstack((self.Pileouts.AFObj.t,self.Pileouts.PileObj.t)),nlegs), \
                       np.tile(self.Legouts.LegObj.t,nlegs), np.tile(self.Xbrcouts.LLURObj.t,nlegs),np.tile(self.Xbrcouts.LLURObj.t,nlegs),\
                       np.tile(self.Mbrcouts.brcObj.t,nlegs),np.tile(self.Hbrcouts.brcObj.t,nlegs), np.tile(self.TPouts.stmpObj.t,nlegs),\
                       self.TPouts.stemObj.t, np.tile(self.TPouts.girObj.t,nlegs),   np.tile(self.TPouts.brcObj.t,nlegs), \
                       np.tile(self.TPouts.strtObj.t,nlegs),\
                       self.Twrouts.TwrObj.t,self.Twrouts.Twr2RNAObj.t))

        Lgths=np.hstack(( np.tile(np.hstack((self.Pileouts.AFObj.L,self.Pileouts.PileObj.L)),nlegs), \
                       np.tile(self.Legouts.LegObj.L,nlegs), np.tile(self.Xbrcouts.LLURObj.L,nlegs),np.tile(self.Xbrcouts.LLURObj.L,nlegs),\
                       np.tile(self.Mbrcouts.brcObj.L,nlegs),np.tile(self.Hbrcouts.brcObj.L,nlegs), np.tile(self.TPouts.stmpObj.L,nlegs),\
                       self.TPouts.stemObj.L, np.tile(self.TPouts.girObj.L,nlegs),   np.tile(self.TPouts.brcObj.L,nlegs), \
                       np.tile(self.TPouts.strtObj.L,nlegs),\
                       self.Twrouts.TwrObj.L,self.Twrouts.Twr2RNAObj.L))

        Kbucks=np.hstack(( np.tile(np.hstack((self.Pileouts.AFObj.Kbuck,self.Pileouts.PileObj.Kbuck)),nlegs), \
                       np.tile(self.Legouts.LegObj.Kbuck,nlegs), np.tile(self.Xbrcouts.LLURObj.Kbuck,nlegs),np.tile(self.Xbrcouts.LLURObj.Kbuck,nlegs),\
                       np.tile(self.Mbrcouts.brcObj.Kbuck,nlegs),np.tile(self.Hbrcouts.brcObj.Kbuck,nlegs), np.tile(self.TPouts.stmpObj.Kbuck,nlegs),\
                       self.TPouts.stemObj.Kbuck, np.tile(self.TPouts.girObj.Kbuck,nlegs),   np.tile(self.TPouts.brcObj.Kbuck,nlegs), \
                       np.tile(self.TPouts.strtObj.Kbuck,nlegs),\
                       self.Twrouts.TwrObj.Kbuck,self.Twrouts.Twr2RNAObj.Kbuck))

        mats=np.hstack(( np.tile(np.hstack((self.Pileouts.AFObj.mat,self.Pileouts.PileObj.mat)),nlegs), \
                       np.tile(self.Legouts.LegObj.mat,nlegs), np.tile(self.Xbrcouts.LLURObj.mat,nlegs),np.tile(self.Xbrcouts.LLURObj.mat,nlegs),\
                       np.tile(self.Mbrcouts.brcObj.mat,nlegs),np.tile(self.Hbrcouts.brcObj.mat,nlegs), np.tile(self.TPouts.stmpObj.mat,nlegs),\
                       self.TPouts.stemObj.mat, np.tile(self.TPouts.girObj.mat,nlegs),   np.tile(self.TPouts.brcObj.mat,nlegs), \
                       np.tile(self.TPouts.strtObj.mat,nlegs),\
                       self.Twrouts.TwrObj.mat,self.Twrouts.Twr2RNAObj.mat))

        #CJB TODO: This technique creates redundant property sets. A better method would be to search the PropSetSub array and only make a new property set for unique values.
        i=0 #CJB+
        EArray=[] #CJB+
        GArray=[] #CJB+
        rhoArray=[] #CJB+
        PropSetSub=[] #CJB+
        SDPropSet=[1,1,1,1,1,1] #CJB+
        for i in range(len(mats)): #CJB+
            EArray=np.append(EArray,mats[i].E) #CJB+
            GArray=np.append(GArray,mats[i].G) #CJB+
            rhoArray=np.append(rhoArray,mats[i].rho) #CJB+
            PropSetSub=np.hstack(((i+1),EArray[i],GArray[i],rhoArray[i],Ds[i],ts[i])) #CJB+
            SDPropSet=np.vstack((SDPropSet,PropSetSub)) #CJB+

        self.SDPropSet=np.delete(SDPropSet, (0), axis=0) #CJB+

        self.JcktGeoOut.TubeObjs=Tube(Ds,ts,Lgths,Kbucks,mats)

        #store matrix of [nx2] with the nodes for each member and direction cosine matrix for each member
        X1=self.JcktGeoOut.nodes[self.JcktGeoOut.mems[:,0]-1,:].T  #-1 for python
        X2=self.JcktGeoOut.nodes[self.JcktGeoOut.mems[:,1]-1,:].T
        #X1=nds[props[:,1].astype(int)-1,1:4].transpose()  #coordinates of node i's of the members #note '-1' due to frame3DD starting at 1 python at 0
        #X2=nds[props[:,2].astype(int)-1,1:4].transpose()  #coordinates of node j's of the members (Note transpose to make it faster index last)
        #Ls=CalcDist(X1,X2)  #member lengths
        #roll=nds[props[:,1].astype(int)-1,4] #roll angle of element
        rolls=self.JcktGeoOut.props[:,-2]
        self.JcktGeoOut.XnsfM=Melm(X1,X2,rolls)  #[3,3,nmems]

        #_____Now take care of Extra Joint Inertias____#
        # TP extra mass data
        TPjnt_inertia=self.TPouts.TPlumpedMass
        # RNA extra mass data
        RNAjnt_inertia=self.Twrouts.TopMass
        #Tower concentrated masses
        Twrjnt_inertia=self.Twrouts.TwrlumpedMass  #(n,7)

        #joint IDs for lumped masses
        TP_jnt_in_no  =out1[idxSm]*(TPjnt_inertia[0] !=0)
        RNA_jnt_in_no =out1[-1]*(RNAjnt_inertia[0] !=0)

        #For the tower I need to find the nearest nodes to the assigned heights;
        Twr_jnt_in_no =0 #initialize
        nTwrmas=np.nonzero(Twrjnt_inertia[:,1])[-1]
        if nTwrmas.size:
            nTwrmas=nTwrmas[-1]+1 #+1 for python
              #making sure the array is not all zeros
            ztwr=self.JcktGeoOut.nodes[out1[idxtwr:]-1,2]  #-1 for python
            Twr_jnt_in_no=np.zeros(nTwrmas)
            for ii in range(0,nTwrmas):
               Twr_jnt_in_no[ii]=out1[idxtwr+np.argmin(np.abs(ztwr-Twrjnt_inertia[ii,0]))]

            #I now need to make sure I am not trumping the RNA lumped mass in case there is a yaw collar or concentrated mass up there
            #This works as coded only if the last inertia
            idx=np.in1d(Twr_jnt_in_no,RNA_jnt_in_no).nonzero()[0]
            if idx.size:
                #sum to RNA mass
                RNAjnt_inertia += Twrjnt_inertia[idx[0],1:]
                #Remove from twr mass list
                nTwrmas -= 1
                Twrjnt_inertia=np.delete(Twrjnt_inertia,idx[0],0)
                Twr_jnt_in_no=np.delete(Twr_jnt_in_no,idx[0],0)
        else:
            nTwrmas=0


        nj_int=int(TP_jnt_in_no !=0) + int(RNA_jnt_in_no !=0) + nTwrmas

        junk=np.hstack(([TP_jnt_in_no,RNA_jnt_in_no],Twr_jnt_in_no))  #need to do it with hstack as Twr_jnt_in_no may be an array
        jnt_masses_no=junk[junk.nonzero()].reshape(-1,1)
        inertias=np.vstack((TPjnt_inertia,RNAjnt_inertia,Twrjnt_inertia[:,1:]))[junk.nonzero()[0],:]

        self.JcktGeoOut.jnt_masses=np.hstack((jnt_masses_no,inertias)).reshape(-1,11)
        #Account for potential yaw
        RNAjnt_inertia2=self.Twrouts.TopMass_yaw
        inertias2=np.vstack((TPjnt_inertia,RNAjnt_inertia2,Twrjnt_inertia[:,1:]))[junk.nonzero()[0],:]
        self.JcktGeoOut.jnt_masses_yaw=np.hstack((jnt_masses_no,inertias2)).reshape(-1,11)


        #Also make these joints heavy for self-weight along z
        #junk1=np.array([TPjnt_inertia[0],RNAjnt_inertia[0]])
        #junk1=junk1[junk.flatten().nonzero()] #in case some concentrated masses are not used
        #junk=np.dot(junk1.reshape(nj_int,1), np.array([aux['g_x'],aux['g_y'],aux['g_z']]).reshape(1,3))
        #Now add zeros at the end since we are not adding Mxx,Myy,Mzz concentrated moments
        #junk1=np.zeros([nj_int,6])
        #junk1[:,:-3]=junk
        #jnt_frc=np.hstack((jnt_in_no,junk1)).reshape(-1,7)

        #SET REACTION NODE IDS
        ReactNIDs=out1[0:totnNodeleg:nNodesPileg].reshape([-1,1])  #IDs of reaction nodes
        self.JcktGeoOut.Reacts=np.hstack((ReactNIDs,np.ones([nlegs,1]).repeat(6,axis=1)))

        #Store pile and leg nodes (1 leg only) z-coordinate for use later in loading
        self.pillegZs=pilegnodes[0:nNodesPileg,2]  #
        self.pillegDs=np.concatenate((self.Pileouts.PileObj.D,self.Legouts.LegObj.D))
        self.twrDs=self.Twrouts.TwrObj.D

        #Store Xjoint IDs, Xbraces + TP Xbrace at the top, but that only if we want to include the TP in the verification
        #self.XjntIDs=np.hstack((XjntIDs.flatten(), out1[idxSm]))
        self.XjntIDs=XjntIDs.flatten()
        self.KjntIDs=np.unique(np.hstack((LLUR1jntIDs.flatten(),LLUR2jntIDs.flatten(),mudjntIDs)))  #Top hbrace joints included regardless of whether or not H brace is present;mudjoints are added in case the mudbrace is offset from bottom xbrace

#______________________________________________________________________________#

class PileEmbdL(Component):
    #This component calculates pile embedment length to sustain axial load, no lateral capacity yet

    #inputs
    SoilObj    =Instance(Klass=SoilC,     iotype='in', desc='Object of Class Soil used for Embedment Length')
    Lp         =Float(        units='m',  iotype='in', desc='Input Embedment Length')  #User defined or optimizer assigned
    PileObjout = Instance(Klass=Tube,     iotype='in', desc='Object of Class Tube for Non-AF portion of Pile, valid also when ndiv=0')

    MbrFrcs   =   Array(dtype=np.float,   iotype='in', desc='Forces and Moments at Element Nodes from Frame3DD')
    JcktGeoOut = VarTree(JcktGeoOutputs(),iotype='in', desc='Geometry of the Jacket -Node Coordinates and Member Connectivity') #From Build, used only for Reacts, should modify this

    gravity    =Float(units='m/s**2',iotype='in',  desc='Acceleration Gravity.')
        #outputs
    Lp0rat=Float(        units='m',  iotype='out', desc='(Lp-Lp0)/Lp0 Normalized difference between current and needed embedment length')  #necessary to sustain load
    Mpiles  =   Float(units='kg',    iotype='out', desc='Embedded Pile Mass (all piles) with the given Lp, not necessarily = to needed embedment length')

    def execute(self):
        #Reactions are in GLOBAL coordinate system, they should be rotated in the local to account for axial force properly
        #I could make sure of the member forces which are in local coordinate system (Fx for frame3dd is along the member)
        #So find the members that contain the reaction points
        ReactNds=self.JcktGeoOut.Reacts[:,0]  #reaction nodes
        npiles=ReactNds.size #(There should be as many reactions as the legs, i do so I do not have to pass nlegs)
        MbrFrcs=self.MbrFrcs.squeeze().T
        Reacts=MbrFrcs[np.in1d(MbrFrcs[:,1],ReactNds),2]  #Fx [N]
        Nhead=np.max(np.abs(Reacts))
        #pile properties pre-calculated elsewhere, unless no pile whatsoever is used in that case use leg properties

        Dp=self.PileObjout.D
        tp=self.PileObjout.t
        rho_p=self.PileObjout.mat[0].rho  # density

        Lp0=-EmbedLength(Dp,tp,rho_p,Nhead,self.SoilObj,gravity=abs(self.gravity))#*self.SoilObj.SoilSF passed to friction calculation instead

        print('>>>>>>>>>>  needed embedment Lp0=',Lp0) #s

        self.Mpiles=npiles*rho_p*self.PileObjout.Area*self.Lp  #Mass of all piles with input Lp

        self.Lp0rat=(self.Lp-Lp0)/Lp0  #This is for the optimizer must be >0 constraint

#_____________________________________________________#
class SPIstiffness(Component):
    """This component calculates the soil-pile equivalent stiffness matrix at mudline.\n
       NOTES:1. It uses either Matlock&Reese(1960) or Pender's (1993) approximations. \n
               In either case, Es is assumed varying linearly with depth. \n
             2. Use only if AF is not enabled!!! \n
    """

    #inputs
    SoilObj   =Instance(Klass=SoilC,iotype='in', desc='Object of Class Soil.')

    Lp        =Float(               iotype='in', desc='Pile embedment length.')

    PileObjout=Instance(Klass=Tube, iotype='in', desc='Object of Class Tube for Non-AF portion of Pile, valid also when ndiv=0.')
    VPFlag=    Bool(units=None,     iotype='in', desc='Vertical Pile Flag [T/F].')
    batter=    Float(units=None,    iotype='in', desc='2D batter [e.g., 10] for jacket legs hence piles.')
    innr_ang  =Float(units='rad',    iotype='in',desc='Angle between radial direction and base side, in rad.')
    H         =Float(0.,units='N',      iotype='in',desc='Shear at pile head (possibly above mudline), in N. Needed if Pender''s approx is used')
    M         =Float(0.,units='N*m',    iotype='in',desc='Moment at pile head (possibly above mudline), in N. Needed if Pender''s approx is used')
    loadZ     =Float(0.,units='m',      iotype='in',desc='Load level about mudline in meters. Default=0.=mudline.')

    #outputs
    SPI_Kmat  =   Array(np.array([]), dtype=np.float, iotype='out', desc='Soil-pile Equivalent Stiffness Matrix at mudline (1pile)')

    def execute(self):
        #pile properties from pile object
        Ep=self.PileObjout.mat[0].E
        Gp=self.PileObjout.mat[0].G
        Dp=self.PileObjout.D
        tp=self.PileObjout.t
        Jxx_p=self.PileObjout.Jxx
        #Start by calculating a coefficient of sibgrade reaction [N/m3] along theh pile length
        ks_avg=SubgrReact(self.SoilObj,self.Lp, sndflg=self.SoilObj.sndflg, bwtable=True)  #offshore pile is always below water table

        #Here is the Stiffness matrix
        self.SPI_Kmat=SoilPileStiffness(ks_avg,Dp,self.Lp,Ep,Gp,Jxx_p,loadZ=self.loadZ,PenderSwtch=self.SoilObj.PenderSwtch,sndflg=self.SoilObj.sndflg,H=self.H,M=self.M,batter=self.batter*(not(self.VPFlag)),psi=self.innr_ang)

#_____________________________________________________#


#_____________________________________________________#
#     BRACE CRITERIA/CONSTRAINTS TO BE SATISFIED      #
#_____________________________________________________#

class BrcCrOutputs(VariableTree):
    """Braces Criteria to be satisfied. These are also in function BrcFind, not used in optimization."""

    brc_crit01 = Array(dtype=np.float, desc='Cross-sectional Area-Ab>=0)')
    brc_crit02 = Array(dtype=np.float, desc='Dbrc/Dleg>=0.3')
    brc_crit03 = Array(dtype=np.float, desc='Dbrc/tbrc>=31.')
    brc_crit04 = Array(dtype=np.float, desc='brc/tbrc<=hydroconst')
    brc_crit05 = Array(dtype=np.float, desc='brc_Klr<=Klr_mx=70')

class BrcCriteria(Component):
    """This component check brace criteria. They used to be in function BrcFind, not used in optimization anymore."""

    #inputs
    nbays    = Int(iotype='in', units=None, desc='Number of Bays') #Optional, needed only for Xbrace criteria, not Mudbrace
    wdepth   = Float(           units='m',  iotype='in', desc='Water Depth')
    XBrcObj   = Instance(Klass=Tube,        iotype='in',desc='Object of tube class as coming from Component Xbraces')
    MudBrcObj = Instance(Klass=Tube,        iotype='in',desc='Object of tube class as coming from Component Mudbraces')
    LegObj    = Instance(Klass=Tube,        iotype='in',desc='Object of tube class as coming from Component Legs')

    #outputs
    MudBrcCriteria=VarTree(BrcCrOutputs(), iotype='out', desc='Basic Criteria Braces need to satisfy.')
    XBrcCriteria  =VarTree(BrcCrOutputs(), iotype='out', desc='Basic Criteria Braces need to satisfy.')

    def execute(self):  #Should add Andrew's stiffness calculation program
        #simp[lify nomenclature
        nbays=self.nbays
        wdepth=self.wdepth
        XBrcObj=self.XBrcObj
        MudBrcObj=self.MudBrcObj
        LegObj=self.LegObj

        Klr_mx=70.  #Maximum recommended Klr for braces

        hydrocnst=250./(wdepth*(1e3/25.4/12))**(1./3.) #to verify hydrostatic collapse

       #take care of leg properties first
        Dleg=LegObj.D[-nbays:] #ODs at various levels, excluding bottom stumps
        tleg=LegObj.t[-nbays:] #ts at various levels, excluding bottom stumps
        Aleg=LegObj.Area[-nbays:]

        #then take care of Xbrc and MudProperties: we may have ndiv>1 as well
        Xidx=range(0,XBrcObj.Area.size,XBrcObj.Area.size/nbays)
        #for mudbraces, just use the first value, i.e. index=0

        #Criteria #1: Area-Ab>=0
        Ab=0.1*Aleg
        self.XBrcCriteria.brc_crit01=(XBrcObj.Area[Xidx]-Ab)/Ab  #Non-Dimensionalized with Ab
        self.MudBrcCriteria.brc_crit01=np.asarray([(MudBrcObj.Area[0]-Ab[0])/Ab[0]])  #Non-Dimensionalized with Ab; make it into an array since the criteria are delcalred as arrays

        #Criteria #2: Dbrc/Dleg>=0.3
        self.XBrcCriteria.brc_crit02=(XBrcObj.D[Xidx]/Dleg)-0.3  #Non-Dimensional
        self.MudBrcCriteria.brc_crit02=np.asarray([(MudBrcObj.D[0]/Dleg[0])-0.3])  #Non-Dimensional

        #Criteria #3: Dbrc/tbrc>=31.
        self.XBrcCriteria.brc_crit03=(XBrcObj.D[Xidx]/XBrcObj.t[Xidx])-31.  #Non-Dimensional
        self.MudBrcCriteria.brc_crit03=np.asarray([(MudBrcObj.D[0]/MudBrcObj.t[0])-31.])  #Non-Dimensional

        #Criteria #4: Dbrc/tbrc<=hydroconst
        self.XBrcCriteria.brc_crit04=-(XBrcObj.D[Xidx]/XBrcObj.t[Xidx])+hydrocnst  #Non-Dimensional
        self.MudBrcCriteria.brc_crit04=np.asarray([-(MudBrcObj.D[0]/MudBrcObj.t[0])+hydrocnst])  #Non-Dimensional

        #Criteria #5: brc_Klr<=Klr_mx
        self.XBrcCriteria.brc_crit05  =-(XBrcObj.Klr[Xidx]-Klr_mx)/Klr_mx  #Non-Dimensionalized with Klr_mx
        self.MudBrcCriteria.brc_crit05=-np.array([(MudBrcObj.Klr[0]-Klr_mx)/Klr_mx ])#Non-Dimensionalized with Klr_mx

#        Lbrc2=bay_bs[0]-Dleg[1-VPFlag]  #Mud brace effective length
 #       nbrcs=1  #Initialize as if we were calculating only 1 brace (mudbrace or TopHbrace)

##    if not(mudbrc): #If it is a X-brace we are after
##       #1st bay height
##        h1=bay_hs[0]
##        #2nd bay bay width (top of bay 1 width)
##        wbas1=bay_bs[1]
##        #First bay leg segment, it works with 3 and 4 legs
##        AB=h1*np.sqrt(1.+np.tan(al_bat2D)**2 *(1.+np.tan(innr_ang)**2))
##        #Brace Length of 1st bay, and then the longest portion of it due to intersection: AO/AC=wbas0/(wbas0+wbas1)
##        Lbrc2=np.sqrt( h1**2 + (bay_bs[0]-h1*np.tan(al_bat2D))**2 +(h1*np.tan(al_bat2D)*np.tan(innr_ang))**2 ) * bay_bs[0]/(bay_bs[0]+wbas1)
##        nbrcs=nbays
##
##
##
##    #brace function of D,t; tried to make it work with x an array of nbays*2 length, but slsqp bombs out, so I need to use a loop
##    braces= lambda x: Tube(x[0],x[1],Lbrc2,Kbuck)
##
##    #Next is the main function to minimize (in our case all constraints are inequalities)
##    if minmass:
##        func=lambda x: braces(x).Area  #we are trying to minimize the mass=rho*Area*L, but rho and L are fixed
##        cons.append({'type':'ineq','fun': lambda x: Klr_mx-braces(x).Klr})
##        ieqcons.append(cons[4].get('fun'))
##    else:
##        func=lambda x: abs(braces(x).Klr-Klr_mx)  #we are trying to enforce a Klr=70 in case we do not use mass minimization

class SubDynReturns(Component): #CJB+
    SubDynOuts=VarTree(SubDynAOutputs(), iotype='in', desc='Selected outputs from SubDyn') #CJB+
    def execute(self): #CJB+
        print "SubDynReturns has run"
        self.SDFEMEigenvalues=self.SubDynOuts.FEMeigenvalues #CJB+
        self.SDStructureMass=self.SubDynOuts.StructureMass #CJB+ #Structural mass, does not include RNA lumped mass
        self.SDStructureCM=self.SubDynOuts.StructureCM #CJB+

#______________________________________________________________________________#
#                             RUN FRAME 3DD
#______________________________________________________________________________#

class Frame3DDOutputs(VariableTree):
    """General Frame3DD Output Data"""
    Freqs    = Array(  units='Hz',  desc='1st and 2nd EigenFrequency')
    mass     = Array(dtype=np.float,units='kg',  desc='Structural Mass and Total mass in the Jacket-Tower Structure: float(2)')
    forces   = Array(dtype=np.float,desc='Node Forces and Moments')
    reactions= Array(dtype=np.float,desc='Reaction Forces and Moments')
    top_deflection = Array(np.array([-9999.,-9999.,-9999.]),dtype=np.float,units='m', desc='Deflection of tower top in Global Coordinate System')
    Pileloads= Array(np.array([0.,0.,0.]),dtype=np.float,desc='(highest axial load-)Pile head (bottom of analyzed structure) Axial Force, Shear, and Bending. In the future this will be revised to account for correct geoemtry and orientation, or now it works with vertical pile')
    ElemL    = Array(                     dtype=np.float,desc='Lengths of all elements in the Frame3DD model.')
    ElemMass = Array(                     dtype=np.float,desc='Masses of all elements in the Frame3DD model.')

class RunFrame3DDstatic(Component):
    """This component runs Frame3DD for the jacket-tower assembly in static and dynamic modes. The modal analysis may be redone with stiffness constants at the base.
    """

    #inputs
    SPI_Kmat =Array(np.array([]),dtype=np.float,iotype='in', desc='Equivalent Stiffness Matrix at the pile head as calculated by pile SoilPileStiffness')

    JcktGeoOut = VarTree(JcktGeoOutputs(), iotype='in', desc='Geometry of the Jacket -Node Coordinates and Member Connectivity')
    Loadins    = VarTree(LoadOutputs(),    iotype='in', desc='Jacket Loading - Basic Outputs from Loads Component.')
    FrameAuxIn = VarTree(Frame3DDaux(),    iotype='in', desc='Auxiliary Data for Frame3DD')

    #outputs
    Frameouts=VarTree(Frame3DDOutputs(),   iotype='out', desc='Basic output data from Frame3DD')

    def execute(self):
        nds_frcs=self.Loadins.nds_frc
        nodes = NodeData(np.arange(1,self.JcktGeoOut.nNodesJckt+1), self.JcktGeoOut.nodes[:,0],self.JcktGeoOut.nodes[:,1], self.JcktGeoOut.nodes[:,2], self.JcktGeoOut.radii)

        reactions = ReactionData(self.JcktGeoOut.Reacts[:,0], self.JcktGeoOut.Reacts[:,1], self.JcktGeoOut.Reacts[:,2], self.JcktGeoOut.Reacts[:,3], \
                              self.JcktGeoOut.Reacts[:,4], self.JcktGeoOut.Reacts[:,5], self.JcktGeoOut.Reacts[:,6],1)
        if self.SPI_Kmat.size:  #Take care of the BCs if we are running it for the second time (old runframe3dd modal basically)
            nReacts=self.JcktGeoOut.Reacts[:,0].size

            reactions = ReactionData(self.JcktGeoOut.Reacts[:,0], self.SPI_Kmat[0,0].repeat(nReacts), self.SPI_Kmat[1,1].repeat(nReacts), self.SPI_Kmat[2,2].repeat(nReacts), \
                              self.SPI_Kmat[3,3].repeat(nReacts), self.SPI_Kmat[4,4].repeat(nReacts), self.SPI_Kmat[5,5].repeat(nReacts),-999)


        elements = ElementData(np.arange(1,self.JcktGeoOut.nmems+1), self.JcktGeoOut.mems[:,0], self.JcktGeoOut.mems[:,1],self.JcktGeoOut.props[:,0],\
                           self.JcktGeoOut.props[:,1], self.JcktGeoOut.props[:,2], self.JcktGeoOut.props[:,3],\
                           self.JcktGeoOut.props[:,4], self.JcktGeoOut.props[:,5], self.JcktGeoOut.props[:,6],\
                           self.JcktGeoOut.props[:,7] ,self.JcktGeoOut.props[:,8], self.JcktGeoOut.props[:,9],\
                           )

        options = Options(self.FrameAuxIn.sh_fg, self.FrameAuxIn.geo_fg, self.FrameAuxIn.deltaz)

        frame = Frame(nodes, reactions, elements, options)

        load = StaticLoadCase(self.FrameAuxIn.gvector[0], self.FrameAuxIn.gvector[1], self.FrameAuxIn.gvector[2])  #normally g must be given as along negative z, thus -9.81

        #Now add wind and hydro loads that have been calculated
        load.changePointLoads(nds_frcs[:,0], nds_frcs[:,1], nds_frcs[:,2], nds_frcs[:,3], nds_frcs[:,4], nds_frcs[:,5], nds_frcs[:,6])

        frame.addLoadCase(load)


        #Add modal analysis here
        addGravityLoad=True #This means it will add self-weight of extra point masses

        frame.enableDynamics(self.FrameAuxIn.nModes, self.FrameAuxIn.Mmethod, self.FrameAuxIn.lump, self.FrameAuxIn.tol, self.FrameAuxIn.shift)

                                            #joint ID,                    mass,                              IXX
        frame.changeExtraNodeMass(self.JcktGeoOut.jnt_masses_yaw[:,0], self.JcktGeoOut.jnt_masses_yaw[:,1], self.JcktGeoOut.jnt_masses_yaw[:,2],\
                                            #IXX,                         IYY                                 #IXY
                                   self.JcktGeoOut.jnt_masses_yaw[:,3] , self.JcktGeoOut.jnt_masses_yaw[:,4], self.JcktGeoOut.jnt_masses_yaw[:,5], \
                                           #IXZ                           IYZ                         CMxoff
                                   self.JcktGeoOut.jnt_masses_yaw[:,6],  self.JcktGeoOut.jnt_masses_yaw[:,7], self.JcktGeoOut.jnt_masses_yaw[:,8], \
                                           #CMyoff                                CMzoff
                                   self.JcktGeoOut.jnt_masses_yaw[:,9],  self.JcktGeoOut.jnt_masses_yaw[:,10],addGravityLoad)

        displacements, self.Frameouts.forces, self.Frameouts.reactions, internalForces, mass, modal = frame.run()
#        displacements, self.Frameouts.forces, self.Frameouts.reactions, internalForces, mass,junk = frame.run()

        self.Frameouts.mass=np.array([mass.struct_mass,mass.total_mass])
        self.Frameouts.Freqs=modal.freq[0:2]  #First 2 eigenfreqs
        self.Frameouts.top_deflection=np.array([displacements.dx[0,-1],displacements.dy[0,-1],displacements.dz[0,-1]])

        #Pile head loads: NOTE THIS ASSUME VERTICAL PILES, FOR BATTERED I WOULD HAVE TO DECOMPOSE WITH COS MATRIX
        reactions=self.Frameouts.reactions.squeeze().T #[nreacts,X,Y,Z,XX,YY,ZZ]
        idx=np.argmax(np.abs(reactions[:,3]))#Maximum axial load on pile
        Nhead=np.abs(reactions[idx,3])  #We will assume this as a compressive load to do check on pile
        Thead=np.sqrt(reactions[idx,1]**2+reactions[idx,2]**2) #Shear for highest loaded pile, assumed later along the unfavorable direction as wind and wave
        Mhead=np.sqrt(reactions[idx,4]**2+reactions[idx,5]**2) #Bending moment for highest loaded pile, we are ignoring torsion and direction cosines for now ; we assume it to generate same directional displacement as Thead

        #ReactNds=self.JcktGeoOut.Reacts[:,0]  #reaction nodes
        #MbrFrcs=self.Frameouts.forces.squeeze().T
        #Reacts=MbrFrcs[np.in1d(MbrFrcs[:,1],ReactNds),2]  #Fx [N] aligned with member
        #Nhead=np.max(np.abs(Reacts))

        self.Frameouts.Pileloads=np.array([Nhead,Thead,Mhead])

        #Store mass for each element
        DeltaX2=(frame.nx[elements.N2-1]-frame.nx[elements.N1-1])**2
        DeltaY2=(frame.ny[elements.N2-1]-frame.ny[elements.N1-1])**2
        DeltaZ2=(frame.nz[elements.N2-1]-frame.nz[elements.N1-1])**2
        self.Frameouts.ElemL = np.sqrt( DeltaX2+DeltaY2+DeltaZ2)
        self.Frameouts.ElemMass =elements[-1]* frame.eAx * self.Frameouts.ElemL

        print(str(self.Frameouts.mass)+' self.Frameouts.mass (Structural Mass and Total mass in the Jacket-Tower Structure)') #CJB+
#______________________________________________________________________________#
#        Auxiliary Component to Calculate ABS
class ABSaux(Component):
    """This is an auxiliary component that hopefully we will be able to remove once OpenMDAO gets fixed to get the correct connections when operations are specified"""
    #ins
    varin  = Float(  iotype='in',desc='Variable to be ABSed.')
    #outs
    varout = Float(  iotype='in',desc='ABSed Variable')

    def execute(self):
       self.varout=np.abs(self.varin)
#______________________________________________________________________________#
#        Auxiliary Component to handle Frameouts appropriately
class FrameOutsaux(Component):
    """This is an auxiliary component that hopefully we will be able to remove once OpenMDAO gets fixed to properly handle boundary variables"""
    #ins
    Frameouts_ins  = VarTree(Frame3DDOutputs(),   iotype='in',  desc='Basic output data from Frame3DD')
    #ins
    Frameouts_outs = VarTree(Frame3DDOutputs(),   iotype='out', desc='Basic output data from Frame3DD')

    def execute(self):
       self.Frameouts_outs=self.Frameouts_ins


#______________________________________________________________________________#




#______________________________________________________________________________#
#                                 Assembly
#______________________________________________________________________________#
class LoadFrameOuts(Assembly):
    """This assembly wraps the components that are needed to carry out static and modal \n
        analyses and returns needed outputs. It is created so that 2 instances can be used \n
        for 2 separate DLCs."""
    from Utilization import JacketUtilOutputs,TowerUtilOutputs
    #ins
    Waterinputs  = VarTree(WaterInputs(),  iotype='in', desc='Water Data')
    Windinputs   = VarTree(WindInputs(),   iotype='in', desc='Wind Data')

    RNA_F       = Array(np.zeros(6), iotype='in', desc='DLC 1.6: Unfactored Rotor Forces and Moments, excluding weight. Array(6).')
    RNAinputs   = VarTree(RNAprops(),iotype='in', desc='Basic Inertial Properties of RNA.')
    VPFlag      = Bool(units=None,   iotype='in', desc='Vertical Pile Flag [T/F].')

    nlegs       = Int(  iotype='in', units=None, desc='Number of Legs.')  # comes from JcktGeoIn.nlegs
    nodes       = Array(iotype='in',             desc='Node coordinates.')  # comes from JcktGeoOut.nodes

    al_bat3D    = Float(iotype='in', units='rad',desc='Batter Angle in 3D.')
    TwrRigidTop = Bool( iotype='in', units=None, desc='Rigid Member used in tower top or not.')

    pillegDs    = Array(iotype='in',dtype=np.float, units='m',  desc='ODs of pile and leg #1.')
    twrDs       = Array(iotype='in',                units='m',  desc='ODs of Tower.')

    Twrmems  = Array(dtype=int,  iotype='in', desc='Connectivity Array for all members of the tower portion only (including rigid member at the top if requested).')
    Legmems  = Array(dtype=int,  iotype='in', desc='Connectivity Array for all members of the Leg 1 portion only.')
    Pilemems = Array(dtype=int,  iotype='in', desc='Connectivity Array for all members of the Pile 1 portion only. ')
    TPmems   = Array(dtype=int,  iotype='in', desc='Connectivity Array for all members of the TP portion only.')
    XjntIDs  = Array(dtype=int,  iotype='in', desc='X-Joint IDs, as for Frame3DD, used to do checks later.')
    KjntIDs  = Array(dtype=int,  iotype='in', desc='K-Joint IDs, as for Frame3DD, used to do checks later.')

    gravity  = Float(units='m/s**2', iotype='in',desc='Gravity Acceleration (ABSOLUTE VALUE!).')

    FrameAuxIns= VarTree(Frame3DDaux(),    iotype='in', desc='Auxiliary Data for Frame3DD.')

    JcktGeoOut = VarTree(JcktGeoOutputs(), iotype='in', desc='Geometry of the Jacket -Node Coordinates and Member Connectivity.')
    SoilObj    = Instance(Klass=SoilC,     iotype='in', desc='Object of Class Soil.')
    batter     = Float(units=None,         iotype='in', desc='2D batter [e.g., 10] for jacket legs hence piles.')
    Lp         = Float(units='m',          iotype='in', desc='Pile embedment length.')

    PileObjout = Instance(Klass=Tube, iotype='in', desc='Object of Class Tube for Non-AF portion of Pile, valid also when ndiv=0.')

    innr_ang   = Float(units='rad',   iotype='in',desc='Angle between radial direction and base side, in rad.')

    loadZ     = Float(0.,units='m',      iotype='in',desc='Load level about mudline in meters. Default=0.=mudline.')

    Twr_data     = VarTree(TwrGeoOutputs(), iotype='in', desc='Tower Node data.')
    Dt           = Float( units='m',        iotype='in', desc='TowerTop OD from Tower.')
    tt           = Float( units='m',        iotype='in', desc='TowerTop wall thickness from Tower.')
    L_reinforced = Float( units='m',        iotype='in', desc='reinforcement length.')

    IECpsfIns = VarTree(IEC_PSFS(),      iotype='in', desc='Basic IEC psfs.')

    legjnts = Array(iotype='in', desc='Leg joint coordinates.')  #CJB+ Comes from legs.legjnts

    #outs
    Frameouts = VarTree(Frame3DDOutputs(),           iotype='out', desc='Basic output data from Frame3DD.')
    SPI_Kmat  = Array(np.array([]), dtype=np.float,  iotype='out', desc='Soil-pile Equivalent Stiffness Matrix at mudline (1pile).')
    Lp0rat    = Float(units=None,                    iotype='out', desc='(Lp-Lp0)/Lp0 Normalized difference between current and needed embedment length.')  #necessary to sustain load
    Mpiles  =   Float(units='kg',                    iotype='out', desc='Embedded Pile Mass (all piles) with the given Lp, not necessarily = to needed embedment length.')
    tower_utilization  = VarTree(TowerUtilOutputs(), iotype='out', desc='Tower Utilization Basic Outputs.')
    jacket_utilization = VarTree(JacketUtilOutputs(),iotype='out', desc='Jacket Utilization Basic Outputs.')

    def __init__(self, clamped):

        self.clamped = clamped
        super(LoadFrameOuts,self).__init__()

    def configure(self):

        self.add('Loads', JcktLoad())
        self.add('Frame3DD', RunFrame3DDstatic())
        self.add('SPIstiffness', SPIstiffness())
        self.add('Frame3DD2', RunFrame3DDstatic())  #same component as before, we just use it again in case we have to account for non-clamped conditions
        self.add('Utilization', UtilAssembly())
        self.add('Embedment', PileEmbdL())

        self.driver.workflow.add(['Loads','Frame3DD'])

        if not (self.clamped):
            self.driver.workflow.add(['SPIstiffness','Frame3DD2'])

        self.driver.workflow.add(['Utilization','Embedment'])

        # Loads
        self.connect('Waterinputs',                'Loads.waterIns')
        self.connect('Windinputs',                 'Loads.windIns')
        self.connect('RNA_F',                      'Loads.RNA_F')
        self.connect('RNAinputs',                  'Loads.RNAinputs' )

        self.connect('VPFlag',                     'Loads.VPFlag')
        self.connect('pillegDs',                   'Loads.pillegDs')
        self.connect('twrDs',                      'Loads.twrDs')
        self.connect('Pilemems',                   'Loads.Pilemems')
        self.connect('Legmems',                    'Loads.Legmems')
        self.connect('Twrmems',                    'Loads.Twrmems')

        self.connect('al_bat3D',                   'Loads.al_bat3D')
        self.connect('nlegs' ,                     'Loads.nlegs')
        self.connect('nodes',                      'Loads.nodes')

        self.connect('gravity',                    'Loads.gravity')

        self.connect('TwrRigidTop',                'Loads.TwrRigidTop')

        # Run Frame3DD (Clamped)
        self.connect('FrameAuxIns',          'Frame3DD.FrameAuxIn')
        self.connect('JcktGeoOut',           'Frame3DD.JcktGeoOut')
        self.connect('Loads.Loadouts',       'Frame3DD.Loadins')

        # Calculate Pile Stiffness if SPI is activated
        if not self.clamped:
            self.connect('SoilObj',                        'SPIstiffness.SoilObj')
            self.connect('VPFlag',                         'SPIstiffness.VPFlag')
            self.connect('batter',                         'SPIstiffness.batter')
            self.connect('Lp',                             'SPIstiffness.Lp')
            self.connect('PileObjout',                     'SPIstiffness.PileObjout')
            self.connect('innr_ang',                       'SPIstiffness.innr_ang')
            self.connect('loadZ',                          'SPIstiffness.loadZ')
            self.connect('Frame3DD.Frameouts.Pileloads[1]','SPIstiffness.H')
            self.connect('Frame3DD.Frameouts.Pileloads[2]','SPIstiffness.M')

            #Run Frame3DD 2
            self.connect('FrameAuxIns',          'Frame3DD2.FrameAuxIn')
            self.connect('JcktGeoOut',           'Frame3DD2.JcktGeoOut')
            self.connect('Loads.Loadouts',       'Frame3DD2.Loadins')
            self.connect('SPIstiffness.SPI_Kmat','Frame3DD2.SPI_Kmat')

        # Utilization
            #jacket
        self.connect('Pilemems',  'Utilization.Pilemems')
        self.connect('Legmems',   'Utilization.Legmems')
        self.connect('Twrmems',   'Utilization.Twrmems')
        self.connect('TPmems',    'Utilization.TPmems')
        self.connect('nlegs',     'Utilization.nlegs')
        self.connect('JcktGeoOut','Utilization.JcktGeoOut')

        self.connect('XjntIDs',                    'Utilization.XjntIDs')
        self.connect('KjntIDs',                    'Utilization.KjntIDs')

        self.connect('RNA_F[0:3]',                 'Utilization.RNA_F')
        self.connect('RNA_F[3:6]',                 'Utilization.RNA_M')
        self.connect('RNAinputs.rna_weightM',      'Utilization.rna_weightM')
        self.connect('gravity',                    'Utilization.g')

        self.connect('Loads.windLoads',            'Utilization.towerWindLoads')
        self.connect('Loads.waveLoads',            'Utilization.towerWaveLoads')

        self.connect('Twr_data',                   'Utilization.Twr_data')
        self.connect('Dt',                         'Utilization.Dt')
        self.connect('tt',                         'Utilization.tt')
        self.connect('L_reinforced',               'Utilization.L_reinforced' )

        self.connect('IECpsfIns',                  'Utilization.IECpsfIns')
        #self.connect('Tower.Twrouts.rna_yawedcm','Utilization.r_cm')  #Yawed coordinated of CM
        self.connect('RNAinputs.CMoff',            'Utilization.r_cm')
        self.connect('RNAinputs.Thoff',            'Utilization.r_hub')
        #self.connect('Tower.Twrouts.Thoff_yaw','Utilization.r_hub')
        self.Utilization.tilt = 0.0

        # Embedment calculation
        self.connect('PileObjout',              'Embedment.PileObjout')
        self.connect('SoilObj',                 'Embedment.SoilObj')
        self.connect('JcktGeoOut',              'Embedment.JcktGeoOut')

        if self.clamped:
            self.connect('Frame3DD.Frameouts.forces','Utilization.MbrFrcs')
            self.connect('Frame3DD.Frameouts.forces', 'Embedment.MbrFrcs')
            #output
            self.connect('Frame3DD.Frameouts','Frameouts')
        else:
            self.connect('Frame3DD2.Frameouts.forces','Utilization.MbrFrcs')
            self.connect('Frame3DD2.Frameouts.forces','Embedment.MbrFrcs')
            #output
            self.connect('Frame3DD2.Frameouts','Frameouts')
            self.connect('SPIstiffness.SPI_Kmat','SPI_Kmat')


        self.connect('gravity',            'Embedment.gravity')  #Need absolute value here, Embedment takes care of it anyway, but somehow this abs does not trigger errors in optimzer as opposed to others
        self.connect('Lp',                 'Embedment.Lp')


        #more outputs
        self.connect('Utilization.tower_utilization','tower_utilization')
        self.connect('Utilization.jacket_utilization','jacket_utilization')
        self.connect('Embedment.Mpiles','Mpiles')
        self.connect('Embedment.Lp0rat','Lp0rat')

class JacketSE(Assembly):
    gravity  = Float(units='m/s**2', iotype='in',desc='Gravity Acceleration (ABSOLUTE VALUE!).') #CJB+ test
    dck_botz=Float(units='m', iotype='in', desc='Distance from bottom of TP to MSL') #CJB+ test
    nodes=Array(iotype='in', desc='Nodes from JacketSE') #CJB+ test
    mems=Array(iotype='in',desc='Member connections between') #CJB+ test
    Reacts=Array(iotype='in', desc='Base joint IDs and fixity') #CJB+ test
    SDPropSet=Array(iotype='in', desc='Properties to be used in SubDyn') #CJB+ test

    #ins
    JcktGeoIn=   VarTree(JcktGeoInputs(),  iotype='in', desc='Jacket Geometry Basic Inputs')
    leginputs=   VarTree(LegGeoInputs(),   iotype='in', desc="Leg Input Data")
    Xbrcinputs=  VarTree(XBrcGeoInputs(),  iotype='in', desc="Xbrace Input Data")
    Mbrcinputs=  VarTree(MudBrcGeoInputs(),iotype='in',desc="Mud Brace Input Data")
    Hbrcinputs=  VarTree(HBrcGeoInputs(),  iotype='in', desc="Top H-Brace Input Data")
    Pileinputs=  VarTree(PileGeoInputs(),  iotype='in', desc="Pile Input Data")
    TPinputs=    VarTree(TPGeoInputs(),    iotype='in', desc='Basic input data for Transition Piece')
    FrameAuxIns= VarTree(Frame3DDaux(),    iotype='in', desc='Auxiliary Data for Frame3DD')

    Twrinputs  =   VarTree(TwrGeoInputs(), iotype='in', desc='Basic Input data for Tower')
    TPlumpinputs=  VarTree(TPlumpMass(),   iotype='in', desc='Basic Inertial Properties of TP Lumped Mass')

    RNAinputs=     VarTree(RNAprops(),     iotype='in', desc='Basic Inertial Properties of RNA, in case modified for DLC 1.6')
    RNAinputs2=    VarTree(RNAprops(),     iotype='in', desc='Basic Inertial Properties of RNA, in case modified for DLC 6.1')

    Waterinputs  = VarTree(WaterInputs(),  iotype='in', desc='Water Data for DLC 1.6.')
    Windinputs   = VarTree(WindInputs(),   iotype='in', desc='Wind Data for DLC 1.6.')
    Soilinputs   = VarTree(SoilGeoInputs(),iotype='in', desc='Soil Data for DLC 1.6.')
    Waterinputs2  = VarTree(WaterInputs(),  iotype='in', desc='Water Data for DLC 6.1.')
    Windinputs2   = VarTree(WindInputs(),   iotype='in', desc='Wind Data for DLC 6.1.')
    Soilinputs2   = VarTree(SoilGeoInputs(),iotype='in', desc='Soil Data for DLC 6.1.')

    RNA_F       =Array(np.zeros(6), iotype='in', desc='DLC 1.6: Unfactored Rotor Forces and Moments, excluding weight. Array(6)')
    RNA_F2      =Array(np.zeros(6), iotype='in', desc='DLC 6.1: Unfactored Rotor Forces and Moments, excluding weight. Array(6)')

    IECpsfIns = VarTree(IEC_PSFS(), iotype='in', desc='Basic IEC psfs')

    TwrRigidTop =Bool( units=None,iotype='in',desc='Rigid Member used in tower top or not')

    SubDynOuts=VarTree(SubDynAOutputs(), iotype='in', desc='Selected outputs from SubDyn') #CJB+

    #--------------------------------------------------------------------------------------------------------------------

    #INPUT FILE INPUTS----------------------------------------------------------

    InputFile_name=Str(iotype='in') #CJB+

    #Simulation Control
    Echo=Array(iotype='in') #CJB+
    SDdeltaT=Array(iotype='in') #CJB+
    IntMethod=Array(iotype='in') #CJB+
    SttcSolve=Array(iotype='in') #CJB+

    #FEA and CRAIG-BAMPTON PARAMETERS
    FEMmod=Array(iotype='in') #CJB+
    NDiv=Array(iotype='in') #CJB+
    CBMod=Array(iotype='in') #CJB+
    Nmodes=Array(iotype='in') #CJB+
    JDampings=Array(iotype='in') #CJB+

    #Structure Joints
    SDjointsHeader=Array(iotype='in') #CJB+
    SDjoints=Array(iotype='in', desc='Leg joints from jacket.py. In SD joints are used to build the structure, not nodes.') #CJB+

    #Base Reaction Joints
    BaseRxnJointsHeader=Array(iotype='in') #CJB+
    BaseRxnJoints=Array(iotype='in', desc='Base joints of the structure.') #CJB+

    #Interface Joints
    InterfaceRxnJointsHeader=Array(iotype='in') #CJB+
    InterfaceJointsFlags=Array(iotype='in', desc='Interface joints.') #CJB+

    #Members
    MembersHeader=Array(iotype='in') #CJB+
    Members=Array(iotype='in', desc='Build members from joints, not nodes.') #CJB+

    #MEMBER X-SECTION PROPERTY data 1/2
    NPropSets=Array(iotype='in') #CJB+
    PropSet1Header=Array(iotype='in') #CJB+
    PropSet1=Array(iotype='in') #CJB+

    #MEMBER X-SECTION PROPERTY data 2/2
    PropSet2Header=Array(iotype='in') #CJB+
    PropSet2=Array(iotype='in') #CJB+

    #MEMBER COSINE MATRICES COSM(i,j)
    COSMHeader=Array(iotype='in') #CJB+
    COSMs=Array(iotype='in') #CJB+

    #JOINT ADDITIONAL CONCENTRATED MASSES
    CmassHeader=Array(iotype='in') #CJB+
    Cmass=Array(iotype='in') #CJB+

    #OUTPUT: SUMMARY & OUTFILE
    SSSum=Array(iotype='in') #CJB+
    OutCOSM=Array(iotype='in') #CJB+
    OutAll=Array(iotype='in') #CJB+
    OutSwtch=Array(iotype='in') #CJB+
    TabDelim=Array(iotype='in') #CJB+
    OutDec=Array(iotype='in') #CJB+
    OutFmt=Array(iotype='in') #CJB+
    OutSFmt=Array(iotype='in') #CJB+

    #MEMBER OUTPUT LIST
    MemOutListHeader=Array(iotype='in') #CJB+
    MemOutList=Array(iotype='in', desc='Members whose loads and dynamics will be output.') #CJB+

    #SSOutline
    SSOutlist=Array(iotype='in') #CJB+

    #DRIVER INPUTS--------------------------------------------------------------

    InputandDriverpath=Str(iotype='in') #CJB+
    Driver_path=Str(iotype='in') #CJB+

    EchoD=Array(iotype='in') #CJB+

    #Environmental Conditions
    Gravity=Array(iotype='in') #CJB+
    WtrDpth=Array(iotype='in') #CJB+

    #SubDyn
    SDInputFile=Array(iotype='in') #CJB+
    OutRootName=Array(iotype='in') #CJB+
    NSteps=Array(iotype='in') #CJB+
    TimeInterval=Array(iotype='in') #CJB+
    TP_RefPoint=Array(iotype='in') #CJB+
    SubRotateZ=Array(iotype='in') #CJB+

    #INPUTS
    InputsMod=Array(iotype='in') #CJB+
    InputsFile=Array(iotype='in') #CJB+

    #STEADY INPUTS
    uTPInSteady=Array(iotype='in') #CJB+
    uDotTPInSteady=Array(iotype='in') #CJB+
    uDotDotTPInSteady=Array(iotype='in') #CJB+

    #INPUTS TO RUN SUBDYN-------------------------------------------------------

    SDpath=Str(iotype='in') #CJB+

    #INPUTS TO READ OUTPUT------------------------------------------------------

    Readpath_out=Str(iotype='in') #CJB+
    Readpath_sum=Str(iotype='in') #CJB+
    Delete_file=Bool(iotype='in') #CJB+
    InputFile_path=Str(iotype='in') #CJB+
    Driver_path=Str(iotype='in') #CJB+

#--------------------------------------------------------------------------------------------------------------------

    #outs
    #JcktGeoOut = VarTree(JcktGeoOutputs(), iotype='out', desc='Geometry of the Jacket -Node Coordinates and Member Connectivity')
    Twrouts  = VarTree(TwrGeoOutputs(), iotype='out', desc='Basic Output data for Tower')
    Legouts  = VarTree(LegGeoOutputs(), iotype='out', desc='Leg Geometry Output')

    def __init__(self, clamped=True, AFflag=False, twodlcs=True):  #PrebuildTP=True,

        self.clamped = clamped  #initialize to true
        self.AFflag = AFflag
        ##self.PrebuildTP = PrebuildTP
        #PrebuildTP=  Bool(True, iotype='in', desc='If TRUE then the TP will have dimensions connected to those of legs and braces, and some TP inputs will be overwritten.')
        self.twodlcs=twodlcs
        super(JacketSE,self).__init__()

        if not (self.clamped or self.AFflag):
            self.clamped=False #use this boolean to reduce checks below

    def configure(self):

        #self.add('Pileinputs', VarTree(PileGeoInputs(),iotype='in', desc="Pile Input Data"))
        self.add('PreLeg', PreLegBuild())
        self.add('PreBuild', PreJcktBuild())
        self.add('Soil', Soil())
        self.add('Legs', legs())
        self.add('Piles', Piles())
        self.add('Xbraces', Xbraces())
        self.add('Mudbraces', MudBraces())
        self.add('Hbraces', HBraces())
        self.add('TP', TP())
##        self.add('PreBuildTP', PreBuildTP())
##        self.add('PreBuildTP2', PreBuildTP2())
        self.add('Tower', Tower())
        self.add('Build', BuildGeometry())
        self.add('LoadFrameOuts', LoadFrameOuts(self.clamped))
        self.add('LoadFrameOuts2', LoadFrameOuts(self.clamped))
        self.add('TotMass',ComponentMass())
        self.add('BrcCriteria',BrcCriteria())

        self.add('ABS1',ABSaux())  #aux component
        self.add('ABS2',ABSaux())  #aux component
        self.add('FrameOut',FrameOutsaux()) #auxiliary component

        #pySubDyn CJB+
        self.add('SDpySubDynA', pySubDynA())              #CJB+
        self.add('SDReturns', SubDynReturns()) #CJB+

        # BUILD UP THE DRIVER
        casey = False
        if casey:
           self.driver.workflow.add(['PreLeg','PreBuild', 'Soil', 'Legs', 'Piles', 'Xbraces', 'Mudbraces', 'Hbraces', 'SDpySubDynA', 'SDReturns']) #CJBe
        else:
           self.driver.workflow.add(['PreLeg','PreBuild', 'Soil', 'Legs', 'Piles', 'Xbraces', 'Mudbraces', 'Hbraces'])

##        if self.PrebuildTP:
##            self.driver.workflow.add(['PreBuildTP','PreBuildTP2'])

        self.driver.workflow.add(['TP','Tower','Build','ABS1','LoadFrameOuts'])

        if self.twodlcs:
            self.driver.workflow.add('LoadFrameOuts2')

        self.driver.workflow.add(['FrameOut','BrcCriteria','TotMass' ])
        #________________________#

        #Take care of Connections

        # prelegs
        self.connect('JcktGeoIn.nbays',           'PreLeg.nbays' )
        self.connect('leginputs',             'PreLeg.leginputs' )

        # PreJckBuild
        self.create_passthrough('PreBuild.legbot_stmphin')

        self.connect('Waterinputs', 'PreBuild.WaterNew') #CJB+ water flag

        self.connect('JcktGeoIn',          'PreBuild.JcktGeoIn')
        self.connect('PreLeg.prelegouts.legZbot',  'PreBuild.legZbot')
        self.connect('PreLeg.prelegouts.Dleg[0]' , 'PreBuild.LegbD')
        self.connect('PreLeg.prelegouts.Dleg[-1]', 'PreBuild.LegtD')
        self.connect('PreLeg.prelegouts.tleg[-1]', 'PreBuild.Legtt')

        self.connect('Xbrcinputs',         'PreBuild.Xbrcinputs' )
        self.connect('JcktGeoIn.PreBuildTPLvl', 'PreBuild.PreBuildTPLvl')
        self.connect('Twrinputs.Db'      ,      'PreBuild.TwrDb')
        self.connect('Twrinputs.DTRb'      ,    'PreBuild.TwrDTRb')
        self.connect('TPinputs',                'PreBuild.TPinputs' )
        self.connect('PreBuild.BuildTPouts',    'TP.TPinputs' ) #this provides upgraded and overwritten inputs to TP

        #Soil
        self.connect('Soilinputs',   'Soil.SoilIns' )

        # legs
        self.connect('PreLeg.prelegouts',     'Legs.leginputs' )
        self.connect('JcktGeoIn',             'Legs.JcktPrms')
        self.connect('PreBuild.bay_hs',       'Legs.bay_hs' )
        self.connect('PreBuild.JcktH',        'Legs.JcktH' )
        self.connect('PreBuild.innr_ang',     'Legs.innr_ang')
        self.connect('PreBuild.dck_width',    'Legs.dck_width')  #updates deck_width with what PreBuild Calculates
        self.connect('PreBuild.legbot_stmph', 'Legs.legbot_stmph')
        self.connect('PreBuild.WaterNew.wdepth',    'Legs.wdepth') #CJB+ water flag
        self.connect('PreBuild.WaterNew.wlevel',    'Legs.wlevel') #CJB+ water flag
        self.connect('PreBuild.WaterNew.z_floor',    'Legs.z_floor') #CJB+ water flag

        # Xbraces
        ##self.connect('Xbrcinputs',         'Xbraces.Xbrcinputs' )
        self.connect('PreBuild.Xbrcouts',  'Xbraces.Xbrcinputs' )

        self.connect('JcktGeoIn.nbays',    'Xbraces.nbays')
        self.connect('Legs.legouts.joints','Xbraces.legjnts')
        self.connect('PreBuild.WaterNew.wdepth', 'Xbraces.wdepth') #CJB+ water flag
        self.connect('PreBuild.bay_bs',    'Xbraces.bay_bs')
        self.connect('PreBuild.bay_hs',    'Xbraces.bay_hs')
        self.connect('PreBuild.al_bat2D',  'Xbraces.al_bat2D')
        self.connect('PreBuild.innr_ang',  'Xbraces.innr_ang')
        self.connect('Legs.legouts.LegObj','Xbraces.LegObj')

        # MudBraces
        self.connect('Mbrcinputs',         'Mudbraces.Mbrcinputs' )
        self.connect('Legs.legouts.joints','Mudbraces.legjnts')
        self.connect('JcktGeoIn.VPFlag',   'Mudbraces.VPFlag')
        self.connect('Piles.PileFlag',     'Mudbraces.PileFlag')
        self.connect('PreBuild.WaterNew.wdepth', 'Mudbraces.wdepth') #CJB+ water flag
        self.connect('PreBuild.bay_bs',    'Mudbraces.bay_bs')
        self.connect('Legs.legouts.LegObj','Mudbraces.LegObj')

        # HBraces
        self.connect('Hbrcinputs',         'Hbraces.Hbrcinputs' )
        self.connect('Legs.legouts.joints','Hbraces.legjnts')
        self.connect('Xbraces.HbrcD',      'Hbraces.HbrcD')
        self.connect('Xbraces.Hbrct',      'Hbraces.Hbrct')

        # TP and prebuildTP
        self.connect('Legs.legouts.joints','TP.legjnts')
        self.connect('TPlumpinputs',       'TP.TPlumpinputs')


###        if self.PrebuildTP:#For a version with TP struts =leg OD, also TP stem OD =Tower Db, then diagonal braces and horizontal braces fixed with the top Xbraces
##            #self.connect('JcktGeoIn.PreBuildTPLvl', 'PreBuildTP.PreBuildTPLvl')
##            #self.connect('Legs.legouts.LegObj',     'PreBuildTP.LegObj')
##            #self.connect('Xbraces.Xbrcouts.LLURObj','PreBuildTP.XBrcObj')
##            #self.connect('Twrinputs.Db'      ,      'PreBuildTP.TwrDb')
##            #self.connect('Twrinputs.DTRb'      ,    'PreBuildTP.TwrDTRb')
##            #self.connect('TPinputs',                'PreBuildTP.TPinputs' )
##            #self.connect('PreBuildTP.BuildTPouts',  'TP.TPinputs' ) #this provides upgraded and overwritten inputs to TP
##
###        else:
###            self.connect('TPinputs',                'TP.TPinputs' )

        # Piles
        self.connect('Pileinputs',                 'Piles.Pileinputs' )
        self.connect('JcktGeoIn.AFflag',           'Piles.AFflag')
        self.connect('Legs.legouts.joints',        'Piles.legjnts')
        self.connect('PreLeg.prelegouts.Dleg[0]',  'Piles.LegbD')
        self.connect('PreLeg.prelegouts.tleg[0]',  'Piles.Legbt')
        self.connect('JcktGeoIn.VPFlag',           'Piles.VPFlag')
        self.connect('JcktGeoIn.nlegs' ,           'Piles.nlegs')

        # Tower
        self.connect('TwrRigidTop',        'Tower.RigidTop')
        self.connect('PreBuild.WaterNew.wdepth', 'Tower.wdepth') #CJB+ water flag
        self.connect('Twrinputs',          'Tower.Twrins' )
        self.connect('RNAinputs',          'Tower.RNAinputs' )
        self.connect('TP.TwrBsZ',          'Tower.TwrBsZ')
        self.connect('Windinputs.HH',      'Tower.HH' )   #HUBHEIGHT FROM INPUTS FOR WIND

        # BuildGeometry
        self.connect('TP.TPouts',         'Build.TPouts')
        self.connect('Tower.Twrouts',     'Build.Twrouts')
        self.connect('Piles.Pileouts',    'Build.Pileouts')
        self.connect('Legs.legouts',      'Build.Legouts')
        self.connect('Xbraces.Xbrcouts',  'Build.Xbrcouts')
        self.connect('Mudbraces.Mbrcouts','Build.Mbrcouts')
        self.connect('Hbraces.Hbrcouts',  'Build.Hbrcouts')
        self.connect('JcktGeoIn',         'Build.JcktGeoIn')

        #BrcCriteria
        self.connect('JcktGeoIn.nbays' ,         'BrcCriteria.nbays')
        self.connect('Legs.legouts.LegObj',      'BrcCriteria.LegObj')
        self.connect('Mudbraces.Mbrcouts.brcObj','BrcCriteria.MudBrcObj')
        self.connect('Xbraces.Xbrcouts.LLURObj', 'BrcCriteria.XBrcObj')
        self.connect('PreBuild.WaterNew.wdepth',       'BrcCriteria.wdepth') #CJB+ water flag

        #SubDyn #CJB+
        if casey:
            self.connect('InputFile_name',          'SDpySubDynA.InputFile_name') #CJB+
            self.connect('Echo',                    'SDpySubDynA.Echo') #CJB+
            self.connect('SDdeltaT',                'SDpySubDynA.SDdeltaT') #CJB+
            self.connect('IntMethod',               'SDpySubDynA.IntMethod') #CJB+
            self.connect('SttcSolve',               'SDpySubDynA.SttcSolve') #CJB+
            self.connect('FEMmod',                  'SDpySubDynA.FEMmod') #CJB+
            self.connect('NDiv',                    'SDpySubDynA.NDiv') #CJB+
            self.connect('CBMod',                   'SDpySubDynA.CBMod') #CJB+
            self.connect('Nmodes',                  'SDpySubDynA.Nmodes') #CJB+
            self.connect('JDampings',               'SDpySubDynA.JDampings') #CJB+
            self.connect('SDjointsHeader',          'SDpySubDynA.SDjointsHeader') #CJB+
            self.connect('SDjoints',                'SDpySubDynA.SDjoints') #CJB+
            self.connect('BaseRxnJointsHeader',     'SDpySubDynA.BaseRxnJointsHeader') #CJB+
            self.connect('BaseRxnJoints',           'SDpySubDynA.BaseRxnJoints') #CJB+
            self.connect('InterfaceRxnJointsHeader','SDpySubDynA.InterfaceRxnJointsHeader') #CJB+
            self.connect('InterfaceJointsFlags',    'SDpySubDynA.InterfaceJoints') #CJB+
            self.connect('MembersHeader',           'SDpySubDynA.MembersHeader') #CJB+
            self.connect('Members',                 'SDpySubDynA.Members') #CJB+
            self.connect('NPropSets',               'SDpySubDynA.NPropSets') #CJB+
            self.connect('PropSet1Header',          'SDpySubDynA.PropSet1Header') #CJB+
            self.connect('PropSet1',                'SDpySubDynA.PropSet1') #CJB+
            self.connect('PropSet2Header',          'SDpySubDynA.PropSet2Header') #CJB+
            self.connect('PropSet2',                'SDpySubDynA.PropSet2') #CJB+
            self.connect('COSMHeader',              'SDpySubDynA.COSMHeader') #CJB+
            self.connect('COSMs',                   'SDpySubDynA.COSMs') #CJB+
            self.connect('CmassHeader',             'SDpySubDynA.CmassHeader') #CJB+
            self.connect('Cmass',                   'SDpySubDynA.Cmass') #CJB+
            self.connect('SSSum',                   'SDpySubDynA.SSSum') #CJB+
            self.connect('OutCOSM',                 'SDpySubDynA.OutCOSM') #CJB+
            self.connect('OutAll',                  'SDpySubDynA.OutAll') #CJB+
            self.connect('OutSwtch',                'SDpySubDynA.OutSwtch') #CJB+
            self.connect('TabDelim',                'SDpySubDynA.TabDelim') #CJB+
            self.connect('OutDec',                  'SDpySubDynA.OutDec') #CJB+
            self.connect('OutFmt',                  'SDpySubDynA.OutFmt') #CJB+
            self.connect('OutSFmt',                 'SDpySubDynA.OutSFmt') #CJB+
            self.connect('MemOutListHeader',        'SDpySubDynA.MemOutListHeader') #CJB+
            self.connect('MemOutList',              'SDpySubDynA.MemOutList') #CJB+
            self.connect('SSOutlist',               'SDpySubDynA.SSOutlist') #CJB+
            self.connect('InputandDriverpath',      'SDpySubDynA.InputandDriverpath') #CJB+
            self.connect('EchoD',                   'SDpySubDynA.EchoD') #CJB+
            self.connect('Gravity',                 'SDpySubDynA.Gravity') #CJB+
            self.connect('WtrDpth',                 'SDpySubDynA.WtrDpth') #CJB+
            self.connect('SDInputFile',             'SDpySubDynA.SDInputFile') #CJB+
            self.connect('OutRootName',             'SDpySubDynA.OutRootName') #CJB+
            self.connect('NSteps',                  'SDpySubDynA.NSteps') #CJB+
            self.connect('TimeInterval',            'SDpySubDynA.TimeInterval') #CJB+
            self.connect('TP_RefPoint',             'SDpySubDynA.TP_RefPoint') #CJB+
            self.connect('SubRotateZ',              'SDpySubDynA.SubRotateZ') #CJB+
            self.connect('InputsMod',               'SDpySubDynA.InputsMod') #CJB+
            self.connect('InputsFile',              'SDpySubDynA.InputsFile') #CJB+
            self.connect('uTPInSteady',             'SDpySubDynA.uTPInSteady') #CJB+
            self.connect('uDotTPInSteady',          'SDpySubDynA.uDotTPInSteady') #CJB+
            self.connect('uDotDotTPInSteady',       'SDpySubDynA.uDotDotTPInSteady') #CJB+
            self.connect('SDpath',                  'SDpySubDynA.SDpath') #CJB+
            self.connect('Readpath_out',            'SDpySubDynA.Readpath_out') #CJB+
            self.connect('Readpath_sum',            'SDpySubDynA.Readpath_sum') #CJB+
            self.connect('Delete_file',             'SDpySubDynA.Delete_file') #CJB+
            self.connect('InputFile_path',          'SDpySubDynA.InputFile_path') #CJB+
            self.connect('Driver_path',             'SDpySubDynA.Driver_path') #CJB+
    
            #SubDyn Outputs CJB+
            self.connect('SDpySubDynA.SubDynAOuts', 'SDReturns.SubDynOuts') #CJB+
    
            #SubDyn Inputs from Jacket
            #gravity, nodes (see below)
            self.connect('PreBuild.WaterNew.wdepth', 'SDpySubDynA.wdepth') #CJB+ #CJB+ water flag
            self.connect('JcktGeoIn.dck_botz', 'SDpySubDynA.dck_botz') #CJB+
            self.connect('Build.JcktGeoOut.mems','SDpySubDynA.mems') #CJB+
            self.connect('Build.JcktGeoOut.Reacts','SDpySubDynA.Reacts') #CJB+
            self.connect('Build.SDPropSet', 'SDpySubDynA.SDPropSet') #CJB+
            self.connect('RNAinputs',          'SDpySubDynA.RNAinputs' ) #CJB+ test

        # LoadFrameOuts and LoadFrameOuts2
        self.connect('RNA_F',                      'LoadFrameOuts.RNA_F')
        self.connect('RNAinputs',                  'LoadFrameOuts.RNAinputs' )
        self.connect('PreBuild.WaterNew',                'LoadFrameOuts.Waterinputs') #CJB+ water flag
        self.connect('Windinputs',                 'LoadFrameOuts.Windinputs')

        self.connect('JcktGeoIn.VPFlag',           ['LoadFrameOuts.VPFlag'      ,    'LoadFrameOuts2.VPFlag'])
        self.connect('Build.pillegDs',             ['LoadFrameOuts.pillegDs'    ,    'LoadFrameOuts2.pillegDs'])
        self.connect('Build.twrDs',                ['LoadFrameOuts.twrDs'       ,    'LoadFrameOuts2.twrDs'])
        self.connect('Build.Pilemems',             ['LoadFrameOuts.Pilemems'    ,    'LoadFrameOuts2.Pilemems'])
        self.connect('Build.Legmems',              ['LoadFrameOuts.Legmems'     ,    'LoadFrameOuts2.Legmems'])
        self.connect('Build.Twrmems',              ['LoadFrameOuts.Twrmems'     ,    'LoadFrameOuts2.Twrmems'])
        self.connect('Build.TPmems',               ['LoadFrameOuts.TPmems'      ,    'LoadFrameOuts2.TPmems'])
        self.connect('PreBuild.al_bat3D',          ['LoadFrameOuts.al_bat3D'    ,    'LoadFrameOuts2.al_bat3D'])
        self.connect('JcktGeoIn.nlegs' ,           ['LoadFrameOuts.nlegs'       ,    'LoadFrameOuts2.nlegs'])
        self.connect('Build.JcktGeoOut.nodes',     ['LoadFrameOuts.nodes'       ,    'LoadFrameOuts2.nodes', 'SDpySubDynA.nodes']) #CJBe
        self.connect('Build.JcktGeoOut',           ['LoadFrameOuts.JcktGeoOut'  ,    'LoadFrameOuts2.JcktGeoOut'])
        self.connect('TwrRigidTop',                ['LoadFrameOuts.TwrRigidTop' ,    'LoadFrameOuts2.TwrRigidTop'])
        self.connect('Build.XjntIDs',              ['LoadFrameOuts.XjntIDs'     ,    'LoadFrameOuts2.XjntIDs'])
        self.connect('Build.KjntIDs',              ['LoadFrameOuts.KjntIDs'     ,    'LoadFrameOuts2.KjntIDs'])
        self.connect('Tower.Twrouts',              ['LoadFrameOuts.Twr_data'    ,    'LoadFrameOuts2.Twr_data'])
        self.connect('Tower.Dt',                   ['LoadFrameOuts.Dt'          ,    'LoadFrameOuts2.Dt'])
        self.connect('Tower.tt',                   ['LoadFrameOuts.tt'          ,    'LoadFrameOuts2.tt'])
        self.connect('Twrinputs.TwrSecH',          ['LoadFrameOuts.L_reinforced',    'LoadFrameOuts2.L_reinforced'])

        self.connect('IECpsfIns',                  ['LoadFrameOuts.IECpsfIns'   ,    'LoadFrameOuts2.IECpsfIns'])
        self.connect('FrameAuxIns',                ['LoadFrameOuts.FrameAuxIns' ,    'LoadFrameOuts2.FrameAuxIns'])

        self.connect('FrameAuxIns.gvector[2]','ABS1.varin')
        self.connect('ABS1.varout',                ['LoadFrameOuts.gravity'     ,    'LoadFrameOuts2.gravity', 'SDpySubDynA.gravity' ]) #CJBe test

        self.connect('Pileinputs.Lp',              ['LoadFrameOuts.Lp'          ,  'LoadFrameOuts2.Lp' ])
        self.connect('Piles.PileObjout',           ['LoadFrameOuts.PileObjout'  ,  'LoadFrameOuts2.PileObjout'])

        self.connect('Soil.SoilObj',               ['LoadFrameOuts.SoilObj'     ,'LoadFrameOuts2.SoilObj'])
        self.connect('JcktGeoIn.batter',           ['LoadFrameOuts.batter'      ,'LoadFrameOuts2.batter'])
        self.connect('PreBuild.innr_ang',          ['LoadFrameOuts.innr_ang'    ,'LoadFrameOuts2.innr_ang'])
        self.connect('PreLeg.prelegouts.legZbot',  ['LoadFrameOuts.loadZ'       ,'LoadFrameOuts2.loadZ'])

        # LoadFrameOuts2
        self.connect('RNA_F2',                     'LoadFrameOuts2.RNA_F')
        self.connect('RNAinputs2',                 'LoadFrameOuts2.RNAinputs' )
        self.connect('PreBuild.WaterNew',               'LoadFrameOuts2.Waterinputs') #CJB+ water flag I think this should work, but I'm not sure
        self.connect('Windinputs2',                'LoadFrameOuts2.Windinputs')

        #Total Mass
        self.connect('JcktGeoIn.nlegs',          'TotMass.nlegs')
        self.connect('Piles.Pileouts.PileObj',   'TotMass.PileObj')
        self.connect('Piles.PileObjout',         'TotMass.PileObjout')
        self.connect('Pileinputs.Lp',            'TotMass.Lp')

        # overall assembly input coonections and passthroughs
        #outs
        self.create_passthrough('Build.JcktGeoOut')
        self.create_passthrough('Tower.Twrouts.mass')
        self.create_passthrough('PreBuild.bay_hs')
        self.create_passthrough('PreBuild.dck_width')

        self.connect('LoadFrameOuts.Frameouts','FrameOut.Frameouts_ins')  #Needed to please the optimizer in OPENmdao

        if not self.clamped:
            self.create_passthrough('LoadFrameOuts.SPI_Kmat')


        self.create_passthrough('TotMass.ToTMass')

        self.create_passthrough('PreBuild.wbase')  #Width at the mudline send it as output
        self.create_passthrough('LoadFrameOuts.Mpiles') #Mass of all piles with given Lp
        self.create_passthrough('BrcCriteria.MudBrcCriteria') #MudBrace constraints
        self.create_passthrough('BrcCriteria.XBrcCriteria') #XBrace constraints
        self.create_passthrough('PreBuild.beta2D')#used by ANSYS

        self.connect('Tower.Twrouts','Twrouts')
        self.connect('Legs.legouts','Legouts')

#_____________________________________________________#

class ComponentMass(Component):  #TO DO
    """This component Calculates masses of the model + embedded piles"""
    #inputs
    nlegs =     Int(units=None,     iotype='in',desc='Number of Legs [3 or 4]??')
    PileObj=    Instance(Klass=Tube,iotype='in',desc='Object of Class Tube for Non-AF portion of Pile')
    PileObjout = Instance(Klass=Tube,     iotype='in', desc='Object of Class Tube for Non-AF portion of Pile, valid also when ndiv=0')
    Lp         =Float(  units='m',  iotype='in', desc='Input Embedment Length')  #User defined or optimizer assigned
   # outputs: internal mass of the jacket model

    ToTMass = Float(iotype='out', units='kg', desc='Jacket Dry Mass from Model')
    def execute(self):

       """This function returns:
          1.mass, dry mass of the Jacket members only (no TP, no Tower, no piles, no non-structural mass)
       """
       #piles' mass
       pileM=self.nlegs*(self.PileObjout.mat[0].rho*self.Lp*self.PileObjout.Area) #Embedded Portion
       if self.PileObj.L.size:  #This means we have pile stumps
            pileM +=self.nlegs*self.PileObj.mat[0].rho*self.PileObj.Area[0]*self.PileObj.L[0] #Note constant material and cross section assumed, also L is the buckling length, but in this case it is also the assumed pile stump length

       self.ToTMass=pileM
       #calculate mass
       #rhos=np.zeros(self.n_mems)  #initialize an array to contain mat density
       #for ii,jj in enumerate(self.mbr_strct.mat[0:self.n_mems]):
       #    rhos[ii]=jj.rho

       #self.mass=np.sum(rhos*self.mbr_strct.Area[0:self.n_mems]*self.mbr_strct.L[0:self.n_mems])
#_____________________________________________________#

#_____________#
#Other Functions
#_____________#



#______________________________________________________________________________#
def FindBrcAng(Hbays,nbays,wbas0,al_bat2D):   #Does it need to be a component???
    """This function calculates the brace angle and the bay heights for a Jacket\n
        with constant brace angle solving a non-linear equation.\n
        Input is an object of class Jacket.\n
        Output is a tuple containing [bay-base widths=array(nbays), \n
        bay-heights=array(nbays), brace angle with leg (2D) =beta2D] """

    def bay_b(beta2D): #bay base-widths  This function could be made recursive
        """This function calculates the bay base widths given the brace angle and the 2D batter"""
        bay_bs=np.zeros(nbays)
        bay_bs[0]=wbas0
        for ii in range(1,nbays):
            ##bay_bs[ii]=bay_bs[ii-1]*\
            ##          (1-2*tan(pi/2-beta2D-al_bat2D)*np.sin(al_bat2D)* \
            ##           np.sin(beta2D+al_bat2D)/np.sin(2*al_bat2D+beta2D))
            bay_bs[ii]=bay_bs[ii-1]*(1.-2.*tan(al_bat2D)*tan(pi/2-beta2D-al_bat2D)/\
                                        (1.+tan(al_bat2D)*tan(pi/2-beta2D-al_bat2D)))
        return bay_bs  #array [nbays]

    def bay_h(beta2D):   #bay heights known if beta2D is known
        """This function calculates the bay heights given the brace angle and the 2D batter and the bay widths"""
        ##return bay_b(beta2D)*tan(pi/2-beta2D-al_bat2D)*\
        ##      (1-np.cos(beta2D+al_bat2D)*np.sin(al_bat2D)/np.sin(2*al_bat2D+beta2D))
        return bay_b(beta2D)*tan(pi/2.-beta2D-al_bat2D)/\
                (1.+tan(al_bat2D)*tan(pi/2.-beta2D-al_bat2D))

    #solve the main equation that gives me both beta2D and bay_hs'
    def htot2solve(beta2D):
        """Main equation that needs to be solved htot2solve=0"""
        return Hbays-bay_h(beta2D).sum()

    guess=np.deg2rad(40.)
    beta_2D=fsolve(htot2solve,guess)
    #Recalculate the final values of bay base widths and heights
    bay_h=bay_h(beta_2D)
    bay_b=bay_b(beta_2D)
    return bay_b, bay_h, beta_2D[0]
#______________________________________________________________________________#

def FindBeta3D(beta_2D, al_bat2D, al_bat3D, innr_ang, wbas0):
    """This function calculates the 3D beta angle based on the bottom bay. This \n
    is the angle you would compare against NORSOK 30 degrees.
    It also returns the angle between X-braces (facing the botom bay width)
    INPUTS:
        Jacket- object of class JacketC"""
    #1st bay height
    ##h1=wbas0*tan(pi/2-beta_2D-al_bat2D)*(1.-np.cos(beta_2D+al_bat2D)*np.sin(al_bat2D)/np.sin(2.*al_bat2D+beta_2D))
    h1=wbas0*tan(pi/2-beta_2D-al_bat2D)/(1.+tan(al_bat2D)*tan(pi/2-beta_2D-al_bat2D))
    #2nd joint width
    ##BC=wbas0*(1.-2.*tan(pi/2-beta_2D-al_bat2D)*np.sin(al_bat2D)*np.sin(beta_2D+al_bat2D)/np.sin(2.*al_bat2D+beta_2D))
    BC=wbas0-2.*h1*tan(al_bat2D)
    #First bay leg segment
    ##AB=h1/np.cos(al_bat2D)*np.sqrt(1.+np.sin(al_bat2D)**2)
    AB=h1/np.cos(al_bat3D)
    #Brace Length
    AC=np.sqrt( h1**2 + (wbas0-h1*tan(al_bat2D))**2 +(h1*tan(al_bat2D)*tan(innr_ang))**2 )

    def BCeqn2solve(beta_3D):
        """Function containing Extended Pythagora's: BC^2=AB^2+AC^2-2*AB*AC*np.cos(beta_3D)"""
        return BC**2-(AB**2+AC**2-2*AB*AC*np.cos(beta_3D))

    guess=beta_2D
    beta_3D=fsolve(BCeqn2solve,guess)

    al_Xbrc=np.pi-2.*np.arccos((AC**2+wbas0**2-AB**2)/(2.*AC*wbas0))
    return beta_3D[0], al_Xbrc

#______________________________________________________________________________#
def Melm(X1,X2,roll):
    """This function finds the transformation matrix M to go from global X,Y,Z to \n
    local element x,y,z.  Note x along member per Frame3DD. \n
    INPUT \n
    X1    -float(3,n), X,Y,Z values for n elements, i-th (1st)node [m] \n
    X2    -float(3,n), X,Y,Z values for n elements, j-th (2nd)node [m] \n
    roll  -float(n),  roll angle as defined by Frame3DD about local x axis of element  [rad] \n
    OUTPUT    \n
    out   -float(3,3,n), direction cosine matrices of coordinate system having
             x-axis along element, y-axis along element principal axis, z-axis along other principal axis, for starting and ending x-sections \n"""
    from numpy import cos,sin,arctan2,arctan
    out=np.zeros([3,3,np.shape(X1)[1]])
    DX=X2-X1 #(3,n)
    L=CalcDist(X1,X2)
    psi=np.arctan2(DX[1,:],DX[0,:]) #psi angle about global Z
    tht=-np.arctan2(DX[2,:],np.sqrt(DX[0,:]**2+DX[1,:]**2)) #tht angle about global Y
    phi=roll  #phi angle about local x
    out[:,:,:]=np.array([[cos(psi)*cos(tht),  cos(tht)*sin(psi),  -sin(tht)],\
       [cos(psi)*sin(phi)*sin(tht)-sin(psi)*cos(phi),  cos(phi)*cos(psi)  + sin(phi)*sin(psi)*sin(tht),  cos(tht)*sin(phi)],\
       [cos(psi)*cos(phi)*sin(tht)+sin(psi)*sin(phi),  -sin(phi)*cos(psi) + cos(phi)*sin(psi)*sin(tht),  cos(tht)*cos(phi)]])
    return out
#______________________________________________________________________________#


def BrcFind(legObj,bay_bs,bay_hs,innr_ang,al_bat2D,wdepth,mudbrc=False,VPFlag=False,Kbuck=0.8,Klr_mx=70,minmass=True,PipeSch=True):
    #JcktInfo=Jacket()
    """This Function calculates a brace size given a leg size and follwong a few empirical constraints as follows:\n
    !!!It needs to be augmented to include horizontal bracing !!!
    INPUT
        legObj      -Tube object as coming from component legs
        innr_ang    -float,  Angle between radial direction and base side [rad] (pi/6 * (nlegs==3) + pi/4 *(nlegs==4))
        al_bat2D    -float, angle between leg and vertical in 2D projection [rad]
        wdepth      -float, Water Depth [m]
        bay_bs      -float(nbays), Widths of bays [m]
        bay_hs      -float(nbays), Heights of bays [m]
        mudbrc   -boolean, flag to see whether we are calculating horizontal (mud-) brace
    OPTIONAL INPUT
        VPFlag -boolean, Vertical Pile Flag (used for mudbrace sizing)
        Kbuck  -float, buckling constant for Klr
        Klr_mx -float, maximum valu for Klr
    OUTPUT
        Dbrc -float, brace diameter [m]
        tbrc -float, brace thickness [m]
    NOTE:
        1.Criteria to Find Braces:
                1.1. A>0.1*Aleg (stiffness criteria)
                1.2. Dbrc/Dleg>0.3  (joint capacity)
                1.3. Db/tb>=31  (between 19 and 60, >=31 floats)
                1.4. Db/tb> [250/ (wdth in ft)^(1/3)] (to avoid hydrostatic issues)
                1.5. Klr<=70  (70<Klr<90 is a good interval to avoid buckling)
        2.For bottom horizontal brace use Klr<=90,Kbuck=0.9
        3.Assumed: minimum Dbrc=0.25m minimum tbrc=0.25 in
   """
    import scipy.optimize

    nbays=bay_bs.size
    #here takce care of the braces, linking them to the Dleg via constraints


    #we will need to find the D,t for the brace given everything else
    #Find function that gets brace diameter and thickness linked to Dleg
    #need to calculate some extra parameters that would be in buildjck() but we do not want to build yet as we are finding braces

    Dleg=legObj.D[-nbays:] #ODs at various levels
    tleg=legObj.t[-nbays:] #ts at various levels

    Ab=0.1*legObj.Area[-nbays:]  #Constraint Ab

    Lbrc2=bay_bs[0]-Dleg[1-VPFlag]  #Mud brace effective length
    nbrcs=1  #Initialize as if we were calculating only 1 brace (mudbrace or TopHbrace)

    if not(mudbrc): #If it is a X-brace we are after
       #1st bay height
        h1=bay_hs[0]
        #2nd bay bay width (top of bay 1 width)
        wbas1=bay_bs[1]
        #First bay leg segment, it works with 3 and 4 legs
        AB=h1*np.sqrt(1.+np.tan(al_bat2D)**2 *(1.+np.tan(innr_ang)**2))
        #Brace Length of 1st bay, and then the longest portion of it due to intersection: AO/AC=wbas0/(wbas0+wbas1)
        Lbrc2=np.sqrt( h1**2 + (bay_bs[0]-h1*np.tan(al_bat2D))**2 +(h1*np.tan(al_bat2D)*np.tan(innr_ang))**2 ) * bay_bs[0]/(bay_bs[0]+wbas1)
        nbrcs=nbays
    else:
        al_bat2D=0.  #Initialize for the print out later on

    #Initialize output
    Dbrcs=np.empty([nbrcs])
    tbrcs=np.empty([nbrcs])

    hydrocnst=250./(wdepth*(1e3/25.4/12))**(1./3.)
    #brace function of D,t; tried to make it work with x an array of nbays*2 length, but slsqp bombs out, so I need to use a loop
    braces= lambda x: Tube(x[0],x[1],Lbrc2,Kbuck)
    #constraints x[0]=Dbrc, x[1]=tbrc, x[2]=lbrc
    cons=[{'type':'ineq',\
          'fun': lambda x: braces(x).Area-Ab},\
          {'type':'ineq',\
          'fun': lambda x: braces(x).D/Dleg-0.3},\
          {'type':'ineq',\
          'fun': lambda x: braces(x).D/braces(x).t-31.},\
          {'type':'ineq',\
          'fun': lambda x: -braces(x).D/braces(x).t+hydrocnst},\
         ]
    ieqcons=[cons[0].get('fun'),cons[1].get('fun'),cons[2].get('fun'),cons[3].get('fun')]

    #Next is the main function to minimize (in our case all constraints are inequalities)
    if minmass:
        func=lambda x: braces(x).Area  #we are trying to minimize the mass=rho*Area*L, but rho and L are fixed
        cons.append({'type':'ineq','fun': lambda x: Klr_mx-braces(x).Klr})
        ieqcons.append(cons[4].get('fun'))
    else:
        func=lambda x: abs(braces(x).Klr-Klr_mx)  #we are trying to enforce a Klr=70 in case we do not use mass minimization

    for ib in range(0,nbrcs): #repeat the optimizer for each bay
        #initial guess
        guess=np.hstack((Dleg[ib]/4.,tleg[ib]*.75))
        tmin=0.0254*0.25    #minimum brace t  (used to be tiled times nbays, but slsqp was not working)
        tmax=tleg[ib]       #max brace t
        Dmin=0.25           #minimum brace D (used to be tiled times nbays, but slsqp was not working)
        Dmax=Dleg[ib]       #max brace D
        bounds=[(Dmin,Dmax),(tmin,tmax)]
    #solve the minimization
    #res=scipy.optimize.fmin_slsqp(func,guess,eqcons=[cons[1].get('fun')],ieqcons=[cons[0].get('fun')],iprint=2, full_output=1)
        res=list(range(5))
        cc=0

        while res[3]!=0 and cc<=10:
        #HERE IS THE OPTIMIZER
            res=scipy.optimize.fmin_slsqp(func,guess,bounds=bounds,ieqcons=ieqcons,iprint=0, full_output=1) #iprint=1 or 2 to have more debugging info
            ##res=scipy.optimize.fmin_slsqp(func,guess,bounds=map(tuple,bounds),\
            ##ieqcons=ieqcons,iprint=3, full_output=True) #iprint=1 or 2 to have more debugging info
            if res[3] == 0:
                Dbrcs[ib]=res[0][0]
                tbrcs[ib]=res[0][1]
            else:
                cc+=1#increase counter for attempts at making it converge
                #change guess and see if it works
                guess[0]=Dleg[ib]/2.+(Dleg[ib]-Dleg[ib]/2.)/9.*(cc-1) #i.e. start with a D at Dleg/2 and go all the way up
                guess[1]=guess[0]/31.
        if res[3] != 0:
            print 'Minimization for Brace Calculation Unsuccessful for batter={:3.2f}, Dleg={:3.2f}, tleg={:6.3f}, nbays={:d}'.format(1./np.tan(al_bat2D),Dleg[ib],tleg[ib],nbays)

        if PipeSch: #ask for scheduling
            Dbrcs[ib]=SchPipe(Dbrcs[ib])
            tbrcs[ib]=np.max(Dbrcs[ib]/hydrocnst,tbrcs[ib])

    return Dbrcs, tbrcs, res[3]

#______________________________________________________________________________#

#______________________________________________________________________________#
def SchPipe(D):
    """This function returns a scheduled pipe Diameter given an initial OD. \n
    It would be correct for pipes up to 18in in OD, but we use it to be conservative also beyond that.
    INPUT
        D  - float(n), OD of pipes in [m]
    OUTPUT
        SchPipeOD   -float(n), scheduled OD  in [m]         """
    Din=D/0.0254  #in inches
    SchPipe=np.ceil(Din) #inches
    SchPipe+=np.mod(SchPipe,2)  #inches
    SchPipe*=0.0254 #in meters
    return SchPipe
#______________________________________________________________________________#


if __name__ == '__main__':
    #PyObject *f = PySys_GetObject("stdout")
    #PyFile_WriteString

    """This 'main' function runs an example of jacket. It can also carry out an optimization, \n
       all you need to do is reviewing input parameters, design variables and constraints.\n
       Main keywords are set at the top (first 2 lines below):\n
        optimize  -boolean, if True it will perform an optimization with a method to be set by OPTswitch
        OPTswitch -string, ['Cobyla'/ 'PyOptSNOPT'/ 'PyOPTCobyla'] for python Cobyla, or PyOPT SNOPT, or PyOPT Cobyla."""

    optimize = False        #Set this one to True if you want a test on optimization
    #OPTswitch= 'Cobyla'     #'Cobyla', 'PyOptSNOPT', 'PyOPTCobyla'

    if len(sys.argv)>1 and len(sys.argv)<5: #This means individual case
        optimize=sys.argv[1]
        OPTswitch=sys.argv[2].lower() #== 'true'


    #--- Example Input Parameters ---#

    #Note: in the following input lines, certain inputs can be sawpped for those following lines " ## if turbine_jacket "
    #       to test the turbineSE inputs as passed to jacket.

    #--- Set Jacket Input Parameters ---#
    Jcktins=JcktGeoInputs()
    Jcktins.nlegs =4
    Jcktins.nbays =5 #CJBe Makes the program run faster for debugging
    Jcktins.batter=12.
    Jcktins.dck_botz =16.
    Jcktins.dck_width=2*6.
    Jcktins.weld2D   =0.5
    Jcktins.VPFlag = True    #vertical pile T/F;  to enable piles in frame3DD set pileinputs.ndiv>0
    Jcktins.clamped= False    #whether or not the bottom of the structure is rigidly connected. Use False when equivalent spring constants are being used.
    Jcktins.AFflag = False  #whether or not to use apparent fixity piles
    Jcktins.PreBuildTPLvl = 5  #if >0, the TP is prebuilt according to rules per PreBuildTP

    #Soil inputs
    Soilinputs=SoilGeoInputs()
    Soilinputs.zbots   =-np.array([3.,5.,7.,15.,30.,50.])
    Soilinputs.gammas  =np.array([10000.,10000.,10000.,10000.,10000.,10000.])
    Soilinputs.cus     =np.array([60000.,60000.,60000.,60000.,60000.,60000.])
    Soilinputs.phis    =np.array([26.,26.,26.,26.,26.,26])#np.array([36.,33.,26.,37.,35.,37.5])#np.array([36.,33.,26.,37.,35.,37.5])
    Soilinputs.delta   =25.
    Soilinputs.sndflg   =True
    Soilinputs.PenderSwtch   =False #True
    Soilinputs.SoilSF   =1.

    Soilinputs2=copy.copy(Soilinputs) #Parked case. We assume same stiffness although this may not be the case under a different load

    #Water and wind inputs
    Waterinputs=WaterInputs()
    Waterinputs.wdepth   =30.
    #Waterinputs.wlevel   =30. #CJB- #Distance from bottom of structure to surface  THIS, I believe is no longer needed as piles may be negative in z, to check and remove in case
        #CJB The above line sets z's based on bottom of pile; initially z=0 is at the bottom of the piles. wlevel is used (exlusively in jacket.py) to shift the z-coordinates
        #so that z=0 is at the MSL. This is used instead of wdepth as the piles can go deeper than the water.


    Waterinputs.T=12.  #Wave Period
    Waterinputs.HW=10. #Wave Height
    Waterinputs.Cd=3.  #Drag Coefficient, enhanced to account for marine growth and other members not calculated
    Waterinputs.Cm=8.#2.  #Added mass Coefficient

    Waterinputs2=copy.copy(Waterinputs)  #PARKED CONDITIONS - still max wave hereyou

    Waterinputs.T=8.  #Wave Period
    Waterinputs.HW=4. #Wave Height

    Windinputs=WindInputs()
    Windinputs.Cdj=4.  #Drag Coefficient for jacket members, enhanced to account for TP drag not calculated otherwise
    Windinputs.Cdt=2  #Drag Coefficient for tower, enhanced to account for TP drag not calculated otherwise
    Windinputs.HH=100. #CHECK HOW THIS COMPLIES....
    Windinputs.U50HH=30. #assumed gust speed

    ## if turbine_jacket
    ##Windinputs.HH=90. #CHECK HOW THIS COMPLIES....
    ##Windinputs.U50HH=11.7373200354 # using rated loads
    ##Windinputs.rho = 1.225
    ##Windinputs.mu = 1.81206e-05

    Windinputs2=copy.copy(Windinputs)
    Windinputs2.U50HH=70. #assumed gust speed

    #Pile data
    Pilematin=MatInputs()
    Pilematin.matname=np.array(['steel'])
    Pilematin.E=np.array([ 25.e9])
    Dpile=2.5#0.75 # 2.0
    tpile=0.01
    Lp=20. #45

    Pileinputs=PileGeoInputs()
    Pileinputs.Pilematins=Pilematin
    Pileinputs.ndiv=0 #3
    Pileinputs.Dpile=Dpile
    Pileinputs.tpile=tpile
    Pileinputs.Lp=Lp #[m] Embedment length

    #Legs data
    legmatin=MatInputs()
    legmatin.matname=(['heavysteel','heavysteel','heavysteel','heavysteel'])
    #legmatin.E=np.array([2.0e11])
    Dleg=np.array([1.5,1.5,1.5,1.5,1.5,1.5])
    tleg=1.5*np.array([0.0254]).repeat(Dleg.size)
    leginputs=LegGeoInputs()
    leginputs.legZbot   = 1.0
    leginputs.ndiv=1
    leginputs.legmatins=legmatin
    leginputs.Dleg0=Dleg[0]
    leginputs.tleg0=tleg[0]

    legbot_stmphin =1.5  #Distance from bottom of leg to second joint along z; must be>0

    #Xbrc data
    Xbrcmatin=MatInputs()
    Xbrcmatin.matname=np.array(['heavysteel']).repeat(Jcktins.nbays)
    #Xbrcmatin.E=np.array([ 2.2e11, 2.0e11,2.0e11,2.0e11,2.0e11])
    Dbrc=np.array([1.,1.,1.0,1.0,1.0])
    tbrc=np.array([1.,1.,1.0,1.0,1.0])*0.0254

    Xbrcinputs=XBrcGeoInputs()
    Xbrcinputs.Dbrc0=Dbrc[0]
    Xbrcinputs.tbrc0=tbrc[0]
    Xbrcinputs.ndiv=2#2
    Xbrcinputs.Xbrcmatins=Xbrcmatin
    Xbrcinputs.precalc=False #True   #This can be set to true if we want Xbraces to be precalculated in D and t, in which case the above set Dbrc and tbrc would be overwritten

    #Mbrc data
    Mbrcmatin=MatInputs()
    Mbrcmatin.matname=np.array(['heavysteel'])
    #Mbrcmatin.E=np.array([ 2.5e11])
    Dbrc_mud=1.5

    Mbrcinputs=MudBrcGeoInputs()
    Mbrcinputs.Dbrc_mud=Dbrc_mud
    Mbrcinputs.ndiv=2
    Mbrcinputs.Mbrcmatins=Mbrcmatin
    Mbrcinputs.precalc=False #True   #This can be set to true if we want Mudbrace to be precalculated in D and t, in which case the above set Dbrc_mud and tbrc_mud would be overwritten

    #Hbrc data
    Hbrcmatin=MatInputs()
    Hbrcmatin.matname=np.array(['heavysteel'])
    Hbrcmatin.E=np.array([ 2.5e11])
    Dbrc_hbrc=1.1

    Hbrcinputs=HBrcGeoInputs()
    Hbrcinputs.Dbrch=Dbrc_hbrc
    Hbrcinputs.ndiv=0#2
    Hbrcinputs.Hbrcmatins=Hbrcmatin
    Hbrcinputs.precalc=True   #This can be set to true if we want Hbrace to be set=Xbrace top D and t, in which case the above set Dbrch and tbrch would be overwritten

    #TP data
    TPlumpinputs=TPlumpMass()
    TPlumpinputs.mass=200.e3 #[kg]

    TPstmpsmatin=MatInputs()
    TPbrcmatin=MatInputs()
    TPstemmatin=MatInputs()
    TPbrcmatin.matname=np.array(['heavysteel'])
    #TPbrcmatin.E=np.array([ 2.5e11])
    TPstemmatin.matname=np.array(['heavysteel']).repeat(2)
    #TPstemmatin.E=np.array([ 2.1e11]).repeat(2)

    TPinputs=TPGeoInputs()
    TPinputs.TPbrcmatins=TPbrcmatin
    TPinputs.TPstemmatins=TPstemmatin
    TPinputs.TPstmpmatins=TPstmpsmatin
    TPinputs.Dstrut=leginputs.Dleg[-1]
    TPinputs.tstrut=leginputs.tleg[-1]
    TPinputs.Dgir=Dbrc_hbrc
    TPinputs.tgir=0.0254
    TPinputs.Dbrc=1.1
    TPinputs.Dbrc=TPinputs.Dgir
    TPinputs.tbrc=TPinputs.tgir

    TPinputs.hstump=1.0#1.0
    TPinputs.Dstump=1.25#1.0
    TPinputs.stumpndiv=1#2
    TPinputs.brcndiv=1#2
    TPinputs.girndiv=1#2
    TPinputs.strutndiv=1#2
    TPinputs.stemndiv=1#2
    TPinputs.nstems=3
    TPinputs.Dstem=np.array([6.]).repeat(TPinputs.nstems)
    TPinputs.tstem=np.array([0.1,0.11,0.11])
    TPinputs.hstem=np.array([6./TPinputs.nstems]).repeat(TPinputs.nstems)

    #Tower data
    Twrmatin=MatInputs()
    Twrmatin.matname=np.array(['heavysteel'])
    #Twrmatin.E=np.array([ 2.77e11])
    Twrinputs=TwrGeoInputs()
    Twrinputs.Twrmatins=Twrmatin
    #Twrinputs.Htwr=70.  #Trumped by HH
    Twrinputs.Htwr2frac=0.2   #fraction of tower height with constant x-section
    Twrinputs.ndiv=np.array([6,12])  #ndiv for uniform and tapered section
    Twrinputs.DeltaZmax= 6. #[m], maximum FE element length allowed in the tower members (i.e. the uniform and the tapered members)
    Twrinputs.Db=5.6
    Twrinputs.DTRb=130.
    Twrinputs.DTRt=150.
    Twrinputs.Dt=0.55*Twrinputs.Db
    ## if turbine_jacket
    ##Twrinputs.Dt = 3.87

    TwrRigidTop= False #False       #False=Account for RNA via math rather than a physical rigidmember

    #RNA data
    RNAins=RNAprops()
    RNAins.mass=3*350.e3
    RNAins.I[0]=86.579E+6
    RNAins.I[1]=53.530E+6
    RNAins.I[2]=58.112E+6
    RNAins.CMoff[2]=2.34
    RNAins.yawangle=45.  #angle with respect to global X, CCW looking from above, wind from left
    RNAins.rna_weightM=True
    ## if turbine_jacket
    ##RNAins.mass=285598.806453
    ##RNAins.I = np.array([1.14930678e8, 2.20354030e7, 1.87597425e7, 0.0, 5.03710467e5, 0.0])
    ##RNAins.CMoff = np.array([-1.13197635, 0.0, 0.50875268])
    ##RNAins.yawangle=0.0  #angle with respect to global X, CCW looking from above, wind from left
    #RNAins.rna_weightM=True

    RNAins2=copy.copy(RNAins)  #PARKED CASE, for now assume the same

    #RNA loads              Fx-z,         Mxx-zz
    RNA_F=np.array([1000.e3,0.,0.,0.,0.,0.])    #operational
    RNA_F2=np.array([500.e3,0.,0.,0.,0.,0.])    #Parked
    ## if turbine_jacket
    ##RNA_F=np.array([1284744.19620519,0.,-2914124.84400512,3963732.76208099,-2275104.79420872,-346781.68192839])

    #Frame3DD parameters
    FrameAuxIns=Frame3DDaux()
    FrameAuxIns.sh_fg=1               #shear flag-->Timoshenko
    FrameAuxIns.deltaz=5.
    FrameAuxIns.geo_fg=0
    FrameAuxIns.nModes = 6             # number of desired dynamic modes of vibration
    FrameAuxIns.Mmethod = 1            # 1: subspace Jacobi     2: Stodola
    FrameAuxIns.lump = 0               # 0: consistent mass ... 1: lumped mass matrix
    FrameAuxIns.tol = 1e-9             # mode shape tolerance
    FrameAuxIns.shift = 0.0            # shift value ... for unrestrained structures
    FrameAuxIns.gvector=np.array([0.,0.,-9.8065])    #GRAVITY
    ## if turbine_jacket
    ##FrameAuxIns.gvector=np.array([0.,0.,-9.81])    #GRAVITY

    #Decide whether or not to consider DLC 6.1 as well
    twodlcs=False

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    myjckt=set_as_top(JacketSE(Jcktins.clamped,Jcktins.AFflag,twodlcs=twodlcs)) ##(Jcktins.PreBuildTPLvl>0),

    #pySubDyn Parameters CJB+
    #SDpySubDynA = pySubDynA()

    #INPUTS TO RUN SUBDYN-------------------------------------------------------
    #(PATH INFORMATION FOR INPUT FILE AND DRIVER)

    Base_name="Jacket_pySubDyn" #CJB+ #Input name here

    #INPUT FILE PATH  #CJB+
    myjckt.InputFile_name=str(Base_name)+".txt"  #CJB+
    #myjckt.SDpySubDynA.InputFile_name=str(Base_name)+".txt"  #CJB+
    myjckt.InputandDriverpath=os.sep.join(['SubDyn', 'CertTest'])
    #myjckt.InputandDriverpath="C:\wisdem\plugins\JacketSE\src\jacketse\SubDyn\CertTest"  #CJB+
    myjckt.InputFile_path=myjckt.InputandDriverpath+os.sep+str(myjckt.InputFile_name)  #CJB+

    #DRIVER PATH  #CJB+
    myjckt.Driver_name=str(Base_name)+"D"+".txt"  #CJB+
    myjckt.Driver_path=myjckt.InputandDriverpath+os.sep+str(myjckt.Driver_name)  #CJB+
    myjckt.Driver_path=myjckt.InputandDriverpath+os.sep+str(myjckt.Driver_name)  #CJB+

    #PATH TO RUN SUBDYN  #CJB+
    SDEXEpath="C:\wisdem\plugins\JacketSE\src\jacketse\SubDyn"+os.sep+"bin\SubDyn_Win32.exe"  #CJB+
    myjckt.SDpath=str(SDEXEpath)+' '+str(myjckt.Driver_path)  #CJB+
    #test.SDpath='r'+"'''"+test.SDEXEpath+' '+test.SDDriverpath+"'''"  #CJB+

    #PATH TO READ OUTPUT (INPUTS TO READ OUTPUT)  #CJB+
    myjckt.Readpath_out=str(myjckt.InputandDriverpath)+os.sep+str(Base_name)+".SD.out"  #CJB+
    myjckt.Readpath_sum=str(myjckt.InputandDriverpath)+os.sep+str(Base_name)+".SD.sum"  #CJB+
    myjckt.Delete_file=False #Deletes driver, input, and output files. Does not delete Echo file.  #CJB+

    #INPUT FILE INPUTS----------------------------------------------------------  #CJB+

    #Simulation Control  #CJB+
    myjckt.Echo=np.array([False, "Echo", "- Echo input data to ""<rootname>.SD.ech"" (flag)"])  #CJB+
    myjckt.SDdeltaT=np.array(["DEFAULT", "SDdeltaT", "- Local Integration Step. If ""default"", the glue-code integration step will be used."])  #CJB+
    myjckt.IntMethod=np.array([4, "IntMethod", "- Integration Method [1/2/3/4 = RK4/AB4/ABM4/AM2]."])  #CJB+
    myjckt.SttcSolve=np.array([False, "SttcSolve", "- Solve dynamics about static equilibrium point"])  #CJB+

    #FEA and CRAIG-BAMPTON PARAMETERS  #CJB+
    myjckt.FEMmod=np.array([3, "- FEM switch: element model in the FEM. [1= Euler-Bernoulli(E-B);  2=Tapered E-B (unavailable);  3= 2-node Timoshenko;  4= 2-node tapered Timoshenko (unavailable)]"])  #CJB+
    myjckt.NDiv=np.array([1, "NDiv", "- Number of sub-elements per member"])#CJB "HARDWIRED" INTO PYSUBDYN AS 1 TO ALLOW JACKETSE'S NODES TO BE USED AS SUBDYN'S JOINTS  #CJB+
    myjckt.CBMod=np.array([False, "- [T/F] If True perform C-B reduction, else full FEM dofs will be retained. If True, select Nmodes to retain in C-B reduced system."])  #CJB+
    myjckt.Nmodes=np.array([75, "- Number of internal modes to retain (ignored if CBMod=False). If Nmodes=0 --> Guyan Reduction."])  #CJB+
    myjckt.JDampings=np.array([2, "JDampings", "- Damping Ratios for each retained mode (% of critical) If Nmodes>0, list Nmodes structural damping ratios for each retained mode (% of critical),\
     or a single damping ratio to be applied to all retained modes. (last entered value will be used for all remaining modes)."])  #CJB+

    #Structure Joints  #CJB+
    myjckt.SDjointsHeader=np.array([['JointID','JointXss','JointYss','JointZss','[Coordinates of Member joints in SS-Coordinate System]'], ['(-)','(m)','(m)','(m)']])  #CJB+

    #Base Reaction Joints  #CJB+
    myjckt.BaseRxnJointsHeader=np.array([['RJointID','RctTDXss','RctTDYss','RctTDZss','RctRDXss','RctRDYss','RctRDZss','[Global Coordinate System]'], ['(-)',('flag'),('flag'),('flag'),('flag'),('flag'),('flag')]])  #CJB+

    #Interface Joints  #CJB+
    myjckt.InterfaceRxnJointsHeader=np.array([['IJointID','ItfTDXss','ItfTDYss','ItfTDZss','ItfRDXss','ItfRDYss','ItfRDZss','[Global Coordinate System]'], ['(-)',('flag'),('flag'),('flag'),('flag'),('flag'),('flag')]])  #CJB+

    #Members  #CJB+
    myjckt.MembersHeader=np.array([['MemberID','MJointID1','MJointID2','MPropSetID1','MPropSetID2','COSMID'], ['(-)','(-)','(-)','(-)','(-)','(-)']])  #CJB+

    #MEMBER X-SECTION PROPERTY data 1/2  #CJB+
    myjckt.NPropSets=np.array([6, 'NPropSets', '- # of structurally unique x-sections (i.e. # of X-sectional props utilized throughout all of the members)'])  #CJB+
    myjckt.PropSet1Header=np.array([['PropSetID', 'YoungE', 'ShearG', 'MatDens', 'XsecD', 'XsecT'],['(-)', '(N/m2)', '(N/m2)', '(kg/m3)', '(m)', '(m)']])  #CJB+

    #MEMBER X-SECTION PROPERTY data 2/2  #CJB+
    myjckt.PropSet2=np.array([])  #CJB+
    myjckt.PropSet2Header=np.array([["PropSetID","YoungE","ShearG","MatDens","XsecA","XsecAsx","XsecAsy","XsecJxx","XsecJyy","XsecJ0"],["(-)","(N/m2)","(N/m2)","(kg/m3)","(m2)","(m2)","(m2)","(m4)","(m4)","(m4)"]])  #CJB+

    #MEMBER COSINE MATRICES COSM(i,j)  #CJB+
    myjckt.COSMHeader=np.array([["COSMID","COSM11","COSMID12","COSMID13","COSMID21","COSMID22","COSMID23","COSMID31","COSMID32","COSMID33"],["(-)","(-)","(-)","(-)","(-)","(-)","(-)","(-)","(-)","(-)"]])  #CJB+
    myjckt.COSMs=np.array([])  #CJB+

    #JOINT ADDITIONAL CONCENTRATED MASSES  #CJB+
    myjckt.CmassHeader=np.array([["CMJointID","JMass","JMXX","JMYY","JMZZ"],["(-)","(kg)" ,"(kg*m^2)","(kg*m^2)","(kg*m^2)"]])  #CJB+
    myjckt.Cmass=np.array([])  #CJB+

    #OUTPUT: SUMMARY & OUTFILE  #CJB+
    myjckt.SSSum=np.array([True, "SSSum", "- Output a Summary File (flag).It contains: matrices K,M  and C-B reduced M_BB, M-BM, K_BB, K_MM(OMG^2), PHI_R, PHI_L. It can also contain COSMs if requested."])  #CJB+
    myjckt.OutCOSM=np.array([True, "OutCOSM", "- Output cosine matrices with the selected output member forces (flag)"])  #CJB+
    myjckt.OutAll=np.array([True, "OutAll", "- [T/F] Output all members' end forces "])  #CJB+
    myjckt.OutSwtch=np.array([1, "OutSwtch", "- [1/2/3] Output requested channels to: 1=<rootname>.SD.out;  2=<rootname>.out (generated by FAST);  3=both files."])  #CJB+
    myjckt.TabDelim=np.array([True, "TabDelim", "- Generate a tab-delimited output in the <rootname>.SD.out file"])  #CJB+
    myjckt.OutDec=np.array([1, "OutDec", "- Decimation of output in the <rootname>.SD.out file"])  #CJB+
    myjckt.OutFmt=np.array(["Es11.4e2", "OutFmt", "- Output format for numerical results in the <rootname>.SD.out file"])  #CJB+
    myjckt.OutSFmt=np.array(["A11", "OutFmt", "- Output format for header strings in the <rootname>.SD.out file"])  #CJB+

    #MEMBER OUTPUT LIST  #CJB+
    myjckt.MemOutListHeader=np.array([['MemberID','NoutCnt','NodeCnt','[NOutCnt=how many nodes to get output for [< 10]; NodeCnt are local ordinal numbers from the start of the member, and must be >=1 and <= NDiv+1] If NMOutputs=0 leave blank as well.]'],\
                                           ['(-)','(-)','(-)']])  #CJB+

    #SSOutline  #CJB+
    myjckt.SSOutlist=np.array([["ReactFXss, ReactFYss, ReactFZss, ReactMXss, ReactMYss, ReactMZss",'-Base reactions (forces onto SS structure)'],\
                    ["IntfFXss,  IntfFYss,  IntfFZss,  IntfMXss, IntfMYss, IntfMZss",'-Interface reactions (forces from SS structure)'],\
                    ["IntfTDXss,  IntfTDYss,  IntfTDZss,  IntfRDXss, IntfRDYss, IntfRDZss",'-Interface deflections '],\
                    ["IntfTAXss,  IntfTAYss,  IntfTAZss,  IntfRAXss, IntfRAYss, IntfRAZss",'Interface accelerations']])


    myjckt.InterfaceJointsFlags=np.array([[10,1,1,1,1,1,1]]) #CJB+ The first entry is just a placeholder; the 6 ones need to be there.

    myjckt.MemOutList=np.array([[1,2,1,2]])  #CJB+

    #DRIVER INPUTS--------------------------------------------------------------  #CJB+

    myjckt.EchoD=np.array([True, "Echo", "- Echo the input file data (flag)"])  #CJB+

    #Environmental Conditions  #CJB+
    myjckt.Gravity=np.array([9.81, "Gravity", "- Gravity (m/s^2)"])  #CJB+
    myjckt.WtrDpth=np.array([43.127, "WtrDpth", "- Water Depth (m) positive value"])  #CJB+

    #SubDyn  #CJB+
    myjckt.SDInputFile=np.array([myjckt.InputFile_path, "SDInputFile"])  #CJB+
    myjckt.OutRootName=np.array([str(myjckt.InputandDriverpath)+os.sep+str(Base_name), "OutRootName"])  #CJB+
    myjckt.NSteps=np.array([600, "NSteps", "- Number of time steps in the simulations (-)"])  #CJB+
    myjckt.TimeInterval=np.array([0.005, "TimeInterval", "- TimeInterval for the simulation (sec)"])  #CJB+
    myjckt.TP_RefPoint=np.array([0.0, 0.0, 18.15, "TP_RefPoint", "- Location of the TP reference point in global coordinates (m)"])  #CJB+
    myjckt.SubRotateZ=np.array([0.0, "SubRotateZ", "- Rotation angle of the structure geometry in degrees about the global Z axis."])  #CJB+

    #INPUTS  #CJB+
    myjckt.InputsMod=np.array([1, "InputsMod", "- Inputs model {0: all inputs are zero for every timestep, 1: steadystate inputs, 2: read inputs from a file (InputsFile)} (switch)"])  #CJB+
    myjckt.InputsFile=np.array(['""', "InputsFile", "- Name of the inputs file if InputsMod = 2"])  #CJB+

    #STEADY INPUTS  #CJB+
    myjckt.uTPInSteady=np.array([0.1, 0.0, 0.0, 0.0, 0.0, 0.0, "uTPInSteady", "- input displacements and rotations ( m, rads )"])  #CJB+
    myjckt.uDotTPInSteady=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "uDotTPInSteady", "- input translational and rotational velocities ( m/s, rads/s)"])  #CJB+
    myjckt.uDotDotTPInSteady=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "uDotDotTPInSteady", "- input translational and rotational accelerations ( m/s^2, rads/s^2)"])  #CJB+

    #myjckt.SDpySubDynA.run()  #CJB+

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    #-----Launch the assembly-----#

    #myjckt=set_as_top(JacketSE(Jcktins.clamped,Jcktins.AFflag,twodlcs=twodlcs)) ##(Jcktins.PreBuildTPLvl>0),

    #Pass all inputs to assembly
    myjckt.JcktGeoIn=Jcktins

    myjckt.Soilinputs=Soilinputs
    myjckt.Soilinputs2=Soilinputs2   #Parked conditions

    myjckt.Waterinputs=Waterinputs
    myjckt.Windinputs=Windinputs
    myjckt.RNA_F=RNA_F

    #myjckt.SDpySubDynA=SDpySubDynA #CJB+

    myjckt.Waterinputs2=Waterinputs2 #Parked conditions
    myjckt.Windinputs2=Windinputs2   #Parked conditions
    myjckt.RNA_F2=RNA_F2            #Parked conditions

    myjckt.Pileinputs=Pileinputs
    myjckt.leginputs=leginputs
    myjckt.legbot_stmphin =legbot_stmphin
    myjckt.Xbrcinputs=Xbrcinputs
    myjckt.Mbrcinputs=Mbrcinputs
    myjckt.Hbrcinputs=Hbrcinputs
    myjckt.TPlumpinputs=TPlumpinputs
    myjckt.TPinputs=TPinputs
    myjckt.RNAinputs=RNAins
    myjckt.RNAinputs2=RNAins2
    myjckt.Twrinputs=Twrinputs
    myjckt.TwrRigidTop=TwrRigidTop
    myjckt.FrameAuxIns=FrameAuxIns

    #--- RUN OPTIMIZATION ---#
    if optimize:

        #Set Optimization Bounds for the various variables:
        #          x=  [ batter,  Dpile,    tpile,   Lp,   Dleg,     tleg,       Dbrc,   tbrc,     Dbrc_mud,   tbrc_mud,    Db,   DTRb   Dt,   DTRt   Htwr2fac ]
        MnCnst=np.array([ 8.,      1.,    1.*0.0254, 30.,   1.,       1.*0.0254,  1.,    1.*0.0254,   1.,     1.*0.0254,    5.,   120.,  3.,   120.,     0. ])
        MxCnst=np.array([ 15.,     2.5,   5.*0.0254, 50.,   2.5,      5.*0.0254,  2.,    5.*0.0254,   2.,     5.*0.0254,    7.,   200.,  4.,   200.,   0.25 ])

        # --- optimizer imports ---
        from pyopt_driver.pyopt_driver import pyOptDriver
        from openmdao.lib.casehandlers.api import DumpCaseRecorder
        # ----------------------

        # --- Setup Optimizer ---
        if OPTswitch == 'cobyla':
            myjckt.replace('driver', COBYLAdriver())

            myjckt.driver.rhobeg=0.01
            myjckt.driver.rhoend=1.e-3
            myjckt.driver.maxfun=2000
            myjckt.driver.iprint=1
        else:
            myjckt.replace('driver', pyOptDriver())
            if  OPTswitch== 'pyoptsnopt':
                myjckt.driver.optimizer = 'SNOPT'
                myjckt.driver.options = {'Major feasibility tolerance': 1e-3,
                                     'Minor feasibility tolerance': 1e-3,
                                     'Major optimality tolerance': 1e-3,
                                     'Function precision': 1e-3}
            elif OPTswitch== 'pyoptcobyla':
                myjckt.driver.optimizer = 'COBYLA'
                myjckt.driver.options = {'RHOEND':1.e-2,'RHOEND':1.e-3,'MAXFUN':2000,'IPRINT':1}
            else:
                sys.exit('Error: OPTSwitch must be set to ''Cobyla'', ''pyoptsnopt'' or ''pyoptcobyla''!!!')
        # ----------------------

        # --- Objective ---
        myjckt.driver.add_objective('(LoadFrameOuts.Frameouts.mass[0]+LoadFrameOuts.Mpiles)/1.e6')
        # ----------------------

        # --- Design Variables ---
        myjckt.driver.add_parameter('JcktGeoIn.batter',   low=MnCnst[0],  high=MxCnst[0])
        myjckt.driver.add_parameter('Pileinputs.Dpile',   low=MnCnst[1],  high=MxCnst[1])
        myjckt.driver.add_parameter('Pileinputs.tpile',   low=MnCnst[2],  high=MxCnst[2])
        myjckt.driver.add_parameter('Pileinputs.Lp',      low=MnCnst[3],  high=MxCnst[3])
        myjckt.driver.add_parameter('leginputs.Dleg0',    low=MnCnst[4],  high=MxCnst[4])
        myjckt.driver.add_parameter('leginputs.tleg0',    low=MnCnst[5],  high=MxCnst[5])
        myjckt.driver.add_parameter('Xbrcinputs.Dbrc0',   low=MnCnst[6],  high=MxCnst[6])
        myjckt.driver.add_parameter('Xbrcinputs.tbrc0',   low=MnCnst[7],  high=MxCnst[7])
        myjckt.driver.add_parameter('Mbrcinputs.Dbrc_mud',low=MnCnst[8],  high=MxCnst[8])
        myjckt.driver.add_parameter('Mbrcinputs.tbrc_mud',low=MnCnst[9],  high=MxCnst[9])
        myjckt.driver.add_parameter('Twrinputs.Db',       low=MnCnst[10], high=MxCnst[10])
        myjckt.driver.add_parameter('Twrinputs.DTRb',     low=MnCnst[11], high=MxCnst[11])
        myjckt.driver.add_parameter('Twrinputs.Dt',       low=MnCnst[12], high=MxCnst[12])
        myjckt.driver.add_parameter('Twrinputs.DTRt',     low=MnCnst[13], high=MxCnst[13])
        myjckt.driver.add_parameter('Twrinputs.Htwr2frac',low=MnCnst[14], high=MxCnst[14])

       #--- Constraints ---#
        myjckt.driver.add_constraint('LoadFrameOuts.Frameouts.Freqs[0] >=0.22')
        myjckt.driver.add_constraint('max(LoadFrameOuts.tower_utilization.GLUtil) <=1.0')
        myjckt.driver.add_constraint('max(LoadFrameOuts.tower_utilization.EUshUtil) <=1.0')
        myjckt.driver.add_constraint('max(LoadFrameOuts2.tower_utilization.GLUtil) <=1.0')
        myjckt.driver.add_constraint('max(LoadFrameOuts2.tower_utilization.EUshUtil) <=1.0')

        myjckt.driver.add_constraint('LoadFrameOuts.jacket_utilization.t_util <=1.0')
        myjckt.driver.add_constraint('LoadFrameOuts.jacket_utilization.cb_util <=1.0')
        myjckt.driver.add_constraint('LoadFrameOuts.jacket_utilization.KjntUtil <= 1.0')
        myjckt.driver.add_constraint('LoadFrameOuts.jacket_utilization.XjntUtil <= 1.0')
        myjckt.driver.add_constraint('LoadFrameOuts2.jacket_utilization.t_util <=1.0')
        myjckt.driver.add_constraint('LoadFrameOuts2.jacket_utilization.cb_util <=1.0')
        myjckt.driver.add_constraint('LoadFrameOuts2.jacket_utilization.KjntUtil <= 1.0')
        myjckt.driver.add_constraint('LoadFrameOuts2.jacket_utilization.XjntUtil <= 1.0')

        myjckt.driver.add_constraint('PreBuild.wbase <= 30.')
        myjckt.driver.add_constraint('leginputs.Dleg0/leginputs.tleg0 >= 22.')
        myjckt.driver.add_constraint('Xbrcinputs.Dbrc0/Xbrcinputs.tbrc0 >= 22.')
        myjckt.driver.add_constraint('Mbrcinputs.Dbrc_mud/Mbrcinputs.tbrc_mud >= 22.')

        myjckt.driver.add_constraint('LoadFrameOuts.Lp0rat >= 0.')
        myjckt.driver.add_constraint('LoadFrameOuts2.Lp0rat >= 0.')
        # ----------------------

        # --- recorder ---
        myjckt.recorders = [DumpCaseRecorder('optimization.dat')]
        # ----------------------

    #--- RUN JACKET ---#
    print myjckt.SDpySubDynA.InputFile_name
    myjckt.run()
    # ---------------

    #_____________________________________#
    #Now show results
    print('Minimum found at: \n')
    printJacket(myjckt)

    #CJB+ Prepare utilization data to be plotted
    #TODO Format arrays in a dedicated component
    #Inputs: nodes=self.JcktGeoOut.nodes, t_util=self.jacket_utilization.t_util, cb_util=self.jacket_utilization.cb_util, XjntUtil=self.jacket_utilization.XjntUtil, KjntUtil=self.jacket_utilization.KjntUtil,
             #StressUtil=self.tower_utilization.StressUtil, GLUtil=self.tower_utilization.GLUtil, EUshUtil=self.tower_utilization.EUshUtil, myjckt.JcktGeoOut.nmems, myjckt.LoadFrameOuts.TPmems, myjckt.LoadFrameOuts.Twrmems
    #Outputs: member_check, KXcoords, KXcolors

    t_util_list=np.array([[1,1]])  #CJB+
    count_t=1  #CJB+
    for i in myjckt.LoadFrameOuts.jacket_utilization.t_util:  #CJB+
        if i>=1:  #CJB+
            step=np.array([count_t,'r'])  #CJB+
        elif .8<=i<1:  #CJB+
            step=np.array([count_t,'m'])  #CJB+
        elif .6<=i<.8:  #CJB+
            step=np.array([count_t,'y'])  #CJB+
        elif .4<=i<.6:  #CJB+
            step=np.array([count_t,'g'])  #CJB+
        elif .2<=i<.4:  #CJB+
            step=np.array([count_t,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_t,'k'])  #CJB+
        count_t+=1  #CJB+
        t_util_list=np.vstack((t_util_list, step))  #CJB+
    t_util_list=np.delete(t_util_list,0,0)  #CJB+
    #print t_util_list  #CJB+

    cb_util_list=np.array([[1,1]])  #CJB+
    count_cb=1  #CJB+
    for j in myjckt.LoadFrameOuts.jacket_utilization.cb_util:  #CJB+
        if j>=1:  #CJB+
            step=np.array([count_cb,'r'])  #CJB+
        elif .8<=j<1:  #CJB+
            step=np.array([count_cb,'m'])  #CJB+
        elif .6<=j<8:  #CJB+
            step=np.array([count_cb,'y'])  #CJB+
        elif .4<=j<6:  #CJB+
            step=np.array([count_cb,'g'])  #CJB+
        elif .2<=j<.4:  #CJB+
            step=np.array([count_cb,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_cb,'k'])  #CJB+
        count_cb+=1  #CJB+
        cb_util_list=np.vstack((cb_util_list, step))  #CJB+
    cb_util_list=np.delete(cb_util_list,0,0)  #CJB+
    #print cb_util_list  #CJB+

    tower_Stressutil_list=np.array([[1,1]])  #CJB+
    count_tower=1  #CJB+
    for k in myjckt.LoadFrameOuts.tower_utilization.StressUtil:  #CJB+
        if k>=1:  #CJB+
            step=np.array([count_tower,'r'])  #CJB+
        elif .8<=k<1:  #CJB+
            step=np.array([count_tower,'m'])  #CJB+
        elif .6<=k<.8:  #CJB+
            step=np.array([count_tower,'y'])  #CJB+
        elif .4<=k<.6:  #CJB+
            step=np.array([count_tower,'g'])  #CJB+
        elif .2<=k<.4:  #CJB+
            step=np.array([count_tower,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_tower,'k'])  #CJB+
        count_tower+=1  #CJB+
        tower_Stressutil_list=np.vstack((tower_Stressutil_list, step))  #CJB+
    tower_Stressutil_list=np.delete(tower_Stressutil_list,0,0)  #CJB+
    #print tower_Stressutil_list  #CJB+

    tower_GLutil_list=np.array([[1,1]])  #CJB+
    count_tower=1  #CJB+
    for k in myjckt.LoadFrameOuts.tower_utilization.GLUtil:  #CJB+
        if k>=1:  #CJB+
            step=np.array([count_tower,'r'])  #CJB+
        elif .8<=k<1:  #CJB+
            step=np.array([count_tower,'m'])  #CJB+
        elif .6<=k<.8:  #CJB+
            step=np.array([count_tower,'y'])  #CJB+
        elif .4<=k<.6:  #CJB+
            step=np.array([count_tower,'g'])  #CJB+
        elif .2<=k<.4:  #CJB+
            step=np.array([count_tower,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_tower,'k'])  #CJB+
        count_tower+=1  #CJB+
        tower_GLutil_list=np.vstack((tower_GLutil_list, step))  #CJB+
    tower_GLutil_list=np.delete(tower_GLutil_list,0,0)  #CJB+
    #print tower_GLutil_list  #CJB+

    tower_EUshutil_list=np.array([[1,1]])  #CJB+
    count_tower=1  #CJB+
    for k in myjckt.LoadFrameOuts.tower_utilization.EUshUtil:  #CJB+
        if k>=1:  #CJB+
            step=np.array([count_tower,'r'])  #CJB+
        elif .8<=k<1:  #CJB+
            step=np.array([count_tower,'m'])  #CJB+
        elif .6<=k<.8:  #CJB+
            step=np.array([count_tower,'y'])  #CJB+
        elif .4<=k<.6:  #CJB+
            step=np.array([count_tower,'g'])  #CJB+
        elif .2<=k<.4:  #CJB+
            step=np.array([count_tower,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_tower,'k'])  #CJB+
        count_tower+=1  #CJB+
        tower_EUshutil_list=np.vstack((tower_EUshutil_list, step))  #CJB+
    tower_EUshutil_list=np.delete(tower_EUshutil_list,0,0)  #CJB+
    #print tower_EUshutil_list  #CJB+

    t_condensed=np.array([[1,1]])  #CJB+
    for m in range(len(t_util_list)):  #CJB+
        counter=(m/8+1)  #CJB+
        if t_util_list[m][1]=='r':  #CJB+
            step=np.array([counter, 'r'])  #CJB+
        elif t_util_list[m][1]=='m':  #CJB+
            step=np.array([counter, 'm'])  #CJB+
        elif t_util_list[m][1]=='y':  #CJB+
            step=np.array([counter, 'y'])  #CJB+
        elif t_util_list[m][1]=='g':  #CJB+
            step=np.array([counter, 'g'])  #CJB+
        elif t_util_list[m][1]=='b':  #CJB+
            step=np.array([counter, 'b'])  #CJB+
        else:  #CJB+
            step=np.array([counter, 'k'])  #CJB+
        t_condensed=np.vstack((t_condensed,step))  #CJB+
    t_condensed=np.delete(t_condensed,0,0)  #CJB+
    t_condensed=t_condensed[0:-1:8]  #CJB+
    #print t_condensed  #CJB+

    cb_condensed=np.array([[1,1]])  #CJB+
    for m in range(len(cb_util_list)):  #CJB+
        counter=(m/4+1)  #CJB+
        if cb_util_list[m][1]=='r':  #CJB+
            step=np.array([counter, 'r'])  #CJB+
        elif cb_util_list[m][1]=='m':  #CJB+
            step=np.array([counter, 'm'])  #CJB+
        elif cb_util_list[m][1]=='y':  #CJB+
            step=np.array([counter, 'y'])  #CJB+
        elif cb_util_list[m][1]=='g':  #CJB+
            step=np.array([counter, 'g'])  #CJB+
        elif cb_util_list[m][1]=='b':  #CJB+
            step=np.array([counter, 'b'])  #CJB+
        else:  #CJB+
            step=np.array([counter, 'k'])  #CJB+
        cb_condensed=np.vstack((cb_condensed,step))  #CJB+
    cb_condensed=np.delete(cb_condensed,0,0)  #CJB+
    cb_condensed=cb_condensed[0:-1:4]  #CJB+
    #print cb_condensed  #CJB+

    jacket_check=np.array([[1,1]])  #CJB+
    for n in range(len(t_condensed)):  #CJB+
        counter=n+1  #CJB+
        if t_condensed[n][1]=='r' or cb_condensed[n][1]=='r':  #CJB+
            step=np.array([counter, 'r'])  #CJB+
        elif t_condensed[n][1]=='m' or cb_condensed[n][1]=='m':  #CJB+
            step=np.array([counter, 'm'])  #CJB+
        elif t_condensed[n][1]=='y' or cb_condensed[n][1]=='y':  #CJB+
            step=np.array([counter, 'y'])  #CJB+
        elif t_condensed[n][1]=='g' or cb_condensed[n][1]=='g':  #CJB+
            step=np.array([counter, 'g'])  #CJB+
        elif t_condensed[n][1]=='b' or cb_condensed[n][1]=='b':  #CJB+
            step=np.array([counter, 'b'])  #CJB+
        else:  #CJB+
            step=np.array([counter, 'k'])  #CJB+
        jacket_check=np.vstack((jacket_check,step))  #CJB+
    jacket_check=np.delete(jacket_check,0,0)  #CJB+
    #print jacket_check  #CJB+

    tower_check=np.array([[1,1]])  #CJB+
    for p in range(len(tower_Stressutil_list)):  #CJB+
        counter=p+len(jacket_check)+1  #CJB+
        if tower_Stressutil_list[p][1]=='r' or tower_GLutil_list[p][1]=='r' or tower_EUshutil_list[p][1]=='r':  #CJB+
            step=np.array([counter, 'r'])  #CJB+
        elif tower_Stressutil_list[p][1]=='m' or tower_GLutil_list[p][1]=='m' or tower_EUshutil_list[p][1]=='m':  #CJB+
            step=np.array([counter, 'm'])  #CJB+
        elif tower_Stressutil_list[p][1]=='y' or tower_GLutil_list[p][1]=='y' or tower_EUshutil_list[p][1]=='y':  #CJB+
            step=np.array([counter, 'y'])  #CJB+
        elif tower_Stressutil_list[p][1]=='g' or tower_GLutil_list[p][1]=='g' or tower_EUshutil_list[p][1]=='g':  #CJB+
            step=np.array([counter, 'g'])  #CJB
        elif tower_Stressutil_list[p][1]=='b' or tower_GLutil_list[p][1]=='b' or tower_EUshutil_list[p][1]=='b':  #CJB+
            step=np.array([counter, 'b'])  #CJB+
        else:  #CJB+
            step=np.array([counter, 'k'])  #CJB+
        tower_check=np.vstack((tower_check,step))  #CJB+
    tower_check=np.delete(tower_check,0,0)  #CJB+
    #print tower_check  #CJB+

    member_check=np.vstack((jacket_check,tower_check))  #CJB+
    #print member_check  #CJB+

    Kjnts_check=np.array([[1,1,1]])  #CJB+
    for i in myjckt.Build.KjntIDs:  #CJB+
        Kjnts_check=np.vstack((Kjnts_check,myjckt.JcktGeoOut.nodes[i-1]))  #CJB+
    Kjnts_check=np.delete(Kjnts_check,0,0)  #CJB+

    K_count=1  #CJB+
    K_joints_color=np.array([[1,1]])  #CJB+
    for i in myjckt.LoadFrameOuts.jacket_utilization.KjntUtil:  #CJB+
        if i>=1:  #CJB+
            step=np.array([K_count,'r^'])  #CJB+
        elif .8<=i<1:  #CJB+
            step=np.array([K_count,'mo'])  #CJB+
        elif .6<=i<.8:  #CJB+
            step=np.array([K_count,'yo'])  #CJB+
        elif .4<=i<.6:  #CJB+
            step=np.array([K_count,'go'])  #CJB+
        elif .2<=i<.4:  #CJB+
            step=np.array([K_count,'bo'])  #CJB+
        else:  #CJB+
            step=np.array([K_count,'ko'])  #CJB+
        K_count+=1  #CJB+
        K_joints_color=np.vstack((K_joints_color,step))  #CJB+
    K_joints_color=np.delete(K_joints_color,0,0)  #CJB+
    #print K_joints_color  #CJB+

    Xjnts_check=np.array([[1,1,1]])  #CJB+
    for i in myjckt.Build.XjntIDs:  #CJB+
        Xjnts_check=np.vstack((Xjnts_check,myjckt.JcktGeoOut.nodes[i-1]))  #CJB+
    Xjnts_check=np.delete(Xjnts_check,0,0)  #CJB+

    X_count=1  #CJB+
    X_joints_color=np.array([[1,1]])  #CJB+
    for i in myjckt.LoadFrameOuts.jacket_utilization.XjntUtil:  #CJB+
        if i>=1:  #CJB+
            step=np.array([X_count,'r^'])  #CJB+
        elif .8<=i<1:  #CJB+
            step=np.array([X_count,'mo'])  #CJB+
        elif .6<=i<.8:  #CJB+
            step=np.array([X_count,'yo'])  #CJB+
        elif .4<=i<.6:  #CJB+
            step=np.array([X_count,'go'])  #CJB+
        elif .2<=i<.4:  #CJB+
            step=np.array([X_count,'bo'])  #CJB+
        else:  #CJB+
            step=np.array([X_count,'ko'])  #CJB+
        X_count+=1  #CJB+
        X_joints_color=np.vstack((X_joints_color,step))  #CJB+
    X_joints_color=np.delete(X_joints_color,0,0)  #CJB+
    #print X_joints_color  #CJB+

    KXcoords=np.vstack((Kjnts_check,Xjnts_check))  #CJB+
    KXcolors=np.vstack((K_joints_color,X_joints_color))  #CJB+
    #print KXcoords  #CJB+
    #print KXcolors  #CJB+

    #Plot geometry

    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = Axes3D(fig)
    nodes=myjckt.JcktGeoOut.nodes  #coordinates
    mems=myjckt.JcktGeoOut.mems
    XYZ1=nodes[mems[:,0]-1,:]
    XYZ2=nodes[mems[:,1]-1,:]
    #ax.scatter(nodes[:,0],nodes[:,1],nodes[:,2])
    x,y,z=zip(*nodes)
    #ax.scatter(x,y,z)
    Xs=np.vstack((XYZ1[:,0],XYZ2[:,0])).T
    Ys=np.vstack((XYZ1[:,1],XYZ2[:,1])).T
    Zs=np.vstack((XYZ1[:,2],XYZ2[:,2])).T
    #ax.plot([XYZ1[1:5,0],XYZ2[1:5,0]],[XYZ1[1:5,1],XYZ2[1:5,1]],[XYZ1[1:5,2],XYZ2[1:5,2]])
    for i in range(0,mems.shape[0]): #mems.shape[0])
        ax.plot([XYZ1[i,0],XYZ2[i,0]],[XYZ1[i,1],XYZ2[i,1]],[XYZ1[i,2],XYZ2[i,2]],member_check[i][1]) #CJBe Member check gives memebrs utilization colors
        #CJB TODO If RigidTop=False, then the code does not make the top member. However, the utilization plot (Figure 2) will still show the utiliazation of the top member.  (There is a different number of tower members being plotted)
    axisEqual3D(ax)
#    ax.set_aspect('equal')
#    ax.auto_scale_xyz([min(XYZ1[:,0]),max(XYZ1[:,0])],[min(XYZ1[:,1]),max(XYZ1[:,1])],[min(XYZ1[:,2]),max(XYZ1[:,2])])
    plt.show()

   #Plot tower utilization
    #twr_zs=myjckt.Tower.Twrouts.nodes[2, :-int(myjckt.TwrRigidTop)] #CJB-
    twr_zs=myjckt.Tower.Twrouts.nodes[2,0:myjckt.Tower.Twrouts.nNodes-int(myjckt.TwrRigidTop)]-myjckt.Tower.Twrouts.nodes[2,0] #CJB+ This line, taken from PlotJacket.py, plots the utilization with or without a RigidMember
    twr_ds=np.hstack((myjckt.Tower.Twrouts.TwrObj.D,myjckt.Tower.Dt))

    fig2= plt.figure(2,figsize=(6.0, 3.5))

    ax1=plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=3)

    ax1.plot([],[],'k',label='Tower Profile') #CJB+
    ax1.plot(myjckt.LoadFrameOuts.tower_utilization.StressUtil,twr_zs , label='VonMises Util')
    ax1.plot(myjckt.LoadFrameOuts.tower_utilization.GLUtil,twr_zs, label='GL Util')
    ax1.plot(myjckt.LoadFrameOuts.tower_utilization.EUshUtil, twr_zs, label='EUsh Util')
    if myjckt.LoadFrameOuts2.tower_utilization.StressUtil[0] != -9999.:
        ax1.plot(myjckt.LoadFrameOuts2.tower_utilization.StressUtil,twr_zs , label='VonMises Util2')
        ax1.plot(myjckt.LoadFrameOuts2.tower_utilization.GLUtil,twr_zs, label='GL Util2')
        ax1.plot(myjckt.LoadFrameOuts2.tower_utilization.EUshUtil, twr_zs, label='EUsh Util2')

    ax1.legend(bbox_to_anchor=(1.05, 1.0), loc=2)

    ax2=plt.twiny(ax1)# .axes(ax1.get_position(),frameon=False)
    #hh2=plt.axes('Position', pos,'NextPlot','Add','XtickLabel','','Xtick',[],frameon=False);
    ax2.plot(twr_ds/2,twr_zs,'k');
    ax2.plot(- twr_ds/2,twr_zs,'k');
    ax2.set_aspect('equal')
    ax2.set_frame_on(False)
    ax2.set_xticklabels('')
    ax2.set_xticks([])
    ax1.set_xlim([0,2]);

    plt.xlabel('utilization')
    plt.ylabel('height along tower (m)')
    plt.show()

    fig3= plt.figure(3) #CJB+
    ax = Axes3D(fig3) #CJB+
    for i in range(myjckt.JcktGeoOut.nmems-len(myjckt.LoadFrameOuts.TPmems)-len(myjckt.LoadFrameOuts.Twrmems)): #CJB+
        plt.plot([XYZ1[i,0],XYZ2[i,0]],[XYZ1[i,1],XYZ2[i,1]],[XYZ1[i,2],XYZ2[i,2]],'.75') #CJB+
    for i in range(0,KXcoords.shape[0]): #CJB+
        plt.plot([KXcoords[i,0]],[KXcoords[i,1]],[KXcoords[i,2]], KXcolors[i][1]) #CJB+
    axisEqual3D(ax) #CJB+
    plt.show() #CJB+

##________________SAMPLE CALL from outside IDE________________##
## python jacket.py True PyOPTSnopt
##___________________________________________##
