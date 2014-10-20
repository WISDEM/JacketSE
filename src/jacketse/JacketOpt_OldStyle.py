#-------------------------------------------------------------------------------
# Name:        JacketOptMDAO_OldStyle.py
# Purpose:     It solves for minimum mass problem with constraints on max utilization.
#              It uses the old style wrapper with the new mdao framework. See JacketOptMDAO.py for
#              the all-in optimizer, all included in the framework.
# Author:      rdamiani
#
# Created:     4/25/2014 based off of JacketOpt.py
# Copyright:   (c) rdamiani 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np
import scipy.optimize
import os
import cPickle as pickle;

from openmdao.main.api import Assembly
from openmdao.main.datatypes.api import Int, Float, Array, VarTree, Bool, Dict,Instance
from openmdao.lib.drivers.api import COBYLAdriver
from openmdao.main.api import enable_console
from openmdao.lib.casehandlers.api import DumpCaseRecorder,ListCaseRecorder,CSVCaseRecorder, DBCaseRecorder

#enable_console()
import logging
logging.getLogger().setLevel(logging.DEBUG)

from jacket import JacketAsmly,JcktGeoInputs,LegGeoInputs,XBrcGeoInputs,MudBrcGeoInputs,HBrcGeoInputs,PileGeoInputs,TPGeoInputs,TwrGeoInputs,RNAprops,TPlumpMass,WaterInputs,WindInputs,SoilGeoInputs,Frame3DDaux,IEC_PSFS

from PlotJacket import main as PlotJacket

import SetJacketInputsUH  as SetJacketInputs      #This is for the specific project!!!!!

    #Build up a wrapepr that forms the Jacket and Runs it, Calculates stresses and utilization
def JcktWrapper(x,myjckt,avgcnst):
    # x=  [ batter,  Dpile,    tpile,   Lp,   Dleg,     tleg,   Db,   DTRb   Dt,   DTRt,   Htwr2fac ,     Dbrc,   tbrc,     Dbrc_mud,   tbrc_mud   ]
    """This function builds up the actual model and calculates Jacket stresses and utilization.
    INPUT
        x      -list(11, or 15), as follows:
            batter - float, 2D batter value (e.g. 8)
            Dpile  - float, pile OD [m]
            tpile  - float, pile t [m]
            Dleg   -float, leg OD [m]
            tleg   -float, leg thickness [m]
            Db     -float, tower base OD [m]
            DTRb    -float, tower D/t
            Dt     -float, tower base OD [m]
            DTRt    -float, tower D/t
            HTwr2frac  -float, tower constant D segment length as fraction of overall length
            Dbrace   -float, brace OD [m]
            tbrace   -float, brace t [m]
            Dbrc_mud  -float, mudbrace OD [m]
            tbrc_mud  -float, mudbrace t [m]

            !!!!ALL NORMALIZED BY THEIR AVERAGES!!!!

                If 11 elements, then Dbrc,tbrc,Dbrc,tbrch are calculated with an internal optimizer (not recommended)
            myjckt    -assmebly of jacketSE
        """
    global  xlast,f1,max_GLUtil,max_EUUtil,max_KjntUtil,max_XjntUtil,max_tutil,max_cbutil,offset,\
            MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat


    myjckt.JcktGeoIn.batter=x[0]*avgcnst[0] #remove this after testing

    if len(x)>8:
        offset=6
        myjckt.Pileinputs.Dpile=x[1]*avgcnst[1]
        myjckt.Pileinputs.tpile=x[2]*avgcnst[2]
        myjckt.Pileinputs.Lp   =x[3]*avgcnst[3]
        myjckt.leginputs.Dleg  =np.asarray([x[4]*avgcnst[4]]).repeat(myjckt.JcktGeoIn.nbays+1)
        myjckt.leginputs.tleg  =np.asarray([x[5]*avgcnst[5]]).repeat(myjckt.JcktGeoIn.nbays+1)


    myjckt.Twrinputs.Db    =x[0+offset]*avgcnst[0+offset]
    myjckt.Twrinputs.DTRb  =x[1+offset]*avgcnst[1+offset]
    myjckt.Twrinputs.Dt    =x[2+offset]*avgcnst[2+offset]
    myjckt.Twrinputs.DTRt  =myjckt.Twrinputs.DTRb#x[3+offset]*avgcnst[3+offset]  #myjckt.Twrinputs.DTRb.copy() #x[3+offset]*avgcnst[3+offset]  #can use DTRt=DTRb here for simplicity

    myjckt.Twrinputs.Htwr2frac =x[4+offset]*avgcnst[4+offset]

    myjckt.TPinputs.Dgir  =x[5+offset]*avgcnst[5+offset]
    myjckt.TPinputs.tgir  =x[6+offset]*avgcnst[6+offset]

    myjckt.Xbrcinputs.precalc=True #Initialize
    myjckt.Mbrcinputs.precalc=True #Initialize

    if len(x)>13:
        myjckt.Xbrcinputs.precalc=False
        myjckt.Xbrcinputs.Dbrc    =np.asarray([x[13]*avgcnst[13]]).repeat(myjckt.JcktGeoIn.nbays)
        myjckt.Xbrcinputs.tbrc    =np.asarray([x[14]*avgcnst[14]]).repeat(myjckt.JcktGeoIn.nbays)
        myjckt.Mbrcinputs.precalc=False
        myjckt.Mbrcinputs.Dbrc_mud =x[15]*avgcnst[15]
        myjckt.Mbrcinputs.tbrc_mud =x[16]*avgcnst[16]

    #Run the assembly and get main output
    myjckt.run()

    #Get Frame3dd mass
    mass=myjckt.Frameouts2.mass[0] +myjckt.Mpiles  #Total structural mass
    #Get Frame3dd-calculated 1st natural frequency
    f1=myjckt.Frameouts2.Freqs[0]

    #Get Model calculated TP mass
    TPmass=myjckt.TP.TPouts.mass

    #Get Utilizations
    max_GLUtil=np.nanmax(myjckt.tower_utilization.GLUtil)
    #max_GLUtil=myjckt.tower_utilization.GLUtil
    max_EUUtil=np.nanmax(myjckt.tower_utilization.EUshUtil)

    #Member checks
    max_tutil=np.nanmax(myjckt.jacket_utilization.t_util)
    max_cbutil=np.nanmax(myjckt.jacket_utilization.cb_util)

    #Joint checks
    max_XjntUtil=np.nanmax(myjckt.jacket_utilization.XjntUtil)
    max_KjntUtil=np.nanmax(myjckt.jacket_utilization.KjntUtil)
    #t_util[t_util<0]=np.NaN
    #cb_util[cb_util<0]=np.NaN

    #Brace Criteria
    MudCrit01=np.nanmax(myjckt.MudBrcCriteria.brc_crit01)
    MudCrit02=np.nanmax(myjckt.MudBrcCriteria.brc_crit02)
    MudCrit03=np.nanmax(myjckt.MudBrcCriteria.brc_crit03)
    MudCrit04=np.nanmax(myjckt.MudBrcCriteria.brc_crit04)
    MudCrit05=np.nanmax(myjckt.MudBrcCriteria.brc_crit05)
    XBrcCrit01=np.nanmax(myjckt.XBrcCriteria.brc_crit01)
    XBrcCrit02=np.nanmax(myjckt.XBrcCriteria.brc_crit02)
    XBrcCrit03=np.nanmax(myjckt.XBrcCriteria.brc_crit03)
    XBrcCrit04=np.nanmax(myjckt.XBrcCriteria.brc_crit04)
    XBrcCrit05=np.nanmax(myjckt.XBrcCriteria.brc_crit05)

    #Pile Embedment Criteria
    Lp0rat=myjckt.Lp0rat
    #__________________________________________#

    #calc width at seabed proportional to stiffness

    xlast=x.copy() #update the latest set of input params to the current


    print('Jwrapper SOLUTION: bat={:5.2f}, Dpile={:5.2f}, tpile={:5.3f}, Lp={:5.1f} Dleg{:5.2f}, tleg{:5.3f}').\
            format(myjckt.JcktGeoIn.batter,myjckt.Piles.Pileinputs.Dpile,myjckt.Piles.Pileinputs.tpile,myjckt.Piles.Pileinputs.Lp,myjckt.Legs.leginputs.Dleg[0],myjckt.Legs.leginputs.tleg[0])
    print('           Dbrc={:5.2f}, tbrc={:5.3f}, Dmudbrc={:5.2f}, tmudbrc={:5.3f}').\
            format(myjckt.Xbrcinputs.Dbrc[0],myjckt.Xbrcinputs.tbrc[0],myjckt.Mbrcinputs.Dbrc_mud,myjckt.Mbrcinputs.tbrc_mud)

    print('from Jwrapper Db={:5.2f}, DTRb={:5.2f}, Dt={:5.2f}, DTRt={:5.2f},H2twrfrac={:5.2f}, Dgir={:5.2f},tgir={:5.3f}, Twrmass={:6.3f}, PilesMass ={:6.3f}, TPmass= {:8.3e}, Frame3DD+Piles Totmass={:10.3f}'.\
            format(myjckt.Twrinputs.Db,myjckt.Twrinputs.DTRb,myjckt.Twrinputs.Dt,\
                   myjckt.Twrinputs.DTRt,myjckt.Twrinputs.Htwr2frac,myjckt.TPinputs.Dgir,myjckt.TPinputs.tgir,myjckt.Tower.Twrouts.mass,myjckt.Mpiles, myjckt.TP.TPouts.mass,mass))
    print(' \n')

    return mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat


#constraints x[0]=batter x[1]=Dbrc, x[2]=tbrc, x[2]=lbrc

def f0Cnstrt1(x,myjckt,avgcnst):  #f1>f0
    global xlast,f1
    #print('from const x={:10.7f}'.format(x[0]))
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    cnstrt=(f1-f0)/f0
    print('f0Cnstrt1=',cnstrt)
    return cnstrt

def f0Cnstrt2(x,myjckt,avgcnst): #f1<(f0*(1+f0eps))
    global xlast,f1,f0eps
    #print('from const x={:10.7f}'.format(x[0]))
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    cnstrt=(-f1+f0*(1+f0eps))/f0
    print('f0Cnstrt2=',cnstrt)
    return cnstrt

def cbCnstrt(x,myjckt,avgcnst): #cb<1
    global xlast,f1,max_cbutil
    #print('from const x={:10.7f}'.format(x[0]))
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    cnstrt=1.-max_cbutil
    print('cbCnstrt=',cnstrt)
    return cnstrt

def tCnstrt(x,myjckt,avgcnst): #tUtil<1
    global xlast,max_tutil
    #print('from const x={:10.7f}'.format(x[0]))
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    cnstrt=1.-max_tutil
    print('tutil constraint=',cnstrt)
    return cnstrt


def KjntCnstrt(x,myjckt,avgcnst): #KUtil<1
    global xlast,max_KjntUtil
    #print('from const x={:10.7f}'.format(x[0]))
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    cnstrt=1.-max_KjntUtil
    print('KjntCnstrt=',cnstrt)
    return cnstrt

def XjntCnstrt(x,myjckt,avgcnst): #XUtil<1
    global xlast,max_XjntUtil
    #print('from const x={:10.7f}'.format(x[0]))
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    cnstrt=1.-max_XjntUtil
    print('XjntCnstrt=',cnstrt)
    return cnstrt

def GLCnstrt(x,myjckt,avgcnst): #GLUtil<1
    global xlast,max_GLUtil
    #print('from const x={:10.7f}'.format(x[0]))
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    cnstrt=1.-max_GLUtil
    print('GL constraint=',cnstrt)
    return cnstrt

def EUCnstrt(x,myjckt,avgcnst): #EUUtil<1
    global xlast,max_EUUtil
    #print('from const x={:10.7f}'.format(x[0]))
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
          MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
          XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    cnstrt=1.-max_EUUtil
    print('EU constraint=',cnstrt)
    return cnstrt

#NOW use constraints for bounds on batter, Dleg, tleg, Dp, Db, DTR
def batCnsrt1(x,myjckt,avgcnst):  #bound of batter for cobyla batter>7
    cnstrt=x[0]*batavg-batmin
    print('batter constraint1=',cnstrt)
    return cnstrt   #9.5=batvg
def batCnsrt2(x,myjckt,avgcnst):  #bound of batter for cobyla batter<max7
    cnstrt=-x[0]*batavg+batmax
    print('batter constraint2=',cnstrt)
    return cnstrt   #9.5=batvg

#Leg Constraints
def DlegCnstrt(x,myjckt,avgcnst):  #bound of Dleg for cobyla Dleg>1
    cnstrt=x[4]*Dlegavg-Dlegmin
    print('Dleg constraint=', cnstrt)
    return cnstrt
def tlegCnstrt(x,myjckt,avgcnst):  #bound of t for cobyla t>1"
    cnstrt=x[5]*tlegavg-tlegmin
    print('tleg constraint=',cnstrt)
    return cnstrt#

#Brace Size Constraints
def DbrcCnstrt(x,myjckt,avgcnst):  #bound of Dbrc for cobyla Dleg>1
    cnstrt=x[13]*Dbrcavg-Dbrcmin
    print('Dbrc constraint=', cnstrt)
    return cnstrt
def tbrcCnstrt(x,myjckt,avgcnst):  #bound of t for cobyla t>1"
    cnstrt=x[14]*tbrcavg-tbrcmin
    print('tbrc constraint=',cnstrt)
    return cnstrt#

#Mudbrace Size Constraints
def DmudCnstrt(x,myjckt,avgcnst):  #bound of Dbrc for cobyla Dleg>1
    cnstrt=x[15]*Dmudavg-Dmudmin
    print('Dmud constraint=', cnstrt)
    return cnstrt
def tmudCnstrt(x,myjckt,avgcnst):  #bound of t for cobyla t>1"
    cnstrt=x[16]*tmudavg-tmudmin
    print('tmud constraint=',cnstrt)
    return cnstrt#

#Pile Constraints
def DpCnstrtbattered(x):  #bound of Dp>1
    Dp=x[1]*Davg-2*x[2]*tavg-3.5*0.0254 #pile diameter for battered-pile (driven through legs) configuration
    return Dp-1.  #this is a proxy to the results in mcad showing that a Dp<1.1m would not pass API checks, It should boost up the Dleg to get towards 1.713

def DpCnstrt(x,myjckt,avgcnst):  #bound of Dp>1
    cnstrt=x[1]*Dpavg-Dpmin
    print('Dp constraint=', cnstrt)
    return cnstrt

def tpCnstrt(x,myjckt,avgcnst):  #tp >0.0254
    cnstrt=x[2]*tpavg-tpmin
    print('tp constraint=', cnstrt)
    return cnstrt


#TP constraints
def girDCnsrt1(x,myjckt,avgcnst):  #bound of girder D>Dmin
    cnstrt=x[5+offset]*girDavg-girDmin
    print('gir D constraint1=',cnstrt)
    return cnstrt
def girDCnsrt2(x,myjckt,avgcnst):  #bound of girder D<Dmax
    cnstrt=-x[5+offset]*girDavg+girDmax
    print('gir D constraint2=',cnstrt)
    return cnstrt
def girtCnsrt1(x,myjckt,avgcnst):  #bound of girder t>tmin
    cnstrt=x[6+offset]*girtavg-girtmin
    print('gir t constraint1=',cnstrt)
    return cnstrt
def girtCnsrt2(x,myjckt,avgcnst):  #bound of girder t<tmax
    cnstrt=-x[6+offset]*girtavg+girtmax
    print('gir t constraint2=',cnstrt)
    return cnstrt


 #Xbrc constraints
def XbCrit01(x,myjckt,avgcnst): #XbrcCrit01
    global xlast,XBrcCrit01
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Xbrc constraint01=', XBrcCrit01)
    return XBrcCrit01

def XbCrit02(x,myjckt,avgcnst): #XbrcCrit02
    global xlast,XBrcCrit02
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Xbrc constraint02=', XBrcCrit02)
    return XBrcCrit02

def XbCrit03(x,myjckt,avgcnst): #XbrcCrit03
    global xlast,XBrcCrit03
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Xbrc constraint03=', XBrcCrit03)
    return XBrcCrit03

def XbCrit04(x,myjckt,avgcnst): #XbrcCrit03
    global xlast,XBrcCrit04
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Xbrc constraint04=', XBrcCrit04)
    return XBrcCrit04

def XbCrit05(x,myjckt,avgcnst): #XbrcCrit05
    global xlast,XBrcCrit05
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Xbrc constraint05=', XBrcCrit05)
    return XBrcCrit05

 #Mudbrc constraints
def MbCrit01(x,myjckt,avgcnst): #MbrcCrit01
    global xlast,MudCrit01
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Mbrc constraint01=', MudCrit01)
    return MudCrit01

def MbCrit02(x,myjckt,avgcnst): #XbrcCrit02
    global xlast,MudCrit02
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Mbrc constraint02=', MudCrit02)
    return MudCrit02

def MbCrit03(x,myjckt,avgcnst): #XbrcCrit03
    global xlast,MudCrit03
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Mbrc constraint03=', MudCrit03)
    return MudCrit03

def MbCrit04(x,myjckt,avgcnst): #XbrcCrit03
    global xlast,MudCrit04
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Mbrc constraint04=', MudCrit04)
    return MudCrit04

def MbCrit05(x,myjckt,avgcnst): #XbrcCrit05
    global xlast,MudCrit05
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05=JcktWrapper(x,myjckt,avgcnst)

        xlast=x.copy()
    print('Mbrc constraint05=', MudCrit05)
    return MudCrit05

 #Tower constraints
def TwrCnstrt01(x,myjckt,avgcnst): #Maximum tower Diameter
    cnstrt=Dbmax-x[0+offset]*Dbavg
    print('Db TwrCnstrt01=',cnstrt )
    return cnstrt
def TwrCnstrt02(x,myjckt,avgcnst): #Maximum tower Diameter
    cnstrt=Dtmax-x[2+offset]*Dtavg
    print('Dt TwrCnstrt02=', cnstrt)
    return cnstrt
def TwrCnstrt03(x,myjckt,avgcnst): #Minimum tower Diameter at the top
    cnstrt=x[2+offset]*Dtavg-Dtmin
    print('Dt TwrCnstrt03=', cnstrt)
    return cnstrt
def TwrCnstrt04(x,myjckt,avgcnst):  #Minimum DTRb>120
    cnstrt=x[1+offset]*DTRavg-DTRmin
    #print('DTR min TwrCnstrt04=', cnstrt)
    return cnstrt
def TwrCnstrt05(x,myjckt,avgcnst):  #Max DTRb<200
    cnstrt=DTRmax-x[1+offset]*DTRavg
    #print('DTR max TwrCnstrt05=',cnstrt)
    return cnstrt
def TwrCnstrt06(x,myjckt,avgcnst):  #Minimum DTRt>120
    cnstrt=x[3+offset]*DTRavg-DTRmin
    #print('DTR min TwrCnstrt06=', cnstrt)
    return cnstrt
def TwrCnstrt07(x,myjckt,avgcnst):  #Max DTRt<200
    cnstrt=DTRmax-x[3+offset]*DTRavg
    #print('DTR max TwrCnstrt07='cnstrt)
    return cnstrt
def TwrCnstrt08(x,myjckt,avgcnst):  #Maximum Htwr2 < Htwr/4
    cnstrt=Htwr2max-x[4+offset]*Htwr2avg
    print('Htwr2 max TwrCnstrt08=',cnstrt)
    return cnstrt


#Embedment length constraint
def LpCnstrt(x,myjckt,avgcnst):  #Maximum Htwr2 < Htwr/4
    global xlast,Lp0rat
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
         MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,\
         XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,avgcnst)
        xlast=x.copy()
    print('Lp0rat constraint=',Lp0rat)
    return Lp0rat

#Embedment length constraint2
def LpCnstrt2(x,myjckt,avgcnst):  #Maximum Lp < Lpmax
    cnstrt=Lpmax-x[3]*Lpavg
    print('Lp constraint2=',cnstrt)
    return cnstrt

    #Function to minimize
def mass(x,myjckt,avgcnst):
    """This function assembles output data from the assembly in terms of mass: \n
        it is mostly frame3dd mass+piles'' mass."""
    return  JcktWrapper(x,myjckt,avgcnst)[0]/1.5e6


#______________________________________________________________________________#


def main():
    global batmax,batmin,batavg,Dpmax,Dpmin,Dpavg,tpmin,tpmax,tpavg,Lpmax,Lpmin,Lpavg,Dlegmax,Dlegmin,Dlegavg,tlegmax,tlegmin,tlegavg,\
           Dbmax,Dbmin,Dbavg,Dtmax,Dtavg,Dtmin,DTRmin,DTRmax,DTRavg,Htwr2avg,Htwr2max,girDmin,girDmax,girDavg,girtmin,girtmax,girtavg,\
           Dbrcmax,Dbrcmin,Dbrcavg,tbrcmax,tbrcmin,tbrcavg,Dmudmax,Dmudmin,Dmudavg,tmudmax,tmudmin,tmudavg,\
           f0,f0eps,offset

    offset=1 #index offset for testing

    filename="C:\PROJECTS\OFFSHORE_WIND\UH_REACT\PYTHON_OPT\JCKdataLumpedMass.dat"

    #1. set the target frequency and the epsilon of acceptance
    f0=0.21  #[Hz]
    f0eps=0.03#(1+f0eps)*f0 will be the maximum allowed frequency

    #2.Set bounds for optimization variables

    #          x=  [ batter,  Dpile,    tpile,   Lp,   Dleg,     tleg,          Db,   DTRb   Dt,   DTRt   Htwr2fac, Dbrc,   tbrc,     Dbrc_mud,   tbrc_mud,  Dgir,   tgir ]
    #                  0        1        2       3       4           5          6      7     8    9          10       11      12           13        14      15       16
    MnCnst=np.array([ 8.,      1.5,   1.5*0.0254, 30.,   1.5,     1.5*0.0254,    6.,   120.,  3.5,   120.,     0. ,      0.8,   1.*0.0254,   1.,     1.*0.0254,  0.8,   0.75*0.0254])
    MxCnst=np.array([ 30.,     2.5,   5.*0.0254, 60.,   2.5,      5.*0.0254,    7.,   190.,  4.,   190.,     0.25,     2.,    5.*0.0254,   2.,     5.*0.0254,  1.2,    5.*0.0254])
    avgcnst=(MnCnst+MxCnst)/2.


    #3. Set up main variables and create assembly with main inputs via SetJacketInputs.py
                                 #batter,   Dpile,tpile,        Lp     ,Dleg,    tleg,    Dbrc,       tbrc,      Dbrc_mud,    tbrc_mud       Dgir,     tgir,         Db,      DTRb,         Dt,       DTRt,     Htwr2frac
    myjckt=SetJacketInputs.main(avgcnst[0],2.0, avgcnst[2],avgcnst[3],   2.2,avgcnst[5], avgcnst[11],avgcnst[12],avgcnst[13],avgcnst[14],avgcnst[15],avgcnst[16], avgcnst[6],avgcnst[7], avgcnst[8],avgcnst[9],avgcnst[10])
    #myjckt.run()
    myjckt.recorders = [CSVCaseRecorder(filename=filename)]

   #4. Initial guess :Make Sure the size is appropriate for the number of design variables to explore


    #MnCnst=MnCnst[6:11]
    #MxCnst=MxCnst[6:11]
    #avgcnst=avgcnst[6:11]
    batmax=MxCnst[0]
    batmin=MnCnst[0]
    batavg=avgcnst[0]
    Dbmax=MxCnst[6]
    Dbmin=MnCnst[6]
    Dbavg=avgcnst[6]
    Dtmax=MxCnst[8]
    Dtmin=MnCnst[8]
    Dtavg=avgcnst[8]
    DTRmax=MxCnst[7]
    DTRmin=MnCnst[7]
    DTRavg=avgcnst[7]

    Htwr2max=MxCnst[10]
    Htwr2avg=avgcnst[10]

    girDmax=MxCnst[15]
    girDmin=MnCnst[15]
    girDavg=avgcnst[15]
    girtmax=MxCnst[16]
    girtmin=MnCnst[16]
    girtavg=avgcnst[16]

    Dpmax=MxCnst[1]
    Dpmin=MnCnst[1]
    Dpavg=avgcnst[1]
    tpmax=MxCnst[2]
    tpmin=MnCnst[2]
    tpavg=avgcnst[2]

    Lpmax=MxCnst[3]
    Lpmin=MnCnst[3]
    Lpavg=avgcnst[3]


    Dlegmax=MxCnst[4]
    Dlegmin=MnCnst[4]
    Dlegavg=avgcnst[4]
    tlegmax=MxCnst[5]
    tlegmin=MnCnst[5]
    tlegavg=avgcnst[5]

    Dbrcmax=MxCnst[11]
    Dbrcmin=MnCnst[11]
    Dbrcavg=avgcnst[11]
    tbrcmax=MxCnst[12]
    tbrcmin=MnCnst[12]
    tbrcavg=avgcnst[12]

    Dmudmax=MxCnst[13]
    Dmudmin=MnCnst[13]
    Dmudavg=avgcnst[13]
    tmudmax=MxCnst[14]
    tmudmin=MnCnst[14]
    tmudavg=avgcnst[14]

    #make sure you get the right size and cosntants here, depending on the design variable subset
    avgcnstin=avgcnst[np.array([0,1,2,3,4,5,6,7,8,9,10,15,16,11,12,13,14])]
    #   x 0     1     2         3     4    5       6  7        8  9        10        11    12        13      14      15      16
    #   batter,Dp,   tp,       Lp,   Dleg,tleg,    Db,DTRb,   Dt,DTRt,     H2frac,   Dgir,tgir       Dbrc   tbrc    Dmdbrc   tmdbrc
    guess=[9., 1.6,1.*0.0254,40.,   1.6,1.*0.0254, 6.6, 120., 3.5, 120., 0.25,   1.0,  1.*0.0254, 1.0,1.*0.0254, 1.0,1.*0.0254]/avgcnstin

    import time
    tt = time.time()

    #Cobyla : If it returns 3 elements , it means it did not converge;  constraint2,constraint3
    res1=scipy.optimize.fmin_cobyla(mass,guess,[f0Cnstrt1,f0Cnstrt2,batCnsrt1, batCnsrt2, \
               cbCnstrt,tCnstrt,KjntCnstrt,XjntCnstrt, GLCnstrt,EUCnstrt,\
                girDCnsrt1,girDCnsrt2,girtCnsrt1,girtCnsrt2,\
               TwrCnstrt01,TwrCnstrt02,TwrCnstrt03,TwrCnstrt04,TwrCnstrt05,TwrCnstrt06,TwrCnstrt07,\
               XbCrit01,XbCrit02,XbCrit03,XbCrit04,XbCrit05,MbCrit01,MbCrit02,MbCrit03,MbCrit04,MbCrit05,\
               LpCnstrt,LpCnstrt2, DlegCnstrt,tlegCnstrt,DpCnstrt,tpCnstrt,\
               DbrcCnstrt,tbrcCnstrt,DmudCnstrt,tmudCnstrt],args=[myjckt,avgcnstin],consargs=None,rhobeg=0.01,maxfun=2000,disp=1)

    print "\n"
    print "Minimum mass Mjacket, MPiles, TPmass = %f %f %f" %(myjckt.Frameouts2.mass[0],myjckt.Mpiles,myjckt.TP.TPouts.mass)
    print "Minimum mass Tower, Jacket(no tower no piles) = %f %f" %(myjckt.Tower.Twrouts.mass,myjckt.Frameouts2.mass[0]-myjckt.Tower.Twrouts.mass)
    print "Minimum found at Dpile=%f, tpile=%f  Lp=%f " % (myjckt.Piles.Pileinputs.Dpile,myjckt.Piles.Pileinputs.tpile,myjckt.Piles.Pileinputs.Lp)
    print "Minimum found at Dbrc=%f, tbrc=%f  " % (myjckt.Xbraces.Xbrcouts.LLURObj.D[0],myjckt.Xbraces.Xbrcouts.LLURObj.t[0])
    print "Minimum found at Dbrcmud=%f, tbrcmud=%f  " % (myjckt.Mudbraces.Mbrcouts.brcObj.D[0],myjckt.Mudbraces.Mbrcouts.brcObj.t[0])
    print "Minimum found at batter=%f, Dleg=%f, tleg=%f,  " % (myjckt.JcktGeoIn.batter,myjckt.leginputs.Dleg[0],myjckt.leginputs.tleg[0])
    print "Minimum found at Dgir=%f, tgir=%f " % (myjckt.TPinputs.Dgir,myjckt.TPinputs.tgir)
    print "Minimum found at Db=%f DTRb=%f Dt=%f DTRt=%f H2frac=%f " % (myjckt.Tower.Twrins.Db,myjckt.Tower.Twrins.DTRb,myjckt.Tower.Twrins.Dt,myjckt.Tower.Twrins.DTRt,myjckt.Tower.Twrins.Htwr2frac)
    print "Minimum found at Freq %f"  % (myjckt.Frameouts2.Freqs[0])
    print "Minimum found at GLutil=%f EUutil=%f"  % (np.nanmax(myjckt.tower_utilization.GLUtil),np.nanmax(myjckt.tower_utilization.EUshUtil))
    print "Minimum found at Mudline Footprint=%f"  % (myjckt.wbase)
    print "Elapsed time: ", time.time()-tt, "seconds"
    print "Execution count: ", myjckt.exec_count

    #plot jacket and tower utilization
    PlotJacket(myjckt,util=True)
    #now Plot utilization with no thrust load whatsoever to see the effect of self-weight alone
##    RNA_F=np.copy(myjckt.RNA_F)
##    myjckt.RNA_F=np.array([0.,0.,0.,0.,0.,0.])
##    myjckt.run()
##    PlotJacket(myjckt,util=True)
##    #Restore the load
##    myjckt.RNA_F=np.copy(RNA_F)
##    myjckt.run()

    #Save results
    filename="C:\PROJECTS\OFFSHORE_WIND\UH_REACT\PYTHON_OPT\JCKdataNewCd.dat"
    fp=open(filename,'w')

    myjckt.recorders = [DumpCaseRecorder(fp)]
    myjckt.run()  # just to write case
    fp.close()

    #filename="D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\UH_REACT\JCKdataNoMass.p"
    filename="C:\PROJECTS\OFFSHORE_WIND\UH_REACT\PYTHON_OPT\JCKdataNewCd.txt"
    with open(filename,'w') as fp:
        pickle.dump([res1,myjckt.Frameouts2.mass,myjckt.Mpiles,myjckt.Tower.Twrouts.mass,myjckt.Piles.Pileinputs.Dpile,myjckt.Piles.Pileinputs.tpile,myjckt.Piles.Pileinputs.Lp,\
                    myjckt.Xbraces.Xbrcouts.LLURObj.D,myjckt.Xbraces.Xbrcouts.LLURObj.t,myjckt.Mudbraces.Mbrcouts.brcObj.D,myjckt.Mudbraces.Mbrcouts.brcObj,\
                    myjckt.JcktGeoIn.batter,myjckt.leginputs.Dleg,myjckt.leginputs.tleg,myjckt.TPinputs.Dgir,myjckt.TPinputs.tgir,myjckt.TP.TPouts.TPlumpedMass,\
                    myjckt.Tower.Twrins.Db,myjckt.Tower.Twrins.DTRb,myjckt.Tower.Twrins.Dt,myjckt.Tower.Twrins.DTRt,myjckt.Tower.Twrins.Htwr2frac,\
                    myjckt.Frameouts2.Freqs,myjckt.tower_utilization.GLUtil,myjckt.tower_utilization.EUshUtil,myjckt.wbase,myjckt.Xbraces.bay_hs,myjckt.Xbraces.bay_bs],\
                    fp)
        #pickle.dump(myjckt,fp) #does not work
    return res1,myjckt
    #5. store data for future use import cPickle as pickle;

    #6. recall data
   # with open(filename,'rb') as fp:
   #     res1=pickle.load(fp)
   #     myjckt=pickle.load(fp)


def retrievedata(filename="C:\PROJECTS\OFFSHORE_WIND\UH_REACT\PYTHON_OPT\JCKdataLumpedMass.p"):
    import cPickle as pickle;
    with open(filename,'r') as fp:
        res1=pickle.load(fp)
        myjckt=pickle.load(fp)
    return res1,myjckt


if __name__ == '__main__':
    #This is how you call this function
    main()
