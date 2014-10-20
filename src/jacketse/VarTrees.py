#-------------------------------------------------------------------------------
# Name:        VarTrees.py
# Purpose:   Stores common VarTrees that are used by multiple modules (files).
#            This avoids circular calls.
#
# Author:      rdamiani
#
# Created:     17/02/2014
# Copyright:   (c) rdamiani 2014
# Licence:     Apache C 2014
#-------------------------------------------------------------------------------
from openmdao.main.api import VariableTree,  Component, Assembly, set_as_top
from openmdao.main.datatypes.api import Int, Float, Array, VarTree, Bool, Dict,Instance
from commonse.Tube import Tube
from commonse.RigidMember import RigidMember
import numpy as np
#______________________________________________________________________________#

class RNAprops(VariableTree):
    """Basic Inertial and Geometric Properties of RNA"""
    mass=Float(   units='kg',    desc='RNA mass')    #RNA mass [kg]
    Ixx=Float(    units='kg*m**2', desc='RNA IXX @tower top flange')    #RNA Ixx
    Iyy=Float(    units='kg*m**2', desc='RNA IYY @tower top flange')    #RNA Iyy
    Izz=Float(    units='kg*m**2', desc='RNA IZZ @tower top flange')    #RNA Izz
    Ixy=Float(0.0,units='kg*m**2', desc='RNA  @tower top flange ')    #RNA Ixy
    Ixz=Float(0.0,units='kg*m**2', desc='RNA  @tower top flange')    #RNA Ixz
    Iyz=Float(0.0,units='kg*m**2', desc='RNA  @tower top flange')    #RNA Iyz

    CMxoff=Float(0.,   units='m', desc='RNA CM x offset from Tower Top Flange')       # RNA CMxoff [m]
    CMyoff=Float(0.,   units='m', desc='RNA CM y offset from Tower Top Flange')       # RNA CMyoff [m]
    CMzoff=Float(0.,   units='m', desc='RNA CM z offset from Tower Top Flange')       # RNA CMzoff [m]

    Thxoff=Float(0.,    units='m', desc='Rotor Hub Center x offset from Tower Top Flange')       # Thrust point of application [m]
    Thyoff=Float(0.,    units='m', desc='Rotor Hub Center y offset from Tower Top Flange')       #
    Thzoff=Float(0.,    units='m', desc='Rotor Hub Center z offset from Tower Top Flange')       #

    rna_weightM = Bool(True, units=None, desc='flag to consider or not the RNA weight effect on Moment')

    yawangle=Float(0.,   units='deg', desc='YAW angle CCW RH rule, to account for possible nacelle weight contribution.')       # RNA Yaw angle [deg]



class JcktGeoOutputs(VariableTree):
    """Node Coordinates and Member Connectivity."""
    nNodesJckt = Int(  units=None,                 desc='Total Number of nodes in the Jacket, no Tower')
    nodes =      Array(units='m',  dtype=np.float, desc='Node''s coordinates in the Substructure')
    radii =      Array(units='m',  dtype=np.float, desc='Node''s Radii')
    Reacts   =   Array(units=None, dtype=int,      desc='Node IDs for reactions + Fixity values (1/0)')  #Fixed for the time being with all fixity (6 dofs per node fixed)
    nmems  =     Int(  units=None,                 desc='Total Number of Members in the Jacket, no Tower')
    mems   =     Array(units=None, dtype=int,      desc='Member Connectivity Node i Node j for every member (i.e.,element)')
    props  =     Array(            dtype=np.float, desc='Jacket Member''s xsectional and material properties')
    XnsfM  =     Array(            dtype=np.float, desc='Jacket Member''s Global to Local DIrection Cosine Matrices [3,3,nelems]')
    TubeObjs=    Instance(Klass=Tube,                  desc='Object of Class Tube for Jacket Elements. Tube Objects, one per element')
    jnt_masses=  Array(            dtype=np.float, desc='Jacket Concentrated Masses at Joints: (joint no, mass, IXX, IYY, IZZ)')
    jnt_masses_yaw=  Array(            dtype=np.float, desc='Jacket Concentrated Masses at Joints: (joint no, mass, IXX, IYY, IZZ) accounting for RNA yaw')

class TwrGeoOutputs(VariableTree):
    """Basic Geometric Outputs needed to build Tower of Jacket."""
    TwrObj   =  Instance(Klass=Tube,                             desc='Object of Class Tube for Tower portion of Jacket')
    Twr2RNAObj =Instance(Klass=RigidMember,                      desc='Object of Class RigidMember for Rigid Tower portion')
    joints   = Array(np.array([]),units='m', dtype=np.float, desc='Pile Joint (with legs) Coordinates (3,nlegs=nfaces)')
    nodes    = Array(np.array([]),units='m', dtype=np.float, desc='Tower ALL Nodes'' Coordinates (3,nNodes)')
    nNodes   = Int(0,           units=None,                  desc='Number of Tower Nodes INCLUDING joint at TP')
    mass =     Float(           units='kg',                  desc='Tower Mass')
    TopMass  = Array(np.zeros([10]),         dtype=np.float, desc='Tower Top mass, Ixx, Iyy, Izz,Ixy,Ixz,Iyz,CMxoff,CMyoff,CMzoff from RNA properties in input')
    TopMass_yaw  = Array(np.zeros([10]),        dtype=np.float, desc='Tower Top mass, Ixx, Iyy, Izz,Ixy,Ixz,Iyz,CMxoff,CMyoff,CMzoff from RNA properties in input including yaw angle w.r.t. Global XYZ')

    HH       = Float(            units='m', desc='Hub-Height')
    Htwr     = Float(            units='m', desc='Tower Length')
    Htwr2     = Float(            units='m', desc='Tower at constant cross section Length')
    Thoff_yaw= Array(np.zeros([3]),units='m',  dtype=np.float, desc='Tower-top to Hub-center vector in yawed coordinate system (TO BE MOVED SOMEWHERE ELSE at one point)')
    rna_yawedcm  = Array(np.zeros([3]),        dtype=np.float, desc='Tower Top mass CMxoff,CMyoff,CMzoff in yawed coordinate system(TO BE MOVED SOMEWHERE ELSE at one point)')

def main():
    pass

if __name__ == '__main__':
    main()
