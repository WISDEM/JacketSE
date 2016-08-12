#-------------------------------------------------------------------------------
# Name:        MyTowerInputs.py
# Purpose:     This module is  a template to set up a Tower Assembly with basic input & optimization parameters (to be set (harwired here))
#              THIS IS JUST A TEMPLATE:
#              in order to modify this: 1. Copy this file and rename as 'yourinputfile.py'
#                                       2. Edit all the tower inputs till the line that prohibits to further edit
#                                       3. Launch TowerSEOpt_Py&MDAOopt.py pointing to 'yourinputfile.py'
# Author:      rdamiani
#
# Created:     11/24/2014
# Copyright:   (c) rdamiani 2014
# Licence:     Apache (2014)
#-------------------------------------------------------------------------------
import numpy as np
from scipy.interpolate import interp1d

from openmdao.main.api import set_as_top
from openmdao.main.datatypes.api import Float

from commonse.environment import PowerWind, TowerSoil, LinearWaves
from commonse.Material import Material

from towerse.tower import TowerMonopileSE,TowerWithFrame3DD, TowerWithpBEAM

from jacketse.VarTrees import Frame3DDaux

from PlotTower import main as PlotTower  #COMMENT THIS ONE OUT FOR PEREGRINE"S SAKE


def main(): #\
    """Function to Instantiate a TowerSE Assembly: \n
       INPUTS \n
             All hardwired, so edit the quantities below all the way to the line "#________________ DO NOT MODIFY THE FOLLOWING ________________#" \n
             -See TowerSEOpt_Py&MDAOopt.py for more information. \n
       OUTPUTS \n
             mytwr -tower assembly instance \n\n

             Optimization parameters:    \n\n

             f0          -float, target frequency [Hz]
             f0epsilon   -float,  f0*(1+f0epsilon) will not be exceeded \n
             guesses     -Float(n), guesses for all design variables check out DesVar class. \n
             bounds      -Float(n,2), bounds for all design variables check out DesVar class. \n\n
             SAMPLE CALLS: \n
             1.OPTIMIZATION: python towerOpt_ExtCobyla.py C:\RRD\PYTHON\WISDEM\towerSE\src\towerse\MyTowerInputs.py \n
             2.OPTIMIZATION: python TowerSEOpt_Py&MDAOopt.py C:\RRD\PYTHON\WISDEM\towerSE\src\towerse\MytowerInputs.py True \n
             3.BUILD Tower: python >>> mytwr=C:\RRD\PYTHON\WISDEM\TowerSE\src\towerse\MyTowerInputs.py \n
        """

    # I need to put this at the top as it is not set up nicely as jacket: do not modify next line and go to inputs below
    mytwr=set_as_top(TowerMonopileSE())

    # __________Frame3DD or PBeam___________#
    mytwr.replace('tower1', TowerWithpBEAM())
    mytwr.replace('tower2', TowerWithpBEAM())
    #mytwr.replace('tower1', TowerWithFrame3DD())
    #mytwr.replace('tower2', TowerWithFrame3DD())

    # __________Material___________#
    mytwr.material=Material(matname='heavysteel',E=2.1e11,G=8.08e10, rho=8500.)

    # __________Geometry/Positioning___________#
    mytwr.sea_depth = 20.0
    mytwr.tower_length = 87.60
    mytwr.tower_to_shaft = 2.0
    mytwr.deck_height = 15.0
    mytwr.monopile_extension = 5.0
    mytwr.d_monopile = 6.0 # positioning for several variables now depends on monopile diameter
    mytwr.t_monopile = 0.06
    mytwr.t_jacket = 0.05
    mytwr.d_tower_base = 6.0
    mytwr.d_tower_top = 3.87
    mytwr.t_tower_base = 1.3*0.027
    mytwr.t_tower_top = 1.3*0.019

    mytwr.yaw = 0.0
    mytwr.tilt = 5.0

    # full geometry now specified by positioning component for monopile application

    # __________Environment___________#

    mytwr.replace('wind1', PowerWind())
    mytwr.replace('wind2', PowerWind())

    # wind
    mytwr.wind1.shearExp = 0.2
    mytwr.wind2.shearExp = 0.2

    # waves
    mytwr.replace('wave1', LinearWaves())
    mytwr.replace('wave2', LinearWaves())

    mytwr.wave1.Uc = 0.0
    mytwr.wave1.hs = 8.0 * 1.86
    mytwr.wave1.T = 10.0
    mytwr.wave1.g = 9.81
    mytwr.wave1.betaWave = 0.0

    mytwr.wave2.Uc = 0.0
    mytwr.wave2.hs = 8.0 * 1.86
    mytwr.wave2.T = 10.0
    mytwr.wave2.g = 9.81
    mytwr.wave2.betaWave = 0.0

    # __________Soil___________#
    mytwr.replace('soil', TowerSoil())

    mytwr.soil.rigid = 6*[True]

    #________RNA mass Properties_________#
    mytwr.top_m = 350000. #Float(iotype='in', units='m', desc='RNA (tower top) mass')
    mytwr.top_I = np.array([114930678.00, 22035403.00,  18759742.50, 0.00,  503710.47,  0.00]) #Array(iotype='in', units='kg*m**2', desc='mass moments of inertia. order: (xx, yy, zz, xy, xz, yz)')
    mytwr.top_cm = np.array([-1.13, 0.00, 0.51]) #Array(iotype='in', units='m', desc='RNA center of mass')

    #_______________Loads________________#

    # max Thrust case
    mytwr.wind_Uref1 = 11.73732
    mytwr.top1_F = np.array([1284744.196, 0.0,  -112400.5527]) #Array(iotype='in', units='N', desc='Aerodynamic forces')
    mytwr.top1_M = np.array([3963732.762, 896380.8464,  -346781.6819]) #Array(iotype='in', units='N*m', desc='Aerodynamic moments')

    # max wind speed case
    mytwr.wind_Uref2 = 70.0
    mytwr.top2_F = np.array([188038.8045, 0,  -16451.2637]) #Array(iotype='in', units='N', desc='Aerodynamic forces') # with all blades feathered
    mytwr.top2_M = np.array([0.0, 131196.8431,  0.0]) #Array(iotype='in', units='N*m', desc='Aerodynamic moments')


    # fatigue
    mytwr.z_DEL = np.array([0.000]) #, 1.327, 3.982, 6.636, 9.291, 11.945, 14.600, 17.255, 19.909, 22.564, 25.218, 27.873, 30.527, 33.182, 35.836, 38.491, 41.145, 43.800, 46.455, 49.109, 51.764, 54.418, 57.073, 59.727, 62.382, 65.036, 67.691, 70.345, 73.000, 75.655, 78.309, 80.964, 83.618, 86.273, 87.600])
    mytwr.M_DEL = 1e3*np.array([1.0]) #8.2940E+003, 8.1518E+003, 7.8831E+003, 7.6099E+003, 7.3359E+003, 7.0577E+003, 6.7821E+003, 6.5119E+003, 6.2391E+003, 5.9707E+003, 5.7070E+003, 5.4500E+003, 5.2015E+003, 4.9588E+003, 4.7202E+003, 4.4884E+003, 4.2577E+003, 4.0246E+003, 3.7942E+003, 3.5664E+003, 3.3406E+003, 3.1184E+003, 2.8977E+003, 2.6811E+003, 2.4719E+003, 2.2663E+003, 2.0673E+003, 1.8769E+003, 1.7017E+003, 1.5479E+003, 1.4207E+003, 1.3304E+003, 1.2780E+003, 1.2673E+003, 1.2761E+003])
    mytwr.gamma_fatigue = 1.35*1.3*1.0
    mytwr.life = 20.0
    mytwr.m_SN = 4

    # Frame3DD parameters
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


    #_____ Safety Factors______#
    mytwr.gamma_f = 1.35
    mytwr.gamma_m = 1.3
    mytwr.gamma_n = 1.0
    mytwr.gamma_b = 1.1

    #______________________________________________#
    #______________________________________________#

    # OTHER AUXILIARY CONSTRAINTS AND TARGETS FOR OPTIMIZATION # !!! NOT USED IN THIS FILE!
    #______________________________________________#
    #______________________________________________#

    # _______Geometric constraints__________#

    mytwr.min_taper = 0.4

    DTRsdiff = True  #Set whether or not DTRt=DTRb

    #______Set target frequency [Hz] and f0epsilon, i.e. fmax=(1+f0eps)*f0_________#
    f0=0.28
    f0epsilon=0.1

    #________Set Optimization Bounds and guesses for the various variables_________#
    #          x=  [      Db,   DTRb   Dt,   DTRt   Htwr2fac  ]
    MnCnst = np.array([  5.,   120.,  3.,   120.,     0.05  ])
    MxCnst = np.array([  7.,   200.,  4.,   200.,     0.25])
    guesses= np.array([  6.,   140.,  3.5,  150.,     0.2 ])


   #_____________________________________________________________#
   #________________ DO NOT MODIFY THE FOLLOWING ________________#
   #_____________________________________________________________#

    mytwr.min_d_to_t = np.min(MnCnst[[1,3]])
    bounds=np.vstack((MnCnst,MxCnst))
    desvarmeans=np.mean(bounds,1)

    mytwr.FrameAuxIns=FrameAuxIns

    return mytwr,f0,f0epsilon,DTRsdiff,guesses,bounds.T

def run_tower(run_results=False):

    mytwr= main()[0]
    #--- RUN JACKET ---#
    mytwr.run()
    # ---------------- #

    #_____________________________________#
    #Now show results of modal analysis
    print('First two Freqs.= {:5.4f} and {:5.4f} Hz'.format(mytwr.tower1.f1,mytwr.tower1.f2))
    #print component masses

    print('tower mass [kg] = {:6.0f}'.format(mytwr.mass))

    #print tower top displacement
    print('Tower Top Displacement (load cases 1 and 2) in Global Coordinate System [m] ={:5.4f} & {:5.4f}'.format(mytwr.top_deflection1,mytwr.top_deflection2))
    #print max Utilizations
    print('MAX GL buckling = {:5.4f}'.format(np.max((mytwr.buckling1,mytwr.buckling2))))
    print('MAX Shell buckling = {:5.4f}'.format(np.max((mytwr.shellBuckling1,mytwr.shellBuckling2))))


    from PlotTower import main as PlotTower  #COMMENT THIS ONE OUT FOR PEREGRINE"S SAKE
    PlotTower(mytwr,util=True)
    print mytwr.jp.z
    print mytwr.jp.towerHeight
    print (mytwr.jp.z*mytwr.jp.towerHeight)-(mytwr.jp.d_monopile*1.5-mytwr.jp.monopile_extension)
    print mytwr.z_nodes

    if run_results:
        sea_depth = np.array([23.50, 12.80, 40.80, 9.70, 65.80, 50.60])
        hmax_50yr = np.array([18.33, 9.98, 17.63, 7.57, 22.60, 14.19])
        tp_50yr = np.array([13.34, 10.91, 12.48, 13.11, 14.13, 11.20])
        deck_height = np.array([13.20, 13.00, 16.00, 7.00, 18.00, 12.70])
        d_monopile = np.array([6.90, 6.24, 7.72, 6.13, 9.81, 7.96])
        t_monopile = np.array([0.0369, 0.0292, 0.0454, 0.0279, 0.0645, 0.0477])
        t_jacket = np.array([0.0207, 0.0228, 0.0190, 0.0232, 0.0172, 0.0186])
        d_tower_base = np.array([5.55, 5.57, 5.66, 5.72, 4.87, 5.522])
        d_tower_top = np.array([4.05, 4.07, 4.08, 4.22, 3.92, 4.05])
        t_tower_base = np.array([0.0260, 0.0259, 0.0246, 0.0266, 0.0282, 0.0263])
        t_tower_top = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
        #t_tower_top = 2 * t_tower_top
        t_monopile = d_monopile / 90.0

        for i in range(len(sea_depth)):
            mytwr.sea_depth = sea_depth[i]
            mytwr.deck_height = deck_height[i]
            mytwr.wave1.hs = hmax_50yr[i]
            mytwr.wave1.T = tp_50yr[i]
            mytwr.wave2.hs = hmax_50yr[i]
            mytwr.wave2.T = tp_50yr[i]

            #reset to initial conditions for geometry
            mytwr.d_monopile = d_monopile[i]
            mytwr.t_monopile = t_monopile[i]
            mytwr.t_jacket = t_jacket[i]
            mytwr.d_tower_base = d_tower_base[i]
            mytwr.d_tower_top = d_tower_top[i]
            mytwr.t_tower_base = t_tower_base[i]
            mytwr.t_tower_top = t_tower_top[i]

            mytwr.run()

            #Now show results of modal analysis
            print('First two Freqs.= {:5.4f} and {:5.4f} Hz'.format(mytwr.tower1.f1,mytwr.tower1.f2))
            #print component masses

            print('tower mass [kg] = {:6.0f}'.format(mytwr.mass))

            #print tower top displacement
            print('Tower Top Displacement (load cases 1 and 2) in Global Coordinate System [m] ={:5.4f} & {:5.4f}'.format(mytwr.top_deflection1,mytwr.top_deflection2))
            #print max Utilizations
            print('MAX GL buckling = {:5.4f}'.format(np.max((mytwr.buckling1,mytwr.buckling2))))
            print('MAX Shell buckling = {:5.4f}'.format(np.max((mytwr.shellBuckling1,mytwr.shellBuckling2))))

            PlotTower(mytwr,util=True)
            print mytwr.jp.z
            print mytwr.jp.towerHeight
            print (mytwr.jp.z*mytwr.jp.towerHeight)-(mytwr.jp.d_monopile*1.5-mytwr.jp.monopile_extension)
            print mytwr.z_nodes



def test_grads():

    mytwr= main()[0]

    mytwr.run()

    for i in range(5):

         print 6-i

         mass = mytwr.mass
         mytwr.d_monopile += (10**(-(6-i)))
         mytwr.run()
         print mytwr.d_monopile, mytwr.mass
         print mytwr.mass/mass

         mass = mytwr.mass
         mytwr.t_monopile += (10**-(6-i))
         mytwr.run()
         print mytwr.t_monopile, mytwr.mass
         print mytwr.mass/mass

         mass = mytwr.mass
         mytwr.t_jacket += (10**-(6-i))
         mytwr.run()
         print mytwr.t_jacket, mytwr.mass
         print mytwr.mass/mass

         mass = mytwr.mass
         mytwr.d_tower_base += (10**-(6-i))
         mytwr.run()
         print mytwr.d_tower_base, mytwr.mass
         print mytwr.mass/mass

         mass = mytwr.mass
         mytwr.t_tower_base += (10**(-(6-i)))
         mytwr.run()
         print mytwr.t_tower_base, mytwr.mass
         print mytwr.mass/mass

         mass = mytwr.mass
         mytwr.d_tower_top += (10**(-(6-i)))
         mytwr.run()
         print mytwr.d_tower_top, mytwr.mass
         print mytwr.mass/mass

         mass = mytwr.mass
         mytwr.t_tower_top += (10**(-(6-i)))
         mytwr.run()
         print mytwr.t_tower_top, mytwr.mass
         print mytwr.mass/mass


def opt_tower(sea_depth=np.array([]), hmax_50yr=None, tp_50yr=None, deck_height=None, d_monopile=None, t_monopile=None, t_jacket=None, \
              d_tower_base=None, d_tower_top=None, t_tower_base=None, t_tower_top=None):

    mytwr= main()[0]

    MDAOswitch = 'md_pysnopt'

    # optimization
    import pyOpt
    from pyopt_driver.pyopt_driver import pyOptDriver
    from openmdao.lib.casehandlers.api import DumpCaseRecorder
    from openmdao.lib.drivers.api import COBYLAdriver

    if MDAOswitch == 'md_cobyla':
        mytwr.replace('driver', COBYLAdriver())

        mytwr.driver.rhobeg=0.01
        mytwr.driver.rhoend=1.e-3
        mytwr.driver.maxfun=2000
        mytwr.driver.iprint=1
    else:
        mytwr.replace('driver', pyOptDriver())
        if  MDAOswitch== 'md_pysnopt':
            mytwr.driver.optimizer = 'SNOPT'
            mytwr.driver.pyopt_diff = True
            #mytwr.driver.sens_step = 1e-8
            mytwr.driver.options = {'Major feasibility tolerance': 1e-4,
                                 'Minor feasibility tolerance': 1e-4,
                                 'Major optimality tolerance': 1e-4,
                                 'Function precision': 1e-6}
        elif MDAOswitch== 'md_pyCobyla':
            mytwr.driver.optimizer = 'COBYLA'
            mytwr.driver.options = {'RHOEND':1.e-2,'RHOEND':1.e-3,'MAXFUN':2000,'IPRINT':1}
        else:
            sys.exit('Error: MDAOswitch must be set to ''pyCobyla'' or ''pySNOPT'' or ''md_cobyla'' or ''md_pySNOPT'' or ''md_pyCobyla''!!!')
    # ----------------------

    # --- Objective ---
    mytwr.driver.add_objective('(tower1.mass) / 500.e3')
    # ----------------------

    #diameter and thickness variables
    mytwr.driver.add_parameter('d_monopile', low=3.0, high=15.0) # this is OD monopile
    mytwr.driver.add_parameter('t_monopile', low=0.01, high=0.1) # this is t of monopile
    mytwr.driver.add_parameter('t_jacket', low=0.01, high = 0.1) # this is t of jacket # d jacket fixed by monopile + t_jacket + t_grout

    mytwr.driver.add_parameter('d_tower_base', low= 3.0,   high= 15.0)   #This is OD at the base
    mytwr.driver.add_parameter('t_tower_base', low= 0.01,  high=0.3)   #This is t at the base

    mytwr.driver.add_parameter('d_tower_top',    low= 3.87,   high= 8.0)   #This is OD at the top # OD at top should be fixed
    mytwr.driver.add_parameter('t_tower_top', low= 0.01,  high=0.1)   #This is t at top

    # node positioning variables
    #mytwr.driver.add_parameter('jp.z[5]', low= 0.0,  high= 0.99) # this is H2 (can't move this - will have to add extra variable to make that happen)

    #--- Constraints ---#
    # frequency
    mytwr.driver.add_constraint('tower1.f1 >= 0.20') # from 0.28 to 0.26
    #mytwr.driver.add_constraint('tower1.f1 <= 0.30')

    # utilization
    mytwr.driver.add_constraint('tower1.stress <= 1.0')
    mytwr.driver.add_constraint('tower2.stress <= 1.0')
    mytwr.driver.add_constraint('tower1.buckling <= 1.0')
    mytwr.driver.add_constraint('tower2.buckling <= 1.0')
    mytwr.driver.add_constraint('tower1.shellBuckling <= 1.0')
    mytwr.driver.add_constraint('tower2.shellBuckling <= 1.0')
    mytwr.driver.add_constraint('tower1.damage <= 1.0')
    #mytwr.driver.add_constraint('gc.weldability <= 0.0') # this goes for entire structure which we don't want
    #mytwr.driver.add_constraint('gc.manufacturability <= 0.0') # just going to use explicit constraints below

    # uniformity of diameter and thickness of tower
    mytwr.driver.add_constraint('d_tower_top <= d_tower_base')
    mytwr.driver.add_constraint('d_tower_base <= d_monopile')
    #mytwr.driver.add_constraint('d_monopile - d_tower_base <= 0.4')
    mytwr.driver.add_constraint('t_tower_top <= t_tower_base')

    # manufacturing and installation
    mytwr.driver.add_constraint('d_tower_top/d_tower_base >= 0.40') # manufacturability - already covered; taper ratio for manufacturing
    mytwr.driver.add_constraint('d_tower_base/t_tower_base >= 120.0') # weldibility - already covered; tower DTR for rolling operation
    mytwr.driver.add_constraint('d_tower_top/t_tower_top >= 120.0')
    #mytwr.driver.add_constraint('d_monopile/t_monopile <= 100.0') # pile DTR for installation
    #mytwr.driver.add_constraint('d_tower_base/t_tower_base <= 200.0') # tower DTR for welding (not covered)
    #mytwr.driver.add_constraint('d_tower_top/t_tower_top <= 200.0')

    ## mytwr.driver.add_constraint('Lp0rat >= 0.') # was for embedment - not in model
    # ----------------------

    # set optimization variables
    if sea_depth.size:

        for i in range(len(sea_depth)):
            mytwr.sea_depth = sea_depth[i]
            mytwr.deck_height = deck_height[i]
            mytwr.wave1.hs = hmax_50yr[i]
            mytwr.wave1.T = tp_50yr[i]
            mytwr.wave2.hs = hmax_50yr[i]
            mytwr.wave2.T = tp_50yr[i]

            #reset to initial conditions for geometry
            mytwr.d_monopile = d_monopile[i]
            mytwr.t_monopile = t_monopile[i]
            mytwr.t_jacket = t_jacket[i]
            mytwr.d_tower_base = d_tower_base[i]
            mytwr.d_tower_top = d_tower_top[i]
            mytwr.t_tower_base = t_tower_base[i]
            mytwr.t_tower_top = t_tower_top[i]

            # print IC
            print 'Initial condition\n'
            print 'sea_depth=', mytwr.sea_depth
            print 'deck_height=', mytwr.deck_height
            print 'significant wave height=', mytwr.wave1.hs
            print 'wave period=', mytwr.wave1.T
            print 'd_monopile=', mytwr.d_monopile
            print 't_monopile=', mytwr.t_monopile
            print 't_jacket=', mytwr.t_jacket
            print 'd_tower_base=', mytwr.d_tower_base
            print 'd_tower_top=', mytwr.d_tower_top
            print 't_tower_base=', mytwr.t_tower_base
            print 't_tower_top=', mytwr.t_tower_top

            mytwr.driver.options.update({'Summary file': 'Case_' + str(i) + '_summary.out',
                                 'Print file': 'Case_' + str(i) + '_print.out'})

            if sea_depth[i] > 60.0: # deepwater case
                mytwr.driver.options['Major optimality tolerance']=1e-3
            else:
                mytwr.driver.options['Major optimality tolerance']=1e-3

            #RUN
            import time
            tt = time.time()
            mytwr.run()

            print "Final conditions"
            print "Minimum found at Db=%f, Dt=%f, tb=%f, tt=%f; mass= (%f)" % (mytwr.tower1.d[0],mytwr.tower1.d[-1],mytwr.tower1.t[0],mytwr.tower1.t[-1],mytwr.tower1.mass)
            print "Minimum found at z1=%f, D1=%f, t1=%f" % (mytwr.tower1.z[1],mytwr.tower1.d[1],mytwr.tower1.t[1])
            print "Minimum found at DTRb DTRt(%f, %f)" % (mytwr.tower1.d[0]/mytwr.tower1.t[0],mytwr.tower1.d[-1]/mytwr.tower1.t[-1])

            print "\n"
            # print FC
            print 'd_monopile=', mytwr.d_monopile
            print 't_monopile=', mytwr.t_monopile
            print 't_jacket=', mytwr.t_jacket
            print 'd_tower_base=', mytwr.d_tower_base
            print 'd_tower_top=', mytwr.d_tower_top
            print 't_tower_base=', mytwr.t_tower_base
            print 't_tower_top=', mytwr.t_tower_top
            print 'zs=', mytwr.tower1.z
            print 'ds=', mytwr.tower1.d
            print 'ts=', mytwr.tower1.t

            print "\n"
            print "Minimum found at Freq %f"  % (mytwr.tower1.f1)
            print "Minimum found at GLutil EUutil %f %f"  % (np.max(np.vstack((mytwr.tower1.buckling,mytwr.tower2.buckling))),np.max(np.vstack((mytwr.tower1.shellBuckling,mytwr.tower2.shellBuckling))))
            print "Minimum found at GLutil 1 and 2"  , mytwr.tower1.buckling,mytwr.tower2.buckling
            print "Minimum found at EUutil 1 and 2"  , mytwr.tower1.shellBuckling,mytwr.tower2.shellBuckling

            print "Elapsed time: ", time.time()-tt, "seconds\n"

            #Plot
            PlotTower(mytwr,util=True,savefileroot=str(i))

    else:

        #RUN
        import time
        tt = time.time()
        mytwr.run()

        print "\n"
        print "Minimum found at Db=%f, Dt=%f, tb=%f, tt=%f; mass= (%f)" % (mytwr.tower1.d[0],mytwr.tower1.d[-1],mytwr.tower1.t[0],mytwr.tower1.t[-1],mytwr.tower1.mass)
        print "Minimum found at z1=%f, D1=%f, t1=%f" % (mytwr.tower1.z[1],mytwr.tower1.d[1],mytwr.tower1.t[1])
        print "Minimum found at DTRb DTRt(%f, %f)" % (mytwr.tower1.d[0]/mytwr.tower1.t[0],mytwr.tower1.d[-1]/mytwr.tower1.t[-1])

        print "\n"
        print 'zs=', mytwr.tower1.z
        print 'ds=', mytwr.tower1.d
        print 'ts=', mytwr.tower1.t

        print "\n"
        print "Minimum found at Freq %f"  % (mytwr.tower1.f1)
        print "Minimum found at GLutil EUutil %f %f"  % (np.max(np.vstack((mytwr.tower1.buckling,mytwr.tower2.buckling))),np.max(np.vstack((mytwr.tower1.shellBuckling,mytwr.tower2.shellBuckling))))
        print "Minimum found at GLutil 1 and 2"  , mytwr.tower1.buckling,mytwr.tower2.buckling
        print "Minimum found at EUutil 1 and 2"  , mytwr.tower1.shellBuckling,mytwr.tower2.shellBuckling

        print "Elapsed time: ", time.time()-tt, "seconds"
        print "Execution count: ", mytwr.exec_count

        #Plot
        PlotTower(mytwr,util=True)

        print mytwr.jp.z
        print mytwr.jp.towerHeight
        print (mytwr.jp.z*mytwr.jp.towerHeight)-(mytwr.d_monopile*1.5-mytwr.monopile_extension)
        print mytwr.z_nodes


if __name__ == '__main__':

    run_tower() # run once
    #test_grads() # check gradients
    #opt_tower() # run single optimization for OC3
    #run_tower(True) # run through results

    # set optimization parameters (if not using default)
    sea_depth = np.array([23.50, 12.80, 40.80, 9.70, 65.80, 50.60])
    hmax_50yr = np.array([18.33, 9.98, 17.63, 7.57, 22.60, 14.19])
    tp_50yr = np.array([13.34, 10.91, 12.48, 13.11, 14.13, 11.20])
    deck_height = np.array([13.20, 13.00, 16.00, 7.00, 18.00, 12.70])
    d_monopile = np.array([6.0, 6.0, 6.0, 6.0, 6.0, 6.0])
    t_monopile = np.array([0.06, 0.06, 0.06, 0.06, 0.06, 0.06])
    t_jacket = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    d_tower_base = np.array([6.0, 6.0, 6.0, 6.0, 6.0, 6.0])
    d_tower_top = np.array([3.87, 3.87, 3.87, 3.87, 3.87, 3.87])
    t_tower_base = np.array([1.3*0.027, 1.3*0.027, 1.3*0.027, 1.3*0.027, 1.3*0.027, 1.3*0.027])
    t_tower_top = np.array([1.3*0.019, 1.3*0.019 , 1.3*0.019 ,1.3*0.019 ,1.3*0.019 ,1.3*0.019 ])

    # run optimization for series of sites
    #opt_tower(sea_depth, hmax_50yr, tp_50yr, deck_height, d_monopile, \
    #          t_monopile, t_jacket, d_tower_base, d_tower_top, t_tower_base, t_tower_top)