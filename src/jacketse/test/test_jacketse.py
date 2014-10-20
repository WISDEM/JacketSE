#!/usr/bin/env python
# encoding: utf-8
"""
test_jacketse.py

Created by RRD on 2013-12-20.
Copyright (c) NREL. All rights reserved.
"""

import unittest
import numpy as np
from commonse.utilities import check_gradient, check_gradient_unit_test
from jacketse.jacket import JcktGeoInputs,TwrGeoInputs,TPlumpMass,TPGeoInputs,XBrcGeoInputs, MudBrcGeoInputs, HBrcGeoInputs, MatInputs, LegGeoInputs, SoilGeoInputs,PileGeoInputs,Frame3DDaux,JacketAsmly
from jacketse.VarTrees import RNAprops
from jacketse.loads import WaterInputs, WindInputs
from openmdao.main.api import  Assembly, set_as_top

# Jacket Components
class TestJacket(unittest.TestCase):

    def setUp(self):

        # simple test of module

        Jcktins=JcktGeoInputs()
        Jcktins.nlegs =4
        Jcktins.nbays =4
        Jcktins.batter=15.
        Jcktins.dck_botz =16.
        Jcktins.VPFlag= True    #vertical pile T/F;  to enable piles in frame3DD set pileinputs.ndiv>0
        Jcktins.clamped= False    #whether or not the bottom of the structure is rigidly connected. Use False when equivalent spring constants are being used.

        self.jacket=set_as_top(JacketAsmly())

        self.jacket.JcktGeoIn=Jcktins

        Soilinputs=SoilGeoInputs()
        Soilinputs.zbots   =-np.array([3.,5.,7.,15.,30.,50.])
        Soilinputs.gammas  =np.array([10000.,10000.,10000.,10000.,10000.,10000.])
        Soilinputs.cus     =np.array([60000.,60000.,60000.,60000.,60000.,60000.])
        Soilinputs.phis    =np.array([26.,26.,26.,26.,26.,26])#np.array([36.,33.,26.,37.,35.,37.5])#np.array([36.,33.,26.,37.,35.,37.5])
        Soilinputs.delta   =25.
        Soilinputs.sndflg   =True
        Soilinputs.PenderSwtch   =False #True
        Soilinputs.SoilSF   =1.

        self.jacket.Soilinputs=Soilinputs

     #Water and wind inputs
        Waterinputs=WaterInputs()
        Waterinputs.wdepth   =30.
        Waterinputs.wlevel   =30. #Distance from bottom of structure to surface  THIS, I believe is no longer needed as piles may be negative in z, to check and remove in case
        Waterinputs.T=12.  #Wave Period
        Waterinputs.HW=10. #Wave Height
        Windinputs=WindInputs()
        Windinputs.HH=100. #hub height set here

        self.jacket.Waterinputs=Waterinputs
        self.jacket.Windinputs=Windinputs

    #RNA loads                          Fx-z,         Mxx-zz
        self.jacket.RNA_F=np.array([1000.e3,0.,0.,0.,0.,0.])


    #Legs Data
        legmatin=MatInputs()
        legmatin.matname=(['steel'])

        Dleg=np.array([2.0]).repeat(Jcktins.nbays+1)
        tleg=np.array([0.0254]).repeat(Dleg.size)

        leginputs=LegGeoInputs()
        leginputs.legZbot   = 1.0

        leginputs.legmatins=legmatin
        leginputs.Dleg=Dleg
        leginputs.tleg=tleg

        self.jacket.leginputs=leginputs

    #The following is a passthrough variables
        self.jacket.legbot_stmphin =1.5  #Distance from bottom of leg to second joint along z; must be>0

    #Xbrc data
        Xbrcmatin=MatInputs()
        Xbrcmatin.matname=np.array(['steel']).repeat(Jcktins.nbays)

        Xbrcinputs=XBrcGeoInputs()
        Xbrcinputs.ndiv=2
        Xbrcinputs.Xbrcmatins=Xbrcmatin
        Xbrcinputs.precalc=True   #This can be set to true if we want Xbraces to be precalculated in D and t, in which case the above set Dbrc and tbrc would be overwritten
        self.jacket.Xbrcinputs=Xbrcinputs

    #Mbrc data
        Mbrcmatin=MatInputs()
        Mbrcmatin.matname=np.array(['steel'])

        Mbrcinputs=MudBrcGeoInputs()
        Mbrcinputs.ndiv=2
        Mbrcinputs.Mbrcmatins=Mbrcmatin
        Mbrcinputs.precalc=True   #This can be set to true if we want Mudbrace to be precalculated in D and t, in which case the above set Dbrc_mud and tbrc_mud would be overwritten

        self.jacket.Mbrcinputs=Mbrcinputs

    #Hbrc data
        Hbrcmatin=MatInputs()
        Hbrcmatin.matname=np.array(['steel'])
        Hbrcinputs=HBrcGeoInputs()
        Hbrcinputs.ndiv=0#2
        Hbrcinputs.Hbrcmatins=Hbrcmatin
        Hbrcinputs.precalc=True   #This can be set to true if we want Hbrace to be set=Xbrace top D and t, in which case the above set Dbrch and tbrch would be overwritten

        self.jacket.Hbrcinputs=Hbrcinputs

    #TP data
        TPlumpinputs=TPlumpMass()
        TPlumpinputs.mass=300.e3 #[kg]
        self.jacket.TPlumpinputs=TPlumpinputs

        TPstmpsmatin=MatInputs()
        TPbrcmatin=MatInputs()
        TPstemmatin=MatInputs()
        TPbrcmatin.matname=np.array(['steel'])
        TPstemmatin.matname=np.array(['steel']).repeat(2)

        TPinputs=TPGeoInputs()
        TPinputs.TPbrcmatins=TPbrcmatin
        TPinputs.TPstemmatins=TPstemmatin
        TPinputs.TPstmpmatins=TPstmpsmatin
        TPinputs.Dstrut=1.6
        TPinputs.hstump=0.0#1.0
        TPinputs.nstems=3
        TPinputs.Dstem=np.array([6.]).repeat(TPinputs.nstems)
        TPinputs.tstem=np.array([0.1,0.11,0.11])
        TPinputs.hstem=np.array([4.,3.,1.])

        self.jacket.TPinputs=TPinputs

        #Pile data
        Pilematin=MatInputs()
        Pilematin.matname=np.array(['steel'])
        Dpile=2.5#0.75 # 2.0
        tpile=0.01
        Lp=30. #45

        Pileinputs=PileGeoInputs()
        Pileinputs.Pilematins=Pilematin
        Pileinputs.ndiv=0 #3
        Pileinputs.Dpile=Dpile
        Pileinputs.tpile=tpile
        Pileinputs.AFflag=False
        Pileinputs.Lp=Lp #[m] Embedment length

        self.jacket.Pileinputs=Pileinputs

    #RNA and Tower data
        RNAins=RNAprops()
        RNAins.mass=3.*350e3
        RNAins.Ixx=86.579E+6
        RNAins.Iyy=53.530E+6
        RNAins.Izz=58.112E+6
        RNAins.CMzoff=2.34
        RNAins.yawangle=45.  #angle with respect to global X, CCW looking from above, wind from left

        self.jacket.TwrRigidTop=False       #False=Account for RNA via math rather than a physical rigidmember
        self.jacket.RNAinputs=RNAins

   # Tower data
        Twrmatin=MatInputs()
        Twrmatin.matname=np.array(['steel'])
        Db=5.9
        tb=0.05
        Dt=Db*0.55

        Twrinputs=TwrGeoInputs()
        Twrinputs.Twrmatins=Twrmatin
        #Twrinputs.Htwr=70.  #Trumped by HH
        Twrinputs.Htwr2frac=0.2   #fraction of tower height with constant x-section
        Twrinputs.ndiv=np.array([6,6])  #ndiv for uniform and tapered section
        Twrinputs.Db=Db
        Twrinputs.DTRb=Db/tb
        Twrinputs.Dt=Dt

        self.jacket.Twrinputs=Twrinputs

    # Frame3DD parameters
        FrameAuxIns=Frame3DDaux()
        FrameAuxIns.deltaz=5.
        FrameAuxIns.nModes = 6             # number of desired dynamic modes of vibration
        FrameAuxIns.tol = 1e-9             # mode shape tolerance

        self.jacket.FrameAuxIns=FrameAuxIns

    def test_functionality(self):

        self.jacket.run()

        self.maxDiff=None
        self.assertEqual(np.round(self.jacket.Tower.Twrouts.mass,2), 365295.51)
        self.assertEqual(np.round(self.jacket.Frameouts2.Freqs,4).tolist(), [ 0.2164,  0.2175])
        #tower utilizations
        GLUtil=[ 0.47028029 , 0.457376 ,   0.44461144 , 0.43198971 , 0.41951374 , 0.40718626, \
            0.39500987 , 0.40627729 , 0.41836944 , 0.4313658  , 0.44534452,  0.46037907,\
            0.47653772 , 0.49389097,  0.51253912 , 0.53268978 , 0.55490678]

        StressUtil=[ 0.36817623 , 0.35509984,  0.34216735,  0.32938195  ,0.31674666,  0.30426433,\
            0.29193768 , 0.29926594 , 0.30596969,  0.3115559  , 0.31528316 , 0.31603233 ,\
            0.31210165 , 0.30087775 , 0.27830493,  0.23806529 , 0.17081804]

        EUshUtil=[ 0.20414757 , 0.19072449 , 0.177874  ,  0.16558686 , 0.15385331 , 0.14266318,\
            0.13200589 , 0.13878835 , 0.14524854,  0.15092411 , 0.15510409 , 0.15671953   ,\
            0.15420726 , 0.1453824  , 0.12746531,  0.0977601  , 0.0566553 ]

        self.assertEqual(np.round(self.jacket.tower_utilization.GLUtil,8).tolist(), GLUtil)
        self.assertEqual(np.round(self.jacket.tower_utilization.EUshUtil,8).tolist(), EUshUtil)
        self.assertEqual(np.round(self.jacket.tower_utilization.StressUtil,8).tolist(), StressUtil)
        #jacket utilizations
        KjntUtil=[ 0.13671728,  0.16566477,  0.12377022,  0.10054838,  1.01250082,\
            0.16081369,  0.27087764,  0.1264021 ,  0.08689998,  0.6435632 ,\
            0.27978504,  0.11972022,  0.15779066,  0.09069759,  1.77188842,\
            0.16113317,  0.26925897,  0.12738972,  0.08688897,  0.64354391]
        XjntUtil= [ 0.34329747,  0.25255367,  0.10399988,  0.15595033,  0.55505745, \
            0.22226466,  0.10929611,  0.13713145,  0.55206134,  0.2202987 ,\
            0.10919169,  0.13712363,  0.34420208,  0.2533687 ,  0.10410395,\
            0.15607568]
        self.assertEqual(np.round(self.jacket.jacket_utilization.KjntUtil,8).tolist(), KjntUtil)
        self.assertEqual(np.round(self.jacket.jacket_utilization.XjntUtil,8).tolist(), XjntUtil)

        cb_util= np.round([  1.73768127e-01,  -7.22967924e-02,   1.35860433e-01,\
            -9.57247251e-02,   9.53173661e-02,  -7.80043161e-02,\
            8.94897729e-02,  -7.72294704e-02,   7.52457685e-02,\
            -2.54906904e-02,  -7.99669349e-02,   1.58896916e-01,\
            -1.03366134e-01,   1.51646700e-01,  -1.04647104e-01,\
             1.31858689e-01,  -1.16257763e-01,   1.45320552e-01,\
            -9.83991439e-02,   1.88393557e-01,  -2.69519879e-01,\
             4.26461255e-01,  -3.04270145e-01,   3.80237758e-01,\
            -3.23154595e-01,   3.61020739e-01,  -3.03221577e-01,\
             3.47096760e-01,  -2.93329683e-01,   3.71726789e-01,\
            -7.82249335e-02,   1.57213216e-01,  -1.01937694e-01,\
             1.50284056e-01,  -1.03802598e-01,   1.30945237e-01,\
            -1.16256588e-01,   1.45206055e-01,  -9.83518263e-02,\
             1.88354013e-01,   7.69435543e-02,  -6.59329206e-02,\
             7.13198227e-02,  -5.35850977e-02,   8.79547277e-02,\
            -6.75331055e-02,   7.53390036e-02,  -5.65404276e-02,\
             4.56538493e-02,  -3.79376551e-02,   4.66699321e-02,\
            -2.74206320e-02,   6.24250105e-02,  -3.85089156e-02,\
             5.06035012e-02,  -3.15357977e-02,   1.30745849e-02,\
            -5.25795777e-03,   1.13220327e-02,   7.60288806e-03,\
             2.66986219e-02,  -1.83558811e-03,   1.94180995e-02,\
             1.08648053e-02,   5.35962627e-02,  -3.86145241e-02,\
             4.39569559e-02,  -2.21412561e-02,   8.52359622e-02,\
            -2.46877177e-02,   6.45273100e-02,   1.86044779e-03,\
             4.08895920e-02,  -1.76296954e-02,   2.85594197e-02,\
            -1.55916368e-02,   8.36869927e-02,  -2.72162308e-02,\
             2.94521710e-02,   3.91864781e-03,   8.42922986e-02,\
            -4.68326324e-02,   5.64891278e-02,  -3.54539269e-02,\
             1.02678141e-01,  -5.14275775e-02,   5.62724780e-02,\
            -1.93720036e-02,   4.34821193e-02,  -1.45153970e-02,\
             2.82637773e-02,  -1.18071213e-02,   8.05928796e-02,\
            -2.29311085e-02,   2.63265016e-02,   1.51460948e-02,\
             7.57040806e-02,  -2.05735474e-02,   3.70839028e-02,\
            -2.98834239e-03,   1.07671225e-01,  -2.12856366e-02,\
             4.40466836e-02,   5.05057371e-02,  -9.97076885e-02,\
             1.29048184e-01,  -1.16175167e-01,   1.28964313e-01,\
            -7.37410008e-02,   1.21596750e-01,  -1.05729889e-01,\
             1.46546884e-01,   9.98635189e-03,   5.06216314e-02,\
            -4.26121940e-02,   8.58562278e-02,  -1.50241042e-02,\
             5.09680257e-02,  -3.10213066e-02,   6.62052252e-02,\
             3.00518843e-02,   1.36078198e-02,  -1.50161062e-03,\
             3.90682295e-02,   3.05597712e-02,   1.32485471e-02,\
             8.96118019e-03,   3.48847751e-02,   7.01886411e-02,\
            -1.01213799e-02,   1.90875047e-02,   3.41492896e-02,\
             7.25117903e-02,  -1.44480737e-03,   3.94956793e-02,\
             5.07343913e-02,  -6.44653593e-02,   8.46727673e-02,\
            -6.94494495e-02,   8.06758105e-02,  -6.48397657e-02,\
             7.85369386e-02,  -6.25669177e-02,   7.47904050e-02,\
            -3.10821154e-02,   6.15680910e-02,  -5.03927507e-02,\
             6.22716631e-02,  -4.67365901e-02,   5.71258565e-02,\
            -4.34634792e-02,   6.03364567e-02,  -5.05938366e-03,\
             2.71254306e-02,  -1.41774993e-02,   2.55928616e-02,\
            -8.22106071e-03,   2.03010392e-02,  -1.20029311e-02,\
             3.26075184e-02,   4.68397222e-02,  -1.76223811e-02,\
             3.65431091e-02,  -2.03977044e-02,   5.12240943e-02,\
            -2.76616471e-02,   3.70813606e-02,  -1.47730481e-02,\
            -6.56924377e-02,   7.98887742e-02,  -6.39212561e-02,\
             7.95195704e-02,  -6.94483666e-02,   8.28454807e-02,\
            -6.75886832e-02,   9.01558776e-02,  -3.93522614e-02,\
             5.78583078e-02,  -4.41788605e-02,   5.62099170e-02,\
            -4.65644690e-02,   6.04029054e-02,  -4.92043190e-02,\
             8.18061959e-02,   1.50163315e-03,   2.05123153e-02,\
            -1.22212060e-02,   2.57535697e-02,  -1.37552556e-02,\
             2.68003404e-02,  -1.38410533e-02,   3.75233269e-02,\
             5.16982175e-02,  -2.80681222e-02,   3.74867928e-02,\
            -1.27117898e-02,   3.48640520e-02,  -1.72907294e-02,\
             3.62121054e-02,  -5.59893948e-03,  -8.25584682e-02,\
             1.24779429e-01,  -1.08897533e-01,   1.59752943e-01,\
            -1.12195541e-01,   1.27053637e-01,  -1.14216750e-01,\
             1.46210501e-01,  -1.53615793e-02,   5.21415801e-02,\
            -3.21821348e-02,   6.99701296e-02,  -5.36634310e-03,\
             5.01775388e-02,  -4.22295087e-02,   1.05247582e-01,\
             3.18800211e-02,   1.34173037e-02,   8.80424308e-03,\
             3.64696463e-02,   2.55418604e-02,   1.36167062e-02,\
            -1.53057169e-03,   4.68703251e-02,   9.33094691e-02,\
            -1.66820938e-03,   3.97207872e-02,   3.25416087e-02,\
             6.46569622e-02,  -9.91175208e-03,   1.88705679e-02,\
             4.25452815e-02,   6.31462546e-02,  -2.76282980e-02,\
             2.98793132e-02,   2.81419521e-02,   3.21746144e-02,\
            -1.69518615e-02,   2.78734574e-02,  -2.60800471e-03,\
             9.11183072e-02,  -5.22028815e-02,   5.70637168e-02,\
            -4.60819392e-03,   6.86761406e-02,  -4.57682530e-02,\
             5.53543954e-02,  -1.63436554e-02,   6.63388115e-02,\
            -2.32745701e-02,   2.66637140e-02,   3.22389442e-02,\
             3.22905557e-02,  -1.41918731e-02,   2.79036320e-02,\
             2.66771360e-03,   1.17802697e-01,  -2.16920596e-02,\
             4.44522337e-02,   4.29510295e-02,   5.58013593e-02,\
            -2.02146911e-02,   3.67254361e-02,   1.97010084e-02,\
             9.07611202e-02,  -6.97662042e-02,   7.75828112e-02,\
            -5.56307789e-02,   8.51922777e-02,  -6.51020968e-02,\
             7.05303852e-02,  -5.74572154e-02,   6.05462792e-02,\
            -3.97317496e-02,   5.18365925e-02,  -2.63990302e-02,\
             5.88499203e-02,  -3.75770449e-02,   4.63563361e-02,\
            -3.67540043e-02,   3.38860550e-02,  -2.11945126e-03,\
             1.97144486e-02,   6.59087151e-03,   2.58379419e-02,\
            -5.25982457e-03,   1.13371278e-02,  -1.86850887e-03,\
             9.30908082e-02,  -2.52850309e-02,   6.51271691e-02,\
            -3.42760511e-03,   6.14749455e-02,  -3.81459711e-02,\
             4.34811595e-02,  -2.71231664e-02,   2.48248156e-02,\
             1.38300322e-02,   1.49019285e-02,   3.26853885e-02,\
             3.91835533e-02,  -5.18802020e-04,   2.69021271e-02,\
             2.06454612e-02,   4.81507286e-02,  -4.85858212e-04,\
             2.68729161e-02,   1.17375898e-02,   3.37609596e-02,\
             1.38629918e-02,   1.48727020e-02,   2.37477691e-02,\
             8.54904385e-03,   4.54968289e-02,   6.38268381e-02,\
             8.52006189e-02,   1.05285260e-01,   9.84352011e-02,\
             1.82729149e-01,  -1.32106016e-01,   2.63925516e-01,\
            -2.18060320e-01,   2.51151559e-01,  -2.05943669e-01,\
             1.95116248e-01,  -1.44045332e-01,  -6.49774047e-01,\
             7.46754354e-01,   3.26644486e-01,  -4.01050573e-02,\
             1.17329030e+00,  -8.36299493e-01,   3.26645663e-01,\
            -4.01027920e-02,   1.45592571e-01,   6.03149217e-02,\
            -7.61726862e-02,   2.79812070e-01,  -3.16583204e-01,\
             6.75893254e-01,  -7.61667358e-02,   2.79814754e-01,\
             1.80001411e-01,  -6.05261588e-02,   1.32462405e-01,\
            -8.42771165e-02,   8.97089007e-02,  -6.84275724e-02,\
             8.50520769e-02,  -6.98541999e-02,   7.05317909e-02,\
            -1.12558317e-02,  -6.42638966e-02,   1.57223250e-01,\
            -8.61942139e-02,   1.44075097e-01,  -8.93201934e-02,\
             1.22308607e-01,  -1.05114414e-01,   1.40041417e-01,\
            -8.45425697e-02,   1.91006857e-01,  -2.33393192e-01,\
             4.18137516e-01,  -2.59760656e-01,   3.50701863e-01,\
            -2.90386758e-01,   3.35813106e-01,  -2.73105537e-01,\
             3.25490418e-01,  -2.64641605e-01,   3.57381346e-01,\
            -6.25954877e-02,   1.55623412e-01,  -8.48973970e-02,\
             1.42855197e-01,  -8.85339350e-02,   1.21441484e-01,\
            -1.05143000e-01,   1.39936403e-01,  -8.44985961e-02,\
             1.90972051e-01,   6.98236475e-02,  -5.66073627e-02,\
             6.29443880e-02,  -4.16043614e-02,   8.31617338e-02,\
            -5.89837519e-02,   6.81666480e-02,  -4.56740127e-02,\
             4.18591794e-02,  -3.25103135e-02,   4.27833791e-02,\
            -1.98054776e-02,   6.15511573e-02,  -3.32464481e-02,\
             4.74751741e-02,  -2.47850755e-02,   1.33609294e-02,\
            -3.94420974e-03,   1.10784120e-02,   1.14210797e-02,\
             2.88553346e-02,   5.74142257e-04,   2.01111549e-02,\
             1.57128522e-02,   5.26012262e-02,  -3.48240209e-02,\
             4.11091869e-02,  -1.52021282e-02,   8.93868977e-02,\
            -1.80573578e-02,   6.49272023e-02,   1.34080525e-02,\
             4.13319576e-02,  -1.36796835e-02,   2.65380665e-02,\
            -1.09373399e-02,   9.04288385e-02,  -2.39083991e-02,\
             2.65388938e-02,   1.31070227e-02,   8.44527367e-02,\
            -4.02915235e-02,   5.16517741e-02,  -2.65541052e-02,\
             1.06214725e-01,  -4.59444640e-02,   5.16442307e-02,\
            -7.83223891e-03,   4.55557819e-02,  -1.12869550e-02,\
             2.74614680e-02,  -7.87011865e-03,   8.85715700e-02,\
            -2.06529780e-02,   2.46475473e-02,   2.44058190e-02,\
             8.17466779e-02,  -1.67892250e-02,   3.62131004e-02,\
             4.12686061e-03,   1.18665960e-01,  -1.70098699e-02,\
             4.37874913e-02,   6.77565591e-02,  -7.91370099e-02,\
             1.14317266e-01,  -9.91764011e-02,   1.14533117e-01,\
            -5.34111177e-02,   1.10445226e-01,  -9.17809711e-02,\
             1.39710273e-01,   2.56896296e-02,   4.61980534e-02,\
            -3.67754103e-02,   8.77178165e-02,  -6.27648519e-03,\
             4.88556366e-02,  -2.53891850e-02,   6.69413334e-02,\
             3.76492624e-02,   1.39646602e-02,   2.77933162e-04,\
             4.41234763e-02,   3.66946708e-02,   1.50347771e-02,\
             1.10943133e-02,   4.06718835e-02,   7.89501036e-02,\
            -8.15451556e-03,   1.87028876e-02,   4.41677211e-02,\
             8.03803945e-02,   3.34126055e-03,   4.14244181e-02,\
             6.49535872e-02,  -5.18561227e-02,   7.60338246e-02,\
            -5.81258417e-02,   7.16987841e-02,  -5.57823142e-02,\
             7.21551630e-02,  -5.33679269e-02,   6.80427387e-02,\
            -2.01066976e-02,   5.63850383e-02,  -4.32380304e-02,\
             5.74725728e-02,  -4.11595673e-02,   5.35924968e-02,\
            -3.75193972e-02,   5.75729035e-02,  -1.18753130e-04,\
             2.63226715e-02,  -1.10898568e-02,   2.47439035e-02,\
            -5.33068873e-03,   1.97332226e-02,  -9.97075380e-03,\
             3.43804464e-02,   4.82998968e-02,  -1.37643393e-02,\
             3.60239463e-02,  -1.68487926e-02,   5.24400186e-02,\
            -2.45941867e-02,   3.56761697e-02,  -9.24846626e-03,\
            -5.56058906e-02,   7.19803809e-02,  -5.31967118e-02,\
             7.12651928e-02,  -6.04604199e-02,   7.58880975e-02,\
            -5.79400358e-02,   8.41240812e-02,  -3.18814488e-02,\
             5.34370346e-02,  -3.73440259e-02,   5.12731984e-02,\
            -4.01943610e-02,   5.62324513e-02,  -4.30579215e-02,\
             8.10490621e-02,   6.01615582e-03,   1.97051509e-02,\
            -9.95092311e-03,   2.56662594e-02,  -1.11679209e-02,\
             2.63053637e-02,  -1.10591718e-02,   3.86955866e-02,\
             5.21982970e-02,  -2.45978617e-02,   3.56785974e-02,\
            -6.65829942e-03,   3.42701632e-02,  -1.37668480e-02,\
             3.60272451e-02,  -1.68566455e-04,  -6.12188107e-02,\
             1.11040845e-01,  -9.23605843e-02,   1.51276791e-01,\
            -9.76959698e-02,   1.14882089e-01,  -9.97822752e-02,\
             1.36842815e-01,  -6.12844630e-03,   4.92391229e-02,\
            -2.57578793e-02,   6.98832186e-02,   6.09270001e-03,\
             4.65355293e-02,  -3.71850458e-02,   1.10821138e-01,\
             3.79365099e-02,   1.51599601e-02,   1.09830351e-02,\
             4.20757134e-02,   3.18016441e-02,   1.40726795e-02,\
             1.46298643e-04,   5.29646338e-02,   1.04185865e-01,\
             3.37384411e-03,   4.13938089e-02,   4.35089676e-02,\
             7.23162713e-02,  -8.12724909e-03,   1.86670262e-02,\
             5.34568039e-02,   6.50704664e-02,  -2.37265399e-02,\
             2.63747559e-02,   4.18242922e-02,   3.11128248e-02,\
            -1.35186756e-02,   2.63675413e-02,   3.08580900e-03,\
             9.10854575e-02,  -4.57697092e-02,   5.14881519e-02,\
             1.02938868e-02,   6.68026700e-02,  -4.01664948e-02,\
             5.14441070e-02,  -5.66105013e-03,   7.10005158e-02,\
            -2.06327061e-02,   2.46199168e-02,   4.46109642e-02,\
             3.23622767e-02,  -1.12844906e-02,   2.74159396e-02,\
             8.36726956e-03,   1.29721509e-01,  -1.70073570e-02,\
             4.37839015e-02,   5.90423422e-02,   5.84445108e-02,\
            -1.67860782e-02,   3.62104396e-02,   3.00671134e-02,\
             8.38294753e-02,  -5.95652356e-02,   6.87603969e-02,\
            -4.30836814e-02,   8.03829192e-02,  -5.71682740e-02,\
             6.35541736e-02,  -4.84253837e-02,   5.78183936e-02,\
            -3.36165691e-02,   4.78572338e-02,  -1.81022079e-02,\
             5.75829220e-02,  -3.28561770e-02,   4.31846271e-02,\
            -3.21394454e-02,   3.66959718e-02,   4.63252869e-04,\
             2.02367267e-02,   1.05209004e-02,   2.80476327e-02,\
            -4.05625464e-03,   1.12060210e-02,  -2.73306552e-04,\
             9.75954406e-02,  -1.80875383e-02,   6.49602145e-02,\
             7.53672383e-03,   6.20863189e-02,  -3.48598365e-02,\
             4.11364991e-02,  -2.20402933e-02,   2.89895387e-02,\
             1.64851502e-02,   1.73171559e-02,   3.86709665e-02,\
             4.05626701e-02,   4.87415077e-03,   2.61646724e-02,\
             2.98607213e-02,   5.10877479e-02,   4.90048198e-03,\
             2.61427344e-02,   1.93313941e-02,   3.95138243e-02,\
             1.65114807e-02,   1.72952172e-02,   2.81421338e-02,\
             8.23694836e-03,   5.56132224e-02,   7.30028564e-02,\
             1.02523746e-01,   1.21577287e-01,   1.18160587e-01,\
             1.67052854e-01,  -1.07255467e-01,   2.41358154e-01,\
            -1.87818456e-01,   2.26724009e-01,  -1.73192111e-01,\
             1.81361691e-01,  -1.21558710e-01,  -5.98129836e-01,\
             7.09049131e-01,   3.35131690e-01,   2.13819038e-03,\
             1.11056225e+00,  -7.08644672e-01,   3.35133351e-01,\
             2.14055979e-03,   1.48285832e-01,   9.45691408e-02,\
            -5.71338133e-02,   2.97211678e-01,  -2.84492229e-01,\
             7.07317875e-01,  -5.71265675e-02,   2.97214595e-01],4)

        self.assertEqual(np.round(self.jacket.jacket_utilization.cb_util,4).tolist(), cb_util.tolist())

if __name__ == "__main__":
    unittest.main()
