#-------------------------------------------------------------------------------
# Name:        MyJacketInputs.py
# Purpose:     This module is  a template to set up a Jacket Assembly with basic input & optimization parameters (to be set (harwired here))
#              THIS IS JUST A TEMPLATE:
#              in order to modify this: 1. Copy this file and rename as 'yourinputfile.py'
#                                       2. Edit all the jacket inputs till the line that prohibits to further edit
#                                       3. Launch JacketOpt_PyOPT pointing to 'yourinputfile.py'
# Author:      rdamiani
#
# Created:     11/04/2014
# Copyright:   (c) rdamiani 2014
# Licence:     Apache (2014)
#-------------------------------------------------------------------------------

import numpy as np
from openmdao.main.api import set_as_top
from jacket import JcktGeoInputs,SoilGeoInputs,WaterInputs,WindInputs,RNAprops,TPlumpMass,Frame3DDaux,\
                    MatInputs,LegGeoInputs,XBrcGeoInputs,MudBrcGeoInputs,HBrcGeoInputs,TPGeoInputs,PileGeoInputs,\
                    TwrGeoInputs,JacketSE
from scipy.interpolate import interp1d

def main(): #\
    """Function to Instantiate a JacketSE Assembly: \n
       INPUTS \n
             All hardwired, so edit the quantities below all the way to the line "#________________ DO NOT MODIFY THE FOLLOWING ________________#" \n
             -See JacketOpt_PyOPT for more information. \n
       OUTPUTS \n
             myjckt -jacket assembly instance \n\n

             Optimization parameters:    \n\n

             f0          -float, target frequency [Hz]
             f0epsilon   -float,  f0*(1+f0epsilon) will not be exceeded \n
             jcktDTRmin  -Float, minimum jacket member DTR allowed, mostly for 60+waterdepths. \n
             mxftprint   -Float, max allowed foot print [m]
             guesses     -Float(n), guesses for all design variables check out DesVar class. \n
             bounds      -Float(n,2), bounds for all design variables check out DesVar class. \n\n
             SAMPLE CALLS: \n
             1.OPTIMIZATION: python JacketOpt_ExtCobyla.py C:\RRD\PYTHON\WISDEM\JacketSE\src\jacketse\MyJacketInputs.py \n
             2.OPTIMIZATION: python JacketOpt_PyOPT.py C:\RRD\PYTHON\WISDEM\JacketSE\src\jacketse\MyJacketInputs.py True \n
             3.BUILD JACKET: python >>> myjacket=C:\RRD\PYTHON\WISDEM\JacketSE\src\jacketse\MyJacketInputs.py \n
        """

    #Set inputs							     ###----ALL USER INPUT----###
    Jcktins=JcktGeoInputs()

    Jcktins.nlegs =4
    Jcktins.nbays =5
    Jcktins.batter=12.
    Jcktins.dck_botz =16.
    Jcktins.dck_width=2*6.
    Jcktins.weld2D   =0.5
    Jcktins.VPFlag = True    #vertical pile T/F;  to enable piles in frame3DD set pileinputs.ndiv>0
    Jcktins.clamped= False    #whether or not the bottom of the structure is rigidly connected. Use False when equivalent spring constants are being used.
    Jcktins.AFflag = False  #whether or not to use apparent fixity piles
    Jcktins.PreBuildTPLvl = 5  #if >0, the TP is prebuilt according to rules per PreBuildTP

    #______________________________________________#

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
    #______________________________________________#

    #Water and wind inputs
    Waterinputs=WaterInputs()
    Waterinputs.wdepth   =30.
    Waterinputs.wlevel   =30. #Distance from bottom of structure to surface  THIS, I believe is no longer needed as piles may be negative in z, to check and remove in case
    Waterinputs.T=12.  #Wave Period
    Waterinputs.HW=10. #Wave Height
    Waterinputs.Cd=3.  #Drag Coefficient, enhanced to account for marine growth and other members not calculated
    Waterinputs.Cm=8.#2.  #ADded mass Coefficient

    Windinputs=WindInputs()
    Windinputs.HH=100. #CHECK HOW THIS COMPLIES....
    Windinputs.U50HH=30. #assumed gust speed
    Windinputs.Cdj=4.  #Drag Coefficient for jacket members, enhanced to account for TP drag not calculated otherwise
    Windinputs.Cdt=2  #Drag Coefficient for tower, enhanced to account for TP drag not calculated otherwise
    #______________________________________________#

    #Pile data
    Pilematin=MatInputs()
    Pilematin.matname=np.array(['heavysteel'])

    Pileinputs=PileGeoInputs()
    Pileinputs.Pilematins=Pilematin
    Pileinputs.ndiv=0 #3			   ###----USER INPUT----###
    Pileinputs.Dpile=2.5
    Pileinputs.tpile=0.0254
    Pileinputs.Lp=20. #[m] Embedment length
    #______________________________________________#

    #Legs data
    legmatin=MatInputs()
    legmatin.matname=(['heavysteel','heavysteel','heavysteel','heavysteel'])
    Dleg=np.asarray([1.5]).repeat(Jcktins.nbays+1)                      #e.g., np.array([2.0,1.8,1.6,1.6,1.6])
    tleg=np.asarray([1.5*0.0254]).repeat(Jcktins.nbays+1)

    leginputs=LegGeoInputs()
    leginputs.legZbot   = 1.0		 ###----USER INPUT----###
    leginputs.ndiv=1     			 ###----USER INPUT----###
    leginputs.legmatins=legmatin
    leginputs.Dleg=Dleg
    leginputs.tleg=tleg


    #The following is a passthrough variables
    legbot_stmphin =1.5  #Distance from bottom of leg to second joint along z; must be>0
    #______________________________________________#

    #Xbrc data
    Xbrcmatin=MatInputs()
    Xbrcmatin.matname=np.array(['heavysteel']).repeat(Jcktins.nbays)
    Dbrc=np.asarray([1.]).repeat(Jcktins.nbays)#np.array([1.,1.,0.8,0.8])
    tbrc=np.asarray([0.0254]).repeat(Jcktins.nbays)

    Xbrcinputs=XBrcGeoInputs()
    Xbrcinputs.Dbrc=Dbrc
    Xbrcinputs.tbrc=tbrc
    Xbrcinputs.ndiv=1				 ###----USER INPUT----###
    Xbrcinputs.Xbrcmatins=Xbrcmatin
    Xbrcinputs.precalc=False   #This can be set to true if we want Xbraces to be precalculated in D and t, in which case the above set Dbrc and tbrc would be overwritten
    #______________________________________________#

    #Mbrc data
    Mbrcmatin=MatInputs()
    Mbrcmatin.matname=np.array(['heavysteel'])

    Mbrcinputs=MudBrcGeoInputs()
    Mbrcinputs.Dbrc_mud=1.5                  ###----USER INPUT----###
    Mbrcinputs.tbrc_mud=1.5*0.0254
    Mbrcinputs.ndiv=2					    ###----USER INPUT----###
    Mbrcinputs.Mbrcmatins=Mbrcmatin
    Mbrcinputs.precalc=False   #This can be set to true if we want Mudbrace to be precalculated in D and t, in which case the above set Dbrc_mud and tbrc_mud would be overwritten
    #______________________________________________#

    #Hbrc data
    Hbrcmatin=MatInputs()
    Hbrcmatin.matname=np.array(['heavysteel'])
    Dbrc_hbrc=1.1				           ###----USER INPUT----###

    Hbrcinputs=HBrcGeoInputs()
    Hbrcinputs.Dbrch=Dbrc_hbrc					   ###----USER INPUT----###
    Hbrcinputs.ndiv=0#2
    Hbrcinputs.Hbrcmatins=Hbrcmatin
    Hbrcinputs.precalc=True   #This can be set to true if we want Hbrace to be set=Xbrace top D and t, in which case the above set Dbrch and tbrch would be overwritten
    #______________________________________________#

    #TP data


    #Note PrebuildTPLvl is set in JacketIns				           ###----USER INPUT----###
    #TP lumped mass data
    TPlumpinputs=TPlumpMass()						      ###----ALL USER INPUT----###
    TPlumpinputs.mass=200.e3 #[kg]
    TPlumpinputs.CMoff= np.zeros(3)
    TPlumpinputs.I=np.zeros(6) #Note that the code will assume the mass distributed at the corners of the TP for Ixx,Iyy,Izz purposes

    TPstrtmatin=MatInputs()
    TPstmpsmatin=MatInputs()
    TPgirdmatin=MatInputs()
    TPbrcmatin=MatInputs()
    TPstemmatin=MatInputs()
    TPstmpsmatin.matname=np.array(['heavysteel'])
    TPstrtmatin.matname=np.array(['heavysteel'])
    TPbrcmatin.matname=np.array(['heavysteel'])
    TPgirdmatin.matname=np.array(['heavysteel'])
    TPstemmatin.matname=np.array(['heavysteel']).repeat(2) ###----ALL USER INPUT----###

    TPinputs=TPGeoInputs()
    TPinputs.TPstrtmatins=TPstrtmatin
    TPinputs.TPbrcmatins=TPbrcmatin
    TPinputs.TPstemmatins=TPstemmatin
    TPinputs.TPstmpmatins=TPstmpsmatin
    TPinputs.TPgirdmatins=TPgirdmatin

    #Set TP dimensions as leg and brace dimensions
    TPinputs.Dstrut=leginputs.Dleg[-1]				    ###----ALL USER INPUT----###
    TPinputs.tstrut=leginputs.tleg[-1]
    TPinputs.Dgir=Dbrc_hbrc
    TPinputs.tgir=0.0254
    TPinputs.Dbrc=TPinputs.Dgir
    TPinputs.tbrc=TPinputs.tgir
							     ###----ALL USER INPUT----###
    TPinputs.hstump=0.0#1.0
    TPinputs.Dstump=1.2
    TPinputs.tstump=0.04
    TPinputs.stumpndiv=1
    TPinputs.brcndiv=1
    TPinputs.girndiv=1
    TPinputs.strutndiv=1
    TPinputs.stemndiv=1
    TPinputs.nstems=3
    TPinputs.Dstem=np.array([6.]).repeat(TPinputs.nstems)
    TPinputs.tstem=np.array([0.1,0.11,0.11])
    TPinputs.hstem=np.array([6./TPinputs.nstems]).repeat(TPinputs.nstems)
    #______________________________________________#

    #Tower data
    Twrmatin=MatInputs()
    Twrmatin.matname=np.array(['heavysteel'])

    Twrinputs=TwrGeoInputs()
    Twrinputs.Twrmatins=Twrmatin
    #Twrinputs.Htwr=0.  #Trumped by HH
    Twrinputs.Htwr2frac=0.2  #fraction of tower height with constant x-section
    Twrinputs.ndiv=np.array([6,12])  #ndiv for uniform and tapered section			   ###----USER INPUT----###
    Twrinputs.DeltaZmax= 6. #[m], maximum FE element length allowed in the tower members (i.e. the uniform and the tapered members)
    Twrinputs.Db=5.6
    Twrinputs.DTRb=130.
    Twrinputs.Dt=0.55*Twrinputs.Db
    Twrinputs.DTRt=150.
        #Set whether or not DTRb and DTRt for the tower are the same. Note if next set to False it will trump DTRt setting above
    Twrinputs.DTRsdiff=False    ##SET THIS TO TRUE IF YOU WANT DTRs to be different between base and top

        #If you use the following 14 lines, The geometry defined above is ignored
    ###ztwr=np.array([  52.   ,   54.522,   57.044,   59.566,   62.088,   64.61 , 67.132,   72.176,   77.22 ,   82.264,   87.308,   92.352, 97.396,  102.44 ,  107.484,  112.528,  117.572,  122.616, 127.66  ]]) -52.
    ###Dtwr=np.array([ 5.6 ,  5.6 ,  5.6 ,  5.6 ,  5.6 ,  5.6 ,  5.6 ,  5.39,  5.18, 4.97,  4.76,  4.55,  4.34,  4.13,  3.92,  3.71,  3.5 ,  3.29, 3.08])
    ###ttwr=np.array([ 0.04307692,  0.04307692,  0.04307692,  0.04307692,  0.04307692, 0.04307692,  0.04307692,  0.04119829,  0.03931966,  0.03744103, 0.03556239,  0.03368376,  0.03180513,  0.0299265 ,  0.02804786, 0.02616923,  0.0242906 ,  0.02241197, 0.02053])
    ###pmtwr=[]  #(3,2) first col z's, second weights
    ####Interpolate data to refine tower discretization
    ###dz=1. #[m] maximum deltaz allowed in discretization
    ###ztwr2=np.linspace(ztwr[0],ztwr[-1],round((ztwr[-1]-ztwr[0])/dz))   #New discretization
    ###Dtwr_interp=interp1d(ztwr,Dtwr)
    ###ttwr_interp=interp1d(ztwr,ttwr)
    ###Twrinputs.ztwr=ztwr2
    ###Twrinputs.Dtwr=Dtwr_interp(ztwr2)
    ###Twrinputs.ttwr=ttwr_interp(ztwr2)
    ###Twrinputs.TwrlumpedMass=np.zeros([pmtwr.shape[0],11])
    ###Twrinputs.TwrlumpedMass[:,0:2]=pmtwr

    TwrRigidTop=True           #False=Account for RNA via math rather than a physical rigidmember
    #______________________________________________#

    #RNA data
    RNAins=RNAprops()
    RNAins.mass=3*350.e3  #[kg]
    RNAins.I[0]=86.579E+6  #[kg m2]
    RNAins.I[1]=53.530E+6  #[kg m2]
    RNAins.I[2]=58.112E+6  #[kg m2]
    RNAins.CMoff[2]=2.34  #[m]
    RNAins.Thoff[2]=RNAins.CMoff[2]  #[m]
    RNAins.yawangle=45.  #angle with respect to global X, CCW looking from above, wind from left
    RNAins.rna_weightM=True
    #______________________________________________#

    #RNA loads              Fx-z,         Mxx-zz
    RNA_F=np.array([1000.e3,0.,0.,0.,0.,0.])   #unfactored thrust, though accounting for gust and dynamic effects (no IEC PSF though)
    #______________________________________________#

    # Frame3DD parameters           									      ###----ALL USER INPUT----###
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
    #______________________________________________#
    #______________________________________________#

# OTHER AUXILIARY CONSTRAINTS AND TARGETS FOR OPTIMIZATION #
    #______________________________________________#
    #______________________________________________#

    #Set Optimization Bounds and guesses for the various variables:
    #          x=  [ batter,  Dpile,    tpile,        Lp,   Dleg,     tleg,       Dbrc,   tbrc,     Dbrc_mud,   tbrc_mud,   Dgit,      tgir,      Db,   DTRb   Dt,   DTRt   Htwr2fac        dck_widthfact]
    MnCnst = np.array([ 8.,      1.,    1.*0.0254,   20.,   1.,       1.*0.0254,  1.,    1.*0.0254,   1.,     1.*0.0254,      1.,     1.*0.0254,  5.,   120.,  3.,   120.,     0.05,         2.])
    MxCnst = np.array([ 15.,     2.5,   5.*0.0254,   50.,   2.5,      5.*0.0254,  2.,    5.*0.0254,   2.,     5.*0.0254,      2.,     5.*0.0254,  7.,   200.,  4.,   200.,     0.25,         3.])
    guesses= np.array([  10.,    1.5,   1.5*0.0254,  26.,   1.8,     1.5*0.0254,  1.2,   1.5*0.0254,  1.2,    1.5*0.0254,     1.2,    1.5*0.0254, 6.,   140.,  3.5,  150.,     0.2,          2.])

    #SET Maximum Footprint [m]
    mxftprint =30.

    #Set target frequency [Hz] and f0epsilon, i.e. fmax=(1+f0eps)*f0
    f0=0.22
    f0epsilon=0.1

    #Set the minminimum DTR allowed for jacket members; mostly for 60+waterdepths
    jcktDTRmin=22.

   #_____________________________________________________________#
   #________________ DO NOT MODIFY THE FOLLOWING ________________#
   #_____________________________________________________________#

    bounds=np.vstack((MnCnst,MxCnst))
    desvarmeans=np.mean(bounds,1)


    # Now Launch the assembly and pass all of the inputs

    myjckt=set_as_top(JacketSE(Jcktins.clamped,Jcktins.AFflag,Jcktins.PreBuildTPLvl>0))
    myjckt.JcktGeoIn=Jcktins
    myjckt.Soilinputs=Soilinputs
    myjckt.Waterinputs=Waterinputs
    myjckt.Windinputs=Windinputs
    myjckt.Pileinputs=Pileinputs
    myjckt.leginputs=leginputs
    myjckt.legbot_stmphin =legbot_stmphin  #Distance from bottom of leg to second joint along z; must be>0
    myjckt.Xbrcinputs=Xbrcinputs
    myjckt.Mbrcinputs=Mbrcinputs
    myjckt.Hbrcinputs=Hbrcinputs
    myjckt.TPlumpinputs=TPlumpinputs
    myjckt.TPinputs=TPinputs

    myjckt.Twrinputs=Twrinputs
    myjckt.TwrRigidTop=TwrRigidTop       #Account for RNA via math rather than a physical rigidmember
    myjckt.RNAinputs=RNAins
    myjckt.RNA_F=RNA_F

    myjckt.FrameAuxIns=FrameAuxIns


    return myjckt,f0,f0epsilon,jcktDTRmin,mxftprint,guesses,bounds.T



if __name__ == '__main__':
    from PlotJacket import main as PlotJacket  #COMMENT THIS ONE OUT FOR PEREGRINE"S SAKE

    myjckt= main()[0]
    #--- RUN JACKET ---#
    myjckt.run()
    # ---------------- #

    #_____________________________________#
    #Now show results of modal analysis
    print('First two Freqs.= {:5.4f} and {:5.4f} Hz'.format(*myjckt.Frameouts.Freqs))
    #print component masses
    print('jacket+TP(structural+lumped) mass (no tower, no piles) [kg] = {:6.0f}'.format(myjckt.Frameouts.mass[0]+myjckt.TP.TPlumpinputs.mass-myjckt.Tower.Twrouts.mass))
    print('tower mass [kg] = {:6.0f}'.format(myjckt.Tower.Twrouts.mass))
    print('TP mass structural + lumped mass [kg] = {:6.0f}'.format(myjckt.TP.TPouts.mass+myjckt.TP.TPlumpinputs.mass))
    print('piles (all) mass (for assigned (not optimum, unless optimization is run) Lp [kg] = {:6.0f}'.format(myjckt.Mpiles))
    print('frame3dd model mass (structural + TP lumped) [kg] = {:6.0f}'.format(myjckt.Frameouts.mass[0]+myjckt.TP.TPlumpinputs.mass))
    print('frame3dd model mass (structural + TP lumped) + Pile Mass [kg] = {:6.0f}'.format(myjckt.Frameouts.mass[0]+myjckt.TP.TPlumpinputs.mass+myjckt.Mpiles))
    print('frame3dd model mass (structural only) no piles no tower [kg] = {:6.0f}'.format(myjckt.Frameouts.mass[0]-myjckt.Tower.Twrouts.mass))
    #print tower top displacement
    print('Tower Top Displacement in Global Coordinate System [m] ={:5.4f}'.format(*myjckt.Frameouts.top_deflection))
    #print max API code checks
    print('MAX member compression-bending utilization at joints = {:5.4f}'.format(np.max(myjckt.jacket_utilization.cb_util)))
    print('MAX member tension utilization at joints = {:5.4f}'.format(np.max(myjckt.jacket_utilization.t_util)))
    print('MAX X-joint  utilization at joints = {:5.4f}'.format(np.max(myjckt.jacket_utilization.XjntUtil)))
    print('MAX K-joint  utilization at joints = {:5.4f}'.format(np.max(myjckt.jacket_utilization.KjntUtil)))

    PlotJacket(myjckt,util=True)
