#-------------------------------------------------------------------------------
# Name:        SetJacketInputsPeregrine.py
# Purpose:  Call this module to set up a Jacket Assembly with basic input parameters (to be set (harwired here)
#           and with optimization parameters to be set via arguments.
#                THIS IS FOR 6 MW LCOE MACHINE 40 m water depth
# Author:      rdamiani
#
# Created:     11/04/2014
# Copyright:   (c) rdamiani 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np
from openmdao.main.api import set_as_top
from jacket import JcktGeoInputs,SoilGeoInputs,WaterInputs,WindInputs,RNAprops,TPlumpMass,Frame3DDaux,\
                    MatInputs,LegGeoInputs,XBrcGeoInputs,MudBrcGeoInputs,HBrcGeoInputs,TPGeoInputs,PileGeoInputs,\
                    TwrGeoInputs,JacketSE

def main(DesPrms,DesVars): #\
    """Function to Instantiate a JacketSE Assembly with Assigned Design Parameters:
       INPUTS
             DesVars: Class (structure) containing: batter,Dpile,tpile,Lp,Dleg,tleg,Dbrc,tbrc,Dbrc_mud,tbrc_mud,Dgir,tgir,Db,DTRb,Dt,DTRt,Htwr2frac. See JacketOpt_Peregrine.py.
             DesPrms: Class (structure) containing: dck_botz,wdepth,wlevel,Tp50,HW50,HH,U50HH,RNA_F,RNAins, TPlumpmass. See JacketOpt_Peregrine.py.
        """

         #wdepth=wdepth,wlevel=wlevel,): #we call this module with the main optimization parameters

    #______________________________________________#

              # NON-OPTIMIZATION VARIABLES #
    #______________________________________________#

    #Set inputs							     ###----ALL USER INPUT----###
    Jcktins=JcktGeoInputs()
    Jcktins.nlegs =4
    Jcktins.nbays =DesPrms.nbays
    Jcktins.batter=DesVars.batter
    Jcktins.dck_width=DesVars.dck_widthfact*DesVars.Db
    Jcktins.dck_botz =DesPrms.dck_botz
    Jcktins.weld2D   =0.5
    Jcktins.VPFlag= True    #vertical pile T/F;  to enable piles in frame3DD set pileinputs.ndiv>0
    Jcktins.clamped= False #False    #whether or not the bottom of the structure is rigidly connected. Use False when equivalent spring constants are being used.

    #The following is a passthrough variables
    legbot_stmphin =1.5  #Distance from bottom of leg to second joint along z; must be>0

    #______________________________________________#
    #Soil inputs
    Soilinputs=SoilGeoInputs()							 ###----ALL USER INPUT----###
    Soilinputs.zbots   =-np.array([3.,5.,7.,15.,30.,50.])
    Soilinputs.gammas  =np.array([10000.,10000.,10000.,10000.,10000.,10000.])
    Soilinputs.cus     =np.array([60000.,60000.,60000.,60000.,60000.,60000.])
    Soilinputs.phis    =np.array([36.,33.,26.,37.,35.,37.5])#np.array([26.,26.,26.,26.,26.,26])#np.array([36.,33.,26.,37.,35.,37.5])#np.array([36.,33.,26.,37.,35.,37.5])
    Soilinputs.delta   =25.
    Soilinputs.sndflg   =True
    Soilinputs.SoilSF   =1.25 #Safety factor for Lp and stiffness from soil
    Soilinputs.PenderSwtch   =False #True

    #______________________________________________#
    #Water and wind inputs
    Waterinputs=WaterInputs()					      ###----ALL USER INPUT----###
    Waterinputs.wdepth   =DesPrms.wdepth
    Waterinputs.wlevel   =DesPrms.wdepth #Distance from bottom of structure to surface  THIS, I believe is no longer needed as piles may be negative in z, to check and remove in case
    Waterinputs.T=        DesPrms.Tp50  #Wave Period
    Waterinputs.HW=       DesPrms.HW50 #Wave Height: peak-to-peak for Andrew's load routine, 0.5*peak-to-peak for mine in case I activate.
    Waterinputs.Cd=3.  #Drag Coefficient, enhanced to account for marine growth and other members not calculated
    Waterinputs.Cm=8.#2.  #ADded mass Coefficient

    Windinputs=WindInputs()
    Windinputs.HH=    DesPrms.HH #CHECK HOW THIS COMPLIES....
    Windinputs.U50HH=DesPrms.U50HH #20.  # Since we are assuming operational conditions, pick a wind speed near rated
    Windinputs.Cdj=4.  #Drag Coefficient for jacket members, enhanced to account for TP drag not calculated otherwise
    Windinputs.Cdt=2  #Drag Coefficient for tower, enhanced to account for TP drag not calculated otherwise
    #______________________________________________#
    #RNA loads              Fx-z,         Mxx-zz
    RNA_F=DesPrms.RNA_F    #unfactored thrust, though accounting for gust and dynamic effects (no IEC PSF though)

    #Torque at rated with 95% efficiency generator: +12564863.93   Nm

    #______________________________________________#
    #RNA Properties
    RNAins=RNAprops()						    ###----ALL USER INPUT----###
    RNAins.mass=DesPrms.RNAins.mass
    RNAins.Ixx=DesPrms.RNAins.Ixx    #249.667E+6           #
    RNAins.Iyy=DesPrms.RNAins.Iyy    #169.667E+6         #
    RNAins.Izz=DesPrms.RNAins.Izz   #162.200E+6          #
    RNAins.CMzoff=DesPrms.RNAins.CMzoff
    RNAins.CMxoff=DesPrms.RNAins.CMxoff#positive means downwind  -4.95#
    RNAins.Thzoff=DesPrms.RNAins.Thzoff  #From UHreact Excel

    RNAins.yawangle=45.
    TwrRigidTop=True #False       #False=Account for RNA via math rather than a physical rigidmember


    #______________________________________________#
    #TP mass data
    TPlumpinputs=TPlumpMass()						      ###----ALL USER INPUT----###
    TPlumpinputs.mass=DesPrms.TPlumpmass#0. #200.e3 #[kg]
    #______________________________________________#
									      ###----ALL USER INPUT----###
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
    FrameAuxIns.gvector=np.array([0.,0.,9.8065])    #GRAVITY



    #______________________________________________#

              # OPTIMIZATION VARIABLES #
              #Note: materials are fixed
    #______________________________________________#

    #Legs Data
    legmatin=MatInputs()
    legmatin.matname=(['heavysteel','heavysteel','heavysteel','heavysteel'])
    Dleg=np.asarray([DesVars.Dleg]).repeat(Jcktins.nbays+1)                      #e.g., np.array([2.0,1.8,1.6,1.6,1.6])
    tleg=np.asarray([DesVars.tleg]).repeat(Jcktins.nbays+1)

    leginputs=LegGeoInputs()
    leginputs.legZbot   = 0.0		 ###----USER INPUT----###
    leginputs.ndiv=DesPrms.legndiv			 ###----USER INPUT----###
    leginputs.legmatins=legmatin
    leginputs.Dleg=Dleg
    leginputs.tleg=tleg


    #Xbrc data
    Xbrcmatin=MatInputs()
    Xbrcmatin.matname=np.array(['heavysteel']).repeat(Jcktins.nbays)
    Dbrc=np.asarray([DesVars.Dbrc]).repeat(Jcktins.nbays)#np.array([1.,1.,0.8,0.8])
    tbrc=np.asarray([DesVars.tbrc]).repeat(Jcktins.nbays)

    Xbrcinputs=XBrcGeoInputs()
    Xbrcinputs.Dbrc=Dbrc
    Xbrcinputs.tbrc=tbrc
    Xbrcinputs.ndiv=1				 ###----USER INPUT----###
    Xbrcinputs.Xbrcmatins=Xbrcmatin
    Xbrcinputs.precalc=False   #This can be set to true if we want Xbraces to be precalculated in D and t, in which case the above set Dbrc and tbrc would be overwritten


    #Mbrc data
    Mbrcmatin=MatInputs()
    Mbrcmatin.matname=np.array(['heavysteel'])

    Mbrcinputs=MudBrcGeoInputs()
    Mbrcinputs.Dbrc_mud=DesVars.Dbrc_mud
    Mbrcinputs.tbrc_mud=DesVars.tbrc_mud
    Mbrcinputs.ndiv=2					    ###----USER INPUT----###
    Mbrcinputs.Mbrcmatins=Mbrcmatin
    Mbrcinputs.precalc=False   #This can be set to true if we want Mudbrace to be precalculated in D and t, in which case the above set Dbrc_mud and tbrc_mud would be overwritten


    #Hbrc data
    Hbrcmatin=MatInputs()
    Hbrcmatin.matname=np.array(['heavysteel'])
    Dbrc_hbrc=1.1				           ###----USER INPUT----###

    Hbrcinputs=HBrcGeoInputs()
    Hbrcinputs.Dbrch=Dbrc_hbrc					   ###----USER INPUT----###
    Hbrcinputs.ndiv=0#2
    Hbrcinputs.Hbrcmatins=Hbrcmatin
    Hbrcinputs.precalc=True   #This can be set to true if we want Hbrace to be set=Xbrace top D and t, in which case the above set Dbrch and tbrch would be overwritten


    #TP data
    TPstrtmatin=MatInputs()
    TPstmpsmatin=MatInputs()
    TPgirdmatin=MatInputs()
    TPbrcmatin=MatInputs()
    TPstemmatin=MatInputs()
    TPstmpsmatin.matname=np.array(['heavysteel'])
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
    TPinputs.Dstrut=DesVars.Dleg				    ###----ALL USER INPUT----###
    TPinputs.tstrut=DesVars.Dleg
    TPinputs.Dgir=DesVars.Dgir
    TPinputs.tgir=DesVars.tgir
    TPinputs.Dbrc=TPinputs.Dgir
    TPinputs.tbrc=TPinputs.tgir
							     ###----ALL USER INPUT----###
    TPinputs.hstump=0.0#1.0
    TPinputs.stumpndiv=1
    TPinputs.brcndiv=1
    TPinputs.girndiv=1
    TPinputs.strutndiv=1
    TPinputs.stemndiv=1
    TPinputs.nstems=3
    TPinputs.Dstem=np.array([DesVars.Db]).repeat(TPinputs.nstems)
    TPinputs.tstem=1.5*np.array([DesVars.Db/DesVars.DTRb]).repeat(TPinputs.nstems)
    TPinputs.hstem=np.array([6./TPinputs.nstems]).repeat(TPinputs.nstems)

    #Pile data
    Pilematin=MatInputs()
    Pilematin.matname=np.array(['heavysteel'])

    Pileinputs=PileGeoInputs()
    Pileinputs.Pilematins=Pilematin
    Pileinputs.ndiv=0 #3			   ###----USER INPUT----###
    Pileinputs.Dpile=DesVars.Dpile
    Pileinputs.tpile=DesVars.tpile
    Pileinputs.AFflag=False			   ###----USER INPUT----###
    Pileinputs.Lp=DesVars.Lp #[m] Embedment length



    #Tower data
    Twrmatin=MatInputs()
    Twrmatin.matname=np.array(['heavysteel'])

    Twrinputs=TwrGeoInputs()
    Twrinputs.Twrmatins=Twrmatin
    #Twrinputs.Htwr=0.  #Trumped by HH
    Twrinputs.Htwr2frac=DesVars.Htwr2frac   #fraction of tower height with constant x-section
    Twrinputs.ndiv=np.array([6,12])  #ndiv for uniform and tapered section			   ###----USER INPUT----###
    Twrinputs.DeltaZmax= 6. #[m], maximum FE element length allowed in the tower members (i.e. the uniform and the tapered members)
    Twrinputs.Db=DesVars.Db
    Twrinputs.DTRb=DesVars.DTRb
    Twrinputs.Dt=DesVars.Dt
    Twrinputs.DTRt=DesVars.DTRt


    # Now Launch the assembly and pass all of the inputs

    myjckt=set_as_top(JacketAsmly())

    myjckt.JcktGeoIn=Jcktins

    myjckt.Soilinputs=Soilinputs
    myjckt.Waterinputs=Waterinputs
    myjckt.Windinputs=Windinputs

    myjckt.RNAinputs=RNAins
    myjckt.TwrRigidTop=TwrRigidTop       #Account for RNA via math rather than a physical rigidmember
    myjckt.RNA_F=RNA_F

    myjckt.leginputs=leginputs
        #The following is a passthrough variables
    myjckt.legbot_stmphin =legbot_stmphin  #Distance from bottom of leg to second joint along z; must be>0

    myjckt.Xbrcinputs=Xbrcinputs
    myjckt.Mbrcinputs=Mbrcinputs

    myjckt.Hbrcinputs=Hbrcinputs
    myjckt.TPlumpinputs=TPlumpinputs

    myjckt.PrebuildTP=True				###----USER INPUT----###
    myjckt.PreBuildTPLvl=5                 ###----USER INPUT----###

    myjckt.TPinputs=TPinputs
    myjckt.Pileinputs=Pileinputs
    myjckt.Twrinputs=Twrinputs
    myjckt.FrameAuxIns=FrameAuxIns

    return myjckt



if __name__ == '__main__':
    main()
