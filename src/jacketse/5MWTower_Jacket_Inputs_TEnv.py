#-------------------------------------------------------------------------------
# Name:        5MW Tower_Jacket_Inputs_JEnv.py
# Purpose:     This module is a template to set up a Jacket Assembly with basic input & optimization parameters (to be set (harwired here))
#              It follows the generic template myjacketinputs.py
#              It is based on the OC4 Jacket design with NREL's 5-MW baseline turbine's tower.
#              It uses the environment from the NREL's 5-MW baseline turbine's tower
#
# Author:      cbroslawski
#              SULI intern under the guidance of Dr. Rick Damiani
#
# Created:     8/2016

#All lines marked CJB+, CJB-, or CJBe have been added, "removed" (commented out), or edited respectively by Casey Broslawski. Summer 2016
#-------------------------------------------------------------------------------

import numpy as np
import os
from openmdao.main.api import set_as_top
from jacket import JcktGeoInputs,SoilGeoInputs,WaterInputs,WindInputs,RNAprops,TPlumpMass,Frame3DDaux,\
                    MatInputs,LegGeoInputs,XBrcGeoInputs,MudBrcGeoInputs,HBrcGeoInputs,TPGeoInputs,PileGeoInputs,\
                    TwrGeoInputs,JacketSE
from Utilization import UtilAssembly #CJBtest test
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
    Jcktins.nbays =4
    Jcktins.batter=30.5495 #=(-43.127+24.614)/(5.939-5.333)
    #Jcktins.dck_botz =15.651
    Jcktins.dck_botz =15 #CJBe
    Jcktins.dck_width= 8.+1.2
    Jcktins.weld2D   = 0.
    Jcktins.VPFlag = True    #vertical pile T/F;  to enable piles in frame3DD set pileinputs.ndiv>0
    Jcktins.clamped= True    #whether or not the bottom of the structure is rigidly connected. Use False when equivalent spring constants are being used.
    Jcktins.AFflag = False  #whether or not to use apparent fixity piles
    Jcktins.PreBuildTPLvl = 0  #if >0, the TP is prebuilt according to rules per PreBuildTP

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

    #Water and wind inputs
    Waterinputs=WaterInputs() #CJBtest Max thrust case
    Waterinputs.wdepth   =20. #CJBtest Max thrust case
    Waterinputs.wlevel   =20. #CJBtest Max thrust case #Distance from bottom of structure to surface  THIS, I believe is no longer needed as piles may be negative in z, to check and remove in case
    Waterinputs.T=10.  #CJBtest Max thrust case #Wave Period
    Waterinputs.HW=8.0*1.86 #CJBtest Max thrust case #Wave Height
    Waterinputs.Cd=3.  #CJBtest Max thrust case #Drag Coefficient, enhanced to account for marine growth and other members not calculated
    Waterinputs.Cm=8.#2.  #CJBtest Max thrust case #ADded mass Coefficient
    WaterInputs.Uc=0.0 #CJBtest Max thrust case test

    Windinputs=WindInputs() #CJBtest Max thrust case
    Windinputs.HH=87.6+.51+Jcktins.dck_botz+4 #CJBtest Max thrust case (4 comes from the height of the TP) #CHECK HOW THIS COMPLIES....
    Windinputs.Cdj=4.   #CJBtest Max thrust case #Drag Coefficient for jacket members, enhanced to account for TP drag not calculated otherwise
    Windinputs.Cdt=2   #CJBtest Max thrust case #Drag Coefficient for tower, enhanced to account for TP drag not calculated otherwise
    Windinputs.al_shear=.2 #CJBtest Max thrust case test
    #Windinputs.U50HH=11.73732 #CJBtest Max thrust case test #assumed gust speed
    Windinputs.U50HH=70.  #CJBtest Max wind speed case test #assumed gust speed

    #______________________________________________#

    #Pile data
    Pilematin=MatInputs()
    Pilematin.matname=np.array(['steel']) #CJB Changed this from RC into steel
    Pilematin.E=np.array([ 2.1e11])
    Pilematin.G=np.array([8.07690e+10])
    Pilematin.rho=np.array([3339.12]) *34274.82/36876.   #From Test04.txt in SD, to check with official FAST certtest

    Pileinputs=PileGeoInputs()
    Pileinputs.Pilematins=Pilematin
    Pileinputs.ndiv=0 #CJB Change this from 1 to 0 #3			   ###----USER INPUT----###
    Pileinputs.Dpile=2.082
    Pileinputs.tpile=0.491
    Pileinputs.Lp=0. #[m] Embedment length
    #______________________________________________#

    #Legs data
    legmatin=MatInputs()
    legmatin.matname=(['steel'])
    legmatin.rho=np.array([7850.])
    Dleg=np.asarray([1.2]).repeat(Jcktins.nbays+1) #e.g., np.array([2.0,1.8,1.6,1.6,1.6])
    tleg=np.asarray([0.05,0.05,0.035,0.035,0.035]) #CJBtest Length must be 1 more than (nbays)

    leginputs=LegGeoInputs()
    leginputs.legZbot   = 4.5		 ###----USER INPUT----###
    leginputs.ndiv=1     			 ###----USER INPUT----###
    leginputs.legmatins=legmatin
    leginputs.Dleg=Dleg
    leginputs.tleg=tleg


    #The following is a passthrough variables
    legbot_stmphin =(45.5-43.127)  #=2.373 Distance from bottom of leg to second joint along z; must be>0
    #______________________________________________#

    #Xbrc data
    Xbrcmatin=MatInputs()
    Xbrcmatin.matname=np.array(['steel']).repeat(Jcktins.nbays)
    Xbrcmatin.rho=np.array([7850.])

    Dbrc=np.asarray([0.8]).repeat(Jcktins.nbays)
    tbrc=np.asarray([0.02]).repeat(Jcktins.nbays)

    Xbrcinputs=XBrcGeoInputs()
    Xbrcinputs.Dbrc=Dbrc
    Xbrcinputs.tbrc=tbrc
    Xbrcinputs.ndiv=1				 ###----USER INPUT----###
    Xbrcinputs.Xbrcmatins=Xbrcmatin
    Xbrcinputs.precalc=False   #This can be set to true if we want Xbraces to be precalculated in D and t, in which case the above set Dbrc and tbrc would be overwritten
    #______________________________________________#

    #Mbrc data
    Mbrcmatin=MatInputs()
    Mbrcmatin.matname=np.array(['steel'])
    Mbrcmatin.rho=np.array([7850.])

    Mbrcinputs=MudBrcGeoInputs()
    Mbrcinputs.Dbrc_mud=0.8                  ###----USER INPUT----###
    Mbrcinputs.tbrc_mud=0.02
    Mbrcinputs.ndiv=2					    ###----USER INPUT----###
    Mbrcinputs.Mbrcmatins=Mbrcmatin
    Mbrcinputs.precalc=False   #This can be set to true if we want Mudbrace to be precalculated in D and t, in which case the above set Dbrc_mud and tbrc_mud would be overwritten
    #______________________________________________#

    #Hbrc data
    Hbrcmatin=MatInputs()
    Hbrcmatin.matname=np.array(['steel'])
    Hbrcmatin.rho=np.array([7850.])
    Dbrc_hbrc=1.1				           ###----USER INPUT----###

    Hbrcinputs=HBrcGeoInputs()
    Hbrcinputs.Dbrch=Dbrc_hbrc					   ###----USER INPUT----###
    Hbrcinputs.ndiv=0
    Hbrcinputs.Hbrcmatins=Hbrcmatin
    Hbrcinputs.precalc=True   #This can be set to true if we want Hbrace to be set=Xbrace top D and t, in which case the above set Dbrch and tbrch would be overwritten
    #______________________________________________#

    #TP data


    #Note PrebuildTPLvl is set in JacketIns				           ###----USER INPUT----###
    #TP lumped mass data
    TPlumpinputs=TPlumpMass()						      ###----ALL USER INPUT----###
    TPlumpinputs.mass = 666.e3#-98385.33+( 7850.*4*np.pi/4.*(1.2**2-(1.2-2.*0.04)**2)*4)  #[kg]  TO MODIFY AFTER WE ASSESS OVERALL TP MASS, to be reduced for the steel weight
    TPlumpinputs.CMoff= np.array([0.,0.,2.]) #CG of concrete block is 2 m above intersection of diagonal braces
    TPlumpinputs.I    = 1./12*TPlumpinputs.mass*np.array([8.**2+4.**2,8.**2+4.**2,8.**2+8.**2,0.,0.,0.])

    TPstrtmatin=MatInputs()
    TPstmpsmatin=MatInputs()
    TPgirdmatin=MatInputs()
    TPbrcmatin=MatInputs()
    TPstemmatin=MatInputs()
    TPstmpsmatin.matname=np.array(['steel'])
    TPstmpsmatin.rho=np.array([7850.])

    TPstrtmatin.matname=np.array(['steel'])
    TPstrtmatin.rho=np.array([1350.])
    TPbrcmatin.matname=np.array(['steel'])
    TPbrcmatin.rho=np.array([1350.])#np.array([7850.]) tpstrucmass=np.pi/4*(1.2**2-(1.2-2*0.04)**2)*7850*4*4 7850.*tpstrucmass/myjckt.TP.TPouts.mass=1462.185
    TPgirdmatin.matname=np.array(['steel'])
    TPgirdmatin.rho=np.array([1350.])
    TPstemmatin.matname=np.array(['steel']).repeat(2) ###----ALL USER INPUT----###
    TPstemmatin.rho=np.array([1350.])

    TPinputs=TPGeoInputs()
    TPinputs.TPstrtmatins=TPstrtmatin
    TPinputs.TPbrcmatins=TPbrcmatin
    TPinputs.TPstemmatins=TPstemmatin
    TPinputs.TPstmpmatins=TPstmpsmatin
    TPinputs.TPgirdmatins=TPgirdmatin

    #Set TP dimensions as leg and brace dimensions
    TPinputs.Dstrut=1.2				    ###----ALL USER INPUT----###
    TPinputs.tstrut=0.04
    TPinputs.Dgir=TPinputs.Dstrut
    TPinputs.tgir=TPinputs.tstrut
    TPinputs.Dbrc=TPinputs.Dstrut
    TPinputs.tbrc=TPinputs.tstrut
							     ###----ALL USER INPUT----###
    TPinputs.hstump=0.499#(16.15-15.651)
    TPinputs.Dstump=1.2
    TPinputs.tstump=0.04

    TPinputs.stumpndiv=1
    TPinputs.brcndiv=1
    TPinputs.girndiv=1
    TPinputs.strutndiv=1
    TPinputs.stemndiv=1
    TPinputs.nstems=3
    TPinputs.Dstem=np.array([5.6]).repeat(TPinputs.nstems)
    TPinputs.tstem=np.array([0.032,0.032,0.032])
    TPinputs.hstem=np.array([(4.)/TPinputs.nstems]).repeat(TPinputs.nstems)
    #______________________________________________#

    #Tower data
    Twrmatin=MatInputs()
    Twrmatin.matname=np.array(['heavysteel']) #CJBtest
    Twrmatin.rho=np.array([8500.]) #CJBtest
    #Twrmatin.matname=np.array(['steel'])
    #Twrmatin.rho=np.array([7850.])

    Twrinputs=TwrGeoInputs()
    Twrinputs.Twrmatins=Twrmatin
    #Twrinputs.Htwr=88.15  #Trumped by HH
    Twrinputs.Htwr2frac=1./87.6 #CJBtest  #fraction of tower height with constant x-section
    Twrinputs.ndiv=np.array([1,1])  #ndiv for uniform and tapered section			   ###----USER INPUT----###
    Twrinputs.DeltaZmax= 5. #[m], maximum FE element length allowed in the tower members (i.e. the uniform and the tapered members)
    Twrinputs.Db=6. #CJBtest
    Twrinputs.DTRb=Twrinputs.Db/(1.3*.027) #CJBtest
    Twrinputs.Dt=3.87 #CJBtest
    Twrinputs.DTRt=Twrinputs.Dt/(1.3*.019) #CJBtest
        #Set whether or not DTRb and DTRt for the tower are the same. Note if next set to False it will trump DTRt setting above
    Twrinputs.DTRsdiff=True    ##SET THIS TO TRUE IF YOU WANT DTRs to be different between base and top

        #If you use the following 12 lines, The geometry defined above is ignored
        #CJBtest Don't use stations
    #ztwr=np.array([20.15,21.15,32.15,42.15,54.15,64.15,74.15,83.15,88.15]) -20.15
    #Dtwr=np.array([5.6,5.577,5.318,5.082,4.8,4.565,4.329,4.118,4])
    #ttwr=np.array([0.032,0.032,0.03,0.028,0.024,0.022,0.02,0.03,0.03])
    #pmtwr=np.array([ np.array([20.15,54.15,88.15])-20.15,[1.9e3,1.4e3,1.0e3]]).T  #(3,2) first col z's, second weights
    #Interpolate data to refine tower discretization
    #dz=1. #[m] maximum deltaz allowed in discretization
    #ztwr2=np.linspace(ztwr[0],ztwr[-1],round((ztwr[-1]-ztwr[0])/dz))   #New discretization
    #Dtwr_interp=interp1d(ztwr,Dtwr)
    #ttwr_interp=interp1d(ztwr,ttwr)
    #Twrinputs.ztwr=ztwr2
    #Twrinputs.Dtwr=Dtwr_interp(ztwr2)
    #Twrinputs.ttwr=ttwr_interp(ztwr2)
    #Twrinputs.TwrlumpedMass=np.zeros([pmtwr.shape[0],11])
    #Twrinputs.TwrlumpedMass[:,0:2]=pmtwr


    TwrRigidTop=False           #False=Account for RNA via math rather than a physical rigidmember
    #______________________________________________#

    #RNA data
    RNAins=RNAprops()
    RNAins.mass=350.e3 #CJBtest #[kg]
    RNAins.I[0]=114930678 #CJBtest #[kg m2]
    RNAins.I[1]=22035403 #CJBtest  #[kg m2]
    RNAins.I[2]=18759742.50 #CJBtest  #[kg m2]
    RNAins.I[4]=503710.47 #CJBtest  #[kg m2]
    RNAins.CMoff[0]=-1.13 #CJBtest    #[m]
    RNAins.CMoff[2]=.51 #CJBtest    #[m]
    #RNAins.Thoff[2]=2.4 #CJBtest     #[m]
    RNAins.yawangle=0. #CJBtest  #angle with respect to global X, CCW looking from above, wind from left
    RNAins.rna_weightM=True

    UtilAssembly.tilt=5.0 #CJBtest test
    #______________________________________________#

    #RNA loads              Fx-z,         Mxx-zz
    #RNA_F=np.array([1284744.196,0.,-112400.5527,3963732.762,896380.8464,-346781.6819]) #CJBtest Max thrust case  #unfactored thrust, though accounting for gust and dynamic effects (no IEC PSF though)
    RNA_F=np.array([188038.8045,0.,-16451.2637,0.0,131196.8431,0.0]) #CJBtest Max wind speed case  #unfactored thrust, though accounting for gust and dynamic effects (no IEC PSF though)
    #RNA_F=np.array([1000.e3,0.,0.,0.,0.,0.])
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

    #Decide whether or not to consider DLC 6.1 as well
    twodlcs=False

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
    #myjckt=set_as_top(JacketSE(Jcktins.clamped,Jcktins.AFflag,twodlcs=twodlcs, wlevel_Ulist=29))#CJB and JQ
    myjckt=set_as_top(JacketSE(Jcktins.clamped,Jcktins.AFflag,twodlcs=twodlcs))
        #pySubDyn Parameters CJB+
    #SDpySubDynA = pySubDynA()

    #INPUTS TO RUN SUBDYN-------------------------------------------------------
    #(PATH INFORMATION FOR INPUT FILE AND DRIVER)

    Base_name="5MW_Tower_TEnv_pySubDyn" #Input name here

    #INPUT FILE PATH
    myjckt.InputFile_name=str(Base_name)+".txt"
    #myjckt.SDpySubDynA.InputFile_name=str(Base_name)+".txt"
    myjckt.InputandDriverpath="C:\wisdem\plugins\JacketSE\src\jacketse\SubDyn\CertTest"
    myjckt.InputFile_path=myjckt.InputandDriverpath+os.sep+str(myjckt.InputFile_name)

    #DRIVER PATH
    myjckt.Driver_name=str(Base_name)+"D"+".txt"
    myjckt.Driver_path=myjckt.InputandDriverpath+os.sep+str(myjckt.Driver_name)
    myjckt.Driver_path=myjckt.InputandDriverpath+os.sep+str(myjckt.Driver_name)

    #PATH TO RUN SUBDYN
    SDEXEpath="C:\wisdem\plugins\JacketSE\src\jacketse\SubDyn"+os.sep+"bin\SubDyn_Win32.exe"
    myjckt.SDpath=str(SDEXEpath)+' '+str(myjckt.Driver_path)
    #test.SDpath='r'+"'''"+test.SDEXEpath+' '+test.SDDriverpath+"'''"

    #PATH TO READ OUTPUT (INPUTS TO READ OUTPUT)
    myjckt.Readpath_out=str(myjckt.InputandDriverpath)+os.sep+str(Base_name)+".SD.out"
    myjckt.Readpath_sum=str(myjckt.InputandDriverpath)+os.sep+str(Base_name)+".SD.sum"
    myjckt.Delete_file=False #Deletes driver, input, and output files. Does not delete Echo file.

    #INPUT FILE INPUTS----------------------------------------------------------

    #Simulation Control
    myjckt.Echo=np.array([False, "Echo", "- Echo input data to ""<rootname>.SD.ech"" (flag)"])
    myjckt.SDdeltaT=np.array(["DEFAULT", "SDdeltaT", "- Local Integration Step. If ""default"", the glue-code integration step will be used."])
    myjckt.IntMethod=np.array([4, "IntMethod", "- Integration Method [1/2/3/4 = RK4/AB4/ABM4/AM2]."])
    myjckt.SttcSolve=np.array([False, "SttcSolve", "- Solve dynamics about static equilibrium point"])

    #FEA and CRAIG-BAMPTON PARAMETERS
    myjckt.FEMmod=np.array([3, "- FEM switch: element model in the FEM. [1= Euler-Bernoulli(E-B);  2=Tapered E-B (unavailable);  3= 2-node Timoshenko;  4= 2-node tapered Timoshenko (unavailable)]"])
    myjckt.NDiv=np.array([1, "NDiv", "- Number of sub-elements per member"])#CJB "HARDWIRED" INTO PYSUBDYN AS 1 TO ALLOW JACKETSE'S NODES TO BE USED AS SUBDYN'S JOINTS
    myjckt.CBMod=np.array([False, "- [T/F] If True perform C-B reduction, else full FEM dofs will be retained. If True, select Nmodes to retain in C-B reduced system."])
    myjckt.Nmodes=np.array([75, "- Number of internal modes to retain (ignored if CBMod=False). If Nmodes=0 --> Guyan Reduction."])
    myjckt.JDampings=np.array([2, "JDampings", "- Damping Ratios for each retained mode (% of critical) If Nmodes>0, list Nmodes structural damping ratios for each retained mode (% of critical), or a single damping ratio to be applied to all retained modes. (last entered value will be used for all remaining modes)."])

    #Structure Joints
    myjckt.SDjointsHeader=np.array([['JointID','JointXss','JointYss','JointZss','[Coordinates of Member joints in SS-Coordinate System]'], ['(-)','(m)','(m)','(m)']])

    #Base Reaction Joints
    myjckt.BaseRxnJointsHeader=np.array([['RJointID','RctTDXss','RctTDYss','RctTDZss','RctRDXss','RctRDYss','RctRDZss','[Global Coordinate System]'], ['(-)',('flag'),('flag'),('flag'),('flag'),('flag'),('flag')]])

    #Interface Joints
    myjckt.InterfaceRxnJointsHeader=np.array([['IJointID','ItfTDXss','ItfTDYss','ItfTDZss','ItfRDXss','ItfRDYss','ItfRDZss','[Global Coordinate System]'], ['(-)',('flag'),('flag'),('flag'),('flag'),('flag'),('flag')]])

    #Members
    myjckt.MembersHeader=np.array([['MemberID','MJointID1','MJointID2','MPropSetID1','MPropSetID2','COSMID'], ['(-)','(-)','(-)','(-)','(-)','(-)']])

    #MEMBER X-SECTION PROPERTY data 1/2
    myjckt.NPropSets=np.array([6, 'NPropSets', '- # of structurally unique x-sections (i.e. # of X-sectional props utilized throughout all of the members)'])
    myjckt.PropSet1Header=np.array([['PropSetID', 'YoungE', 'ShearG', 'MatDens', 'XsecD', 'XsecT'],['(-)', '(N/m2)', '(N/m2)', '(kg/m3)', '(m)', '(m)']])
    myjckt.PropSet1=np.array([[1, 2.10000e+11, 8.07690e+10, 7850.00, 0.800000, 0.020000],\
                   [2, 2.10000e+11, 8.07690e+10, 7850.00, 1.200000, 0.050000],\
                   [3, 2.10000e+11, 8.07690e+10, 7850.00, 1.200000, 0.035000],\
                   [4, 2.10000e+11, 8.07690e+10, 7850.00, 1.200000, 0.040000],\
                   [5, 2.10000e+11, 8.07690e+10, 3339.12, 2.082000, 0.491000],\
                   [6, 2.10000e+11, 8.07690e+10, 7850.00, 2.082000, 0.060000]])

    #MEMBER X-SECTION PROPERTY data 2/2
    myjckt.PropSet2=np.array([])
    myjckt.PropSet2Header=np.array([["PropSetID","YoungE","ShearG","MatDens","XsecA","XsecAsx","XsecAsy","XsecJxx","XsecJyy","XsecJ0"],["(-)","(N/m2)","(N/m2)","(kg/m3)","(m2)","(m2)","(m2)","(m4)","(m4)","(m4)"]])

    #MEMBER COSINE MATRICES COSM(i,j)
    myjckt.COSMHeader=np.array([["COSMID","COSM11","COSMID12","COSMID13","COSMID21","COSMID22","COSMID23","COSMID31","COSMID32","COSMID33"],["(-)","(-)","(-)","(-)","(-)","(-)","(-)","(-)","(-)","(-)"]])
    myjckt.COSMs=np.array([])

    #JOINT ADDITIONAL CONCENTRATED MASSES
    myjckt.CmassHeader=np.array([["CMJointID","JMass","JMXX","JMYY","JMZZ"],["(-)","(kg)" ,"(kg*m^2)","(kg*m^2)","(kg*m^2)"]])
    myjckt.Cmass=np.array([])

    #OUTPUT: SUMMARY & OUTFILE
    myjckt.SSSum=np.array([True, "SSSum", "- Output a Summary File (flag).It contains: matrices K,M  and C-B reduced M_BB, M-BM, K_BB, K_MM(OMG^2), PHI_R, PHI_L. It can also contain COSMs if requested."])
    myjckt.OutCOSM=np.array([True, "OutCOSM", "- Output cosine matrices with the selected output member forces (flag)"])
    myjckt.OutAll=np.array([True, "OutAll", "- [T/F] Output all members' end forces "])
    myjckt.OutSwtch=np.array([1, "OutSwtch", "- [1/2/3] Output requested channels to: 1=<rootname>.SD.out;  2=<rootname>.out (generated by FAST);  3=both files."])
    myjckt.TabDelim=np.array([True, "TabDelim", "- Generate a tab-delimited output in the <rootname>.SD.out file"])
    myjckt.OutDec=np.array([1, "OutDec", "- Decimation of output in the <rootname>.SD.out file"])
    myjckt.OutFmt=np.array(["Es11.4e2", "OutFmt", "- Output format for numerical results in the <rootname>.SD.out file"])
    myjckt.OutSFmt=np.array(["A11", "OutFmt", "- Output format for header strings in the <rootname>.SD.out file"])

    #MEMBER OUTPUT LIST
    myjckt.MemOutListHeader=np.array([['MemberID','NoutCnt','NodeCnt','[NOutCnt=how many nodes to get output for [< 10]; NodeCnt are local ordinal numbers from the start of the member, and must be >=1 and <= NDiv+1] If NMOutputs=0 leave blank as well.]'],\
                                           ['(-)','(-)','(-)']])

    #SSOutline
    myjckt.SSOutlist=np.array([["ReactFXss, ReactFYss, ReactFZss, ReactMXss, ReactMYss, ReactMZss",'-Base reactions (forces onto SS structure)'],\
                    ["IntfFXss,  IntfFYss,  IntfFZss,  IntfMXss, IntfMYss, IntfMZss",'-Interface reactions (forces from SS structure)'],\
                    ["IntfTDXss,  IntfTDYss,  IntfTDZss,  IntfRDXss, IntfRDYss, IntfRDZss",'-Interface deflections '],\
                    ["IntfTAXss,  IntfTAYss,  IntfTAZss,  IntfRAXss, IntfRAYss, IntfRAZss",'Interface accelerations']])

    myjckt.SDjoints=np.array([[1, -5.93900, -5.93900, -43.12700],\
                            [2, 5.93900, -5.93900, -43.12700],\
                            [3, 5.93900, 5.93900, -43.12700],\
                            [4, -5.93900, 5.93900, -43.12700],\
                            [5, -4.01600, -4.01600, 15.65100],\
                            [6, 4.01600, -4.01600, 15.65100],\
                            [7, 4.01600, 4.01600, 15.65100],\
                            [8, -4.01600, 4.01600, 15.65100],\
                            [9, 0.00000, -4.79180, -8.06090],\
                            [10, 4.79180, 0.00000, -8.06090],\
                            [11, 0.00000, 4.79180, -8.06090],\
                            [12, -4.79180, 0.00000, -8.06090]])

    myjckt.BaseRxnJoints=np.array([[1,1,1,1,1,1,1],\
                                 [2,1,1,1,1,1,1],\
                                 [3,1,1,1,1,1,1],\
                                 [4,1,1,1,1,1,1]])

    myjckt.InterfaceJointsFlags=np.array([[5,1,1,1,1,1,1],\
                                   [6,1,1,1,1,1,1],\
                                   [7,1,1,1,1,1,1],\
                                   [8,1,1,1,1,1,1]])

    myjckt.Members=np.array([[1, 1, 5, 2, 2],\
                            [2, 2, 6, 2, 2],\
                            [3, 3, 7, 2, 2],\
                            [4, 4, 8, 2, 2],\
                            [5, 1, 9, 3, 3],\
                            [6, 9, 6, 3, 3],\
                            [7, 5, 9, 3, 3],\
                            [8, 9, 2, 3, 3],\
                            [9, 2, 10, 3, 3],\
                            [10, 10, 7, 3, 3],\
                            [11, 6, 10, 3, 3],\
                            [12, 10, 3, 3, 3],\
                            [13, 3, 11, 3, 3],\
                            [14, 11, 8, 3, 3],\
                            [15, 7, 11, 3, 3],\
                            [16, 11, 4, 3, 3],\
                            [17, 4, 12, 3, 3],\
                            [18, 12, 5, 3, 3],\
                            [19, 8, 12, 3, 3],\
                            [20, 12, 1, 3, 3],\
                            [21, 5, 6, 3, 3],\
                            [22, 6, 7, 3, 3],\
                            [23, 7, 8, 3, 3],\
                            [24, 8, 5, 3, 3]])

    myjckt.MemOutList=np.array([[1,2,1,2]])

    #DRIVER INPUTS--------------------------------------------------------------

    myjckt.EchoD=np.array([True, "Echo", "- Echo the input file data (flag)"])

    #Environmental Conditions
    myjckt.Gravity=np.array([9.81, "Gravity", "- Gravity (m/s^2)"])
    myjckt.WtrDpth=np.array([43.127, "WtrDpth", "- Water Depth (m) positive value"])

    #SubDyn
    myjckt.SDInputFile=np.array([myjckt.InputFile_path, "SDInputFile"])
    myjckt.OutRootName=np.array([str(myjckt.InputandDriverpath)+os.sep+str(Base_name), "OutRootName"])
    myjckt.NSteps=np.array([600, "NSteps", "- Number of time steps in the simulations (-)"])
    myjckt.TimeInterval=np.array([0.005, "TimeInterval", "- TimeInterval for the simulation (sec)"])
    myjckt.TP_RefPoint=np.array([0.0, 0.0, 18.15, "TP_RefPoint", "- Location of the TP reference point in global coordinates (m)"])
    myjckt.SubRotateZ=np.array([0.0, "SubRotateZ", "- Rotation angle of the structure geometry in degrees about the global Z axis."])

    #INPUTS
    myjckt.InputsMod=np.array([1, "InputsMod", "- Inputs model {0: all inputs are zero for every timestep, 1: steadystate inputs, 2: read inputs from a file (InputsFile)} (switch)"])
    myjckt.InputsFile=np.array(['""', "InputsFile", "- Name of the inputs file if InputsMod = 2"])

    #STEADY INPUTS
    myjckt.uTPInSteady=np.array([0.1, 0.0, 0.0, 0.0, 0.0, 0.0, "uTPInSteady", "- input displacements and rotations ( m, rads )"])
    myjckt.uDotTPInSteady=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "uDotTPInSteady", "- input translational and rotational velocities ( m/s, rads/s)"])
    myjckt.uDotDotTPInSteady=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "uDotDotTPInSteady", "- input translational and rotational accelerations ( m/s^2, rads/s^2)"])

    #myjckt.SDpySubDynA.run()
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
    #print 'flag1 ',myjckt.Waterinputs.wlevel #CJB+

    #--- RUN JACKET ---#
    myjckt.run()
    #print 'flag2 ',myjckt.Waterinputs.wlevel #CJB+
    # ---------------

    #_____________________________________#
    #Now show results of modal analysis
    print('First two Freqs.= {:5.4f} and {:5.4f} Hz'.format(*myjckt.LoadFrameOuts.Frameouts.Freqs))
    #print component masses
    print('jacket+TP(structural+lumped) mass (no tower, no piles) [kg] = {:6.0f}'.format(myjckt.LoadFrameOuts.Frameouts.mass[0]+myjckt.TP.TPlumpinputs.mass-myjckt.Tower.Twrouts.mass))
    print('tower mass [kg] = {:6.0f}'.format(myjckt.Tower.Twrouts.mass))
    print('TP mass structural + lumped mass [kg] = {:6.0f}'.format(myjckt.TP.TPouts.mass+myjckt.TP.TPlumpinputs.mass))
    print('piles (all) mass (for assigned (not optimum, unless optimization is run) Lp [kg] = {:6.0f}'.format(myjckt.Mpiles))
    print('frame3dd model mass (structural + TP lumped) [kg] = {:6.0f}'.format(myjckt.LoadFrameOuts.Frameouts.mass[0]+myjckt.TP.TPlumpinputs.mass))
    print('frame3dd model mass (structural + TP lumped) + Pile Mass [kg] = {:6.0f}'.format(myjckt.LoadFrameOuts.Frameouts.mass[0]+myjckt.TP.TPlumpinputs.mass+myjckt.Mpiles))
    print('frame3dd model mass (structural only) no piles no tower [kg] = {:6.0f}'.format(myjckt.LoadFrameOuts.Frameouts.mass[0]-myjckt.Tower.Twrouts.mass))
    #print tower top displacement
    print('Tower Top Displacement in Global Coordinate System [m] ={:5.4f}'.format(*myjckt.LoadFrameOuts.Frameouts.top_deflection))
    #print max API code checks
    print('MAX member compression-bending utilization at joints = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts.jacket_utilization.cb_util)))
    print('MAX member tension utilization at joints = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts.jacket_utilization.t_util)))
    print('MAX X-joint  utilization at joints = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts.jacket_utilization.XjntUtil)))
    print('MAX K-joint  utilization at joints = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts.jacket_utilization.KjntUtil)))

    PlotJacket(myjckt,util=True)
