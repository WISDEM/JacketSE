#-------------------------------------------------------------------------------
# Name:        RunAllJacketTwrs.py
# Purpose:     This program calculates the matrix of experiments for WE14CA01- LCOE_Maureen's Task
#
# Author:      rdamiani
#
# Created:     17/03/2014
# Copyright:   (c) rdamiani 2014
# Licence:     <Apache>
#-------------------------------------------------------------------------------
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import time
import os
import SetJacketInputs
import cProfile
import collections

def main(matfilename,smoothfilename,datasaved=False,):
    """Launch this main after modifying some of the non-optimization parameters (situational parameters) below,\n
        some of the optimization bounds below, setting OTHER PARAMETERS in SetJacketInputs.py!!! (this is for smoothness calc), \n
        and with the following argument Inputs:\n
        matfilename    - string, name of text file that will contain the matrix of experiments. Set to [] to bypass this step and go to smoothness calcs.\n
        smoothfilename - string, name of text file that will contain the smoothness calc data. \n
        datasaved- boolean, set to True if for-loop has already been done for the smoothness calc, and all you want are graphs.\n
        """

    #_____Inputs start_____#
    TurbineMW=np.array([3.,6.,10.])  #MW ratings
    nRs=TurbineMW.size
    baseHH=np.array([80.,103.5,126.5]) #hub-height [m] array for the various turbines
    baseThrust=np.array([507.e3,1014.e3,1687.e3])    #Max unfactored thrusts for the various turbines [N] (contains the 1.5 amplification factor but no IEC PSF yet)
    baseTorque=np.array([2303.e3,6522.e3,13995.e3])  #Max unfactored thrusts for the various turbines [Nm]
    baseRNAmass=np.array([234.e3,500.e3,900.e3])  #RNA masses for the various turbines [kg]
    baseFreq=np.array([0.26,0.21,0.21]) #Target 1st EigenFrequency of the entire system

    #Optimization parameter bounds

    baseBattermin=np.array([8.,8.,7.])    #minimum batter
    baseBattermax=np.array([30.,30.,30.]) #Maximum batter
    batter_guess=np.array([15.,10., 8.])#guess on batter

    baseDpmin=np.array([0.8,0.8,0.8])    #minimum pile OD
    baseDpmax=np.array([2.,3.,4.])   #Maximum pile OD
    Dp_guess=np.array([1.,1.8,1.8])   #guess on Dp

    basetpmin=np.array([1.5,1.5,1.5])*0.0254    #minimum pile t
    basetpmax=np.array([3.,4.,5.])*0.0254   #Maximum pile t
    tp_guess=np.array([1.5,1.5,1.5])*0.0254   #guess on tp

    baseLpmin=np.array([25.,30.,35.])    #minimum pile Length
    baseLpmax=np.array([50.,60.,70.])   #Maximum pile Length
    Lp_guess=np.array([45.,50.,60.])   #guess on Lp

    baseDlegmin=np.array([0.8,0.8,1.0])
    baseDlegmax=np.array([2.0,3.0,4.0])
    Dleg_guess=np.array([0.9,1.1,1.3])   #guess on Dleg
    basetlegmin=np.array([0.5,0.5,0.5])*0.0254
    basetlegmax=np.array([2.,2.5,3.])*0.0254
    tleg_guess=np.array([0.03,0.05,0.07])   #guess on tleg

    baseDbrcmin=np.array([0.4,0.8,0.8])
    baseDbrcmax=baseDlegmax.copy()
    Dbrc_guess=np.array([0.6,0.9,1.0])   #guess on Dbrc
    basetbrcmin=np.array([0.5,0.5,0.5])*0.0254
    basetbrcmax=basetlegmax.copy()
    tbrc_guess=np.array([0.03,0.03,0.05])   #guess on tbrc

    baseDmudmin=np.array([0.6,0.8,1.0])
    baseDmudmax=baseDlegmax.copy()
    Dmud_guess=np.array([1.0,1.0,1.6])   #guess on Dmud
    basetmudmin=np.array([0.5,0.5,0.5])*0.0254
    basetmudmax=basetlegmax.copy()
    tmud_guess=np.array([0.03,0.04,0.05])   #guess on tmud

    baseDgirmin=np.array([0.8,0.8,1.0])
    baseDgirmax=np.array([1.0,1.1,1.3])
    Dgir_guess=np.array([1.0,1.0,1.3])   #guess on Dgir
    basetgirmin=np.array([0.5,0.5,0.5])*0.0254
    basetgirmax=np.array([2.,2.5,2.5])*0.0254
    tgir_guess=np.array([0.03,0.03,0.05])   #guess on tgir

    baseDbmin=np.array([2.5,3.5,4.5]) #Minimum tower base ODs for the various turbines
    baseDbmax=np.array([5.0,6.5,8.]) #Maximum tower base ODs for the various turbines
    Db_guess=np.array([3.0,5.,7.5])   #guess on Db

    baseDtmin=np.array([2.0,3.0,3.5]) #Minimum tower top ODs for the various turbines
    baseDtmax=np.array([3.0,4.,4.5]) #Maximum tower top ODs for the various turbines
    Dt_guess=np.array([2.0,3.5,4.5])   #guess on Dt

    baseDTRbmin=np.array([120.,120.,120.]) #minimum DTR for the tower
    baseDTRbmax=np.array([200.,200.,200.]) #Maximum DTR for the tower
    DTRb_guess=np.array([150.,140.,120.])   #guess on DTRb
    baseDTRtmin=np.array([120.,120.,120.]) #minimum DTR for the tower
    baseDTRtmax=np.array([200.,200.,200.]) #Maximum DTR for the tower
    DTRt_guess=DTRb_guess.copy()   #guess on DTRt

    baseHtwr2min=np.array([0.05,0.05,0.05]) #minimum DTR for the tower  #0.005 for numerics
    baseHtwr2max=np.array([0.25,0.25,0.25]) #Maximum DTR for the tower
    H2_guess=np.array([0.05,0.15,0.25])   #guess on H2twrfrac


    #wave load situations- see excel spreadsheet
    nws=3  #Number of water depths and associated 50-yr RP wave heights, and periods, and deck heights
    wdepths=np.linspace(20.,60.,nws)
    wH50=np.linspace(15.,30.,nws)  #[m]
    Tp50=np.linspace(10.1,14.2,nws)  #[sec]

    dckhgt=np.array([14.25,19.125,24.]) #[m] bottom of deck MSL

    baseftprint=np.array([30.,32.,35.]) #[m] at seabed

    legndiv=np.array([4,4,4]) # number of elements per member in the leg so that MSL level is close to one of them, and always above, no jumping
    nbays=np.array([4,4,5]) #Number of bays function of water depth

    U50HH=np.linspace(30.,30.,nRs)  #[m/s] Speed associated with load case assumed (ULS, 1.6 case)
    TPlumpmass=np.linspace(100.e3,200.e3,nRs)  #[kg] Lump mass at the TP in addition to TP with load case assumed (ULS, 1.6 case)
    CMzoff=np.linspace(2.,4.,nRs)  #[m] RNA CM distance from tower top
    Thzoff=CMzoff.copy()  #[m] Thrust application point distance from tower top

    nHs=3 #Number of hub-height factors
    HHs=np.linspace(1.0,1.3,nHs)

    nTs=1 #number of thrust factors
    maxThrusts=np.linspace(1,1.,nTs)

    nMs=3 #Number of tower top masses to consider including the base one for the size
    RNAmasses=np.linspace(0.8,1.2,nMs)


    #soils
    #Soils=('Soft','Stiff')
    Soils=('Average',)  #leave the comma for 1 element list
    nS=len(Soils)
    #_____Inputs end_____#

    nSims=nws*nHs*nTs*nMs*nS   #Total number of simulations

    if (isinstance(matfilename,str) and (os.path.exists(matfilename) or  os.access(os.path.dirname(matfilename), os.W_OK))):
       outfile=open(matfilename,'w')
       #outfile2=open(os.path.splitext(matfilename)[0]+'guesses.dat','w')

    # write header
       outfile.write('Table of Experiments. Input Generated via RunAllJacketTwrs.py on {:s} \n\n'.format(time.strftime("%c")));
       outfile.write(('{:6d}   {:12d}'+54*' {:>12d}'+18*' {:>12d}'+'\n').format(*range(0,57+18)))
       outfile.write(('{:6s}   {:12s}'+' {:>10s}'+ 53*' {:>11s} '+ 18*' {:>11s} '+'\n').\
           format('Case#','CaseName','Turbine_Rating','Water_Depth','Hmax_50yr','Tp50yr','U50HH','dck_botz','Hub-Height','TPlumpmass',\
           'RNA_mass','Thrust','Torque','CMzoff','Thzoff','EigenFreq.','Soil', 'legndiv','nbays','maxftprint',\
            'batter_min','batter_max','Dp_min','Dp_max','tp_min','tp_max','Lp_min','Lp_max',\
            'Dleg_min','Dleg_max','tleg_min','tleg_max','Dbrc_min','Dbrc_max','tbrc_min','tbrc_max',\
            'Dmud_min','Dmud_max','tmud_min','tmud_max','Dgir_min','Dgir_max','tgir_min','tgir_max',\
            'Db_min','Db_max','DTRb_min','DTRb_max','Dt_min','Dt_max','DTRt_min','DTRt_max','Htwr2min','Htwr2max','dckwdth_factmin','dckwdth_factmax',\
            #guesses
            'batter_g','Dp_g','tp_g','Lp_g',\
            'Dleg_g','tleg_g','Dbrc_g','tbrc_g',\
            'Dmud_g','tmud_g','Dgir_g','tgir_g',\
            'Db_g','DTRb_g','Dt_g','DTRt_g','Htwr2_g','dckwdth_factg'\

            ))

       outfile.write(('{:6s}   {:12s}'+54*' {:>12s}'+18*' {:>12s}'+'\n').format('[-]','[-]','[MW]','[m]','[m]','[s]','[m/s]','[m]','[m]','[kg]',\
           '[kg]',     '[N]',     '[Nm]',     '[m]',      '[m]',   '[Hz]',    '[-]',     '[-]',    '[-]', '[-]',\
            '[-]',     '[-]',     '[m]',      '[m]',      '[m]',   '[m]',     '[m]',     '[m]',\
            '[m]',     '[m]',     '[m]',      '[m]',      '[m]',   '[m]',     '[m]',     '[m]',\
            '[m]',     '[m]',     '[m]',      '[m]',      '[m]',   '[m]',     '[m]',     '[m]',\
            '[m]',     '[m]',     '[-]',      '[-]',      '[m]',   '[m]',     '[-]',     '[-]',    '[-]', '[-]','[-]',     '[-]',\
            #guesses
            '[-]',     '[m]',      '[m]',      '[m]',\
            '[m]',     '[m]',      '[m]',      '[m]',\
            '[m]',     '[m]',      '[m]',      '[m]',\
            '[m]',     '[-]',      '[m]',      '[-]',      '[-]',   '[-]'))

    #write values
       for uu,rating in enumerate(TurbineMW):
            HH0=baseHH[uu]

            TT0=baseThrust[uu]
            Tq0=baseTorque[uu]
            MM0=baseRNAmass[uu]
            F0=baseFreq[uu]
            U50HH0=U50HH[uu]


            battermin=3.   #Hardcode this onebut if water depth is max and turbine is the heavy one let it go to 7.
            battermax=baseBattermax[uu]

            Dpmin=baseDpmin[uu]
            tpmin=basetpmin[uu]
            Lpmin=baseLpmin[uu]

            Dlegmin=baseDlegmin[uu]
            tlegmin=basetlegmin[uu]

            Dbrcmin=baseDbrcmin[uu]
            tbrcmin=basetbrcmin[uu]
            Dmudmin=baseDmudmin[uu]
            tmudmin=basetmudmin[uu]

            Dgirmin=baseDgirmin[uu]
            tgirmin=basetgirmin[uu]


            Dbmin=baseDbmin[uu]
            Dtmin=baseDtmin[uu]
            DTRbmin=baseDTRbmin[uu]
            DTRtmin=baseDTRtmin[uu]

            Dbmax=baseDbmax[uu]

            Dtmax=baseDtmax[uu]
            DTRbmax=baseDTRbmax[uu]
            DTRtmax=baseDTRtmax[uu]

            Htwr2min=baseHtwr2min[uu]
            Htwr2max=baseHtwr2max[uu]




            for ii,WD in enumerate(wdepths):
                if (uu==2):#MAX TURBINE
                    battermin=3.

                maxftprint=min(baseftprint[uu],baseftprint[ii])
                if (uu==1) and (ii==1):  #40m wd and 6 MW
                    maxftprint=30.

                Dpmax=max(baseDpmax[uu],baseDpmax[ii]) #select maximum bounds based on the water depths or turbine capacity
                Dp_g=min(Dp_guess[uu],Dp_guess[ii])
                tpmax=max(basetpmax[uu],basetpmax[ii]) #select maximum bounds based on the water depths or turbine capacity
                Dp_g=Dp_guess[ii]
                tp_g=min(tp_guess[uu],tp_guess[ii])

                Dlegmax=max(baseDlegmax[uu],baseDlegmax[ii]) #select maximum bounds based on the water depths or turbine capacity

                Dleg_g=Dleg_guess[uu]
                tleg_g=min(tleg_guess[uu],tleg_guess[ii])

                if ii==2:
                    Dleg_g=max(Dleg_guess[ii],Dleg_guess[uu])
                    tleg_g=max(tleg_g,tleg_guess[ii])

                tlegmax=basetlegmax[ii] #select maximum bounds based on the water depths


                Dbrcmax=max(baseDbrcmax[uu],baseDbrcmax[ii]) #select maximum bounds based on the water depths or turbine capacity
                Dbrc_g=min(Dbrc_guess[uu],Dbrc_guess[ii])
                tbrcmax=max(basetbrcmax[uu],basetbrcmax[ii]) #select maximum bounds based on the water depths or turbine capacity
                tbrc_g=min(tbrc_guess[uu],tbrc_guess[ii])

                Dmudmax=max(baseDmudmax[uu],baseDmudmax[ii]) #select maximum bounds based on the water depths or turbine capacity
                tmudmax=max(basetmudmax[uu],basetmudmax[ii]) #select maximum bounds based on the water depths or turbine capacity
                Dmud_g=min(Dmud_guess[uu],Dmud_guess[ii])
                tmud_g=min(tmud_guess[uu],tmud_guess[ii])


                Dgirmax=max(baseDgirmax[uu],baseDgirmax[ii]) #select maximum bounds based on the water depths or turbine capacity
                Dgir_g=min(Dgir_guess[uu],Dgir_guess[ii])

                tgirmax=max(basetgirmax[uu],basetgirmax[ii]) #select maximum bounds based on the water depths or turbine capacity
                tgir_g=min(tgir_guess[uu],tgir_guess[ii])

                Lpmax= baseLpmax[ii]#select maximum bounds based on the water depths and turbine size
                if (ii>=1):
                    if (uu==0):
                        Lpmax=50. #leave at 50 for wd<=40m
                        if (ii==2):
                            Lpmax=65.
                if (uu==2): #10MW
                     Lpmax=60.  #set it at 60 for wd<=20
                     if (ii>=1):
                        Lpmax=70.


                Lp_g=Lp_guess[ii]

                Tplump=TPlumpmass[ii]
                if uu>0:
                    Tplump=max(Tplump,TPlumpmass[1])  #For larger turbine sizes I start from the mid TP mass


                for jj,HH in enumerate(HHs):
                    if (uu==2):#MAX TURBINE AND HUBHEIGHT
                       Dbmax=9.
                       if jj==2:
                            Dbmax=9. #set to 9 for 2nd round of 10 MW calculations
                            if ii==2:
                                maxftprint=40.
                    #batter_g=min(batter_guess[uu],batter_guess[jj])
                    #batter_g=min(batter_g,batter_guess[ii]) #if water depth is greater than first value start with second value of batter_guess
                    batter_g=batter_guess[jj]
                    if ii==2: #deepest water
                        batter_g=min(12,batter_g)

                    Db_g=min(Db_guess[uu],Db_guess[jj])
                    #Deck width
                    dckwdth_factmin=2.
                    dckwdth_factmax=3.
                    dckwdth_factg=2.

                    Dt_g=min(Dt_guess[uu],Dt_guess[jj])
                    DTRb_g=min(DTRb_guess[uu],DTRb_guess[jj])
                    DTRt_g=min(DTRt_guess[uu],DTRt_guess[jj])
                    Htwr2_g=H2_guess[jj]

                    for kk,TT in enumerate(maxThrusts):
                        for ll,MM in enumerate(RNAmasses):
                            for mm,SS in enumerate(Soils):
                                counter=mm+1+nS*(ll+nMs*(kk+nTs*(jj+nHs*(ii+nws*uu))))
                                outfile.write(('{:6d} R{:d}D{:d}H{:d}M{:d}T{:d}S{:d} '+14*'{:>12.2f} '+' {:>12s} '+ 2*' {:>12d} '+' {:12.4f}'+36*' {:12.4f}'+18*' {:12.4f}'+'\n').\
                                        format(counter,uu,ii,jj,ll,kk,mm,rating,WD,wH50[ii],Tp50[ii],U50HH0,dckhgt[ii], HH*HH0, Tplump,\
                                                MM*MM0, TT*TT0,Tq0, CMzoff[uu], Thzoff[uu], F0, SS, legndiv[ii], nbays[ii],maxftprint,\
                                                battermin,battermax,Dpmin,Dpmax,tpmin,tpmax,Lpmin,Lpmax,\
                                                Dlegmin,Dlegmax,tlegmin,tlegmax,Dbrcmin,Dbrcmax,tbrcmin,tbrcmax,\
                                                Dmudmin,Dmudmax,tmudmin,tmudmax,Dgirmin,Dgirmax,tgirmin,tgirmax,\
                                                Dbmin,Dbmax,DTRbmin,DTRbmax,Dtmin,Dtmax,DTRtmin,DTRtmax,Htwr2min,Htwr2max,dckwdth_factmin,dckwdth_factmax,\
                                                batter_g,Dp_g,tp_g,Lp_g,\
                                                Dleg_g,tleg_g,Dbrc_g,tbrc_g,Dmud_g,tmud_g,Dgir_g,tgir_g,\
                                                Db_g,DTRb_g,Dt_g,DTRt_g,Htwr2_g,dckwdth_factg) )








       outfile.close()
       print('File {:s} written'.format(matfilename))






    else:
        #____________________SMOOTHNESS CALCULATION STARTS HERE______________________#
       print('No File written, but test of smooth behavior in progress')

        #Now let us check how well behaved this code is for the optimizer, pick
        #one configuration and span the optimdi0ization parameters

       cfg_ID=2 #pick 1 2 or 3 as we have 3 basic configurations
        #length of the various optimiztion spaces
       ngeneric=20 #size of each variable
       mg=10  #characteristic value index (to pick something out of each variable)
       nbatters=ngeneric
       nDbs=ngeneric
       nDts=ngeneric
       nDTRs=ngeneric
       nH2s=ngeneric
       nDps=ngeneric  #pile ODs
       ntps=ngeneric  #pile ts
       nDgs=ngeneric  #leg ODs
       ntgs=ngeneric  #leg ts
       nLps=ngeneric  #Pile Lps

       totEx=nDTRs*nDbs*nH2s*nDts*nbatters*nDps*ntps*nDgs*ntgs #total number of simulations

       batter=np.linspace(baseBattermin[cfg_ID],baseBattermax[cfg_ID],nbatters)

       # DEBUG: batter=np.linspace(9.8,10.3,20.)  #checking what is going on with cbutil

       Dpile= np.linspace(1.0,2.,nDps)
       tpile= np.linspace(0.025,0.075,ntps)
       Lp= np.linspace(35.,45,nLps)
       Dleg= np.linspace(1.4,2.,nDgs)
       tleg= np.linspace(0.025,0.075,ntgs)
       Dbrc= Dleg.copy()
       tbrc= tleg.copy()
       Dbrc_mud=Dbrc.copy()
       tbrc_mud=tbrc.copy()

       Db= np.linspace(baseDbmin[cfg_ID],baseDbmax[cfg_ID],nDbs)
       Dt=np.linspace(baseDtmin[cfg_ID],baseDtmax[cfg_ID],nDts)
       DTR=np.linspace(baseDTRbmin[cfg_ID],baseDTRbmax[cfg_ID],nDTRs)
       DTRb=DTR.copy()
       DTRt=DTR.copy()
       Htwr2frac=np.linspace(0.,0.25,nH2s)

       #Put all the optimization variables here

       allvars=collections.OrderedDict([('batter',batter),('Dleg',Dleg),('tleg',tleg),('Dbrc_mud',Dbrc_mud),('tbrc_mud',tbrc_mud),('Dbrc',Dbrc),('tbrc',tbrc),\
                                        ('Db',Db),('Dt',Dt),('DTRb',DTRb),('DTRt',DTRt),('Htwr2frac',Htwr2frac),('Dpile',Dpile),('tpile',tpile)])
       nvars=allvars.__len__()

       #Initialization
       #max_GLUtil=np.zeros([ngeneric,nDbs,nDts,nDTRs,nH2s,nDps,ntps,nDgs,ntgs])*np.nan
       max_GLUtil=np.zeros([ngeneric,13,nvars])
       max_EUUtil=max_GLUtil.copy()
       Lp0rat=np.zeros([ngeneric,nvars])
       max_cbutil=np.zeros([ngeneric,428,nvars])
       max_tutil=np.zeros([ngeneric,856,nvars])
       max_KjntUtil=np.zeros([ngeneric,20,nvars])
       max_XjntUtil=np.zeros([ngeneric,16,nvars])

       mass_str=Lp0rat.copy()
       mass_tot=Lp0rat.copy()
       f1=Lp0rat.copy()
       f2=Lp0rat.copy()
       wbase=Lp0rat.copy()
       Mpiles=Lp0rat.copy()

       cex=0 #execution counter

       if not(datasaved):
            pathname=os.path.dirname(smoothfilename)
           #return myjckt  #uncomment this for profiling

            keyno=0   #counter for variables to span
            for key,value in allvars.iteritems():
             myjckt=SetJacketInputs.main(15.,Dpile[mg],tpile[mg],Lp[mg],Dleg[mg],tleg[mg],Dbrc[mg],tbrc[mg],\
                                            Dbrc_mud[mg],tbrc_mud,Db[mg],DTRb[mg],Dt[mg],DTRt[mg],Htwr2frac[mg])

             for ib in range(value.size):
                    if key=='batter':
                        myjckt.JcktGeoIn.batter=value[ib]
                    elif key =='Dleg' or key=='tleg':
                        setattr(myjckt.leginputs,key,np.array([value[ib]]).repeat(myjckt.JcktGeoIn.nbays+1))
                    elif key in ['Db','Dt','DTRb','DTRt','Htwr2frac']:
                        setattr(myjckt.Twrinputs,key,value[ib])
                    elif key in ['Dpile','tpile','Lp']:
                        setattr(myjckt.Pileinputs,key,value[ib])
                    elif key in ['Dbrc_mud','tbrc_mud']:
                        setattr(myjckt.Mbrcinputs,key,value[ib])
                    elif key in ['Dbrc','tbrc']:
                        setattr(myjckt.Xbrcinputs,key,np.array([value[ib]]).repeat(myjckt.JcktGeoIn.nbays))

                            #
                            #myjckt=SetJacketInputs.main(batter[ib],Dp[nDp],tp[otp],Lp,Dleg[pDg],tleg[qtg],Db[jD],DTRb[lDTR],Dt[kDt],DTRt[lDTR],Htwr2frac[mH2])
                    myjckt.run()


                    f1[ib,keyno],f2[ib,keyno]=myjckt.Frameouts2.Freqs
                    max_GLUtil[ib,0:myjckt.tower_utilization.GLUtil.size,keyno]=myjckt.tower_utilization.GLUtil
                    max_EUUtil[ib,0:myjckt.tower_utilization.EUshUtil.size,keyno]=myjckt.tower_utilization.EUshUtil
                    max_cbutil[ib,:,keyno]=myjckt.jacket_utilization.cb_util
                    max_tutil[ib,:,keyno]= myjckt.jacket_utilization.t_util
                    max_KjntUtil[ib,:,keyno]=myjckt.jacket_utilization.KjntUtil
                    max_XjntUtil[ib,:,keyno]=myjckt.jacket_utilization.XjntUtil

                    Lp0rat[ib,keyno]=myjckt.Lp0rat
                    wbase[ib,keyno]=myjckt.wbase
                    mass_str[ib,keyno],mass_tot[ib,keyno]=myjckt.Frameouts2.mass
                    Mpiles[ib,keyno]=myjckt.Mpiles

                    cex+=1
                    print 'case no. {:d} of {:d} for var {:s}, this is var number {:d} of {:d}'.format(ib,value.size,key,keyno,nvars)
                    #print "Case no. {:d} of {:d}, batter ={:5.2f}, Db = {:5.2f} m, Dt= {:5.4f} m, DTR={:5.2f}, H2frac={:5.2f} \
                    #              ".format(cex,totEx,batter[ib],Db[jD],Dt[kDt],DTR[lDTR],Htwr2frac[mH2])

             plotting(key,value,pathname,mass_str[:,keyno],f1[:,keyno],max_cbutil[:,:,keyno],max_tutil[:,:,keyno],max_KjntUtil[:,:,keyno],max_XjntUtil[:,:,keyno],Lp0rat[:,keyno])
             keyno +=1

                #sio.savemat(os.path.splitext(OutFrame3DD)[0]+'.mat',{'batter':batter, 'Dleg':Dleg,\
            sio.savemat(smoothfilename,{'batter':batter,\
                            'Dpile':Dpile,'tpile':tpile,'Lp':Lp,'Lp0rat':Lp0rat,'Mpiles':Mpiles,\
                            'Dleg':Dleg,'tleg':tleg,'Db':Db,'Dt':Dt,'DTR':DTR,'Htwr2frac':Htwr2frac,\
                            'f1':f1,'f2':f2,'mass_str':mass_str,'mass_tot':mass_tot,\
                            'max_cbutil':max_cbutil,'max_tutil':max_tutil,\
                            'max_XjntUtil':max_XjntUtil,'max_KjntUtil':max_KjntUtil,\
                            'max_GLUtil':max_GLUtil,'max_EUUtil':max_EUUtil,'wbase':wbase})

       else:       #retrieve data
            datfile=smoothfilename
            a=sio.loadmat(datfile)
            Db=a.get('Db').squeeze()
            Dt=a.get('Dt').squeeze()
            DTR=a.get('DTR').squeeze()
            Dpile=a.get('Dpile').squeeze()
            tpile=a.get('tpile').squeeze()
            Dleg=a.get('Dleg').squeeze()
            tleg=a.get('tleg').squeeze()

            f1=a.get('f1').squeeze()
            f2=a.get('f2').squeeze()
            mass_str=a.get('mass_str').squeeze()
            mass_tot=a.get('mass_tot').squeeze()
            Mpiles=a.get('Mpiles').squeeze()

            max_GLUtil=a.get('max_GLUtil').squeeze()
            max_EUUtil=a.get('max_EUUtil').squeeze()
            max_XjntUtil=a.get('max_XjntUtil').squeeze()
            max_KjntUtil=a.get('max_KjntUtil').squeeze()
            max_cbutil=a.get('max_cbutil').squeeze()
            max_tutil=a.get('max_tutil').squeeze()

            Lp=a.get('Lp').squeeze()
            Lp0rat=a.get('Lp0rat').squeeze()

       #To plot max util vs Db,tb,Dt,tt,Dp,Lp,tp,batter,


def plotting(key,value,pathname,mass_str,f1,max_cbutil,max_tutil,max_KjntUtil,max_XjntUtil,Lp0rat):
    """sub to plot and save figs"""
    fig=plt.figure();


    #Plot mass
    ax1 = fig.add_subplot(711)
    ax1.plot(value,mass_str,'.-')
    ax1.set_ylabel('MASS')
    #Plot mass
    ax2 = fig.add_subplot(712)
    ax2.plot(value,f1,'.-')
    ax2.set_ylabel('f1')
    #cbutil
    ax3 = fig.add_subplot(713)
    ax3.plot(value,max_cbutil[:,40:80],'.-')
    ax3.set_ylabel('max_cbutil')
    #cbutil, tutil, Xutil,Kutil
    ax4 = fig.add_subplot(714)
    ax4.plot(value,max_tutil[:,10:15],'.-')
    ax4.set_ylabel('max_tutil')

    ax5 = fig.add_subplot(715)
    ax5.plot(value,max_KjntUtil[:,10:12],'.-')
    ax5.set_ylabel('max_KjntUtil')

    ax6= fig.add_subplot(716)
    ax6.plot(value,max_XjntUtil[:,10:15],'.-')
    ax6.set_ylabel('max_XjntUtil')

    ax7= fig.add_subplot(717)
    ax7.plot(value,Lp0rat,'.-')
    ax7.set_ylabel('Lp0rat')

    #plt.savefig(pathname+'/'+key+'VsXutil.png',format='png')
    plt.savefig(pathname+'/'+key+'VsMass&util.png',format='png')
    #plt.show()


def test():
    print 'here'
    smoothfilename=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\smoothcheck2.dat'
    datasaved=False
    myjckt=main([],smoothfilename,datasaved=datasaved)
    myjckt.JcktGeoIn.batter=8.0
    myjckt.run()


if __name__ == '__main__':
    #__________INPUT STARTS____________#

    matfilename=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\testmatrix.dat'

    matfilename=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\COBYLA\testmatrix9mDb.dat'

    smoothfilename=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\smoothcheck3.dat'
    #smoothfilename=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\smoothcheck2.dat'
    datasaved=False

    #________INPUT ENDS_________#
    #If looking for matrix of input parameters luanch it with a filename
    main(matfilename,smoothfilename,datasaved=datasaved)

    #else test smoothness via
    ##myjckt=main([],smoothfilename,datasaved=datasaved)

   # cProfile.run('test()','stats.txt')
   # import pstats
   # p = pstats.Stats('stats.txt')
   # p.strip_dirs().sort_stats('time').print_stats(10)



##script to span crucial zones
##batter=np.linspace(15.,16.,20.)
##value=batter
##key='batter'
##keyno=0
#then relaunch the for loop with ctrl+F7
#then plot

##fig=plt.figure()
##plt.plot(batter,max_KjntUtil[:,:,0],'.-')
##plt.show()
