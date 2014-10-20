#-------------------------------------------------------------------------------
# Name:        JacketOpt_Peregrine.py
# Purpose:     It Expects to read data from input files and then calls the optimizer based on that.
#              In particular it uses the table of cases for the LCOE analysis project.
#              It started from JacketOpt_OldStyle.py
# Author:      rdamiani
#
# Created:     4/25/2014 based off of JacketOpt.py
# Copyright:   (c) rdamiani 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np
import scipy.optimize

import cPickle as pickle;

import time
import sys
import os

#from PlotJacket import main as PlotJacket  #COMMENT THIS ONE OUT FOR PEREGRINE"S SAKE

from ReadSaveOpt1Line import ReadTab1Line,SaveOpt1Line,ReadOptFile,DesVars

##class Map(dict):
##    """This class is used later to make Structures"""
##    def __init__(self, **kwargs):
##        super(Map, self).__init__(**kwargs)
##        self.__dict__ = self
##
##DesVars=Map()

#__________________________________________________________#


              #  OPTIMIZATION FUNCTION IS HERE #

#__________________________________________________________#

def JcktOpt1line(casefile,caseno,xlsfilename,f0epsilon=0.05,guessfromoutfile=[]):
    """Optimizer which reads one line caseno out of table file casefile:\n
    INPUT \n
        casefile -string, complete path+filename to table file of cases. \n
        caseno   -int, case number (also row out of the table to be read). \n
        xlsfilename -string, complete path+filename to excel file, whose sheet for the current caseno case will be added/updated. \n
    OPTIONALS \n
        f0epsilon -float, (1+f0eps)*f0 will be the upper bound for the f0 constraint. \n
        guessfromoutfile -string, complete path+name to an optimization output file where to read guesses from. \n

    OUTPUT \n
        It will create a final configuration optimized for the design parameters in the caseno row of casefile table. \n
        Also, it will save text and excel file for design variables, and summary for report, respectively. \n
        Figures will also be made of configuration and tower utilization. (not in peregrine where you can then use reconall.py \n
        casename -string, case name for hte current case usable by other programs
        myjckt   -OPenMdao assembly of JacketSE, final configuration after optimization
        \n
        """
    #This is the optimizer that will act on 1 case from a table
    global f0,f0eps, DTRsdiff, mxftprint  #DTRsdiff   -flag, if trues both DTRb and DTRt are design vars, else only DTRb, and DTRt=DTRb
    import SetJacketInputsPeregrine
    import ntpath

    desvars=DesVars() #instanct of design variables
    #First read design parameters and des var bounds
    Desprms,_,desvarbds,desvarmeans,guesses,casename=ReadTab1Line(casefile,caseno,desvars.dvars.keys())
    #USING bds as desvarbds

    #Assign mxftprint
    mxftprint=Desprms.mxftprint

    #Target frequency
    f0=Desprms.f0
    f0eps=f0epsilon

    #Then set up an initial guess at the Design variables and use it also to instantiate the assembly

    for ii,key in enumerate(desvars.dvars):
        #print key  #debug
        setattr(desvars,key,desvarmeans[ii])

    #Then build the initial assembly
    myjckt=SetJacketInputsPeregrine.main(Desprms,desvars)
#   x       0     1            2              3     4    5       6       7        8        9        10    11    12   13      14      15      16      17
#         batter,Dp,           tp,             Lp,   Dleg,tleg,    Dbrc   tbrc    Dmdbrc   tmdbrc     Dgir,tgir    Db, DTRb,   Dt,         H2frac,   dck_width       DTRt=DTRb,
    #guess=[12., desvarmeans[1], desvarbds[2,0],50.,  1.6, 0.0254,  0.8, 0.0254,     1.5,0.03,         1.0, 0.03,  6.5, 120., 3.5,      120., 0.2,      13.]/desvarmeans
    #guess=[10.88,    2.73, 0.037,50.,                    1.5, 0.054,  0.83, 0.027,     1.54,0.036,        1.0, 0.03,   6.8, 124., 3.5,      124., 0.2,  14.]/desvarmeans
    #guess=desvarmeans/desvarmeans
    guess=guesses/desvarmeans  #GUESS FROM TESTMATRIX TABLE
    #OR guess from previous output file
    if guessfromoutfile:
        desvarguess=ReadOptFile(guessfromoutfile)
        for ii,key in enumerate(desvarguess.dvars):
            guess[ii]=getattr(desvarguess,key)/desvarmeans[ii]

    #NOW MAKE SURE YOU RESIZE x for OPTIMIZER IN CASE DTRt=DTRb

    if not(DTRsdiff):#DTRs FIXED TO EACH OTHER
        idx_DTRt=desvars.dvars.keys().index('DTRt')
        guess=np.delete(guess,idx_DTRt)
        desvarmeans=np.delete(desvarmeans,idx_DTRt)
        desvarbds=np.delete(desvarbds,idx_DTRt,0)

    #Finally call the optimizer
##    res1=guess.copy() #debug
##    ms=mass(guess,myjckt,desvarmeans,desvarbds) #debug

        tt = time.time()
        res1=scipy.optimize.fmin_cobyla(mass,guess,[f0Cnstrt1,f0Cnstrt2,batCnsrt1, batCnsrt2, \
               cbCnstrt,tCnstrt,KjntCnstrt,XjntCnstrt, GLCnstrt,EUCnstrt,\
                girDCnsrt1,girDCnsrt2,girtCnsrt1,girtCnsrt2,\
               TwrCnstrt01,TwrCnstrt02,TwrCnstrt03,TwrCnstrt04,TwrCnstrt05,TwrCnstrt08,TwrCnstrt09,\
               XbCrit01,XbCrit02,XbCrit03,XbCrit04,XbCrit05,MbCrit01,MbCrit02,MbCrit03,MbCrit04,MbCrit05,\
               LpCnstrt,LpCnstrt2, Dleg2BrcCnstrt,Dleg2MudCnstrt, DlegCnstrt,tlegCnstrt,tlegCnstrt2,Dleg2tlegCnstrt,DpCnstrt,tpCnstrt,\
               DbrcCnstrt,tbrcCnstrt,tbrcCnstrt2,Dbrc2tbrcCnstrt,DmudCnstrt,tmudCnstrt,tmudCnstrt2,Dmud2mudCnstrt,\
               dckwidthCnstrt1,dckwidthCnstrt2,ftprintCnstrt,NorsokCnstrt],args=[myjckt,desvarmeans,desvarbds],consargs=None,rhobeg=0.01,maxfun=3000,disp=1)

    #Cobyla : If it returns 3 elements , it means it did not converge;  constraint2,constraint3
    else:  #DTRt variable
        tt = time.time()
        res1=scipy.optimize.fmin_cobyla(mass,guess,[f0Cnstrt1,f0Cnstrt2,batCnsrt1, batCnsrt2, \
               cbCnstrt,tCnstrt,KjntCnstrt,XjntCnstrt, GLCnstrt,EUCnstrt,\
                girDCnsrt1,girDCnsrt2,girtCnsrt1,girtCnsrt2,\
               TwrCnstrt01,TwrCnstrt02,TwrCnstrt03,TwrCnstrt04,TwrCnstrt05,TwrCnstrt06,TwrCnstrt07,TwrCnstrt08,TwrCnstrt09,\
               XbCrit01,XbCrit02,XbCrit03,XbCrit04,XbCrit05,MbCrit01,MbCrit02,MbCrit03,MbCrit04,MbCrit05,\
               LpCnstrt,LpCnstrt2, Dleg2BrcCnstrt,Dleg2MudCnstrt, DlegCnstrt,tlegCnstrt,tlegCnstrt2,Dleg2tlegCnstrt,DpCnstrt,tpCnstrt,\
               DbrcCnstrt,tbrcCnstrt,tbrcCnstrt2,Dbrc2tbrcCnstrt,DmudCnstrt,tmudCnstrt,tmudCnstrt2,Dmud2mudCnstrt,\
               dckwidthCnstrt1,dckwidthCnstrt2,ftprintCnstrt,NorsokCnstrt],args=[myjckt,desvarmeans,desvarbds],consargs=None,rhobeg=0.01,maxfun=3000,disp=1)

    print "\n"
    print "Minimum mass Mjacket, MPiles, TPmass = %f %f %f" %(myjckt.Frameouts2.mass[0],myjckt.Mpiles,myjckt.TP.TPouts.mass)
    print "Minimum mass Tower, Jacket(no tower no piles) = %f %f" %(myjckt.Tower.Twrouts.mass,myjckt.Frameouts2.mass[0]-myjckt.Tower.Twrouts.mass)
    print "Minimum found at Dpile=%f, tpile=%f  Lp=%f " % (myjckt.Piles.Pileinputs.Dpile,myjckt.Piles.Pileinputs.tpile,myjckt.Piles.Pileinputs.Lp)
    print "Minimum found at Dbrc=%f, tbrc=%f  " % (myjckt.Xbraces.Xbrcouts.LLURObj.D[0],myjckt.Xbraces.Xbrcouts.LLURObj.t[0])
    print "Minimum found at Dbrcmud=%f, tbrcmud=%f  " % (myjckt.Mudbraces.Mbrcouts.brcObj.D[0],myjckt.Mudbraces.Mbrcouts.brcObj.t[0])
    print "Minimum found at batter=%f, dckwidth=%f, Dleg=%f, tleg=%f,  " % (myjckt.JcktGeoIn.batter,myjckt.dck_width,myjckt.leginputs.Dleg[0],myjckt.leginputs.tleg[0])
    print "Minimum found at Dgir=%f, tgir=%f " % (myjckt.TPinputs.Dgir,myjckt.TPinputs.tgir)
    print "Minimum found at Db=%f DTRb=%f Dt=%f DTRt=%f H2frac=%f " % (myjckt.Tower.Twrins.Db,myjckt.Tower.Twrins.DTRb,myjckt.Tower.Twrins.Dt,myjckt.Tower.Twrins.DTRt,myjckt.Tower.Twrins.Htwr2frac)
    print "Minimum found at Freq %f"  % (myjckt.Frameouts2.Freqs[0])
    print "Minimum found at GLutil=%f EUutil=%f"  % (np.nanmax(myjckt.tower_utilization.GLUtil),np.nanmax(myjckt.tower_utilization.EUshUtil))
    print "Minimum found at Mudline Footprint=%f beta3D=%f"  % (myjckt.wbase,myjckt.PreBuild.beta3D*180./np.pi)
    print "Elapsed time: ", time.time()-tt, "seconds"
    print "Execution count: ", myjckt.exec_count

    #STORE RESULTS
    rescobyla=res1*desvarmeans #dimensional again

    if not(DTRsdiff):#DTRs FIXED TO EACH OTHER, reexpand array
        idx_DTRb=desvars.dvars.keys().index('DTRb')
        rescobyla=np.insert(rescobyla,idx_DTRt,rescobyla[idx_DTRb])

    outdir=ntpath.dirname(casefile)
    SaveOpt1Line(outdir,caseno,casename,desvars,rescobyla,myjckt,xlsfilename,Desprms)

    #Plot
    #PlotJacket(myjckt,util=True,savefileroot=outdir+'\\'+casename)

    return myjckt,casename
#__________________________________________________________#



def JcktWrapper(x,myjckt,desvarmeans,desvarbds):

    """This function builds up the actual model and calculates Jacket stresses and utilization.
    INPUT
        x         -list(N), as in DesVars with that order, but ALL NORMALIZED BY THEIR AVERAGES!!!! \n
        myjckt    -assmebly of jacketSE.\n
        desvarmeans -array(N), average values for each design variable.\n
        """
    global  DTRsdiff
    global  xlast,f1,max_GLUtil,max_EUUtil,max_KjntUtil,max_XjntUtil,max_tutil,max_cbutil,offset,\
            MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat
#   x 0     1     2         3     4    5       6       7        8        9        10    11    12   13      14      15      16
#   batter,Dp,   tp,       Lp,   Dleg,tleg,    Dbrc   tbrc    Dmdbrc   tmdbrc     Dgir,tgir   Db,DTRb,   Dt,DTRt,     H2frac,


    myjckt.JcktGeoIn.batter=x[0]*desvarmeans[0] #remove this after testing


    myjckt.Pileinputs.Dpile=x[1]*desvarmeans[1]
    myjckt.Pileinputs.tpile=x[2]*desvarmeans[2]
    myjckt.Pileinputs.Lp   =x[3]*desvarmeans[3]

    myjckt.leginputs.Dleg  =np.asarray([x[4]*desvarmeans[4]]).repeat(myjckt.JcktGeoIn.nbays+1)
    myjckt.leginputs.tleg  =np.asarray([x[5]*desvarmeans[5]]).repeat(myjckt.JcktGeoIn.nbays+1)

    myjckt.Xbrcinputs.precalc=False
    myjckt.Xbrcinputs.Dbrc    =np.asarray([x[6]*desvarmeans[6]]).repeat(myjckt.JcktGeoIn.nbays)
    myjckt.Xbrcinputs.tbrc    =np.asarray([x[7]*desvarmeans[7]]).repeat(myjckt.JcktGeoIn.nbays)
    myjckt.Mbrcinputs.precalc=False
    myjckt.Mbrcinputs.Dbrc_mud =x[8]*desvarmeans[8]
    myjckt.Mbrcinputs.tbrc_mud =x[9]*desvarmeans[9]

    myjckt.TPinputs.Dgir  =x[10]*desvarmeans[10]
    myjckt.TPinputs.tgir  =x[11]*desvarmeans[11]

    myjckt.Twrinputs.Db    =x[12]*desvarmeans[12]
    myjckt.Twrinputs.DTRb  =x[13]*desvarmeans[13]
    myjckt.Twrinputs.Dt    =x[14]*desvarmeans[14]
    myjckt.Twrinputs.DTRt  =myjckt.Twrinputs.DTRb#x[15]*desvarmeans[15]  #myjckt.Twrinputs.DTRb.copy()   #can use DTRt=DTRb here for simplicity

    myjckt.Twrinputs.Htwr2frac =x[15+int(DTRsdiff)]*desvarmeans[15+int(DTRsdiff)]

    myjckt.JcktGeoIn.dck_width =x[16+int(DTRsdiff)]*desvarmeans[16+int(DTRsdiff)]*myjckt.Twrinputs.Db  #Deck WIdth

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
    print('            dck_width ={:5.2f}, Dbrc={:5.2f}, tbrc={:5.3f}, Dmudbrc={:5.2f}, tmudbrc={:5.3f}').\
            format(myjckt.dck_width, myjckt.Xbrcinputs.Dbrc[0],myjckt.Xbrcinputs.tbrc[0],myjckt.Mbrcinputs.Dbrc_mud,myjckt.Mbrcinputs.tbrc_mud)

    print('from Jwrapper Db={:5.2f}, DTRb={:5.2f}, Dt={:5.2f}, DTRt={:5.2f},H2twrfrac={:5.2f}, Dgir={:5.2f},tgir={:5.3f}, Twrmass={:6.3f}, PilesMass ={:6.3f}, TPmass= {:8.3e}, Frame3DD+Piles Totmass={:10.3f}'.\
            format(myjckt.Twrinputs.Db,myjckt.Twrinputs.DTRb,myjckt.Twrinputs.Dt,\
                   myjckt.Twrinputs.DTRt,myjckt.Twrinputs.Htwr2frac,myjckt.TPinputs.Dgir,myjckt.TPinputs.tgir,myjckt.Tower.Twrouts.mass,myjckt.Mpiles, myjckt.TP.TPouts.mass,mass))
    print(' \n')

    sys.stdout.flush()  #This for peregrine
    return mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat




def f0Cnstrt1(x,myjckt,desvarmeans,desvarbds):  #f1>f0
    global xlast,f1

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=(f1-f0)/f0
    print('f0Cnstrt1=',cnstrt)
    return cnstrt

def f0Cnstrt2(x,myjckt,desvarmeans,desvarbds): #f1<(f0*(1+f0eps))
    global xlast,f1,f0eps

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=(-f1+f0*(1+f0eps))/f0
    print('f0Cnstrt2=',cnstrt)
    return cnstrt

def cbCnstrt(x,myjckt,desvarmeans,desvarbds): #cb<1
    global xlast,f1,max_cbutil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_cbutil
    print('cbCnstrt=',cnstrt)
    return cnstrt

def tCnstrt(x,myjckt,desvarmeans,desvarbds): #tUtil<1
    global xlast,max_tutil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_tutil
    print('tutil constraint=',cnstrt)
    return cnstrt


def KjntCnstrt(x,myjckt,desvarmeans,desvarbds): #KUtil<1
    global xlast,max_KjntUtil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_KjntUtil
    print('KjntCnstrt=',cnstrt)
    return cnstrt

def XjntCnstrt(x,myjckt,desvarmeans,desvarbds): #XUtil<1
    global xlast,max_XjntUtil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_XjntUtil
    print('XjntCnstrt=',cnstrt)
    return cnstrt

def GLCnstrt(x,myjckt,desvarmeans,desvarbds): #GLUtil<1
    global xlast,max_GLUtil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_GLUtil
    print('GL constraint=',cnstrt)
    return cnstrt

def EUCnstrt(x,myjckt,desvarmeans,desvarbds): #EUUtil<1
    global xlast,max_EUUtil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_EUUtil
    print('EU constraint=',cnstrt)
    return cnstrt

#Ensure Dleg>Dmud and Dbrc
def Dleg2BrcCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of Dleg>Dmud and Dbrc
    idx0=4
    idx1=6
    cnstrt=x[idx0]/desvarmeans[idx0]-x[idx1]/desvarmeans[idx1] #this to work when the bounds are the same for the 2 variables
    if (desvarmeans[idx0]-desvarmeans[idx1]) >0.:
        cnstrt=(x[idx0]-x[idx1])/(desvarmeans[idx0]-desvarmeans[idx1])

    print('Dleg2brc constraint=', cnstrt)
    return cnstrt
def Dleg2MudCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of Dleg>Dmud and Dbrc
    idx0=4
    idx1=8
    cnstrt=x[idx0]/desvarmeans[idx0]-x[idx1]/desvarmeans[idx1]  #this to work when the bounds are the same for the 2 variables
    if (desvarmeans[idx0]-desvarmeans[idx1]) >0.:
        cnstrt=(x[idx0]-x[idx1])/(desvarmeans[idx0]-desvarmeans[idx1])

    print('Dleg2mud constraint=', cnstrt)
    return cnstrt

#NOW use constraints for bounds on batter, Dleg, tleg, Dp, Db, DTR
def batCnsrt1(x,myjckt,desvarmeans,desvarbds):  #bound of batter for cobyla batter>7
    idx=0
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('batter constraint1=',cnstrt)
    return cnstrt   #9.5=batvg
def batCnsrt2(x,myjckt,desvarmeans,desvarbds):  #bound of batter for cobyla batter<max7
    idx=0
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('batter constraint2=',cnstrt)
    return cnstrt   #9.5=batvg

#Leg Constraints
def DlegCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of Dleg for cobyla Dleg>1
    idx=4
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('Dleg constraint=', cnstrt)
    return cnstrt
def tlegCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of t for cobyla t>1"
    idx=5
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('tleg constraint=',cnstrt)
    return cnstrt#

def tlegCnstrt2(x,myjckt,desvarmeans,desvarbds):  #bound of t for cobyla t<max
    idx=5
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('tleg constraint2=',cnstrt)
    return cnstrt#

def Dleg2tlegCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of DTRleg>1
    global jcktDTRmin
    idx0=4
    idx1=5
    cnstrt=cnstrt=(x[idx0]*desvarmeans[idx0]/(x[idx1]*desvarmeans[idx1]) - jcktDTRmin)/jcktDTRmin  #this to work when the bounds are the same for the 2 variables
    print('DTRleg constraint=', cnstrt)
    return cnstrt

#Brace Size Constraints
def DbrcCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of Dbrc for cobyla Dleg>1
    idx=6
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('Dbrc constraint=', cnstrt)
    return cnstrt
def tbrcCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of t for cobyla t>1"
    idx=7
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('tbrc constraint=',cnstrt)
    return cnstrt#
def tbrcCnstrt2(x,myjckt,desvarmeans,desvarbds):  #bound of t for cobyla t<max
    idx=7
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('tbrc constraint2=',cnstrt)
    return cnstrt#

def Dbrc2tbrcCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of DTRbrc>1
    global jcktDTRmin
    idx0=6
    idx1=7
    cnstrt=cnstrt=(x[idx0]*desvarmeans[idx0]/(x[idx1]*desvarmeans[idx1]) - jcktDTRmin)/jcktDTRmin  #this to work when the bounds are the same for the 2 variables
    print('DTRbrc constraint=', cnstrt)
    return cnstrt

#Mudbrace Size Constraints
def DmudCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of Dbrc for cobyla Dleg>1
    idx=8
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('Dmud constraint=', cnstrt)
    return cnstrt
def tmudCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of t for cobyla t>1"
    idx=9
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('tmud constraint=',cnstrt)
    return cnstrt#
def tmudCnstrt2(x,myjckt,desvarmeans,desvarbds):  #bound of t for cobyla t<max
    idx=9
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('tmud constraint2=',cnstrt)
    return cnstrt#

def Dmud2mudCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of DTRmud>1
    global jcktDTRmin
    idx0=8
    idx1=9
    cnstrt=cnstrt=(x[idx0]*desvarmeans[idx0]/(x[idx1]*desvarmeans[idx1]) - jcktDTRmin)/jcktDTRmin  #this to work when the bounds are the same for the 2 variables
    print('DTRmud constraint=', cnstrt)
    return cnstrt

#Pile Constraints

def DpCnstrt(x,myjckt,desvarmeans,desvarbds):  #bound of Dp>1
    idx=1
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('Dp constraint=', cnstrt)
    return cnstrt

def tpCnstrt(x,myjckt,desvarmeans,desvarbds):  #tp >0.0254
    idx=2
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('tp constraint=', cnstrt)
    return cnstrt


#TP constraints
def girDCnsrt1(x,myjckt,desvarmeans,desvarbds):  #bound of girder D>Dmin
    idx=10
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('gir D constraint1=',cnstrt)
    return cnstrt
def girDCnsrt2(x,myjckt,desvarmeans,desvarbds):  #bound of girder D<Dmax
    idx=10
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('gir D constraint2=',cnstrt)
    return cnstrt
def girtCnsrt1(x,myjckt,desvarmeans,desvarbds):  #bound of girder t>tmin
    idx=11
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('gir t constraint1=',cnstrt)
    return cnstrt
def girtCnsrt2(x,myjckt,desvarmeans,desvarbds):  #bound of girder t<tmax
    idx=11
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('gir t constraint2=',cnstrt)
    return cnstrt


 #Xbrc constraints
def XbCrit01(x,myjckt,desvarmeans,desvarbds): #XbrcCrit01
    global xlast,XBrcCrit01
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Xbrc constraint01=', XBrcCrit01)
    return XBrcCrit01

def XbCrit02(x,myjckt,desvarmeans,desvarbds): #XbrcCrit02
    global xlast,XBrcCrit02
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Xbrc constraint02=', XBrcCrit02)
    return XBrcCrit02

def XbCrit03(x,myjckt,desvarmeans,desvarbds): #XbrcCrit03
    global xlast,XBrcCrit03
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Xbrc constraint03=', XBrcCrit03)
    return XBrcCrit03

def XbCrit04(x,myjckt,desvarmeans,desvarbds): #XbrcCrit03
    global xlast,XBrcCrit04
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Xbrc constraint04=', XBrcCrit04)
    return XBrcCrit04

def XbCrit05(x,myjckt,desvarmeans,desvarbds): #XbrcCrit05
    global xlast,XBrcCrit05
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Xbrc constraint05=', XBrcCrit05)
    return XBrcCrit05

def NorsokCnstrt(x,myjckt,desvarmeans,desvarbds): # -Norsok angle leg-brace >30 deg
    NorsokMin=30.*np.pi/180.
    cnstrt=(myjckt.PreBuild.beta3D- NorsokMin)/NorsokMin
    print('Xbrc Norsok constraint=', cnstrt)
    return cnstrt


 #Mudbrc constraints
def MbCrit01(x,myjckt,desvarmeans,desvarbds): #MbrcCrit01
    global xlast,MudCrit01
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Mbrc constraint01=', MudCrit01)
    return MudCrit01

def MbCrit02(x,myjckt,desvarmeans,desvarbds): #XbrcCrit02
    global xlast,MudCrit02
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Mbrc constraint02=', MudCrit02)
    return MudCrit02

def MbCrit03(x,myjckt,desvarmeans,desvarbds): #XbrcCrit03
    global xlast,MudCrit03
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Mbrc constraint03=', MudCrit03)
    return MudCrit03

def MbCrit04(x,myjckt,desvarmeans,desvarbds): #XbrcCrit03
    global xlast,MudCrit04
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Mbrc constraint04=', MudCrit04)
    return MudCrit04

def MbCrit05(x,myjckt,desvarmeans,desvarbds): #XbrcCrit05
    global xlast,MudCrit05
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Mbrc constraint05=', MudCrit05)
    return MudCrit05

 #Tower constraints
def TwrCnstrt01(x,myjckt,desvarmeans,desvarbds): #Maximum tower Diameter
    idx=12
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('Db TwrCnstrt01=',cnstrt )
    return cnstrt
def TwrCnstrt02(x,myjckt,desvarmeans,desvarbds): #Maximum tower Diameter at the top
    idx=14
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('Dt TwrCnstrt02=', cnstrt)
    return cnstrt
def TwrCnstrt03(x,myjckt,desvarmeans,desvarbds): #Minimum tower Diameter at the top
    idx=14
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('Dt TwrCnstrt03=', cnstrt)
    return cnstrt
def TwrCnstrt04(x,myjckt,desvarmeans,desvarbds):  #Minimum DTRb>120
    idx=13
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    #print('DTR min TwrCnstrt04=', cnstrt)
    return cnstrt
def TwrCnstrt05(x,myjckt,desvarmeans,desvarbds):  #Max DTRb<200
    idx=13
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    #print('DTR max TwrCnstrt05=',cnstrt)
    return cnstrt
def TwrCnstrt06(x,myjckt,desvarmeans,desvarbds):  #Minimum DTRt>120
    idx=15
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    #print('DTR min TwrCnstrt06=', cnstrt)
    return cnstrt
def TwrCnstrt07(x,myjckt,desvarmeans,desvarbds):  #Max DTRt<200
    idx=15
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    #print('DTR max TwrCnstrt07='cnstrt)
    return cnstrt
def TwrCnstrt08(x,myjckt,desvarmeans,desvarbds):  #Maximum Htwr2 < Htwr/4
    global    DTRsdiff
    idx=15+int(DTRsdiff)
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('Htwr2 max TwrCnstrt08=',cnstrt)
    return cnstrt
def TwrCnstrt09(x,myjckt,desvarmeans,desvarbds):  #Minimum Htwr2 >0.005
    global    DTRsdiff
    idx=15+int(DTRsdiff)
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('Htwr2 min TwrCnstrt09=',cnstrt)
    return cnstrt
#deck width factor constraints

def dckwidthCnstrt1(x,myjckt,desvarmeans,desvarbds): #Minimum deck width=2*Db
    global    DTRsdiff
    idx=16+int(DTRsdiff)   #deck width factor
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('dck_width cnstrt1=',cnstrt )
    return cnstrt
def dckwidthCnstrt2(x,myjckt,desvarmeans,desvarbds): #Max deck width=3*Db
    global    DTRsdiff
    idx=16+int(DTRsdiff)   #deck width factor
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('dck_width cnstrt2=',cnstrt )
    return cnstrt

def ftprintCnstrt(x,myjckt,desvarmeans,desvarbds): #Max footprint
    global    mxftprint
    cnstrt=(mxftprint-myjckt.wbase)/mxftprint
    print('footprint cnstrt=',cnstrt )
    return cnstrt

##def dckwidthCnstrt2(x,myjckt,desvarmeans,desvarbds): #Max deck width<=factmax*Db
##    global    DTRsdiff
##    idx0=16+int(DTRsdiff)   #deck width factor
##    idx1=12                 #Db
##    Dbcurrent=x[idx1]*desvarmeans[idx1]
##    dckcurrent=x[idx0]*desvarmeans[idx0]*Dbcurrent
##    mx_dckwidth=desvarbds[idx0,1]*Dbcurrent
##    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
##    print('dck_width cnstrt2=',cnstrt )
##    return cnstrt
##def dckwidthCnstrt1(x,myjckt,desvarmeans,desvarbds): #Minimum deck width=2*Db
##    global    DTRsdiff
##    idx0=16+int(DTRsdiff)   #deck width
##    idx1=12                 #Db
##    mn_dckwidth=2.*x[idx1]*desvarmeans[idx1]
##    cnstrt=(x[idx0]*desvarmeans[idx0]- mn_dckwidth)/mn_dckwidth
##    print('dck_width cnstrt1=',cnstrt )
##    return cnstrt
##def dckwidthCnstrt2(x,myjckt,desvarmeans,desvarbds): #Max deck width=2.x*Db
##    global    DTRsdiff,dck_maxfact
##    idx0=16+int(DTRsdiff)   #deck width
##    idx1=12                 #Db
##    mx_dckwidth=dck_maxfact*x[idx1]*desvarmeans[idx1]
##    cnstrt=(-x[idx0]*desvarmeans[idx0]+ mx_dckwidth)/mx_dckwidth
##    print('dck_width cnstrt2=',cnstrt )
##    return cnstrt


#Embedment length constraint
def LpCnstrt(x,myjckt,desvarmeans,desvarbds):  #Maximum Htwr2 < Htwr/4
    global xlast,Lp0rat
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_cbutil,max_KjntUtil,max_XjntUtil,max_GLUtil,max_EUUtil,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat=JcktWrapper(x,myjckt,desvarmeans,desvarbds)
        xlast=x.copy()
    print('Lp0rat constraint=',Lp0rat)
    return Lp0rat

#Embedment length constraint2
def LpCnstrt2(x,myjckt,desvarmeans,desvarbds):  #Maximum Lp < Lpmax
    idx=3
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('Lp constraint2=',cnstrt)
    return cnstrt

    #Function to minimize
def mass(x,myjckt,desvarmeans,desvarbds):
    """This function assembles output data from the assembly in terms of mass: \n
        it is mostly frame3dd mass+piles'' mass."""
    return  JcktWrapper(x,myjckt,desvarmeans,desvarbds)[0]/1.5e6

def maxcnstrt(x,idx,desvarmeans,desvarbds):
    return -x[idx]+desvarbds[idx,1]/desvarmeans[idx]
def mincnstrt(x,idx,desvarmeans,desvarbds):
    return  x[idx]-desvarbds[idx,0]/desvarmeans[idx]


#______________________________________________________________________________#

def main():
        #SET THIS GLOBAL VARIABLE HERE
    global    DTRsdiff ,jcktDTRmin#,maxftprint #,dck_maxfact
    DTRsdiff=False    ##SET THIS TO TRUE IF YOU WANT DTRs to be different between base and top
    jcktDTRmin=22.  ##Set this one to the minimum DTR allowed, mostly for 60+waterdepths
    #dck_maxfact=2.3  ## Set this one to either 2.1 or 2.3 if 6 MW >=40m
    #maxftprint=35.
    guessfromoutfile=[] #initialize

    if len(sys.argv)>1:
        casefile=sys.argv[1]
        casenostart=int(sys.argv[2])
        casenoends=int(sys.argv[3])
        xlsfilename=sys.argv[4]
        if len(sys.argv)>5:
            guessfromoutfile=sys.argv[5]

    else:
        casefile=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\testmatrixspec.dat'
        casefile=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\testmatrixspec.dat'
        casenostart=82  #UHREACT
        casenoends=82

        xlsfilename=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\outputs.xls'
        xlsfilename=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\outputs.xls'


        #guessfromoutfile=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\COBYLA\10MW_55_\58_R2D0H1M0T0S0_out.dat'
        #guessfromoutfile=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\COBYLA\76_R2D2H1M0T0S0_specialin.dat'
        #guessfromoutfile=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\COBYLA\16_R0D1H2M0T0S0_out.dat'

    f0epsilon=0.05 #upper f0 allowed

    for caseno in range(casenostart,casenoends+1):
        if caseno>18:  #For cases at 6 and 10 MW and for 3MW with deep water relax the f0 constraint
            f0epsilon=0.35
        print('#____________________________________________#\n')
        print(('# JacketOpt_Peregrine NOW PROCESSING CASE No. {:d} #\n').format(caseno) )
        print('#____________________________________________#\n')

        myjckt,casename=JcktOpt1line(casefile,caseno,xlsfilename,f0epsilon=f0epsilon,guessfromoutfile=guessfromoutfile)
        guessfromoutfile=os.path.join(os.path.dirname(casefile),casename + '_out.dat')
        sys.stdout.flush()  #This for peregrine
    ##SAMPLE CALL FROM OUTSIDE ENVIRONMENT
    ##python JacketOpt_Peregrine.py D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\testmatrix.dat 58 58 D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\output.xls D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\COBYLA\10MW_55_\58_R2D0H1M0T0S0_out.dat





if __name__ == '__main__':
    #This is how you call this function
    main()
