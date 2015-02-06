#-------------------------------------------------------------------------------
# Name:        JacketOpt_Py_MDAOopt.py
# Purpose:     It Expects Input DATA to be edited directly from 1 input file and then calls the optimizer based on that.
#              !!!!!!!!!!!In particular it used to use a table of cases for the LCOE analysis project.!!!!!!!!!!!!
#              This can also be a non-multicase variant that started from JacketOpt_PyOPT_multicaseTable.py
#              See MyJacketInputs.py for an example of individual input file.
#              This Program will use either OpenMDAO internal pyOPT, or OpenMDAO internal CobylaDriver, or external pyOPT. No external, non-pyOPT Cobyla.
#              NOTE: HERE DESVARS ARE NOT MADE DIMENSIOLESS!!! JUST objfunc and cosntraints! vs JacketOPt_Peregrine
# Author:      rdamiani
#
# Created:     10/2014 based off of JacketOpt_PyOPT_multicaseTable.py
#              11/2014: merged with JacketOpt_Peregrine.py to be able to run external pyOPT optimizations for 1 or multiple cases
#              11/2014: extended to be able to also use internal MDAO optimization drivers (internal pyOPTdriver(SNOPT and cobyla) or non-pyopt cobyladriver)
#
# Copyright:   (c) rdamiani 2014
# License:     Apache (2014)
#-------------------------------------------------------------------------------
import numpy as np
import scipy.optimize

import pyOpt

MPIFlag=True  #INitialize
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
except:
    MPIFlag=False
    #raise ImportError('mpi4py is required for parallelization')

import time
import sys
import os
import ntpath
import copy
import imp

from openmdao.lib.casehandlers.api import DumpCaseRecorder,ListCaseRecorder,CSVCaseRecorder, DBCaseRecorder

from PlotJacket import main as PlotJacket  #COMMENT THIS ONE OUT FOR PEREGRINE"S SAKE
from printJacketres import main as printJacket

from ReadSaveOpt1LineSysEng import ReadTab1Line,SaveOpt1Line,ReadOptFile,DesVars


#__________________________________________________________#


              #  OPTIMIZATION FUNCTION IS HERE #

#__________________________________________________________#

def JcktOpt(prmsfile, SNOPTflag=False, MDAOswitch=[], tablefile=[], caseno=[],xlsfilename=[],guessfromoutfile=[],f0epsset=[],jcktDTRminset=[],towerfixed=False):

    """Optimizer which either reads one line caseno out of table file tablefile,\n
       or an individual input file prmsfile. Note you also need to provide a prmfile, in the multicase table case. \n

    INPUT \n
        prmsfile  -string, complete path+filename to main parameter setting file (e.g.: MyJacketInputFile.py for an individual case, or SetJacketInputsPeregrine.py for multicase table case). \n

    OPTIONALS \n
        SNOPTflag   -boolean, if True it will use SNOPT else COBYLA both from pyOPT. \n
        MDAOswitch  -string, options for optimizer: \n
			 'extCobyla'   : This will use the external (to OpenMDAO) pyOPT with Cobyla. Only in this case variables are made dimensionless.   \n
			 'pyCobyla'    : This will use the external (to OpenMDAO) pyOPT with Cobyla.   \n
             'pySNOPT'     : This will use the external (to OpenMDAO) pyOPT with SNOPT.    \n
             'md_pyCobyla' : This will use the internal (to OpenMDAO) pyOPTdriver with Cobyla.  \n
             'md_pySNOPT'  : This will use the internal (to OpenMDAO) pyOPTdriver with SNOPT.  \n
             'md_Cobyla'   : This will use the internal (to OpenMDAO) Cobyla driver (non-pyOPT).  \n
             \n
             \n
        towerfixed      -Boolean, Set it to True if tower does not need be optimized. \n
               The following are all needed in case of mullticase table input file.\n
        tablefile   -string, complete path+filename to table file of cases, or \n
        caseno      -int, case number (also row out of the table to be read). Only for multicase table case.\n
        xlsfilename -string, complete path+filename to excel file, whose sheet for the current caseno case will be added/updated.Only for multicase table case. \n
        guessfromoutfile -string, complete path+name to an optimization output file where to read guesses from. Only for multicase table case. \n
        f0epsset       -Float,       Set this to a value such that the f0*(1+f0epsilon) will not be exceeded: Only for multicase table case.\n
        jcktDTRminset  -Float,       Set this one to the minimum jacket member DTR allowed, mostly for 60+waterdepths: Only for multicase table case.\n

    OUTPUT \n
        It will create a final configuration optimized for the design parameters in the prmsfile file. \n
        Also, it will save a dump-recorder file, and summary for report, respectively. \n
        Figures will also be made of configuration and tower utilization.  \n
        casename -string, case name for hte current case usable by other programs. \n
        myjckt   -OPenMdao assembly of JacketSE, final configuration after optimization.
        \n
        """
    global DTRsdiff,f0,f0eps,jcktDTRmin,mxftprint,MDAOswitch2,towerfix  #global vars to pass a few parameters around through the constraints

    towerfix=towerfixed

    desvars=DesVars() #instance of design variables


    directory, module_name = os.path.split(prmsfile)
    module_name = os.path.splitext(module_name)[0]

    MyJacketInputs = imp.load_source(module_name,prmsfile)  #This can be either MyJacketInputs or SetJacketInputsPeregrine
    #imp module has been used rather than the __import__ which needed path to be adjusted and that was creating issues as it would always load the module from the current dir first
    #path = list(sys.path)
    #sys.path.insert(0, directory)
    #try:
        #MyJacketInputs =        __import__(module_name)  #This can be either MyJacketInputs or SetJacketInputsPeregrine
    #finally:
    #    sys.path[:] = path # restore

    casename='JacketOpt' #initialize
    outdir=ntpath.dirname(prmsfile)

    if caseno: #This is the case of a table of multiple cases file

        #Set global variables that are passed through keywords
        f0epsilon=f0epsset
        jcktDTRmin=jcktDTRminset

            #First read design parameters and des var bounds
        Desprms,_,desvarbds,desvarmeans,guesses,casename=ReadTab1Line(tablefile,caseno,desvars.dvars.keys())
        #USING bds as desvarbds

        #Assign mxftprint
        mxftprint=Desprms.mxftprint

        #Target frequency
        f0=Desprms.f0
        #Then set up an initial guess at the Design variables and use it also to instantiate the assembly

        for ii,key in enumerate(desvars.dvars):
            #print key  #debug
            #setattr(desvars,key,desvarmeans[ii])  #from mean
            setattr(desvars,key,guesses[ii])       #from user's input NOTE IT IS IMPORTANT TO CHECK ORDER OF variables

        #Then build the initial assembly
        myjckt,f0epsilon,jcktDTRmin=MyJacketInputs.main(Desprms,desvars)

    else:
        #Read Input & Build the initial assembly and get guesses, for the x design variable array check out objfunc
        myjckt,f0,f0epsilon,jcktDTRmin,mxftprint,guesses,desvarbds=MyJacketInputs.main()
        desvarmeans=np.mean(desvarbds,1)

    # Set other global var here
    DTRsdiff=myjckt.Twrinputs.DTRsdiff

    f0eps=f0epsilon

    guess=guesses

    if guessfromoutfile:  #This only used when multiple case file is used
        desvarguess=ReadOptFile(guessfromoutfile)
        for ii,key in enumerate(desvarguess.dvars):
            guess[ii]=getattr(desvarguess,key)#/desvarmeans[ii]

    varlist=desvars.dvars.keys()
    idx_DTRt=varlist.index('DTRt')
    if not(DTRsdiff):#DTRs FIXED TO EACH OTHER
        varlist.pop(idx_DTRt)
        guess=np.delete(guess,idx_DTRt)
        desvarmeans=np.delete(desvarmeans,idx_DTRt)
        desvarbds=np.delete(desvarbds,idx_DTRt,0)

    #SET UP THE OPTIMIZATION PROBLEM
    MDAOswitch=MDAOswitch.lower()

    MDAOswitch2=False #Initialize
    if MDAOswitch=='extcobyla':
    #_____________________________________#
           #  FIRST EXT OPT CASE   #
    #_____________________________________#
        MDAOswitch2=True
        guess=guess/desvarmeans
        if towerfix:
            junk=int(not DTRsdiff)
            guess=np.hstack((guess[0:-6+junk],guess[-1]))
            desvarmeans=np.hstack((desvarmeans[0:-6+junk],desvarmeans[-1]))
            desvarbds=np.vstack((desvarbds[0:-6+junk,:],desvarbds[-1,:]))

        constrfuncs=[f0Cnstrt1,f0Cnstrt2,batCnsrt1, batCnsrt2,cbCnstrt,tCnstrt,KjntCnstrt,XjntCnstrt,\
                     girDCnsrt1,girDCnsrt2,girtCnsrt1,girtCnsrt2,\
                     XbCrit01,XbCrit02,XbCrit03,XbCrit04,XbCrit05,MbCrit01,MbCrit02,MbCrit03,MbCrit04,MbCrit05,\
                     Dleg2BrcCnstrt,Dleg2MudCnstrt, DlegCnstrt,tlegCnstrt,tlegCnstrt2,Dleg2tlegCnstrt,DpCnstrt,tpCnstrt,\
                     DbrcCnstrt,tbrcCnstrt,tbrcCnstrt2,Dbrc2tbrcCnstrt,DmudCnstrt,tmudCnstrt,tmudCnstrt2,Dmud2mudCnstrt,\
                     dckwidthCnstrt1,dckwidthCnstrt2,ftprintCnstrt,NorsokCnstrt,LpCnstrt,LpCnstrt3]  #Initialize
        if not(towerfix):
            constrfuncs.extend([GLCnstrt,EUCnstrt,TwrCnstrt01,TwrCnstrt02,TwrCnstrt03,TwrCnstrt04,TwrCnstrt05,TwrCnstrt08,TwrCnstrt09])
            if DTRsdiff:
                constrfuncs.extend([TwrCnstrt06,TwrCnstrt07])
        if myjckt.twodlcs:
            constrfuncs.extend([cbCnstrt2,tCnstrt2,KjntCnstrt2,XjntCnstrt2,GLCnstrt2,EUCnstrt2,LpCnstrt2])

        #Cobyla : If it returns 3 elements , it means it did not converge;  constraint2,constraint3
        tt = time.time()
        res1=scipy.optimize.fmin_cobyla(mass,guess,constrfuncs,args=[myjckt,desvarmeans,desvarbds],consargs=None,rhobeg=0.01,maxfun=3000,disp=1)

        #Store results for later
        xstr=res1*desvarmeans #dimensional again

    elif (MDAOswitch=='pycobyla') or (MDAOswitch=='pysnopt'):
    #_____________________________________#
           #  THEN PYOPT CASES   #
    #_____________________________________#
        if towerfix:
            junk=int(not DTRsdiff)
            guess=np.hstack((guess[0:-6+junk],guess[-1]))
            desvarmeans=np.hstack((desvarmeans[0:-6+junk],desvarmeans[-1]))
            desvarbds=np.vstack((desvarbds[0:-6+junk,:],desvarbds[-1,:]))
            varlist2=(varlist[0:-6+junk])
            varlist2.append(varlist[-1])
            varlist=varlist2

        opt_prob=pyOpt.Optimization('Jacket Optimization via External PyOPT with {!r:^} '.format(MDAOswitch), objfunc)
        if caseno:
            opt_prob=pyOpt.Optimization('LCOE 1 Line Case no. {:d} Optimization'.format(caseno), objfunc)

        opt_prob.addObj('mass')
        for ii,key in enumerate(varlist):
            opt_prob.addVar(key,'c',lower=desvarbds[ii,0],upper=desvarbds[ii,1],value=guess[ii])

        cnt= 6+ 11+7 #counter for number of constraints- initial value for fixed tower and 1 dlc
        if not(towerfix):
            cnt += 2
        if myjckt.twodlcs:
            cnt += 4 +2+1

        opt_prob.addConGroup('cnstrts',cnt,type='i')
        print opt_prob

        #Finally call the optimizer
        args=(myjckt,desvarmeans,desvarbds)

        strname2='' #initialize
        if MDAOswitch=='pysnopt':
            strname='pyopt_snopt_'+'.hst'

            if caseno:
                strname2=str(caseno).zfill(2)
                strname='pyopt_snopt_'+strname2+'.hst'

            opt_prob.write2file(outfile=os.path.join(os.path.dirname(prmsfile),strname), disp_sols=False, solutions=[])

            #set some strings for MPI incase
            mpistr=''
            printfstr='pyopt_snopt_print_'
            summfstr='pyopt_snopt_summary_'
            if MPIFlag:
                mpistr='pgc'
                printfstr='pyopt_mpisnopt_print_'
                summfstr='pyopt_mpisnopt_summary_'

            opt =pyOpt.SNOPT()  #non-MPI here always
            opt.setOption('Minor print level',1)
            opt.setOption('Major feasibility tolerance',1.e-3)
            opt.setOption('Major optimality tolerance',1.e-3)
            opt.setOption('Minor feasibility tolerance',1.e-3)
            opt.setOption('Print file',os.path.join(os.path.dirname(prmsfile),printfstr+strname2+'.out'))
            opt.setOption('Summary file',os.path.join(os.path.dirname(prmsfile),summfstr+strname2+'.out'))

            opt.setOption('Solution','Yes')
            #Solve
            tt = time.time()

            [fstr, xstr, inform]=opt(opt_prob,'FD',True,True,True,False,mpistr,{},*args)  #for parallel gradient calculations
            #
            #opt_problem={}, sens_type='FD', store_sol=True, disp_opts=False, store_hst=False, hot_start=False, sens_mode='', sens_step={}, *args, **kwargs)

        #COBYLA
        else:
            mpistr=None
            ifilestr='pyopt_cobyla_'
            if MPIFlag:
                mpistr='POA'
                ifilestr='pyopt_mpicobyla_'

            opt =pyOpt.COBYLA(pll_type=mpistr)

            opt.setOption('RHOBEG',0.001)
            opt.setOption('RHOEND',1.e-3)
            opt.setOption('MAXFUN',2000)
            opt.setOption('IPRINT',1)
            opt.setOption('IFILE',os.path.join(os.path.dirname(prmsfile),ifilestr+strname2+'.hst') ) #store cobyla output

            #Solve
            tt = time.time()

            [fstr, xstr, inform]=opt(opt_prob,  True,           False,             False,          False, *args)
        ###opt_problem={}, store_sol=True, disp_opts=False, store_hst=False, hot_start=False

    else:

        #_____________________________________#
                     #THEN MDAO CASES
        #_____________________________________#
        # --- optimizer imports ---
        from pyopt_driver.pyopt_driver import pyOptDriver
        from openmdao.lib.casehandlers.api import DumpCaseRecorder
        from openmdao.lib.drivers.api import COBYLAdriver

        # ----------------------

        # --- Setup Optimizer ---
        if MDAOswitch == 'md_cobyla':
            myjckt.replace('driver', COBYLAdriver())

            myjckt.driver.rhobeg=0.01
            myjckt.driver.rhoend=1.e-3
            myjckt.driver.maxfun=2000
            myjckt.driver.iprint=1
        else:
            myjckt.replace('driver', pyOptDriver())
            myjckt.driver.pyopt_diff=True #This makes pyopt calculate finite differences

            if  MDAOswitch== 'md_pysnopt':
                myjckt.driver.optimizer = 'SNOPT'
                myjckt.driver.options = {'Major feasibility tolerance': 1e-3,\
                                     'Minor feasibility tolerance': 1e-3,\
                                     'Major optimality tolerance': 1e-3,\
                                     'Function precision': 1e-3}

            elif MDAOswitch== 'md_pycobyla':
                myjckt.driver.optimizer = 'COBYLA'
                myjckt.driver.options = {'RHOEND':1.e-2,'RHOEND':1.e-3,'MAXFUN':2000,'IPRINT':1}

            else:
                sys.exit('Error: MDAOswitch must be set to ''pyCobyla'' or ''pySNOPT'' or ''md_Cobyla'' or ''md_pySNOPT'' or ''md_pyCobyla'' or ''extCobyla'' !!!')
        # ----------------------

        # --- Objective ---
        myjckt.driver.add_objective('(LoadFrameOuts.Frameouts.mass[0]+LoadFrameOuts.Mpiles)/1.e6')
        # ----------------------

        # --- Design Variables ---
        #   x 0     1     2         3     4    5       6       7        8        9        10    11    12   13   14   15        16          17
        #   batter,Dp,   tp,       Lp,   Dleg,tleg,    Dbrc   tbrc    Dmdbrc   tmdbrc     Dgir,tgir   Db,DTRb,   Dt,DTRt,     H2frac, dck_withfrac

        myjckt.driver.add_parameter('JcktGeoIn.batter',      low=desvarbds[0,0],  high=desvarbds[0,1])
        myjckt.driver.add_parameter('Pileinputs.Dpile',      low=desvarbds[1,0],  high=desvarbds[1,1])
        myjckt.driver.add_parameter('Pileinputs.tpile',      low=desvarbds[2,0],  high=desvarbds[2,1])
        myjckt.driver.add_parameter('Pileinputs.Lp',         low=desvarbds[3,0],  high=desvarbds[3,1])
        myjckt.driver.add_parameter('leginputs.Dleg0',       low=desvarbds[4,0],  high=desvarbds[4,1])
        myjckt.driver.add_parameter('leginputs.tleg0',       low=desvarbds[5,0],  high=desvarbds[5,1])
        myjckt.driver.add_parameter('Xbrcinputs.Dbrc0',      low=desvarbds[6,0],  high=desvarbds[6,1])
        myjckt.driver.add_parameter('Xbrcinputs.tbrc0',      low=desvarbds[7,0],  high=desvarbds[7,1])
        myjckt.driver.add_parameter('Mbrcinputs.Dbrc_mud',   low=desvarbds[8,0],  high=desvarbds[8,1])
        myjckt.driver.add_parameter('Mbrcinputs.tbrc_mud',   low=desvarbds[9,0],  high=desvarbds[9,1])
        myjckt.driver.add_parameter('TPinputs.Dgir',         low=desvarbds[10,0], high=desvarbds[10,1])
        myjckt.driver.add_parameter('TPinputs.tgir',         low=desvarbds[11,0], high=desvarbds[11,1])
        varcnt=12 #To take care of indices
        #tower  fixed case
        if not(towerfix):
            myjckt.driver.add_parameter('Twrinputs.Db',          low=desvarbds[12,0], high=desvarbds[12,1])
            myjckt.driver.add_parameter('Twrinputs.DTRb',        low=desvarbds[13,0], high=desvarbds[13,1])
            myjckt.driver.add_parameter('Twrinputs.Dt',          low=desvarbds[14,0], high=desvarbds[14,1])

            if DTRsdiff:
                myjckt.driver.add_parameter('Twrinputs.DTRt',        low=desvarbds[15,0], high=desvarbds[15,1])
            ##else:  The else case is already set from input file in case
            ##    myjckt.Twrinputs.DTRsdiff  = False       #Set and Link DTRt=DTRb if needed, within the assembly

            myjckt.driver.add_parameter('Twrinputs.Htwr2frac',   low=desvarbds[15+int(DTRsdiff),0], high=desvarbds[15+int(DTRsdiff),1])
            varcnt=16
       #Note: because we want to use dck_withfrac as a parameter, we need to 0-out the dck_width, else it would trump the dck_withfrac
        myjckt.JcktGeoIn.dck_width=0.
        myjckt.driver.add_parameter('JcktGeoIn.dck_widthfrac',low=desvarbds[varcnt+int(DTRsdiff),0], high=desvarbds[varcnt+int(DTRsdiff),1])

        #--- Constraints ---#
        myjckt.driver.add_constraint('LoadFrameOuts.Frameouts.Freqs[0] >= {:f}'.format(f0))
        myjckt.driver.add_constraint('LoadFrameOuts.Frameouts.Freqs[0] <= {:f}'.format(f0*(1+f0eps)))

        if not(towerfix):
            myjckt.driver.add_constraint('max(LoadFrameOuts.tower_utilization.GLUtil) <=1.0')
            myjckt.driver.add_constraint('max(LoadFrameOuts.tower_utilization.EUshUtil) <=1.0')
            if myjckt.twodlcs:
                myjckt.driver.add_constraint('max(LoadFrameOuts2.tower_utilization.GLUtil) <=1.0')
                myjckt.driver.add_constraint('max(LoadFrameOuts2.tower_utilization.EUshUtil) <=1.0')

        myjckt.driver.add_constraint('numpy.nanmax(LoadFrameOuts.jacket_utilization.t_util) <=1.0')
        myjckt.driver.add_constraint('numpy.nanmax(LoadFrameOuts.jacket_utilization.cb_util) <=1.0')
        myjckt.driver.add_constraint('numpy.nanmax(LoadFrameOuts.jacket_utilization.KjntUtil) <= 1.0')
        myjckt.driver.add_constraint('numpy.nanmax(LoadFrameOuts.jacket_utilization.XjntUtil) <= 1.0')
        if myjckt.twodlcs:
            myjckt.driver.add_constraint('numpy.nanmax(LoadFrameOuts2.jacket_utilization.t_util) <=1.0')
            myjckt.driver.add_constraint('numpy.nanmax(LoadFrameOuts2.jacket_utilization.cb_util) <=1.0')
            myjckt.driver.add_constraint('numpy.nanmax(LoadFrameOuts2.jacket_utilization.KjntUtil) <= 1.0')
            myjckt.driver.add_constraint('numpy.nanmax(LoadFrameOuts2.jacket_utilization.XjntUtil) <= 1.0')

        myjckt.driver.add_constraint('PreBuild.wbase <= {:f}'.format(mxftprint))

        myjckt.driver.add_constraint('leginputs.Dleg0 >= Mbrcinputs.Dbrc_mud')
        myjckt.driver.add_constraint('leginputs.Dleg0 >= Xbrcinputs.Dbrc0')

        myjckt.driver.add_constraint('leginputs.Dleg0/leginputs.tleg0     >= {:f}'.format(jcktDTRmin))
        myjckt.driver.add_constraint('Xbrcinputs.Dbrc0/Xbrcinputs.tbrc0         >= {:f}'.format(jcktDTRmin))
        myjckt.driver.add_constraint('Mbrcinputs.Dbrc_mud/Mbrcinputs.tbrc_mud >= {:f}'.format(jcktDTRmin))

        myjckt.driver.add_constraint('BrcCriteria.XBrcCriteria.brc_crit01 >=0.')
        myjckt.driver.add_constraint('BrcCriteria.XBrcCriteria.brc_crit02 >=0.')
        myjckt.driver.add_constraint('BrcCriteria.XBrcCriteria.brc_crit03 >=0.')
        myjckt.driver.add_constraint('BrcCriteria.XBrcCriteria.brc_crit04 >=0.')
        myjckt.driver.add_constraint('BrcCriteria.XBrcCriteria.brc_crit05 >=0.')
        myjckt.driver.add_constraint('BrcCriteria.MudBrcCriteria.brc_crit01 >=0.')
        myjckt.driver.add_constraint('BrcCriteria.MudBrcCriteria.brc_crit02 >=0.')
        myjckt.driver.add_constraint('BrcCriteria.MudBrcCriteria.brc_crit03 >=0.')
        myjckt.driver.add_constraint('BrcCriteria.MudBrcCriteria.brc_crit04 >=0.')
        myjckt.driver.add_constraint('BrcCriteria.MudBrcCriteria.brc_crit05 >=0.')

        NorsokMin=30.*np.pi/180
        myjckt.driver.add_constraint('PreBuild.beta3D >= {:f}'.format(NorsokMin))

        myjckt.driver.add_constraint('LoadFrameOuts.Lp0rat >= 0.')
        if myjckt.twodlcs:
            myjckt.driver.add_constraint('LoadFrameOuts2.Lp0rat >= 0.')
        # ----------------------

        # --- recorder ---
        myjckt.recorders = [DumpCaseRecorder(ntpath.join(outdir,casename+'_'+MDAOswitch+'.rcd'))]
        # ----------------------

        #RUN
        tt = time.time()
        myjckt.run()
        #Set results into xstr to be written later
        xstr=myjckt.driver.eval_parameters()

    if (MDAOswitch=='pycobyla') or (MDAOswitch=='pysnopt'):
        print opt_prob.solution(0)

    # Now show results
    print "\n\n"
    print "Minimum found with:\n"
    print "Elapsed time: ", time.time()-tt, "seconds"
    print "Execution count: ", myjckt.exec_count
    print "\n"

    #print geometry
    printJacket(myjckt)

    #STORE RESULTS

    if caseno:
        if towerfix: #reexpand array
            xstr=np.hstack((xstr[0:-1],myjckt.Twrinputs.Db,myjckt.Twrinputs.DTRb,myjckt.Twrinputs.Dt,myjckt.Twrinputs.DTRt,myjckt.Twrinputs.Htwr2frac,xstr[-1]))
        if not(DTRsdiff):#DTRs FIXED TO EACH OTHER, reexpand array
            idx_DTRb=desvars.dvars.keys().index('DTRb')
            xstr=np.insert(xstr,idx_DTRt,xstr[idx_DTRb])

        SaveOpt1Line(outdir,caseno,casename,desvars,xstr,myjckt,xlsfilename,Desprms)
    else: #This does not work
        casename=ntpath.splitext(prmsfile)[0]
        caseno=1
        xlsfilename=ntpath.join(outdir,'singleoutput.xlsx')
        filename=casename+'.optidat'

        #SaveOpt1Line(outdir,caseno,casename,desvars,xstr,myjckt,xlsfilename,Desprms)
        #fp=open(filename,'w')

       # myjckt.recorders = [DumpCaseRecorder(filename)]
       # myjckt.run()  # just to write case
        #fp.close()
        import pickle
        filename=filename+'.pik'  #E.g.: "C:\PROJECTS\OFFSHORE_WIND\UH_REACT\PYTHON_OPT\JCKdataNewCd.txt" #filename="D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\UH_REACT\JCKdataNoMass.p"
        with open(filename,'w') as fp:
            pickle.dump([xstr,myjckt.LoadFrameOuts.Frameouts.mass,myjckt.Mpiles,myjckt.Tower.Twrouts.mass,myjckt.Piles.Pileinputs.Dpile,myjckt.Piles.Pileinputs.tpile,myjckt.Piles.Pileinputs.Lp,\
                        myjckt.Xbraces.Xbrcouts.LLURObj.D,myjckt.Xbraces.Xbrcouts.LLURObj.t,myjckt.Mudbraces.Mbrcouts.brcObj.D,myjckt.Mudbraces.Mbrcouts.brcObj,\
                        myjckt.JcktGeoIn.batter,myjckt.leginputs.Dleg,myjckt.leginputs.tleg,myjckt.TPinputs.Dgir,myjckt.TPinputs.tgir,myjckt.TP.TPouts.TPlumpedMass,\
                        myjckt.Tower.Twrins.Db,myjckt.Tower.Twrins.DTRb,myjckt.Tower.Twrins.Dt,myjckt.Tower.Twrins.DTRt,myjckt.Tower.Twrins.Htwr2frac,\
                        myjckt.LoadFrameOuts.Frameouts.Freqs,myjckt.LoadFrameOuts.tower_utilization.GLUtil,myjckt.LoadFrameOuts.tower_utilization.EUshUtil,\
                        myjckt.LoadFrameOuts2.tower_utilization.GLUtil,myjckt.LoadFrameOuts2.tower_utilization.EUshUtil,\
                        myjckt.PreBuild.wbase,myjckt.Xbraces.bay_hs,myjckt.Xbraces.bay_bs],\
                        fp)

    #Plot
    PlotJacket(myjckt,util=True,savefileroot=ntpath.join(outdir,casename))

    return myjckt,casename
#__________________________________________________________#



def JcktWrapper(x,myjckt,desvarmeans,desvarbds):

    """This function builds up the actual model and calculates Jacket stresses and utilization.
    INPUT
        x         -list(N), as in DesVars with that order \n
        myjckt    -assmebly of jacketSE.\n
        desvarmeans -array(N), average values for each design variable.\n
        #via global
        towerfix -boolean, if true the tower is fixed, no optimization on it.
        """
    global DTRsdiff,MDAOswitch2,towerfix, offset
    global  xlast,f1,max_GLUtil,max_EUUtil,max_KjntUtil,max_XjntUtil,max_tutil,max_cbutil,\
            max_GLUtil2,max_EUUtil2,max_KjntUtil2,max_XjntUtil2,max_tutil2,max_cbutil2,\
            MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2
#   x 0     1     2         3     4    5       6       7        8        9        10    11    12   13   14   15        16          17
#   batter,Dp,   tp,       Lp,   Dleg,tleg,    Dbrc   tbrc    Dmdbrc   tmdbrc     Dgir,tgir   Db,DTRb,   Dt,DTRt,     H2frac, dck_wdth_fac

    desvarmeans2=np.ones(desvarmeans.size)
    if MDAOswitch2:
       desvarmeans2=desvarmeans

    myjckt.JcktGeoIn.batter=x[0]*desvarmeans2[0]


    myjckt.Pileinputs.Dpile=x[1]*desvarmeans2[1]
    myjckt.Pileinputs.tpile=x[2]*desvarmeans2[2]
    myjckt.Pileinputs.Lp   =x[3]*desvarmeans2[3]

    myjckt.leginputs.Dleg  =np.asarray([x[4]*desvarmeans2[4]]).repeat(myjckt.JcktGeoIn.nbays+1)
    myjckt.leginputs.tleg  =np.asarray([x[5]*desvarmeans2[5]]).repeat(myjckt.JcktGeoIn.nbays+1)
    myjckt.leginputs.Dleg0 =0.  #Set to 0 for external optimizer
    myjckt.leginputs.tleg0 =0.  #Set to 0 for external optimizer

    myjckt.Xbrcinputs.precalc=False
    myjckt.Xbrcinputs.Dbrc    =np.asarray([x[6]*desvarmeans2[6]]).repeat(myjckt.JcktGeoIn.nbays)
    myjckt.Xbrcinputs.tbrc    =np.asarray([x[7]*desvarmeans2[7]]).repeat(myjckt.JcktGeoIn.nbays)
    myjckt.Xbrcinputs.Dbrc0=0#Set to 0 for external optimizer
    myjckt.Xbrcinputs.tbrc0=0#Set to 0 for external optimizer

    myjckt.Mbrcinputs.precalc=False
    myjckt.Mbrcinputs.Dbrc_mud =x[8]*desvarmeans2[8]
    myjckt.Mbrcinputs.tbrc_mud =x[9]*desvarmeans2[9]

    myjckt.TPinputs.Dgir  =x[10]*desvarmeans2[10]
    myjckt.TPinputs.tgir  =x[11]*desvarmeans2[11]
    varcnt=12 #for indices later
    cnt2=0 #for indices later, DTRsdiff does not matter if towerfixed, so this goes to 0
    if not(towerfix):
        cnt2=int(DTRsdiff)
        myjckt.Twrinputs.Db    =x[12]*desvarmeans2[12]
        myjckt.Twrinputs.DTRb  =x[13]*desvarmeans2[13]
        myjckt.Twrinputs.Dt    =x[14]*desvarmeans2[14]
        myjckt.Twrinputs.DTRt  =myjckt.Twrinputs.DTRb#myjckt.Twrinputs.DTRb.copy()   #can use DTRt=DTRb here for simplicity
        if DTRsdiff:
            myjckt.Twrinputs.DTRt  =x[15]*desvarmeans2[15]

        myjckt.Twrinputs.Htwr2frac =x[15+cnt2]*desvarmeans2[15+cnt2]
        varcnt=16


    myjckt.JcktGeoIn.dck_width =x[varcnt+cnt2]*desvarmeans2[varcnt+cnt2]*myjckt.Twrinputs.Db  #Deck WIdth

    #Run the assembly and get main output
    myjckt.run()

    #Get Frame3dd mass
    mass=myjckt.LoadFrameOuts.Frameouts.mass[0] +myjckt.LoadFrameOuts.Mpiles  #Total structural mass
    #Get Frame3dd-calculated 1st natural frequency
    f1=myjckt.LoadFrameOuts.Frameouts.Freqs[0]

    #Get Model calculated TP mass
    TPmass=myjckt.TP.TPouts.mass

    #Get Utilizations
    max_GLUtil=np.nanmax(myjckt.LoadFrameOuts.tower_utilization.GLUtil)
    max_GLUtil2=np.nanmax(myjckt.LoadFrameOuts2.tower_utilization.GLUtil)
    #max_GLUtil=myjckt.tower_utilization.GLUtil
    max_EUUtil=np.nanmax(myjckt.LoadFrameOuts.tower_utilization.EUshUtil)
    max_EUUtil2=np.nanmax(myjckt.LoadFrameOuts.tower_utilization.EUshUtil)
    #Member checks
    max_tutil=np.nanmax(myjckt.LoadFrameOuts.jacket_utilization.t_util)
    max_cbutil=np.nanmax(myjckt.LoadFrameOuts.jacket_utilization.cb_util)
    max_tutil2=np.nanmax(myjckt.LoadFrameOuts2.jacket_utilization.t_util)
    max_cbutil2=np.nanmax(myjckt.LoadFrameOuts2.jacket_utilization.cb_util)

    #Joint checks
    max_XjntUtil=np.nanmax(myjckt.LoadFrameOuts.jacket_utilization.XjntUtil)
    max_KjntUtil=np.nanmax(myjckt.LoadFrameOuts.jacket_utilization.KjntUtil)
    max_XjntUtil2=np.nanmax(myjckt.LoadFrameOuts2.jacket_utilization.XjntUtil)
    max_KjntUtil2=np.nanmax(myjckt.LoadFrameOuts2.jacket_utilization.KjntUtil)

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
    Lp0rat=myjckt.LoadFrameOuts.Lp0rat
    Lp0rat2=myjckt.LoadFrameOuts2.Lp0rat
    #__________________________________________#

    #calc width at seabed proportional to stiffness

    xlast=x.copy() #update the latest set of input params to the current


    print('Jwrapper SOLUTION: bat={:5.2f}, Dpile={:5.2f}, tpile={:5.3f}, Lp={:5.1f} Dleg{:5.2f}, tleg{:5.3f}').\
            format(myjckt.JcktGeoIn.batter,myjckt.Piles.Pileinputs.Dpile,myjckt.Piles.Pileinputs.tpile,myjckt.Piles.Pileinputs.Lp,myjckt.Legs.leginputs.Dleg[0],myjckt.Legs.leginputs.tleg[0])
    print('        dck_width ={:5.2f},    Dbrc={:5.2f}, tbrc={:5.3f}, Dmudbrc={:5.2f}, tmudbrc={:5.3f} \n').\
            format(myjckt.dck_width, myjckt.Xbrcinputs.Dbrc[0],myjckt.Xbrcinputs.tbrc[0],myjckt.Mbrcinputs.Dbrc_mud,myjckt.Mbrcinputs.tbrc_mud)
    print('        footprint ={:5.2f} \n'.format(myjckt.PreBuild.wbase))
    print('from Jwrapper Db={:5.2f}, DTRb={:5.2f}, Dt={:5.2f}, DTRt={:5.2f},H2twrfrac={:5.2f}, Dgir={:5.2f},tgir={:5.3f}, Twrmass={:6.3f}, PilesMass ={:6.3f}, TPmass= {:8.3e}, Frame3DD+Piles Totmass={:10.3f}'.\
            format(myjckt.Twrinputs.Db,myjckt.Twrinputs.DTRb,myjckt.Twrinputs.Dt,\
                   myjckt.Twrinputs.DTRt,myjckt.Twrinputs.Htwr2frac,myjckt.TPinputs.Dgir,myjckt.TPinputs.tgir,myjckt.Tower.Twrouts.mass,myjckt.Mpiles, myjckt.TP.TPouts.mass,mass))
    print(' \n')

    sys.stdout.flush()  #This for peregrine
    return mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2



def objfunc(x,myjckt,desvarmeans,desvarbds):
    global towerfix
    mass = JcktWrapper(x,myjckt,desvarmeans,desvarbds)[0]/1.5e6
    #          x=  [ batter,  Dpile,    tpile,        Lp,   Dleg,     tleg,       Dbrc,   tbrc,     Dbrc_mud,   tbrc_mud,   Dgir,      tgir,      Db,   DTRb   Dt,   DTRt   Htwr2fac        dck_widthfact]

    cnt= 6+ 11+7 #counter for number of constraints- initial value for fixed tower and 1 dlc

    if not(towerfix):
        cnt += 2

    if myjckt.twodlcs:
        cnt += 4 +2+1

    cnstrts=[0.0]*cnt #given as negatives, since PYOPT wants <0
    #Note the minus signs here for pyOPT's sake (it wants constraints <0, as opposed to regular python cobyla)

    cnstrts[0]=-f0Cnstrt1(x,myjckt,desvarmeans,desvarbds)
    cnstrts[1]=-f0Cnstrt2(x,myjckt,desvarmeans,desvarbds)

    cnstrts[2]=-cbCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnstrts[3]=-tCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnstrts[4]=-KjntCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnstrts[5]=-XjntCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnt=6

    if myjckt.twodlcs:
        cnstrts[6]=-cbCnstrt2(x,myjckt,desvarmeans,desvarbds)
        cnstrts[7]=-tCnstrt2(x,myjckt,desvarmeans,desvarbds)
        cnstrts[8]=-KjntCnstrt2(x,myjckt,desvarmeans,desvarbds)
        cnstrts[9]=-XjntCnstrt2(x,myjckt,desvarmeans,desvarbds)
        cnt=10

    if not(towerfix):
        cnstrts[cnt]  =-GLCnstrt(x,myjckt,desvarmeans,desvarbds)
        cnstrts[cnt+1]=-EUCnstrt(x,myjckt,desvarmeans,desvarbds)
        cnt += 2
        if myjckt.twodlcs:
            cnstrts[cnt]  =-GLCnstrt2(x,myjckt,desvarmeans,desvarbds)
            cnstrts[cnt+1]=-EUCnstrt2(x,myjckt,desvarmeans,desvarbds)
            cnt += 2

    cnstrts[cnt]=-XbCrit01(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+1]=-XbCrit02(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+2]=-XbCrit03(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+3]=-XbCrit04(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+4]=-XbCrit05(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+5]=-MbCrit01(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+6]=-MbCrit02(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+7]=-MbCrit03(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+8]=-MbCrit04(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+9]=-MbCrit05(x,myjckt,desvarmeans,desvarbds)

    cnstrts[cnt+10]=-LpCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnt += 11
    if myjckt.twodlcs:
        cnstrts[cnt]=-LpCnstrt2(x,myjckt,desvarmeans,desvarbds)
        cnt +=1

    cnstrts[cnt]=-Dleg2BrcCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+1]=-Dleg2MudCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+2]=-Dleg2tlegCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+3]=-Dbrc2tbrcCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+4]=-Dmud2mudCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+5]=-ftprintCnstrt(x,myjckt,desvarmeans,desvarbds)
    cnstrts[cnt+6]=-NorsokCnstrt(x,myjckt,desvarmeans,desvarbds)

    #cnstrts[26]=-dckwidthCnstrt1(x,myjckt,desvarmeans,desvarbds)
    #cnstrts[27]=-dckwidthCnstrt2(x,myjckt,desvarmeans,desvarbds)

    fail=0
    return mass,cnstrts,fail

def f0Cnstrt1(x,myjckt,desvarmeans,desvarbds):  #f1>f0
    global xlast,f1

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=(f1-f0)/f0
    print('f0Cnstrt1=',cnstrt)
    return cnstrt

def f0Cnstrt2(x,myjckt,desvarmeans,desvarbds): #f1<(f0*(1+f0eps))
    global xlast,f1,f0eps

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=(-f1+f0*(1+f0eps))/f0
    print('f0Cnstrt2=',cnstrt)
    return cnstrt

def cbCnstrt(x,myjckt,desvarmeans,desvarbds): #cb<1
    global xlast,f1,max_cbutil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_cbutil
    print('cbCnstrt=',cnstrt)
    return cnstrt
def cbCnstrt2(x,myjckt,desvarmeans,desvarbds): #cb<1
    global xlast,f1,max_cbutil2

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_cbutil2
    print('cbCnstrt2=',cnstrt)
    return cnstrt

def tCnstrt(x,myjckt,desvarmeans,desvarbds): #tUtil<1
    global xlast,max_tutil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_tutil
    print('tutil constraint=',cnstrt)
    return cnstrt
def tCnstrt2(x,myjckt,desvarmeans,desvarbds): #tUtil<1
    global xlast,max_tutil2

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_tutil2
    print('tutil constraint2=',cnstrt)
    return cnstrt


def KjntCnstrt(x,myjckt,desvarmeans,desvarbds): #KUtil<1
    global xlast,max_KjntUtil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_KjntUtil
    print('KjntCnstrt=',cnstrt)
    return cnstrt
def KjntCnstrt2(x,myjckt,desvarmeans,desvarbds): #KUtil<1
    global xlast,max_KjntUtil2

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_KjntUtil2
    print('KjntCnstrt2=',cnstrt)
    return cnstrt

def XjntCnstrt(x,myjckt,desvarmeans,desvarbds): #XUtil<1
    global xlast,max_XjntUtil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_XjntUtil
    print('XjntCnstrt=',cnstrt)
    return cnstrt
def XjntCnstrt2(x,myjckt,desvarmeans,desvarbds): #XUtil<1
    global xlast,max_XjntUtil2

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_XjntUtil2
    print('XjntCnstrt2=',cnstrt)
    return cnstrt

def GLCnstrt(x,myjckt,desvarmeans,desvarbds): #GLUtil<1
    global xlast,max_GLUtil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_GLUtil
    print('GL constraint=',cnstrt)
    return cnstrt
def GLCnstrt2(x,myjckt,desvarmeans,desvarbds): #GLUtil<1
    global xlast,max_GLUtil2

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_GLUtil2
    print('GL constraint2=',cnstrt)
    return cnstrt

def EUCnstrt(x,myjckt,desvarmeans,desvarbds): #EUUtil<1
    global xlast,max_EUUtil

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_EUUtil
    print('EU constraint=',cnstrt)
    return cnstrt
def EUCnstrt2(x,myjckt,desvarmeans,desvarbds): #EUUtil<1
    global xlast,max_EUUtil2

    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    cnstrt=1.-max_EUUtil2
    print('EU constraint2=',cnstrt)
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

#Leg Bound Constraints
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
    global jcktDTRmin,MDAOswitch2
    idx0=4
    idx1=5
    factor0=1
    factor1=1
    if MDAOswitch2:
        factor0=desvarmeans[idx0]
        factor1=desvarmeans[idx1]
    cnstrt=(x[idx0]*factor0/(x[idx1]*factor1) - jcktDTRmin)/jcktDTRmin  #this to work when the bounds are the same for the 2 variables
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
    global jcktDTRmin, MDAOswitch2
    idx0=6
    idx1=7
    factor0=1
    factor1=1
    if MDAOswitch2:
        factor0=desvarmeans[idx0]
        factor1=desvarmeans[idx1]
    cnstrt=(x[idx0]*factor0/(x[idx1]*factor1) - jcktDTRmin)/jcktDTRmin  #this to work when the bounds are the same for the 2 variables
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
    global jcktDTRmin,MDAOswitch2
    idx0=8
    idx1=9
    factor0=1
    factor1=1
    if MDAOswitch2:
        factor0=desvarmeans[idx0]
        factor1=desvarmeans[idx1]
    cnstrt=(x[idx0]*factor0/(x[idx1]*factor1) - jcktDTRmin)/jcktDTRmin#this to work when the bounds are the same for the 2 variables
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
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Xbrc constraint01=', XBrcCrit01)
    return XBrcCrit01

def XbCrit02(x,myjckt,desvarmeans,desvarbds): #XbrcCrit02
    global xlast,XBrcCrit02
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Xbrc constraint02=', XBrcCrit02)
    return XBrcCrit02

def XbCrit03(x,myjckt,desvarmeans,desvarbds): #XbrcCrit03
    global xlast,XBrcCrit03
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Xbrc constraint03=', XBrcCrit03)
    return XBrcCrit03

def XbCrit04(x,myjckt,desvarmeans,desvarbds): #XbrcCrit03
    global xlast,XBrcCrit04
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Xbrc constraint04=', XBrcCrit04)
    return XBrcCrit04

def XbCrit05(x,myjckt,desvarmeans,desvarbds): #XbrcCrit05
    global xlast,XBrcCrit05
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

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
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Mbrc constraint01=', MudCrit01)
    return MudCrit01

def MbCrit02(x,myjckt,desvarmeans,desvarbds): #XbrcCrit02
    global xlast,MudCrit02
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Mbrc constraint02=', MudCrit02)
    return MudCrit02

def MbCrit03(x,myjckt,desvarmeans,desvarbds): #XbrcCrit03
    global xlast,MudCrit03
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Mbrc constraint03=', MudCrit03)
    return MudCrit03

def MbCrit04(x,myjckt,desvarmeans,desvarbds): #XbrcCrit03
    global xlast,MudCrit04
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

        xlast=x.copy()
    print('Mbrc constraint04=', MudCrit04)
    return MudCrit04

def MbCrit05(x,myjckt,desvarmeans,desvarbds): #XbrcCrit05
    global xlast,MudCrit05
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)

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

#deck width
def dckwidthCnstrt1(x,myjckt,desvarmeans,desvarbds): #Minimum deck width=2*Db
    global    DTRsdiff,towerfix
    idx=16+int(DTRsdiff)   #deck width factor
    if towerfix:
        idx=12   #deck width factor
    cnstrt=mincnstrt(x,idx,desvarmeans,desvarbds)
    print('dck_width cnstrt1=',cnstrt )
    return cnstrt
def dckwidthCnstrt2(x,myjckt,desvarmeans,desvarbds): #Max deck width=3*Db
    global    DTRsdiff,towerfix
    idx=16+int(DTRsdiff)   #deck width factor
    if towerfix:
        idx=12   #deck width factor
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('dck_width cnstrt2=',cnstrt )
    return cnstrt

def ftprintCnstrt(x,myjckt,desvarmeans,desvarbds): #Max footprint
    global    mxftprint
    cnstrt=(mxftprint-myjckt.PreBuild.wbase)/mxftprint
    print('footprint cnstrt=',cnstrt )
    return cnstrt

#Embedment length constraint
def LpCnstrt(x,myjckt,desvarmeans,desvarbds):  #Maximum Htwr2 < Htwr/4
    global xlast,Lp0rat
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)
        xlast=x.copy()
    print('Lp0rat constraint=',Lp0rat)
    return Lp0rat
def LpCnstrt2(x,myjckt,desvarmeans,desvarbds):  #Maximum Htwr2 < Htwr/4
    global xlast,Lp0rat2
    if xlast==None or np.any([x != xlast]):
        #print('call Jwrapper from const')
        mass,f1,max_tutil,max_tutil2,max_cbutil,max_cbutil2,max_KjntUtil,max_KjntUtil2,max_XjntUtil,max_XjntUtil2,max_GLUtil,max_GLUtil2,max_EUUtil,max_EUUtil2,\
        MudCrit01,MudCrit02,MudCrit03,MudCrit04,MudCrit05,\
        XBrcCrit01,XBrcCrit02,XBrcCrit03,XBrcCrit04,XBrcCrit05,Lp0rat,Lp0rat2=JcktWrapper(x,myjckt,desvarmeans,desvarbds)
        xlast=x.copy()
    print('Lp0rat2 constraint=',Lp0rat2)
    return Lp0rat2

#Embedment length constraint2
def LpCnstrt3(x,myjckt,desvarmeans,desvarbds):  #Maximum Lp < Lpmax
    idx=3
    cnstrt=maxcnstrt(x,idx,desvarmeans,desvarbds)
    print('Lp constraint3=',cnstrt)
    return cnstrt

    #Function to minimize
def mass(x,myjckt,desvarmeans,desvarbds):
    """This function assembles output data from the assembly in terms of mass: \n
        it is mostly frame3dd mass+piles'' mass."""
    return  JcktWrapper(x,myjckt,desvarmeans,desvarbds)[0]/1.5e6

def maxcnstrt(x,idx,desvarmeans,desvarbds):
    global MDAOswitch2
    factor=1
    if MDAOswitch2:
        factor=desvarmeans[idx]
    return (-x[idx]*factor+desvarbds[idx,1])/desvarmeans[idx]
def mincnstrt(x,idx,desvarmeans,desvarbds):
    global MDAOswitch2
    factor=1
    if MDAOswitch2:
        factor=desvarmeans[idx]
    return  (x[idx]*factor-desvarbds[idx,0])/desvarmeans[idx]


#______________________________________________________________________________#

def main(prmsfile='MyJacketInputs.py',MDAOswitch='pyCobyla', tablefile=[],f0epsilon=0.1,jcktDTRmin_set=22., multi=False,towerfixed=False):

    #NOTE SET THE UPPER KEYWORDS TOO FOR MULTICASE TABLE

    """

    INPUT THROUGH KEYWORDS INTO MAIN:
    towerfixed      -Boolean, Set it to True if tower does not need be optimized. \n
        The following are needed here, as keywords above ^^^,  only if multiplecase table is going to be run, else ignore.\n
    multi           -Boolean,     Set it to True if multiple case file, else False.
    The following two are generally not used, but they can be used to override the settings in prmsfile for some cases. Left here as legacy.
    f0epsilon       -Float,       Set this to a value such that the f0*(1+f0epsilon) will not be exceeded: This is to be fixed here only in case of multiple-case tables. Note it is overriden belo, so in case modify, this was for LCOE project. \n
    jcktDTRmin_set  -Float,       Set this one to the minimum jacket member DTR allowed, mostly for 60+waterdepths: This is to be fixed here only in case of multiple-case tables.\n

    """


#take care of possible external (outside of IDE) run

    if len(sys.argv)>1 and len(sys.argv)<5: #This means individual case
        prmsfile=sys.argv[1]
        MDAOswitch=sys.argv[2].lower()
        if len(sys.argv)==4:
            towerfixed=(sys.argv[3].lower()=='true')
    elif len(sys.argv)>1:       #This means multicase table
        multi=True
        guessfromoutfile=[]

        prmsfile=sys.argv[1]
        tablefile=sys.argv[2]
        casenostart=int(sys.argv[3])
        casenoends=int(sys.argv[4])
        xlsfilename=sys.argv[5]
        MDAOswitch='pyCobyla'
        if len(sys.argv)>6:
            MDAOswitch= sys.argv[6].lower()
            if len(sys.argv)>7:
                towerfixed=(sys.argv[7].lower()=='true')
                if len(sys.argv)>8:
                    guessfromoutfile=sys.argv[8]



    #Now check whether we are not passing stuff from command line and we want to assign defaults

    elif not multi:  # INDIVIDUAL FILE
        prmsfile=r'MyJacketInputs.py'
        MDAOswitch='pyCobyla'

    else:           #MULTICASE TABLE FILE
        prmsfile=r'SetJacketInputsPeregrine.py'
        tablefile=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\testmatrix.dat'
      #  tablefile=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\testmatrixspec.dat'
        casenostart=82  #UHREACT
        casenoends=82

        MDAOswitch='pyCobyla'
        xlsfilename=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\outputs.xls'
     #   xlsfilename=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\outputs.xls'

        #guessfromoutfile='D:\\RRD_ENGINEERING\\PROJECTS\\NREL\\OFFSHOREWIND\\LCOE_ANALYSIS\\CASES_55_63_tobereplaced\\55_R2D0H0M0T0S0_out.dat'
        #guessfromoutfile='D:\\RRD_ENGINEERING\\PROJECTS\\NREL\\OFFSHOREWIND\\LCOE_ANALYSIS\\CASES_10_18_tobereplaced\\10_R0D1H0M0T0S0_out.dat'
        #guessfromoutfile=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\COBYLA\08_R0D0H2M1T0S0_out.dat'


#________________RUN NOW________________#
    if not multi: # INDIVIDUAL FILE
        myjckt,casename=JcktOpt(prmsfile,MDAOswitch=MDAOswitch,towerfixed=towerfixed)
        sys.stdout.flush()  #This for peregrine

    else: #MULTICASE TABLE FILE
        ###f0epsilon=0.05 #upper f0 allowed  THIS WAS FOR LCOE's proje
        for caseno in range(casenostart,casenoends+1):
            ###if caseno>27:  THIS IS FOR LCOE PROJECT FY14
            ###    f0epsilon=0.35
            print('#____________________________________________#\n')
            print(('#JacketOpt  NOW PROCESSING CASE No. {:d} #\n').format(caseno) )
            print('#____________________________________________#\n')
            myjckt,casename=JcktOpt(prmsfile,MDAOswitch=MDAOswitch,tablefile=tablefile, caseno=caseno,xlsfilename=xlsfilename,guessfromoutfile=guessfromoutfile, f0epsset=f0epsilon, jcktDTRminset=jcktDTRmin_set,towerfixed=towerfixed)
            ###guessfromoutfile=os.path.join(os.path.dirname(prmsfile),casename + '_out.dat')  #THIS IS TO BE USED FOR SIMILAR CASES IN A ROW
            sys.stdout.flush()  #This for peregrine



##___________________________________________________________##
    ##SAMPLE CALL FROM OUTSIDE IDE with 1 individual file
    ##python JacketOpt_Py_MDAOopt.py C:\RRD\PYTHON\WISDEM\JacketSE\src\jacketse\MyJacketInputs.py  pyCobyla
    # where                              main input file                                           MDAOswitch    : in this case pyOPTcobyla will be used (No OMDAO)
    ##SAMPLE CALL FROM OUTSIDE ENVIRONMENT
    ##python JacketOpt_Py_MDAOopt.py D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\SetJacketInputsPeregrine.py D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\testmatrix.dat 55 55 D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\output.xls pySNOPT false D:\\RRD_ENGINEERING\\PROJECTS\\NREL\\OFFSHOREWIND\\LCOE_ANALYSIS\\CASES_55_63_tobereplaced\\55_R2D0H0M0T0S0_out.dat
    #OR ANOTHER ONE
    ##python JacketOpt_Py_MDAOopt.py C:\PROJECTS\OFFSHORE_WIND\SEJacketTower\SITEdata\SetJacketInputsPeregrine.py C:\PROJECTS\OFFSHORE_WIND\SEJacketTower\SITEdata\SiteData_PythonInput.xlsx 3 3 C:\PROJECTS\OFFSHORE_WIND\SEJacketTower\SITEdata\output.xls md_pySNOPT False

    ##python JacketOpt_Py_MDAOopt.py D:\RRD_ENGINEERING\PROJECTS\BOOK_CHAPTER\CALCS\LongIsland\SetJacketInputsPeregrine.py D:\RRD_ENGINEERING\PROJECTS\BOOK_CHAPTER\CALCS\LongIsland\SiteData_LI.xlsx 3 3 D:\RRD_ENGINEERING\PROJECTS\BOOK_CHAPTER\CALCS\LongIsland\LIoutput.xls md_pySNOPT False

    ##python JacketOpt_Py_MDAOopt.py D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\SEjacketTower\SITEdata\SetJacketInputsPeregrine.py D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\SEjacketTower\SITEdata\SiteData_PythonInput.xlsx 3 3 D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\SEjacketTower\SITEdata\output.xls md_pySNOPT
    ##python JacketOpt_Py_MDAOopt.py D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\SEjacketTower\SITEdata\SetJacketInputsPeregrine.py D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\SEjacketTower\SITEdata\SiteData_PythonInput.xlsx 3 3 D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\SEjacketTower\SITEdata\output.xls extcobyla
##___________________________________________________________##




if __name__ == '__main__':
    #This is how you call this function
    main()
