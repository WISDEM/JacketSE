#-------------------------------------------------------------------------------
# Name:        JacketOptMDAO.py
# Purpose:     It solves for minimum mass problem with constraints on max utilization.
#
# Author:      rdamiani
#
# Created:     4/25/2014 based off of JacketOpt.py
# Copyright:   (c) rdamiani 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np
import scipy.optimize
import os

from openmdao.main.api import Assembly
from openmdao.main.datatypes.api import Int, Float, Array, VarTree, Bool, Dict,Instance
from openmdao.lib.drivers.api import COBYLAdriver
from openmdao.lib.casehandlers.api import DumpCaseRecorder
from openmdao.main.api import enable_console
#enable_console()
import logging
logging.getLogger().setLevel(logging.DEBUG)

#from jacket import JacketAsmly,JcktGeoInputs,LegGeoInputs,XBrcGeoInputs,MudBrcGeoInputs,HBrcGeoInputs,PileGeoInputs,TPGeoInputs,TwrGeoInputs,RNAprops,TPlumpMass,WaterInputs,WindInputs,SoilGeoInputs,Frame3DDaux,IEC_PSFS
import SetJacketInputs

def main(recordfile):

    #1.Set bounds for optimization variables

    #          x=  [ batter,  Dpile,    tpile,   Lp,   Dleg,     tleg,       Dbrc,   tbrc,     Dbrc_mud,   tbrc_mud,    Db,   DTRb   Dt,   DTRt   Htwr2fac ]
    MnCnst=np.array([ 7.,      1.,    1.*0.0254, 30.,   1.,       1.*0.0254,  1.,    1.*0.0254,   1.,     1.*0.0254,    5.,   120.,  3.,   120.,     0. ])
    MxCnst=np.array([ 30.,     2.5,   5.*0.0254, 50.,   2.5,      5.*0.0254,  2.,    5.*0.0254,   2.,     5.*0.0254,    7.,   200.,  4.,   200.,   0.25 ])
    avgcnst=(MnCnst+MxCnst)/2.
    factors=1./avgcnst

    #2. Set up main variables and create assembly with main inputs via SetJacketInputs.py
    myjckt=SetJacketInputs.main(avgcnst[0],avgcnst[1],avgcnst[2],avgcnst[3],avgcnst[4],avgcnst[5],avgcnst[6],avgcnst[7], avgcnst[8],avgcnst[9],avgcnst[10],avgcnst[11],avgcnst[12],avgcnst[13],avgcnst[14])


   #3. Replace driver in main assembly and specify optimizer parameters
    myjckt.replace('driver',COBYLAdriver())

    myjckt.driver.iprint = 1
    #myjckt.driver.rhoend = 0.1
    myjckt.driver.rhobeg=0.1
    myjckt.driver.disp=1

    #4. Objective #TO DO ADD PILE MASS
    myjckt.driver.add_objective('(FrameOut.Frameouts_outs.mass[0]+Embedment.Mpiles)/1.e6')

    #5. Design Variables

 #   myjckt.driver.add_parameter('JcktGeoIn.batter',   low=MnCnst[0]/factors[0], high=MxCnst[0]/factors[0],scaler=factors[0])
    #myjckt.driver.add_parameter('Pileinputs.Dpile',   low=MnCnst[1]/factors[1], high=MxCnst[1]/factors[1],scaler=factors[1])
    #myjckt.driver.add_parameter('Pileinputs.tpile',   low=MnCnst[2]/factors[2], high=MxCnst[2]/factors[2],scaler=factors[2])
    #myjckt.driver.add_parameter('Pileinputs.Lp',      low=MnCnst[3]/factors[3], high=MxCnst[3]/factors[3],scaler=factors[3])
   # myjckt.driver.add_parameter('leginputs.Dleg',     low=MnCnst[4]/factors[4], high=MxCnst[4]/factors[4],scaler=factors[4])
   # myjckt.driver.add_parameter('leginputs.tleg',     low=MnCnst[5]/factors[5], high=MxCnst[5]/factors[5],scaler=factors[5])
   # myjckt.driver.add_parameter('Xbrcinputs.Dbrc',    low=MnCnst[6]/factors[6], high=MxCnst[6]/factors[6],scaler=factors[6])
   # myjckt.driver.add_parameter('Xbrcinputs.tbrc',    low=MnCnst[7]/factors[7], high=MxCnst[7]/factors[7],scaler=factors[7])
    #myjckt.driver.add_parameter('Mbrcinputs.Dbrc_mud',low=MnCnst[8]/factors[8], high=MxCnst[8]/factors[8],scaler=factors[8])
    #myjckt.driver.add_parameter('Mbrcinputs.tbrc_mud',low=MnCnst[9]/factors[9], high=MxCnst[9]/factors[9],scaler=factors[9])
    myjckt.driver.add_parameter('Twrinputs.Db',        low=MnCnst[10], high=MxCnst[10])#/factors[10],scaler=factors[10]
    myjckt.driver.add_parameter('Twrinputs.DTRb',      low=MnCnst[11], high=MxCnst[11])# /factors[11],scaler=factors[11]
    myjckt.driver.add_parameter('Twrinputs.Dt',        low=MnCnst[12], high=MxCnst[12])# /factors[12],scaler=factors[12]
    myjckt.driver.add_parameter('Twrinputs.DTRt',      low=MnCnst[13], high=MxCnst[13])# /factors[13],scaler=factors[13]
    myjckt.driver.add_parameter('Twrinputs.Htwr2frac', low=MnCnst[14], high=MxCnst[14])# /factors[14],scaler=factors[14]

    #6. Constraitns
    myjckt.driver.add_constraint('FrameOut.Frameouts_outs.Freqs[0] >=0.25')
    myjckt.driver.add_constraint('max(Utilization.tower_utilization.GLUtil) <=1.0')

    myjckt.driver.add_constraint('max(Utilization.tower_utilization.EUshUtil) <=1.0')
   # myjckt.driver.add_constraint('Utilization.jacket_utilization.t_util <=1.0')
   # myjckt.driver.add_constraint('Utilization.jacket_utilization.cb_util <=1.0')
   # myjckt.driver.add_constraint('Utilization.jacket_utilization.KjntUtil <= 1.0')
   # myjckt.driver.add_constraint('Utilization.jacket_utilization.XjntUtil <= 1.0')

    #7 recorder

    #if (isinstance(recordfile,str) and (os.path.exists(recordfile) or  os.access(os.path.dirname(recordfile), os.W_OK))):
    fileID=open(recordfile,'a')
    myjckt.driver.recorders = [DumpCaseRecorder(fileID),DumpCaseRecorder()]

    import time
    tt = time.time()

    myjckt.run()

    print "\n"
    print "Minimum found at (%f, %f, %f)" % (myjckt.Tower.Twrins.Db,myjckt.Tower.Twrins.Dt,myjckt.Tower.Twrins.Htwr2frac)
    print "Minimum found at DTRb DTRt(%f, %f)" % (myjckt.Tower.Twrins.DTRb,myjckt.Tower.Twrins.DTRt)
    print "Minimum found at Freq %f"  % (myjckt.FrameOut.Frameouts_outs.Freqs[0])
    print "Minimum found at GLutil EUutil %f %f"  % (np.max(myjckt.tower_utilization.GLUtil),np.max(myjckt.tower_utilization.EUshUtil))
    print "Elapsed time: ", time.time()-tt, "seconds"
    print "Execution count: ", myjckt.exec_count
#    res=JacketOpt(f0,f0eps,MnCnst,MxCnst,infile,InFrame3DD,wwinputfile,SoilInfo) #depending on whether or not AF is needed

    #5. store data for future use import cPickle as pickle; pickle.dump(res,open("C:\PROJECTS\OFFSHORE_WIND\NREL5MW_30mWD_Jacket\PYTHON\JCKdata.p","w"))

    #6. recall data               res=pickle.load(open("C:\PROJECTS\OFFSHORE_WIND\NREL5MW_30mWD_Jacket\PYTHON\JCKdata.p","r"))
    fileID.close()

if __name__ == '__main__':
    #This is how you call this function
    main('jacketopt_test.txt')
