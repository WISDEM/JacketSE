#-------------------------------------------------------------------------------
# Name:        ReadSaveOpt1LineSysEng.py
# Purpose:    It contains functions that can:
#               1. Read from a Table of cases for the design parameters and design var. bounds.
#               2. Save results to a text file in terms of design parameters,
#                  and to an excel file in terms of general data needed to assess the jacket tower.
#
#         Modified version of ReadSaveOpt1Line.py to work with Systems ENgineering Jacket vs Tower Study
# Author:      rdamiani
#
# Created:     11/21/2014
# Copyright:   (c) rdamiani 2014
# Licence:     <APACHE 2014>
#-------------------------------------------------------------------------------
import numpy as np
import os,ntpath
import time
import collections

from xlwt import Workbook,easyxf
from xlrd import open_workbook
from xlutils.copy import copy as xlscopy

from VarTrees import RNAprops
from commonse.GetSheetByName import get_sheet_idx
import SetJacketInputsPeregrine
from PlotJacket import main as PlotJacket
from towerse.DesVarsAux import TwrDesPrms, TwrDesVarBounds
#__________________________________________________________#

class DesVars(object): #Design Variables Grouped in a structure
    """Class containing the values of the optimization design variables."""
#   x 0     1     2         3     4    5       6       7        8        9        10    11    12   13      14      15      16      17
#   batter,Dp,   tp,       Lp,   Dleg,tleg,    Dbrc   tbrc    Dmdbrc   tmdbrc     Dgir,tgir   Db,DTRb,   Dt,    DTRt,     H2frac,  dck_widthfact
    def __init__(self,**kwargs):

        dvars=collections.OrderedDict([ ('batter',7),('Dpile',1.5),('tpile',0.035), ('Lp',40.),    \
              ('Dleg',1.5),('tleg',0.0254),\
              ('Dbrc',0.8,),('tbrc',0.0254),\
              ('Dbrc_mud',1.0),('tbrc_mud',0.0254),\
              ('Dgir',1.),('tgir',0.0254),\
              ('Db',5.),('DTRb',120.),('Dt',3.),('DTRt',120.),('Htwr2frac',0.25),('dck_widthfact',2.)\
               ]) #SI Units
        dvarunits=collections.OrderedDict([('batter','[-]'),('Dpile','[m]'),('tpile','[m]'), ('Lp','[m]'),    \
              ('Dleg','[m]'), ('tleg','[m]'),\
              ('Dbrc','[m]',),('tbrc','[m]'),\
              ('Dbrc_mud','[m]'), ('tbrc_mud','[m]'),\
              ('Dgir','[m]'), ('tgir','[m]'),\
              ('Db','[m]'),   ('DTRb','[-]'),('Dt','[m]'),('DTRt','[-]'),('Htwr2frac','[-]'),('dck_widthfact','[-]')\
               ]) #SI Units

        dvars.update(kwargs) #update in case user put some new params in

        self.dvars=dvars  #This allows me to save the orer of the design variables for loops elsewhere
        self.dvarunits=dvarunits  #This allows me to save the orer of the design variables units for loops elsewhere

        for key in dvars:
            setattr(self,key,dvars[key])

#__________________________________________________________#

class DesVarBounds(DesVars): #Design Variable Bounds
    """Class containing the bounds to the optimization design variables."""
    def __init__(self,**kwargs):

        prms={'batter':np.array([7.,30]),'Dpile':np.array([1.5,2.5]),'tpile':np.array([1.0,2.5]), 'Lp':np.array([25.,60.]),    \
              'Dleg':np.array([1.5,2.5]),'tleg':np.array([0.0254,2.5*0.0254]),\
              'Dbrc':np.array([0.8,2.5]),'tbrc':np.array([0.0254,2.5*0.0254]),\
              'Dbrc_mud':np.array([1.0,2.5]),'tbrc_mud':np.array([0.0254,2.5*0.0254]),\
              'Dgir':np.array([1.,1.2]),'tgir':np.array([0.0254,2.5*0.0254]),\
              'Db':np.array([5.,7.]),'DTRb':np.array([120.,200.]),'Dt':np.array([3.,4.5]),'DTRt':np.array([120.,200.]),'Htwr2frac':np.array([0.,0.25]),'dck_widthfact':np.array([2.0,3.0])\
               } #SI Units

        prms.update(kwargs) #update in case user put some new params in
        for key in prms:  #Initialize material object with these parameters, possibly updated by user' stuff
            setattr(self,key,prms[key])

#__________________________________________________________#

class DesPrms(object): #Design Parameters
    """Class containing the Main Design Parameters (won't be optimized) that may change from case to case."""

    def __init__(self,**kwargs):


        prms={'TurbRating':3., 'wdepth':20.,'HW50':30,'Tp50':12.,'HW50_2':30,'Tp50_2':12.,'dck_botz':16.,'HH':100.,'U50HH':30.,'U50HH_2':70.,'TPlumpmass':200.e3, \
              'RNA_F':np.array([1.5*1700.e3,0.,0.,12564863.93,0.,0.]),'RNA_F2':np.array([700.e3,0.,0.,0.,0.,0.]),\
              'RNAins':RNAprops(),'f0':0.35,'legndiv':3,'nbays':4,'mxftprint':30.} #SI Units excpet for MW for the rating and legndiv

        prms.update(kwargs) #update in case user put some new params in
        for key in prms:  #Initialize material object with these parameters, possibly updated by user' stuff
            setattr(self,key,prms[key])

#__________________________________________________________#

def Recon1LineJckt(casefile,caseno,xlsfilename=[],optfile=[],titlines=3,hdrlines=1,untlines=1,opttitlines=3,opthdrlines=1,optuntlines=1):
    """This function reads results from a fcaseile, in terms of case number and case name, and design parameters,\n
        then it reads optimized design variables from a file "casename+'_out.dat'" and reconstructs the jacket assembly.\n
        Note: the optimization output filename is assumed to be "caseno_casename_out.dat". \n
    INPUT \n
        casefile -string, complete path+filename to table file of cases. \n
        caseno   -int, case number (also row out of the table to be read). \n
    OPTIONALS:
        xlsfilename -string, xcel file name (path included) where to write summary of results for the current case. \n
        optfile     -string, optimization results' file name (path included) where to read configuration geometry for the current case. If left blank it will be grabbed from the casefile directory automatically.\n
        titlines - int, number of lines for the title in the table file \n
        hdrlines - int, number of lines for the variable headers in the table file\n
        untlines - int, number of lines for the units in the table file\n
        opttitlines - int, number of lines for the title in the optimal var file \n
        opthdrlines - int, number of lines for the variable headers in the optimal var file\n
        optuntlines - int, number of lines for the units in the optimal var file\n
    OUTPUT \n
        It creates a plot of the configuration, utilization of tower. \n
        myjckt  -OpenMdao object assembly of JacketSE with the design parameters and optimzation variables found.\n
        """

    desvars=DesVars() #instance of design variables

    #First read design parameters from table file
    Desprms,_,desvarbds,desvarmeans,guesses,casename=ReadTab1Line(casefile,caseno,desvars.dvars.keys(),\
                                                          titlines=titlines,hdrlines=hdrlines,untlines=untlines)

    #Then read optimized design variables from file
    if not(optfile):
        optfile=ntpath.join(ntpath.dirname(casefile),casename+'_out.dat')
    desvars=ReadOptFile(optfile,opttitlines=opttitlines,opthdrlines=opthdrlines,optuntlines=optuntlines)

    #Then build the initial assembly and run it to instantiate everything
    myjckt=SetJacketInputsPeregrine.main(Desprms,desvars)
    myjckt.run()
    #Edit and add a sheet to the xlsfile if requested
    optresults=np.zeros(len(desvars.dvars))
    for ii,key in enumerate(desvars.dvars):
        #print key  #debug
        optresults[ii]= getattr(desvars,key)


    if xlsfilename:
        SaveOpt1Line(ntpath.dirname(casefile),caseno,casename,desvars,optresults,myjckt,xlsfilename,Desprms,xlsfile_only=True)
    #Plot
    PlotJacket(myjckt,util=True,savefileroot=os.path.join(ntpath.dirname(casefile),casename))

    return myjckt
#__________________________________________________________#
#__________________________________________________________#

def ReadOptFile(optfile,opttitlines=3,opthdrlines=1,optuntlines=1):
    """This function reads results from an optimization results' file.\n
        then it reads optimized design variables from a file "casename+'_out.dat'" and reconstructs the jacket assembly.\n
        Note: the optimization output filename is assumed to be "caseno_casename_out.dat". \n
    INPUT \n
        optfile  -string, complete path+filename to optimization file. \n
    OPTIONALS:
        opttitlines - int, number of lines for the title in the optimal var file \n
        opthdrlines - int, number of lines for the variable headers in the optimal var file\n
        optuntlines - int, number of lines for the units in the optimal var file\n
    OUTPUT \n
        desvars -object of class DesVars with optimed values for each variable
        """

    desvars=DesVars() #instance of design variables

    #Then read optimized design variables from file

    fileID=open(optfile,'r')
    #skip title, header, unit lines
    for ii in range(0,opthdrlines+opthdrlines+optuntlines+1):
        _=fileID.readline()

    #read the actual data
    dat=fileID.readline().split()

    for ii,key in enumerate(desvars.dvars):
        #print key  #debug
        setattr(desvars,key,float(dat[ii]))

    fileID.close()

    return desvars
#__________________________________________________________#
def ReadTab1Line(casefile,caseno,desvarnames,towerdata=False,titlines=3,hdrlines=1,untlines=1):
    """This function reads in a text file with case definitions and gets the variable bounds and the case parameters.
    Then it sets adn returns three objects that contain the design variables, design bounds, and design parameters.
    INPUTS \n
        casefile - rstring, complete 'r'+path+name to text or xls file containing cases to run. Format of file is to be found at RunAllJacketTwrs.py. \n
                   If xls, then a sheet called 'TestMatrix' must exist.\n
        caseno   - int, line to read out of the file, corresponding to case number.\n
        desvarnames - list of official names of design variables in the exact and fixed order, usually got from an ordered dictionary from DesVars class.\n
    OPTIONALS:
        towerdata- boolean, If True, then a tower table of experiment is expected.
        titlines - int, number of lines for the title in the text file\n
        hdrlines - int, number of lines for the variable headers in the text file\n
        untlines - int, number of lines for the units in the text file\n
    OUTPUTS \n
        Desprms  -object of class DesPrms, containing design parameter values (they do not change within each case, only from case to case).    \n
        Desbds   -object of class DesVarBounds, containing bounds of design variables (they do not change within each case, only from case to case). \n
        desvarmeans -float(17), mean of the bounds for each design varsariable."""

    ext=ntpath.splitext(casefile)[1].lower()
    if (ext=='.xls') or (ext=='.xlsx'):
        wb=open_workbook(casefile)
        if not towerdata:
            sheet=wb.sheet_by_name('JcktMatrix')
        else:
            sheet=wb.sheet_by_name('TwrMatrix')
        row=titlines+hdrlines+untlines+caseno-1 #-1 for python
        line=[]
        for col in range(sheet.ncols):
            line.append(sheet.cell(row,col).value)
    else:
        fileID=open(casefile,'r')
        lines=fileID.readlines()
        fileID.close()
        #get the case of interest
        line=lines[titlines+hdrlines+untlines+caseno-1].split()#-1 for the usual python stuff

    #Assign design parameters
    if not towerdata:
        Desprms=DesPrms()
        Desvarbds=DesVarBounds()
    else:
        Desprms=TwrDesPrms()
        Desvarbds=TwrDesVarBounds()

    Desprms.TurbRating=float(line[2]) #used only for summary table

    Desprms.wdepth=float(line[3])
    Desprms.HW50=float(line[4])
    Desprms.Tp50=float(line[5])
    Desprms.U50HH=float(line[6])
    Desprms.HW50_2=float(line[7])
    Desprms.Tp50_2=float(line[8])
    Desprms.U50HH_2=float(line[9])

    Desprms.dck_botz=float(line[10])
    Desprms.HH=float(line[11])
    Desprms.TPlumpmass=float(line[12])

    Desprms.RNAins.mass=float(line[13])
    Desprms.RNAins.I=np.array(line[14:20])
    Desprms.RNAins.CMoff=np.array(line[20:23])
    Desprms.RNA_F=np.array(line[23:29])
    Desprms.RNA_F2=np.array(line[29:35])
    Desprms.RNAins.Thoff=np.array(line[35:38])

    Desprms.f0=float(line[38])
    if not(towerdata):
        Desprms.legndiv=int(line[40])
        Desprms.nbays=int(line[41])
        Desprms.mxftprint=float(line[42])
        offset=42 #this is an attempt to miniize issues if the table changes
        mxrange=37
    else:
        offset=40
        mxrange=19
    #Assign bounds

    idx=(np.arange(1,mxrange)+offset)
    idx=idx[::2]
    bds=np.zeros([len(idx),2])
    for ii,indx in enumerate(idx):
        bds[ii,:]=np.array([float(x) for x in line[indx:indx+2]])
        setattr(Desvarbds,desvarnames[ii],bds[ii,:])

    desvarmeans=np.mean(bds,1)

    #Get guesses
    offset2=offset+mxrange
    nvars=len(desvarnames)
    guesses=np.zeros(nvars)
    for ii in range(0,nvars):
        guesses[ii]=float(line[offset2+ii])

    #capture the case name
    casename=str(caseno).zfill(2)+'_'+line[1]

    return Desprms,Desvarbds,bds,desvarmeans,guesses,casename

#__________________________________________________________#
def SaveOpt1Line(outdir,caseno,casename,desvars,rescobyla,myjckt,xlsfilename,Desprms,xlsfile_only=False):
    """This function saves results into a file, in terms of case number and case name, and design variables.\n
    INPUT \n
        outdir      -string, path to where to save a text file \n
        caseno      -int, case number. \n
        casename    -string, case name. \n
        desvars     -object of class DesVars, which contains the ordered dictionary of labels. \n
        rescobyla   -float(nvars), dimensional output from cobyla. \n
        myjckt      -Openmdao assembly, for jacketSE. \n
        xlsfilename -string, xcel file name (path included) where to write summary of results for the current case. \n
        Desprms     -object of class DesPrms, which contains the ordered dictionary of labels and values for current case design parameters. \n
        xlsfile_only -flag, if True the xlsfilename file will be edited, and no textfile will be touched, useful during reconstruction. \n
        """
    if not(xlsfile_only):
        filename=os.path.join(outdir,casename + '_out.dat')

        fileID=open(filename,'w')
        fileID.write('Output of Case {:d} Name {:s}, Generated by JcktOpt1Line/SaveOpt1Line.py on {:s} \n\n'.format(caseno,casename,time.strftime("%c")));
        #Write  design variable optimal values
        ndvars=len(desvars.dvars.keys())

        fileID.write((ndvars*'{:>13s}'+'\n').format(*desvars.dvars.keys())) #Headers
        fileID.write((ndvars*'{:>13s}'+'\n').format(*desvars.dvarunits.values())) #Units
        fileID.write((ndvars*'{:>13f}'+'\n').format(*rescobyla)) #Values

        #Write other stuff of interest that will let me make a good table

    ##    with open(filename,'w') as fileID:
    ##        pickle.dump([res1,myjckt.Frameouts2.mass,myjckt.Tower.Twrouts.mass,,,,\
    ##                    ,
    ##                    myjckt.TPinputs.Dgir,myjckt.TPinputs.tgir,myjckt.TP.TPouts.TPlumpedMass,\
    ##                    ,myjckt.Tower.Twrins.Htwr2frac,\
    ##                    myjckt.Frameouts2.Freqs,myjckt.tower_utilization.GLUtil,myjckt.tower_utilization.EUshUtil,,myjckt.Xbraces.bay_hs,myjckt.Xbraces.bay_bs],\
    ##                    fp)

        fileID.close()

    #Now prepare a sheet in the excel file with things we want to see
    if (isinstance(xlsfilename,str) and (os.path.exists(xlsfilename) or  os.access(os.path.dirname(xlsfilename), os.W_OK))):
        sheetname='caseNO'+str(caseno).zfill(2)
        if os.path.isfile(xlsfilename): #this means already existing file
            book = xlscopy(open_workbook(xlsfilename))
            idx=get_sheet_idx(book,sheetname)
            if idx is not None:
                sheet=book.get_sheet(idx)
            else:
                sheet= book.add_sheet(sheetname)
        else:
            book=  Workbook()
            sheet= book.add_sheet(sheetname)


        #Now Write Title
        sheet.write(0,0,('Optimized Jacket/Tower for Case No. {:d} Name {:s} created on {:s} by Save1OptLine.').format(caseno,casename,time.strftime("%c")),\
                       easyxf('font: name Arial, bold True, underline True'))
        sheet.write(2,1,'PARAMETER',easyxf('font: name Arial, bold True'))
        sheet.write(2,2,'VALUE',easyxf('font: name Arial, bold True'))
        sheet.write(2,3,'UNITS',easyxf('font: name Arial, bold True'))

        #Design and Environmental Parameters
        turbrow=4
        sheet.write(turbrow,0,'Turbine Data',easyxf('font: name Arial, bold True'))
        sheet.write(turbrow+1,1,'Turbine Rating '); sheet.write(turbrow+1,2,Desprms.TurbRating) ;sheet.write(turbrow+1,3,'[MW]')
        sheet.write(turbrow+2,1,'HubHeight '); sheet.write(turbrow+2,2,Desprms.HH) ;sheet.write(turbrow+2,3,'[m]')
        sheet.write(turbrow+3,1,'Unfactored Thrust DLC 1.6'); sheet.write(turbrow+3,2,Desprms.RNA_F[0]/1.e3) ;sheet.write(turbrow+3,3,'[kN]')
        sheet.write(turbrow+4,1,'Unfactored Torque DLC 1.6'); sheet.write(turbrow+4,2,Desprms.RNA_F[3]/1.e3) ;sheet.write(turbrow+4,3,'[kNm]')
        sheet.write(turbrow+5,1,'Unfactored Thrust DLC 6.1'); sheet.write(turbrow+5,2,Desprms.RNA_F2[0]/1.e3) ;sheet.write(turbrow+5,3,'[kN]')
        sheet.write(turbrow+6,1,'Unfactored Torque DLC 6.1'); sheet.write(turbrow+6,2,Desprms.RNA_F2[3]/1.e3) ;sheet.write(turbrow+6,3,'[kNm]')

        sheet.write(turbrow+7,1,'RNA mass'); sheet.write(turbrow+7,2,Desprms.RNAins.mass/1.e3) ;sheet.write(turbrow+7,3,'[tonnes]')
        sheet.write(turbrow+8,1,'RNA CMzoff'); sheet.write(turbrow+8,2,Desprms.RNAins.CMoff[2]) ;sheet.write(turbrow+8,3,'[m]')
        sheet.write(turbrow+9,1,'RNA Thzoff'); sheet.write(turbrow+9,2,Desprms.RNAins.Thoff[2]) ;sheet.write(turbrow+9,3,'[m]')
        sheet.write(turbrow+10,1,'Target f0'); sheet.write(turbrow+10,2,Desprms.f0) ;sheet.write(turbrow+10,3,'[Hz]')

        envrow=turbrow+11
        sheet.write(envrow,0,'Water/Loading Data',easyxf('font: name Arial, bold True'))
        sheet.write(envrow+1,1,'Water Depth'); sheet.write(envrow+1,2,Desprms.wdepth) ;sheet.write(envrow+1,3,'[m]')
        sheet.write(envrow+2,1,'HW50 DLC 1.6'); sheet.write(envrow+2,2,Desprms.HW50) ;sheet.write(envrow+2,3,'[m]')
        sheet.write(envrow+3,1,'Tp50 DLC 1.6'); sheet.write(envrow+3,2,Desprms.Tp50) ;sheet.write(envrow+3,3,'[sec]')
        sheet.write(envrow+4,1,'U50HH DLC 1.6'); sheet.write(envrow+4,2,Desprms.U50HH) ;sheet.write(envrow+4,3,'[m/sec]')
        sheet.write(envrow+5,1,'HW50 DLC 6.1'); sheet.write(envrow+5,2,Desprms.HW50_2) ;sheet.write(envrow+5,3,'[m]')
        sheet.write(envrow+6,1,'Tp50 DLC 6.1'); sheet.write(envrow+6,2,Desprms.Tp50_2) ;sheet.write(envrow+6,3,'[sec]')
        sheet.write(envrow+7,1,'U50HH DLC 6.1'); sheet.write(envrow+7,2,Desprms.U50HH_2) ;sheet.write(envrow+7,3,'[m/sec]')

        #
        pilerow=envrow+9
        sheet.write(pilerow,0,'PILES',easyxf('font: name Arial, bold True'))
        sheet.write(pilerow+1,1,'Pile OD')
        sheet.write(pilerow+2,1,'Pile t')
        sheet.write(pilerow+3,1,'Pile length')
        sheet.write(pilerow+4,1,'Pile Material')
        sheet.write(pilerow+5,1,'Pile Material rho')
        sheet.write(pilerow+6,1,'Pile Material E')
        sheet.write(pilerow+7,1,'Pile Material nu')
        sheet.write(pilerow+8,1,'Pile Material fy')
        sheet.write(pilerow+9,1,'Piles'' mass')
        #values
        sheet.write(pilerow+1,2,myjckt.Piles.Pileinputs.Dpile)
        sheet.write(pilerow+2,2,myjckt.Piles.Pileinputs.tpile)
        sheet.write(pilerow+3,2,myjckt.Piles.Pileinputs.Lp)
        sheet.write(pilerow+4,2,myjckt.Piles.PileObjout.mat[0].matname)
        sheet.write(pilerow+5,2,myjckt.Piles.PileObjout.mat[0].rho)
        sheet.write(pilerow+6,2,myjckt.Piles.PileObjout.mat[0].E)
        sheet.write(pilerow+7,2,myjckt.Piles.PileObjout.mat[0].nu)
        sheet.write(pilerow+8,2,myjckt.Piles.PileObjout.mat[0].fy)
        sheet.write(pilerow+9,2,myjckt.Mpiles/1.e3)
        #units
        sheet.write(pilerow+1,3,'[m]')
        sheet.write(pilerow+2,3,'[m]')
        sheet.write(pilerow+3,3,'[m]')
        sheet.write(pilerow+4,3,'[-]')
        sheet.write(pilerow+5,3,'[kg/m3]')
        sheet.write(pilerow+6,3,'[N/m2]')
        sheet.write(pilerow+7,3,'[-]')
        sheet.write(pilerow+8,3,'[N/m2]')
        sheet.write(pilerow+9,3,'[tonnes]')

        #MAIN LATTICE
        jcktrow=pilerow+11
        sheet.write(jcktrow,0,'MAIN LATTICE',easyxf('font: name Arial, bold True'))
        sheet.write(jcktrow+1,1,'Batter')
        sheet.write(jcktrow+2,1,'No. of legs')
        sheet.write(jcktrow+3,1,'No. of bays')
        sheet.write(jcktrow+4,1,'Footprint')
        sheet.write(jcktrow+5,1,'Leg OD')
        sheet.write(jcktrow+6,1,'Leg t')
        sheet.write(jcktrow+7,1,'X-brace OD')
        sheet.write(jcktrow+8,1,'X-brace t')
        sheet.write(jcktrow+9,1,'Mud-brace OD')
        sheet.write(jcktrow+10,1,'Mud-brace t')
        sheet.write(jcktrow+11,1,'Xbrace-Leg Angle')
        #values
        sheet.write(jcktrow+1,2,myjckt.JcktGeoIn.batter)
        sheet.write(jcktrow+2,2,myjckt.JcktGeoIn.nlegs)
        sheet.write(jcktrow+3,2,myjckt.JcktGeoIn.nbays)
        sheet.write(jcktrow+4,2,myjckt.wbase)
        sheet.write(jcktrow+5,2,myjckt.Legs.legouts.LegObj.D[0])
        sheet.write(jcktrow+6,2,myjckt.Legs.legouts.LegObj.t[0])
        sheet.write(jcktrow+7,2,myjckt.Xbraces.Xbrcouts.LLURObj.D[0])
        sheet.write(jcktrow+8,2,myjckt.Xbraces.Xbrcouts.LLURObj.t[0])
        sheet.write(jcktrow+9,2,myjckt.Mudbraces.Mbrcouts.brcObj.D[0])
        sheet.write(jcktrow+10,2,myjckt.Mudbraces.Mbrcouts.brcObj.t[0])
        sheet.write(jcktrow+11,2,myjckt.PreBuild.beta3D*180./np.pi)
        #units
        sheet.write(jcktrow+1,3,'[-]')
        sheet.write(jcktrow+2,3,'[-]')
        sheet.write(jcktrow+3,3,'[-]')
        sheet.write(jcktrow+4,3,'[m]')
        sheet.write(jcktrow+5,3,'[m]')
        sheet.write(jcktrow+6,3,'[m]')
        sheet.write(jcktrow+7,3,'[m]')
        sheet.write(jcktrow+8,3,'[m]')
        sheet.write(jcktrow+9,3,'[m]')
        sheet.write(jcktrow+10,3,'[m]')
        sheet.write(jcktrow+11,3,'[deg]')

        for ii in range(myjckt.JcktGeoIn.nbays):
            sheet.write(jcktrow+12+ii,1,('Bay {:d} Length').format(ii))
            sheet.write(jcktrow+12+ii,2,(myjckt.Xbraces.bay_hs[ii]))
            sheet.write(jcktrow+12+ii,3,'[m]')

        rowno=jcktrow+12+myjckt.JcktGeoIn.nbays
        sheet.write(rowno,1,'Total Leg Height from mudline')
        sheet.write(rowno+1,1,'Leg bottom z')
        sheet.write(rowno+2,1,'Leg bottom stump')
        sheet.write(rowno+3,1,'Leg top stump')
        sheet.write(rowno+4,1,'Leg Material')
        sheet.write(rowno+5,1,'Leg Material rho')
        sheet.write(rowno+6,1,'Leg Material E')
        sheet.write(rowno+7,1,'Leg Material nu')
        sheet.write(rowno+8,1,'Leg Material fy')
        sheet.write(rowno+9,1,'Xbrace Material ')
        sheet.write(rowno+10,1,'Xbrace Material rho')
        sheet.write(rowno+11,1,'Xbrace Material E')
        sheet.write(rowno+12,1,'Xbrace Material nu')
        sheet.write(rowno+13,1,'Xbrace Material fy')
        sheet.write(rowno+14,1,'Mudbrace Material')
        sheet.write(rowno+15,1,'Mudbrace Material rho')
        sheet.write(rowno+16,1,'Mudbrace Material E')
        sheet.write(rowno+17,1,'Mudbrace Material nu')
        sheet.write(rowno+18,1,'Mudbrace Material fy')
        sheet.write(rowno+19,1,'Mass')
        #values
        sheet.write(rowno,2,myjckt.PreBuild.JcktH)
        sheet.write(rowno+1,2,myjckt.PreBuild.legZbot)
        sheet.write(rowno+2,2,myjckt.PreBuild.legbot_stmph)
        sheet.write(rowno+3,2,myjckt.TPinputs.hstump)
        sheet.write(rowno+4,2,myjckt.Legs.legouts.LegObj.mat[0].matname)
        sheet.write(rowno+5,2,myjckt.Legs.legouts.LegObj.mat[0].rho)
        sheet.write(rowno+6,2,myjckt.Legs.legouts.LegObj.mat[0].E)
        sheet.write(rowno+7,2,myjckt.Legs.legouts.LegObj.mat[0].nu)
        sheet.write(rowno+8,2,myjckt.Legs.legouts.LegObj.mat[0].fy)
        sheet.write(rowno+9,2, myjckt.Xbraces.Xbrcouts.LLURObj.mat[0].matname)
        sheet.write(rowno+10,2,myjckt.Xbraces.Xbrcouts.LLURObj.mat[0].rho)
        sheet.write(rowno+11,2,myjckt.Xbraces.Xbrcouts.LLURObj.mat[0].E)
        sheet.write(rowno+12,2,myjckt.Xbraces.Xbrcouts.LLURObj.mat[0].nu)
        sheet.write(rowno+13,2,myjckt.Xbraces.Xbrcouts.LLURObj.mat[0].fy)
        sheet.write(rowno+14,2, myjckt.Mudbraces.Mbrcouts.brcObj.mat[0].matname)
        sheet.write(rowno+15,2, myjckt.Mudbraces.Mbrcouts.brcObj.mat[0].rho)
        sheet.write(rowno+16,2, myjckt.Mudbraces.Mbrcouts.brcObj.mat[0].E)
        sheet.write(rowno+17,2, myjckt.Mudbraces.Mbrcouts.brcObj.mat[0].nu)
        sheet.write(rowno+18,2, myjckt.Mudbraces.Mbrcouts.brcObj.mat[0].fy)
        sheet.write(rowno+19,2,(myjckt.FrameOut.Frameouts_outs.mass[0]-myjckt.Tower.Twrouts.mass-myjckt.TP.TPouts.mass)/1.e3)  #mass of main lattice no TP
        #units
        sheet.write(rowno,  3,'[m]')
        sheet.write(rowno+1,3,'[m]')
        sheet.write(rowno+2,3,'[m]')
        sheet.write(rowno+3,3,'[m]')
        sheet.write(rowno+4,3,'[-]')
        sheet.write(rowno+5,3,'[kg/m3]')
        sheet.write(rowno+6,3,'[N/m2]')
        sheet.write(rowno+7,3,'[-]')
        sheet.write(rowno+8,3,'[N/m2]')
        sheet.write(rowno+9,3,'[-]')
        sheet.write(rowno+10,3,'[kg/m3]')
        sheet.write(rowno+11,3,'[N/m2]')
        sheet.write(rowno+12,3,'[-]')
        sheet.write(rowno+13,3,'[N/m2]')
        sheet.write(rowno+14,3,'[-]')
        sheet.write(rowno+15,3,'[kg/m3]')
        sheet.write(rowno+16,3,'[N/m2]')
        sheet.write(rowno+17,3,'[-]')
        sheet.write(rowno+18,3,'[N/m2]')
        sheet.write(rowno+19,3,'[tonnes]')

        #TRANSITION PIECE
        tprow=rowno+21
        #params
        sheet.write(tprow,  0,'TRANSITION PIECE',easyxf('font: name Arial, bold True'))
        sheet.write(tprow+1,1,'Deck Height')
        sheet.write(tprow+2,1,'Deck Width')
        sheet.write(tprow+3,1,'Stringer width')
        sheet.write(tprow+4,1,'Stringer height')
        sheet.write(tprow+5,1,'Strut OD')
        sheet.write(tprow+6,1,'Strut t')
        sheet.write(tprow+7,1,'Girder OD')
        sheet.write(tprow+8,1,'Girder t')
        sheet.write(tprow+9,1,'Diagonal brace OD')
        sheet.write(tprow+10,1,'Diagonal brace t')
        sheet.write(tprow+11,1,'Shell OD')
        sheet.write(tprow+12,1,'Shell t')
        sheet.write(tprow+13,1,'TP length')

        sheet.write(tprow+14,1,'Strut Material')
        sheet.write(tprow+15,1,'Strut Material rho')
        sheet.write(tprow+16,1,'Strut Material E')
        sheet.write(tprow+17,1,'Strut Material nu')
        sheet.write(tprow+18,1,'Strut Material fy')
        sheet.write(tprow+19,1,'Girder Material')
        sheet.write(tprow+20,1,'Girder Material rho')
        sheet.write(tprow+21,1,'Girder Material E')
        sheet.write(tprow+22,1,'Girder Material nu')
        sheet.write(tprow+23,1,'Girder Material fy')
        sheet.write(tprow+24,1,'Diag. Brace Material')
        sheet.write(tprow+25,1,'Diag. Brace Material rho')
        sheet.write(tprow+26,1,'Diag. Brace Material E')
        sheet.write(tprow+27,1,'Diag. Brace Material nu')
        sheet.write(tprow+28,1,'Diag. Brace Material fy')
        sheet.write(tprow+29,1,'Shell Material')
        sheet.write(tprow+30,1,'Shell Material rho')
        sheet.write(tprow+31,1,'Shell Material E')
        sheet.write(tprow+32,1,'Shell Material nu')
        sheet.write(tprow+33,1,'Shell Material fy')
        sheet.write(tprow+34,1,'Lumped Mass')
        sheet.write(tprow+35,1,'Tot Mass')
        #value
        sheet.write(tprow+1,2,myjckt.JcktGeoIn.dck_botz)
        sheet.write(tprow+2,2,myjckt.PreBuild.dck_width)
        sheet.write(tprow+3,2,0.05)  #arbitrary for now
        sheet.write(tprow+4,2,0.2)   #arbitrary for now
        sheet.write(tprow+5,2,myjckt.TP.TPouts.strtObj.D[0])
        sheet.write(tprow+6,2,myjckt.TP.TPouts.strtObj.t[0])
        sheet.write(tprow+7,2,myjckt.TP.TPouts.girObj.D[0])
        sheet.write(tprow+8,2,myjckt.TP.TPouts.girObj.t[0])
        sheet.write(tprow+9,2,myjckt.TP.TPouts.brcObj.D[0])
        sheet.write(tprow+10,2,myjckt.TP.TPouts.brcObj.t[0])
        sheet.write(tprow+11,2,myjckt.TP.TPouts.stemObj.D[0])
        sheet.write(tprow+12,2,myjckt.TP.TPouts.stemObj.t[0])
        sheet.write(tprow+13,2,myjckt.TPinputs.hstem.sum()+myjckt.TPinputs.hstump)  #replaced  the deck stringer height (.25) with stump height

        sheet.write(tprow+14,2,myjckt.TP.TPouts.strtObj.mat[0].matname)
        sheet.write(tprow+15,2,myjckt.TP.TPouts.strtObj.mat[0].rho)
        sheet.write(tprow+16,2,myjckt.TP.TPouts.strtObj.mat[0].E)
        sheet.write(tprow+17,2,myjckt.TP.TPouts.strtObj.mat[0].nu)
        sheet.write(tprow+18,2,myjckt.TP.TPouts.strtObj.mat[0].fy)
        sheet.write(tprow+19,2,myjckt.TP.TPouts.girObj.mat[0].matname)
        sheet.write(tprow+20,2,myjckt.TP.TPouts.girObj.mat[0].rho)
        sheet.write(tprow+21,2,myjckt.TP.TPouts.girObj.mat[0].E)
        sheet.write(tprow+22,2,myjckt.TP.TPouts.girObj.mat[0].nu)
        sheet.write(tprow+23,2,myjckt.TP.TPouts.girObj.mat[0].fy)
        sheet.write(tprow+24,2,myjckt.TP.TPouts.brcObj.mat[0].matname)
        sheet.write(tprow+25,2,myjckt.TP.TPouts.brcObj.mat[0].rho)
        sheet.write(tprow+26,2,myjckt.TP.TPouts.brcObj.mat[0].E)
        sheet.write(tprow+27,2,myjckt.TP.TPouts.brcObj.mat[0].nu)
        sheet.write(tprow+28,2,myjckt.TP.TPouts.brcObj.mat[0].fy)
        sheet.write(tprow+29,2,myjckt.TP.TPouts.stemObj.mat[0].matname)
        sheet.write(tprow+30,2,myjckt.TP.TPouts.stemObj.mat[0].rho)
        sheet.write(tprow+31,2,myjckt.TP.TPouts.stemObj.mat[0].E)
        sheet.write(tprow+32,2,myjckt.TP.TPouts.stemObj.mat[0].nu)
        sheet.write(tprow+33,2,myjckt.TP.TPouts.stemObj.mat[0].fy)
        sheet.write(tprow+34,2,myjckt.TP.TPouts.TPlumpedMass[0]/1.e3)
        sheet.write(tprow+35,2,(myjckt.TP.TPouts.mass+myjckt.TP.TPouts.TPlumpedMass[0])/1.e3)

        #units
        sheet.write(tprow+1,3,'[m]')
        sheet.write(tprow+2,3,'[m]')
        sheet.write(tprow+3,3,'[m]')
        sheet.write(tprow+4,3,'[m]')
        sheet.write(tprow+5,3,'[m]')
        sheet.write(tprow+6,3,'[m]')
        sheet.write(tprow+7,3,'[m]')
        sheet.write(tprow+8,3,'[m]')
        sheet.write(tprow+9,3,'[m]')
        sheet.write(tprow+10,3,'[m]')
        sheet.write(tprow+11,3,'[m]')
        sheet.write(tprow+12,3,'[m]')
        sheet.write(tprow+13,3,'[m]')
        sheet.write(tprow+14,3,'[-]')
        sheet.write(tprow+15,3,'[kg/m3]')
        sheet.write(tprow+16,3,'[N/m2]')
        sheet.write(tprow+17,3,'[-]')
        sheet.write(tprow+18,3,'[N/m2]')
        sheet.write(tprow+19,3,'[-]')
        sheet.write(tprow+20,3,'[kg/m3]')
        sheet.write(tprow+21,3,'[N/m2]')
        sheet.write(tprow+22,3,'[-]')
        sheet.write(tprow+23,3,'[N/m2]')
        sheet.write(tprow+24,3,'[-]')
        sheet.write(tprow+25,3,'[kg/m3]')
        sheet.write(tprow+26,3,'[N/m2]')
        sheet.write(tprow+27,3,'[-]')
        sheet.write(tprow+28,3,'[N/m2]')
        sheet.write(tprow+29,3,'[-]')
        sheet.write(tprow+30,3,'[kg/m3]')
        sheet.write(tprow+31,3,'[N/m2]')
        sheet.write(tprow+32,3,'[-]')
        sheet.write(tprow+33,3,'[N/m2]')

        sheet.write(tprow+34,3,'[tonnes]')
        sheet.write(tprow+35,3,'[tonnes]')

        #TOWER
        twrrow=tprow+37
        #params
        sheet.write(twrrow,0,'TOWER',easyxf('font: name Arial, bold True'))
        sheet.write(twrrow+1,1,'Base OD')
        sheet.write(twrrow+2,1,'Base t')
        sheet.write(twrrow+3,1,'Base DTR')
        sheet.write(twrrow+4,1,'Top OD')
        sheet.write(twrrow+5,1,'Top t')
        sheet.write(twrrow+6,1,'Top DTR')
        sheet.write(twrrow+7,1,'Tower Length')
        sheet.write(twrrow+8,1,'Tower H2 Length')
        sheet.write(twrrow+9,1,'Tower Material')
        sheet.write(twrrow+10,1,'Tower Material rho')
        sheet.write(twrrow+11,1,'Tower Material E')
        sheet.write(twrrow+12,1,'Tower Material nu')
        sheet.write(twrrow+13,1,'Tower Material fy')

        sheet.write(twrrow+14,1,'mass')
        #values
        sheet.write(twrrow+1,2,myjckt.Tower.Twrouts.TwrObj.D[0])
        sheet.write(twrrow+2,2,myjckt.Tower.Twrouts.TwrObj.t[0])
        sheet.write(twrrow+3,2,myjckt.Tower.Twrins.DTRb)
        sheet.write(twrrow+4,2,myjckt.Tower.Twrins.Dt)
        sheet.write(twrrow+5,2,myjckt.Tower.Twrins.Dt/myjckt.Tower.Twrins.DTRt)
        #sheet.write(twrrow+4,2,myjckt.Tower.Twrouts.TwrObj.D[-1]) This does not work since it does not give me the exact final node, but its constantxsec element tube
        #sheet.write(twrrow+5,2,myjckt.Tower.Twrouts.TwrObj.t[-1])
        sheet.write(twrrow+6,2,myjckt.Tower.Twrins.DTRt)
        sheet.write(twrrow+7,2,myjckt.Tower.Twrouts.Htwr)
        sheet.write(twrrow+8,2,myjckt.Tower.Twrouts.Htwr2)
        sheet.write(twrrow+9, 2,myjckt.Tower.Twrouts.TwrObj.mat[0].matname)
        sheet.write(twrrow+10,2,myjckt.Tower.Twrouts.TwrObj.mat[0].rho)
        sheet.write(twrrow+11,2,myjckt.Tower.Twrouts.TwrObj.mat[0].E)
        sheet.write(twrrow+12,2,myjckt.Tower.Twrouts.TwrObj.mat[0].nu)
        sheet.write(twrrow+13,2,myjckt.Tower.Twrouts.TwrObj.mat[0].fy)

        sheet.write(twrrow+14,2,myjckt.Tower.Twrouts.mass/1.e3)
        #units
        sheet.write(twrrow+1,3,'[m]')
        sheet.write(twrrow+2,3,'[m]')
        sheet.write(twrrow+3,3,'[-]')
        sheet.write(twrrow+4,3,'[m]')
        sheet.write(twrrow+5,3,'[m]')
        sheet.write(twrrow+6,3,'[-]')
        sheet.write(twrrow+7,3,'[m]')
        sheet.write(twrrow+8,3,'[m]')
        sheet.write(twrrow+9, 3,'[-]')
        sheet.write(twrrow+10,3,'[kg/m3]')
        sheet.write(twrrow+11,3,'[N/m2]')
        sheet.write(twrrow+12,3,'[-]')
        sheet.write(twrrow+13,3,'[N/m2]')
        sheet.write(twrrow+14,3,'[tonnes]')

        #TOTAL MASS
        Totrow=twrrow+16
        sheet.write(Totrow,1,'TOTAL MASS',easyxf('font: name Arial, bold True'))
        sheet.write(Totrow,2,(myjckt.TP.TPouts.TPlumpedMass[0]+myjckt.FrameOut.Frameouts_outs.mass[0]+myjckt.Mpiles)/1.e3,\
                    easyxf('font: name Arial, bold True'))
        sheet.write(Totrow,3,'[tonnes]')

        #Calculated Frequencies
        Freqrow=Totrow+5
        sheet.write(Freqrow,0,'Calculated Frequencies ',easyxf('font: name Arial, bold True'))
        sheet.write(Freqrow+1,1,'1st EigenFrequency ')
        sheet.write(Freqrow+2,1,'2nd EigenFrequency ')
        #values
        sheet.write(Freqrow+1,2,myjckt.FrameOut.Frameouts_outs.Freqs[0], easyxf('font: name Arial, bold True'))
        sheet.write(Freqrow+2,2,myjckt.FrameOut.Frameouts_outs.Freqs[1], easyxf('font: name Arial, bold True'))
        #units
        sheet.write(Freqrow+1,3,'[Hz]')
        sheet.write(Freqrow+2,3,'[Hz]')

        book.save(xlsfilename)

#__________________________________________________________#
def ReconAllLinesJckt(casefile,caserange,titlines=3,hdrlines=1,untlines=1,opttitlines=3,opthdrlines=1,optuntlines=1):
    """This function saves figures from the results of all cases.\n
        Note: the optimization output filename is assumed to be "caseno_casename_out.dat". \n
    INPUT \n
        casefile    -string, complete path+filename to table file of cases. \n
        caserange   -int array(2), first and last case number to read in (also rows out of the table to be read). \n
    OPTIONALS:
        titlines - int, number of lines for the title in the table file \n
        hdrlines - int, number of lines for the variable headers in the table file\n
        untlines - int, number of lines for the units in the table file\n
        opttitlines - int, number of lines for the title in the optimal var file \n
        opthdrlines - int, number of lines for the variable headers in the optimal var file\n
        optuntlines - int, number of lines for the units in the optimal var file\n
    OUTPUT \n
        It creates two plots: of the configuration + utilization of tower. \n

        """
    desvars=DesVars() #instance of design variables

    for caseno in (range(caserange[0],caserange[1])-1):

    #First read design parameters from table file
        Desprms,_,desvarbds,desvarmeans,guesses,casename=ReadTab1Line(casefile,caseno,desvars.dvars.keys(),\
                                                              titlines=titlines,hdrlines=hdrlines,untlines=untlines)

        #Then read optimized design variables from file
        optfile=ntpath.join(ntpath.dirname(casefile),+casename+'_out.dat')
        fileID=open(optfile,'r')
        #skip title, header, unit lines
        for ii in range(0,opthdrlines+opthdrlines+optuntlines+1):
            _=fileID.readline()

        #read the actual data
        dat=fileID.readline().split()

        for ii,key in enumerate(desvars.dvars):
            #print key  #debug
            setattr(desvars,key,float(dat[ii]))

        fileID.close()

        #Then build the initial assembly and run it to instantiate everything
        myjckt=SetJacketInputsPeregrine.main(Desprms,desvars)
        myjckt.run()
        #Plot
        PlotJacket(myjckt,util=True,savefileroot=os.path.join(ntpath.dirname(casefile),casename))

#__________________________________________________________#

def SmoothCheck(casefile,caseno,varnames,titlines=3,hdrlines=1,untlines=1,opttitlines=3,opthdrlines=1,optuntlines=1):
    """This function reads design parameters for caseno case from casefile and then returns   results into a file, in terms of case number and case name, and design variables.\n
        Note: the optimization output filename is assumed to be "caseno_casename_out.dat". \n
    INPUT \n
        casefile -string, complete path+filename to table file of cases. \n
        caseno   -int, case number (also row out of the table to be read). \n
        varnames -Ordered Dictionary, with key:value, where key='batter','Dleg', etc., as from DesVars and value is the corresponding range array for each key
    OPTIONALS:
        titlines - int, number of lines for the title in the table file \n
        hdrlines - int, number of lines for the variable headers in the table file\n
        untlines - int, number of lines for the units in the table file\n
        opttitlines - int, number of lines for the title in the optimal var file \n
        opthdrlines - int, number of lines for the variable headers in the optimal var file\n
        optuntlines - int, number of lines for the units in the optimal var file\n
    OUTPUT \n
        It creates a plot of the configuration, utilization of tower. \n
        myjckt  -OpenMdao object assembly of JacketSE with the design parameters and optimzation variables found.\n
        """
    from RunAllJacketTwrs import plotting
    desvars=DesVars() #instance of design variables

    #First read design parameters from table file
    Desprms,_,desvarbds,desvarmeans,guesses,casename=ReadTab1Line(casefile,caseno,desvars.dvars.keys(),\
                                                          titlines=titlines,hdrlines=hdrlines,untlines=untlines)

    #Then read optimized design variables from file
    optfile=ntpath.join(ntpath.dirname(casefile),casename+'_out.dat')
    fileID=open(optfile,'r')
    #skip title, header, unit lines
    for ii in range(0,opthdrlines+opthdrlines+optuntlines+1):
        _=fileID.readline()

    #read the actual data
    dat=fileID.readline().split()

    for ii,key in enumerate(desvars.dvars):
        #print key  #debug
        setattr(desvars,key,float(dat[ii]))

    fileID.close()

    #Then build the initial assembly and run it to instantiate everything
    myjckt=SetJacketInputsPeregrine.main(Desprms,desvars)
    myjckt.run()
    #Plot
    PlotJacket(myjckt,util=True)#,savefileroot=os.path.join(ntpath.dirname(casefile),casename))


    #NOW SPAN THE Variable passed
    nvars=varnames.__len__()

    print "Now spanning variables\n"
    cex=0 #execution counter
    keyno=0   #counter for variables to span
    for key,value in varnames.iteritems():
        max_GLUtil=np.zeros([value.size,myjckt.tower_utilization.GLUtil.size,nvars])
        max_EUUtil=max_GLUtil.copy()
        Lp0rat=np.zeros([value.size,nvars])
        max_cbutil=np.zeros([value.size,myjckt.jacket_utilization.cb_util.size,nvars])
        max_tutil=np.zeros([value.size,myjckt.jacket_utilization.t_util.size,nvars])
        max_KjntUtil=np.zeros([value.size,myjckt.jacket_utilization.KjntUtil.size,nvars])
        max_XjntUtil=np.zeros([value.size,myjckt.jacket_utilization.XjntUtil.size,nvars])

        mass_str=Lp0rat.copy()
        mass_tot=Lp0rat.copy()
        f1=Lp0rat.copy()
        f2=Lp0rat.copy()
        wbase=Lp0rat.copy()
        Mpiles=Lp0rat.copy()

        for ib in range(value.size):
            if key=='batter' or key=='dck_width':
                setattr(myjckt.JcktGeoIn,key,value[ib])
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

        plotting(key,value,ntpath.dirname(casefile),(mass_str[:,keyno]+Mpiles[:,keyno]),f1[:,keyno],max_cbutil[:,:,keyno],max_tutil[:,:,keyno],max_KjntUtil[:,:,keyno],max_XjntUtil[:,:,keyno],Lp0rat[:,keyno])
        keyno +=1


    return myjckt


#__________________________________________________________#
def main():



 #       casefile=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\COBYLA\testmatrix.dat'

        casefile=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\COBYLA\testmatrix9mDb.dat'
        #caseno=82
        caseno=75
#        xlsfilename=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\COBYLA\outputs.xls'
#        xlsfilename=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\COBYLA\outputs_coby.xls'
        xlsfilename=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\COBYLA\outputs9mDb_coby.xls'
        optfile=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\COBYLA\10MW_LooseBatter\9mDbMax_othercases\75_R2D2H0M2T0S0_out.dat'
        myjckt=Recon1LineJckt(casefile,caseno,xlsfilename=xlsfilename,optfile=optfile)
        print "reconstruction terminated"

        #Smooth check
        casefile=r'C:\PROJECTS\OFFSHORE_WIND\LCOE_COSTANALYSIS\COBYLA\testmatrix.dat'
        casefile=r'D:\RRD_ENGINEERING\PROJECTS\NREL\OFFSHOREWIND\LCOE_ANALYSIS\COBYLA\testmatrix.dat'
        caseno=46
        varnames=collections.OrderedDict([ ('dck_width',np.arange(11.,15.,0.25)),('batter',np.arange(8.,9.2,0.2))  ])   # ('Db',np.arange(3.0,4.5,0.05)) ('Db',np.arange(3.75,4.,0.05))('Dpile',np.arange(0.8,1.5,0.05)), ('Lp',np.arange(48,50,0.05))('batter',np.arange(7.,10.,0.05))
        SmoothCheck(casefile,caseno,varnames)


if __name__ == '__main__':
    main()
