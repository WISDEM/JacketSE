#-------------------------------------------------------------------------------
# Name: pySubDyn
# Purpose: Execute SubDyn within JacketSE
#
#Example inputs correspond to Test03 under CertTest, an OC4 'Jacket' SIMPLIFIED, 1-bay
#
# Author:      Casey Broslawski
#              SULI intern under the guidance of Dr. Rick Damiani
#
# Created:     07/07/2016
#-------------------------------------------------------------------------------

import numpy as np
import subprocess, os
from openmdao.main.api import Component, VariableTree, Assembly
from openmdao.main.datatypes.api import Float, Int, Bool, Array, VarTree, Str
from VarTrees import RNAprops

print "pySubDyn initialize"

class SubDynAOutputs(VariableTree): #test
    FEMeigenvalues=Array(dtype=np.float, untis='Hz', desc='FEM Eigenvalues [Hz]')
    StructureMass=Float(dtype=np.float, units='kg', desc="SubDyn's Total Mass (structural and non-structural)")
    StructureCM=Array(dtype=np.float, units='m', desc="SubDyn's Total Mass CM coordinates (Xcm,Ycm,Zcm)")

class WriteInputFile(Component):

    gravity =Float(units='m/s**2',iotype='in',desc='Gravity Acceleration (ABSOLUTE VALUE!)')
    nodes=Array(iotype='in', desc='Nodes from JacketSE')
    mems=Array(iotype='in',desc='Member connections between joints')
    Reacts=Array(iotype='in', desc='Base joint IDs and fixity')
    SDPropSet=Array(iotype='in', desc='Properties to be used in SubDyn')
    RNAinputs=VarTree(RNAprops(), iotype='in', desc='Basic Inertial Properties of RNA') #test

    TEST=Bool(iotype='in', desc='Check to see if pySubDyn is being run locally') #test flag

    #INPUT FILE INPUTS----------------------------------------------------------

    InputFile_name=Str(iotype='in')
    SDInputFile=Array(iotype='in')

    #Simulation Control
    Echo=Array(iotype='in')
    SDdeltaT=Array(iotype='in')
    IntMethod=Array(iotype='in')
    SttcSolve=Array(iotype='in')

    #FEA and CRAIG-BAMPTON PARAMETERS
    FEMmod=Array(iotype='in')
    NDiv=Array(iotype='in')
    CBMod=Array(iotype='in')
    Nmodes=Array(iotype='in')
    JDampings=Array(iotype='in')

    #Structure Joints
    SDjointsHeader=Array(iotype='in')
    SDjoints=Array(iotype='in', desc='Leg joints from jacket.py. In SD joints are used to build the structure, not nodes.')

    #Base Reaction Joints
    BaseRxnJointsHeader=Array(iotype='in')
    BaseRxnJoints=Array(iotype='in', desc='Base joints of the structure.')

    #Interface Joints
    InterfaceRxnJointsHeader=Array(iotype='in')
    InterfaceJoints=Array(iotype='in', desc='Interface joints.')

    #Members
    MembersHeader=Array(iotype='in')
    Members=Array(iotype='in', desc='Build members from joints, not nodes.')

    #MEMBER X-SECTION PROPERTY data 1/2
    NPropSets=Array(iotype='in')
    PropSet1Header=Array(iotype='in')
    PropSet1=Array(iotype='in')

    #MEMBER X-SECTION PROPERTY data 2/2
    PropSet2Header=Array(iotype='in')
    PropSet2=Array(iotype='in')

    #MEMBER COSINE MATRICES COSM(i,j)
    COSMHeader=Array(iotype='in')
    COSMs=Array(iotype='in')

    #JOINT ADDITIONAL CONCENTRATED MASSES
    CmassHeader=Array(iotype='in')
    Cmass=Array(iotype='in')

    #OUTPUT: SUMMARY & OUTFILE
    SSSum=Array(iotype='in')
    OutCOSM=Array(iotype='in')
    OutAll=Array(iotype='in')
    OutSwtch=Array(iotype='in')
    TabDelim=Array(iotype='in')
    OutDec=Array(iotype='in')
    OutFmt=Array(iotype='in')
    OutSFmt=Array(iotype='in')

    #MEMBER OUTPUT LIST
    MemOutListHeader=Array(iotype='in')
    MemOutList=Array(iotype='in', desc='Members whose loads and dynamics will be output.')

    #SSOutline
    SSOutlist=Array(iotype='in')

    print "Write Input File pre"

    def execute(self):
        print "Write input file"

        self.NDiv_set=1 #SET NDIV AS 1 TO USE JACKETSE'S NODES AS SUBDYN'S JOINTS

        row=0
        Holder=[]
        Interface=[1., 2.]
        for j in self.nodes:
            row+=1
            for i in j:
                if i==np.amax(self.nodes):
                    Holder=np.hstack((row, i))
                    Interface=np.vstack((Interface, Holder))
        Interface=np.delete(Interface,(0),axis=0)

        RNA_CMs=[]
        RNA_CMs=np.hstack((len(self.nodes), self.RNAinputs.mass, self.RNAinputs.I[0], self.RNAinputs.I[1], self.RNAinputs.I[2]))

        #Write SubDyn input file
        file=open(str(self.SDInputFile[0]),"w")

        file.write(('{}\n').format("----------- SubDyn v1.01.x MultiMember Support Structure Input File ------------"))

        file.write(('{}\n').format("Customized jacket input file for SubDyn")) #TODO Set as a variable so users can describe their own jackets

        file.write(('{}\n').format("-------------------------- SIMULATION CONTROL  ---------------------------------"))
        file.write(('{}\n').format(self.Echo[0]+' '+self.Echo[1]+' '+self.Echo[2]))
        file.write(('{}\n').format('"'+self.SDdeltaT[0]+'" '+self.SDdeltaT[1]+' '+self.SDdeltaT[2]))
        file.write(('{}\n').format(self.IntMethod[0]+' '+self.IntMethod[1]+' '+self.IntMethod[2]))
        file.write(('{}\n').format(self.SttcSolve[0]+' '+self.SttcSolve[1]+' '+self.SttcSolve[2]))

        file.write(('{}\n').format("-------------------- FEA and CRAIG-BAMPTON PARAMETERS---------------------------"))
        file.write(('{}\n').format(self.FEMmod[0]+' '+self.FEMmod[1]))
        if not self.TEST:
            file.write(('{}\n').format(str(self.NDiv_set)+' '+self.NDiv[1]+' '+self.NDiv[2]))
        else:
            file.write(('{}\n').format(self.NDiv[0]+' '+self.NDiv[1]+' '+self.NDiv[2]))
        file.write(('{}\n').format(self.CBMod[0]+' '+self.CBMod[1]))
        file.write(('{}\n').format(self.Nmodes[0]+' '+self.Nmodes[1]))
        file.write(('{}\n').format(self.JDampings[0]+' '+self.JDampings[1]+' '+self.JDampings[2]))

        file.write(('{}\n').format("--- STRUCTURE JOINTS: joints connect structure members (~Hydrodyn Input File)---"))
        if not self.TEST:
            self.Njoints=len(self.nodes)
        else:
            self.Njoints=len(self.SDjoints)
        self.NjointsHeader=np.array([self.Njoints, 'Njoints', '- Number of joints (-)'])
        file.write(('{}\n').format(self.NjointsHeader[0]+' '+self.NjointsHeader[1]+' '+self.NjointsHeader[2]))
        file.write(('{}\n').format(self.SDjointsHeader[0][0]+' '+self.SDjointsHeader[0][1]+' '+self.SDjointsHeader[0][2]+' '+self.SDjointsHeader[0][3]+' '+self.SDjointsHeader[0][4]))
        file.write(('{}\n').format(self.SDjointsHeader[1][0]+' '+self.SDjointsHeader[1][1]+' '+self.SDjointsHeader[1][2]+' '+self.SDjointsHeader[1][3]))
        if not self.TEST:
            for i in range(len(self.nodes)):
                self.Nnodes=(range(len(self.nodes)+1))
                file.write(('{}\n').format(str(self.Nnodes[i+1])+' '+str(self.nodes[i][0])+' '+str(self.nodes[i][1])+' '+str(self.nodes[i][2])))
        else:
            for i in range(self.Njoints):
                self.SDjointsInt=int(self.SDjoints[i][0])
                file.write(('{}\n').format(str(self.SDjointsInt)+' '+str(self.SDjoints[i][1])+' '+str(self.SDjoints[i][2])+' '+str(self.SDjoints[i][3])))

        file.write(('{}\n').format("------------------- BASE REACTION JOINTS: 1/0 for Locked/Free DOF @ each Reaction Node ---------------------"))
        if not self.TEST:
            self.NReact=len(self.Reacts)
        else:
            self.NReact=len(self.BaseRxnJoints)
        self.NReactHeader=np.array([self.NReact, 'NReact', '- Number of Joints with reaction forces; be sure to remove all rigid motion DOFs of the structure  (else det([K])=[0])'])
        file.write(('{}\n').format(self.NReactHeader[0]+' '+self.NReactHeader[1]+' '+self.NReactHeader[2]))
        file.write(('{}\n').format(self.BaseRxnJointsHeader[0][0]+' '+self.BaseRxnJointsHeader[0][1]+' '+self.BaseRxnJointsHeader[0][2]+' '+self.BaseRxnJointsHeader[0][3]+' '+self.BaseRxnJointsHeader[0][4]+' '+self.BaseRxnJointsHeader[0][5]+' '+self.BaseRxnJointsHeader[0][6]+' '+self.BaseRxnJointsHeader[0][7]))
        file.write(('{}\n').format(self.BaseRxnJointsHeader[1][0]+' '+self.BaseRxnJointsHeader[1][1]+' '+self.BaseRxnJointsHeader[1][2]+' '+self.BaseRxnJointsHeader[1][3]+' '+self.BaseRxnJointsHeader[1][4]+' '+self.BaseRxnJointsHeader[1][5]+' '+self.BaseRxnJointsHeader[1][6]))
        for i in range(self.NReact):
            if not self.TEST:
                file.write(('{}\n').format(str(self.Reacts[i][0])+' '+str(self.Reacts[i][1])+' '+str(self.Reacts[i][2])+' '+str(self.Reacts[i][3])+' '+str(self.Reacts[i][4])+' '+str(self.Reacts[i][5])+' '+str(self.Reacts[i][6])))
            else:
                file.write(('{}\n').format(str(self.BaseRxnJoints[i][0])+' '+str(self.BaseRxnJoints[i][1])+' '+str(self.BaseRxnJoints[i][2])+' '+str(self.BaseRxnJoints[i][3])+' '+str(self.BaseRxnJoints[i][4])+' '+str(self.BaseRxnJoints[i][5])+' '+str(self.BaseRxnJoints[i][6])))

        file.write(('{}\n').format("------- INTERFACE JOINTS: 1/0 for Locked (to the TP)/Free DOF @each Interface Joint (only Locked-to-TP implemented thus far (=rigid TP)) ---------"))
        if not self.TEST:
            self.NInterf=len(Interface)
        else:
            self.NInterf=len(self.InterfaceJoints)
        self.NInterfHeader=np.array([self.NInterf, 'NInterf', '- Number of interface joints locked to the Transition Piece (TP):  be sure to remove all rigid motion dofs'])
        file.write(('{}\n').format(self.NInterfHeader[0]+' '+self.NInterfHeader[1]+' '+self.NInterfHeader[2]))
        file.write(('{}\n').format(self.InterfaceRxnJointsHeader[0][0]+' '+self.InterfaceRxnJointsHeader[0][1]+' '+self.InterfaceRxnJointsHeader[0][2]+' '+self.InterfaceRxnJointsHeader[0][3]+' '+self.InterfaceRxnJointsHeader[0][4]+' '+self.InterfaceRxnJointsHeader[0][5]+' '+self.InterfaceRxnJointsHeader[0][6]+' '+self.InterfaceRxnJointsHeader[0][7]))
        file.write(('{}\n').format(self.InterfaceRxnJointsHeader[1][0]+' '+self.InterfaceRxnJointsHeader[1][1]+' '+self.InterfaceRxnJointsHeader[1][2]+' '+self.InterfaceRxnJointsHeader[1][3]+' '+self.InterfaceRxnJointsHeader[1][4]+' '+self.InterfaceRxnJointsHeader[1][5]+' '+self.InterfaceRxnJointsHeader[1][6]))
        for i in range(self.NInterf):
            if not self.TEST:
                file.write(('{}\n').format(str(int(Interface[i][0]))+' '+str(self.InterfaceJoints[i][1])+' '+str(self.InterfaceJoints[i][2])+' '+str(self.InterfaceJoints[i][3])+' '+str(self.InterfaceJoints[i][4])+' '+str(self.InterfaceJoints[i][5])+' '+str(self.InterfaceJoints[i][6])))
            else:
                file.write(('{}\n').format(str(self.InterfaceJoints[i][0])+' '+str(self.InterfaceJoints[i][1])+' '+str(self.InterfaceJoints[i][2])+' '+str(self.InterfaceJoints[i][3])+' '+str(self.InterfaceJoints[i][4])+' '+str(self.InterfaceJoints[i][5])+' '+str(self.InterfaceJoints[i][6])))

        file.write(('{}\n').format("----------------------------------- MEMBERS --------------------------------------"))
        if not self.TEST:
            self.NMembers=len(self.mems)
        else:
            self.NMembers=len(self.Members)
        self.NMembersHeader=np.array([self.NMembers, 'NMembers', '- Number of frame members'])
        file.write(('{}\n').format(self.NMembersHeader[0]+' '+self.NMembersHeader[1]+' '+self.NMembersHeader[2]))
        file.write(('{}\n').format(self.MembersHeader[0][0]+' '+self.MembersHeader[0][1]+' '+self.MembersHeader[0][2]+' '+self.MembersHeader[0][3]+' '+self.MembersHeader[0][4]+' '+self.MembersHeader[0][5]))
        file.write(('{}\n').format(self.MembersHeader[1][0]+' '+self.MembersHeader[1][1]+' '+self.MembersHeader[1][2]+' '+self.MembersHeader[1][3]+' '+self.MembersHeader[1][4]+' '+self.MembersHeader[1][5]))
        if not self.TEST:
            if len(self.COSMs)==0: #COSMID column is empty
                for i in range(self.NMembers):
                    file.write(('{}\n').format(str(i+1)+' '+str(self.mems[i][0])+' '+str(self.mems[i][1])+' '+str(int(self.SDPropSet[i][0]))+' '+str(int(self.SDPropSet[i][0]))))
        else:
            if len(self.COSMs)==0: #COSMID column is empty
                for i in range(self.NMembers):
                    file.write(('{}\n').format(str(self.Members[i][0])+' '+str(self.Members[i][1])+' '+str(self.Members[i][2])+' '+str(self.Members[i][3])+' '+str(self.Members[i][4])))
            else:   #COSMID column is not empty
                for i in range(self.NMembers):
                    self.file.write(('{}\n').format(str(self.Members[i][0])+' '+str(self.Members[i][1])+' '+str(self.Members[i][2])+' '+str(self.Members[i][3])+' '+str(self.Members[i][4]+' '+str(self.Members[i][5]))))

        file.write(('{}\n').format("------------------ MEMBER X-SECTION PROPERTY data 1/2 [isotropic material for now: use this table for circular-tubular elements] ------------------------"))
        if not self.TEST:
            self.NPropSets1=len(self.SDPropSet)
        else:
            self.NPropSets1=len(self.PropSet1)
        self.NPropSets1Header=np.array([self.NPropSets1, 'NPropSets', '- Number of structurally unique x-sections (i.e. how many groups of X-sectional properties are utilized throughout all of the members)'])
        file.write(('{}\n').format(self.NPropSets1Header[0]+' '+self.NPropSets1Header[1]+' '+self.NPropSets1Header[2]))
        file.write(('{}\n').format(self.PropSet1Header[0][0]+' '+self.PropSet1Header[0][1]+' '+self.PropSet1Header[0][2]+' '+self.PropSet1Header[0][3]+' '+self.PropSet1Header[0][4]+' '+self.PropSet1Header[0][5]))
        file.write(('{}\n').format(self.PropSet1Header[1][0]+' '+self.PropSet1Header[1][1]+' '+self.PropSet1Header[1][2]+' '+self.PropSet1Header[1][3]+' '+self.PropSet1Header[1][4]+' '+self.PropSet1Header[1][5]))
        if not self.TEST:
            for i in range(self.NPropSets1):
                if i<(self.NPropSets1-1):
                    file.write(('{}\n').format(str(int(self.SDPropSet[i][0]))+' '+str(self.SDPropSet[i][1])+' '+str(self.SDPropSet[i][2])+' '+str(self.SDPropSet[i][3])+' '+str(self.SDPropSet[i][4])+' '+str(self.SDPropSet[i][5])))
                else:
                    file.write(('{}\n').format(str(int(self.SDPropSet[i][0]))+' '+str(int(self.SDPropSet[i][1]))+' '+str(int(self.SDPropSet[i][2]))+' '+str(self.SDPropSet[i][3])+' '+str(self.SDPropSet[i][4])+' '+str(self.SDPropSet[i][5])))
        else:
            for i in range(self.NPropSets1):
                self.PropSet1Int=int(self.PropSet1[i][0])
                file.write(('{}\n').format(str(self.PropSet1Int)+' '+str(self.PropSet1[i][1])+' '+str(self.PropSet1[i][2])+' '+str(self.PropSet1[i][3])+' '+str(self.PropSet1[i][4])+' '+str(self.PropSet1[i][5])))

        file.write(('{}\n').format("------------------ MEMBER X-SECTION PROPERTY data 2/2 [isotropic material for now: use this table if any section other than circular, however provide COSM(i,j) below] ------------------"))
        self.NXPropSets2=len(self.PropSet2)
        self.NXPropSets2Header=np.array([self.NXPropSets2, 'NXPropSets', '- Number of structurally unique non-circular x-sections (if 0 the following table is ignored)'])
        file.write(('{}\n').format(self.NXPropSets2Header[0]+' '+self.NXPropSets2Header[1]+' '+self.NXPropSets2Header[2]))
        file.write(('{}\n').format(self.PropSet2Header[0][0]+' '+self.PropSet2Header[0][1]+' '+self.PropSet2Header[0][2]+' '+self.PropSet2Header[0][3]+' '+self.PropSet2Header[0][4]+' '+self.PropSet2Header[0][5]+' '+self.PropSet2Header[0][6]+' '+self.PropSet2Header[0][7]+' '+self.PropSet2Header[0][8]+' '+self.PropSet2Header[0][9]))
        file.write(('{}\n').format(self.PropSet2Header[1][0]+' '+self.PropSet2Header[1][1]+' '+self.PropSet2Header[1][2]+' '+self.PropSet2Header[1][3]+' '+self.PropSet2Header[1][4]+' '+self.PropSet2Header[1][5]+' '+self.PropSet2Header[1][6]+' '+self.PropSet2Header[1][7]+' '+self.PropSet2Header[1][8]+' '+self.PropSet2Header[1][9]))
        if self.NXPropSets2 != 0:
            for i in range(NXPropSets2):
                self.PropSet2Int=int(self.PropSet2[i][0])
                file.write(('{}\n').format(str(self.PropSet2Int)+' '+str(self.PropSet2[i][1])+' '+str(self.PropSet2[i][2])+' '+str(self.PropSet2[i][3])+' '+str(self.PropSet2[i][4])+' '+str(self.PropSet2[i][5])+' '+str(self.PropSet2[i][6])+' '+str(self.PropSet2[i][7])+' '+str(self.PropSet2[i][8])+' '+str(self.PropSet2[i][9])))

        file.write(('{}\n').format("---------------------- MEMBER COSINE MATRICES COSM(i,j) -------------------------"))
        self.NCOSMs=len(self.COSMs)
        self.NCOSMsHeader=np.array([self.NCOSMs, 'NCOSMs', '- Number of unique cosine matrices (i.e., of unique member alignments including principal axis rotations); ignored if NXPropSets=0   or 9999 in any element below'])
        file.write(('{}\n').format(self.NCOSMsHeader[0]+' '+self.NCOSMsHeader[1]+' '+self.NCOSMsHeader[2]))
        file.write(('{}\n').format(self.COSMHeader[0][0]+' '+self.COSMHeader[0][1]+' '+self.COSMHeader[0][2]+' '+self.COSMHeader[0][3]+' '+self.COSMHeader[0][4]+' '+self.COSMHeader[0][5]+' '+self.COSMHeader[0][6]+' '+self.COSMHeader[0][7]+' '+self.COSMHeader[0][8]+' '+self.COSMHeader[0][9]))
        file.write(('{}\n').format(self.COSMHeader[1][0]+' '+self.COSMHeader[1][1]+' '+self.COSMHeader[1][2]+' '+self.COSMHeader[1][3]+' '+self.COSMHeader[1][4]+' '+self.COSMHeader[1][5]+' '+self.COSMHeader[1][6]+' '+self.COSMHeader[1][7]+' '+self.COSMHeader[1][8]+' '+self.COSMHeader[1][9]))
        if self.NCOSMs != 0:
            for i in range(NCOSMs):
                self.COSMsInt=int(self.COSMs[i][0])
                file.write(('{}\n').format(str(self.COSMsInt)+' '+str(self.COSMs[i][1])+' '+str(self.COSMs[i][2])+' '+str(self.COSMs[i][3])+' '+str(self.COSMs[i][4])+' '+str(self.COSMs[i][5])+' '+str(self.COSMs[i][6])+' '+str(self.COSMs[i][7])+' '+str(self.COSMs[i][8])+' '+str(self.COSMs[i][9])))

        file.write(('{}\n').format("------------------------ JOINT ADDITIONAL CONCENTRATED MASSES--------------------------"))

        if not self.TEST:
            if len(self.Cmass)==0 and self.RNAinputs.mass !=0:
                self.NCmass=1
            else:
                self.NCmass=len(self.Cmass)
        else:
            self.NCmass=len(self.Cmass)

        self.NCmassHeader=np.array([self.NCmass, 'NCmass', '- Number of joints with concentrated masses; Global Coordinate System'])
        file.write(('{}\n').format(' '+self.NCmassHeader[0]+' '+self.NCmassHeader[1]+' '+self.NCmassHeader[2]))
        file.write(('{}\n').format(self.CmassHeader[0][0]+' '+self.CmassHeader[0][1]+' '+self.CmassHeader[0][2]+' '+self.CmassHeader[0][3]+' '+self.CmassHeader[0][4]))
        file.write(('{}\n').format(self.CmassHeader[1][0]+' '+self.CmassHeader[1][1]+' '+self.CmassHeader[1][2]+' '+self.CmassHeader[1][3]+' '+self.CmassHeader[1][4]))

        if not self.TEST:
            if len(self.Cmass)==0 and self.RNAinputs.mass !=0:
                self.CmassInt=int(RNA_CMs[0])
                file.write(('{}\n').format(str(self.CmassInt)+' '+str(RNA_CMs[1])+' '+str(RNA_CMs[2])+' '+str(RNA_CMs[3])+' '+str(RNA_CMs[4])))
        else:
            if self.NCmass != 0:
                for i in range(self.NCmass):
                    self.CmassInt=int(self.Cmass[i][0])
                    file.write(('{}\n').format(str(self.CmassInt)+' '+str(self.Cmass[i][1])+' '+str(self.Cmass[i][2])+' '+str(self.Cmass[i][3])+' '+str(self.Cmass[i][4])))

        file.write(('{}\n').format("---------------------------- OUTPUT: SUMMARY & OUTFILE ------------------------------"))
        file.write(('{}\n').format(self.SSSum[0]+' '+self.SSSum[1]+' '+self.SSSum[2]))
        file.write(('{}\n').format(self.OutCOSM[0]+' '+self.OutCOSM[1]+' '+self.OutCOSM[2]))
        file.write(('{}\n').format(self.OutAll[0]+' '+self.OutAll[1]+' '+self.OutAll[2]))
        file.write(('{}\n').format(self.OutSwtch[0]+' '+self.OutSwtch[1]+' '+self.OutSwtch[2]))
        file.write(('{}\n').format(self.TabDelim[0]+' '+self.TabDelim[1]+' '+self.TabDelim[2]))
        file.write(('{}\n').format(self.OutDec[0]+' '+self.OutDec[1]+' '+self.OutDec[2]))
        file.write(('{}\n').format('"'+self.OutFmt[0]+'" '+self.OutFmt[1]+' '+self.OutFmt[2]))
        file.write(('{}\n').format('"'+self.OutSFmt[0]+'" '+self.OutSFmt[1]+' '+self.OutSFmt[2]))

        file.write(('{}\n').format("------------------------- MEMBER OUTPUT LIST ------------------------------------------"))
        self.NMemOutList=len(self.MemOutList)
        self.NMemOutListHeader=np.array([self.NMemOutList, 'NMOutputs', '- Number of members whose forces/displacements/velocities/accelerations will be output (-) [Must be <= 9].'])
        file.write(('{}\n').format(self.NMemOutListHeader[0]+' '+self.NMemOutListHeader[1]+' '+self.NMemOutListHeader[2]))
        file.write(('{}\n').format(self.MemOutListHeader[0][0]+' '+self.MemOutListHeader[0][1]+' '+self.MemOutListHeader[0][2]+' '+self.MemOutListHeader[0][3]))
        file.write(('{}\n').format(self.MemOutListHeader[1][0]+' '+self.MemOutListHeader[1][1]+' '+self.MemOutListHeader[1][2]))
        if self.NMemOutList!=0:
            for i in range(len(self.MemOutList)):
                if len(self.MemOutList[i])==3:
                    file.write(('{}\n').format(str(self.MemOutList[i][0])+' '+str(self.MemOutList[i][1])+' '+str(self.MemOutList[i][2])))
                if len(self.MemOutList[i])==4:
                    file.write(('{}\n').format(str(self.MemOutList[i][0])+' '+str(self.MemOutList[i][1])+' '+str(self.MemOutList[i][2])+' '+str(self.MemOutList[i][3])))

        file.write(('{}\n').format("------------------------- SSOutList: The next line(s) contains a list of output parameters that will be output in <rootname>.SD.out or <rootname>.out. ------"))
        for i in range(len(self.SSOutlist)):
            file.write(('{}\n').format('"'+str(self.SSOutlist[i][0])+'" '+str(self.SSOutlist[i][1])))

        file.write(('{}\n').format("END of output channels and end of file. (the word ""END"" must appear in the first 3 columns of this line)"))
        file.close()

        print "Write Input File function"

class WriteDriver(Component):

    gravity =Float(  units='m/s**2',iotype='in',desc='Gravity Acceleration (ABSOLUTE VALUE!)')
    wdepth =Float(  units='m',iotype='in',desc='Water depth')
    dck_botz=Float(units='m', iotype='in', desc='Distance from bottom of TP to MSL')
    nodes=Array(iotype='in', desc='Nodes from JacketSE')

    TEST=Bool(iotype='in', desc='Check to see if pySubDyn is being run locally')

    #DRIVER INPUTS--------------------------------------------------------------

    InputandDriverpath=Str(iotype='in')
    Driver_path=Str(iotype='in')

    #Inputs
    EchoD=Array(iotype='in')

    #Environmental Conditions
    Gravity=Array(iotype='in')
    WtrDpth=Array(iotype='in')

    #SubDyn
    SDInputFile=Array(iotype='in')
    OutRootName=Array(iotype='in')
    NSteps=Array(iotype='in')
    TimeInterval=Array(iotype='in')
    TP_RefPoint=Array(iotype='in')
    SubRotateZ=Array(iotype='in')

    #INPUTS
    InputsMod=Array(iotype='in')
    InputsFile=Array(iotype='in')

    #STEADY INPUTS
    uTPInSteady=Array(iotype='in')
    uDotTPInSteady=Array(iotype='in')
    uDotDotTPInSteady=Array(iotype='in')

    print "Write Driver pre"

    def execute(self):

        #Write driver
        file=open(str(self.Driver_path), 'w')
        file.write(('{}\n').format("Customized jacket driver file for SubDyn")) #TODO Set as a variable so users can describe their own jackets
        file.write(('{}\n').format("Compatible with SubDyn v1.01.xx"))

        file.write(('{}\n').format(self.EchoD[0]+' '+self.EchoD[1]+' '+self.EchoD[2]))

        file.write(('{}\n').format("---------------------- ENVIRONMENTAL CONDITIONS ------------------------------"))
        if not self.TEST:
            file.write(('{}\n').format(str(self.gravity)+' '+self.Gravity[1]+' '+self.Gravity[2]))
            file.write(('{}\n').format(str(self.wdepth)+' '+self.WtrDpth[1]+' '+self.WtrDpth[2]))
        else:
            file.write(('{}\n').format(self.Gravity[0]+' '+self.Gravity[1]+' '+self.Gravity[2]))
            file.write(('{}\n').format(self.WtrDpth[0]+' '+self.WtrDpth[1]+' '+self.WtrDpth[2]))

        file.write(('{}\n').format("---------------------- SubDyn ------------------------------------------------"))
        file.write(('{}\n').format('"'+self.SDInputFile[0]+'" '+self.SDInputFile[1]))
        file.write(('{}\n').format('"'+self.OutRootName[0]+'" '+self.OutRootName[1]))
        file.write(('{}\n').format(self.NSteps[0]+' '+self.NSteps[1]+' '+self.NSteps[2]))
        file.write(('{}\n').format(self.TimeInterval[0]+' '+self.TimeInterval[1]+' '+self.TimeInterval[2]))
        if not self.TEST:
            file.write(('{}\n').format(self.TP_RefPoint[0]+' '+self.TP_RefPoint[1]+' '+str(np.amax(self.nodes))+' '+self.TP_RefPoint[3]+' '+self.TP_RefPoint[4]))
        else:
            file.write(('{}\n').format(self.TP_RefPoint[0]+' '+self.TP_RefPoint[1]+' '+str(self.dck_botz)+' '+self.TP_RefPoint[3]+' '+self.TP_RefPoint[4]))
        file.write(('{}\n').format(self.SubRotateZ[0]+' '+self.SubRotateZ[1]+' '+self.SubRotateZ[2]))

        file.write(('{}\n').format("---------------------- INPUTS ------------------------------------------------"))
        file.write(('{}\n').format(self.InputsMod[0]+' '+self.InputsMod[1]))
        file.write(('{}\n').format(self.InputsFile[0]+' '+self.InputsFile[1]))

        file.write(('{}\n').format("---------------------- STEADY INPUTS  ----------------------------------------"))
        file.write(('{}\n').format(self.uTPInSteady[0]+' '+self.uTPInSteady[1]+' '+self.uTPInSteady[2]+' '+self.uTPInSteady[3]+' '+self.uTPInSteady[4]+' '+self.uTPInSteady[5]+' '+self.uTPInSteady[6]+' '+self.uTPInSteady[7]))
        file.write(('{}\n').format(self.uDotTPInSteady[0]+' '+self.uDotTPInSteady[1]+' '+self.uDotTPInSteady[2]+' '+self.uDotTPInSteady[3]+' '+self.uDotTPInSteady[4]+' '+self.uDotTPInSteady[5]+' '+self.uDotTPInSteady[6]+' '+self.uDotTPInSteady[7]))
        file.write(('{}\n').format(self.uDotDotTPInSteady[0]+' '+self.uDotDotTPInSteady[1]+' '+self.uDotDotTPInSteady[2]+' '+self.uDotDotTPInSteady[3]+' '+self.uDotDotTPInSteady[4]+' '+self.uDotDotTPInSteady[5]+' '+self.uDotDotTPInSteady[6]+' '+self.uDotDotTPInSteady[7]))

        file.write(('{}\n').format("END of driver input file"))
        file.close()

        print "Write Driver function"

class RunSubDyn(Component):

    #INPUTS TO RUN SUBDYN

    SDpath=Str(iotype='in')

    print "Run SubDyn pre"

    def execute(self):

        #Way 1
        #print os.popen(self.SDpath).read()
        #way 2 (Preferred)

        def run_command(command):
            print 'running ',command
            p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.PIPE, shell=False)
            return iter(p.stdout.readline, b'')

        check=[]
        for line in run_command(self.SDpath):
           check=np.hstack((check, line))
           print(line) #Prints Command Prompt's output lines

        checkline=check[4].split(' ')
        if str(checkline[len(checkline)-2]) != 'Modal':
            print 'ERROR: SubDyn did not run properly'
            quit()

        print "Run SubDyn function"

class ReadOutput(Component):

    #INPUTS TO READ OUTPUT

    Readpath_out=Str(iotype='in')
    Readpath_sum=Str(iotype='in')
    Delete_file=Bool(iotype='in')
    InputFile_path=Str(iotype='in')
    Driver_path=Str(iotype='in')

    SubDynAOuts=VarTree(SubDynAOutputs(), iotype='out', desc='Selected outputs from SubDyn')
    RNAinputs=VarTree(RNAprops(), iotype='in', desc='Basic Inertial Properties of RNA') #test

    TEST=Bool(iotype='in', desc='Check to see if pySubDyn is being run locally') #test flag

    print "Read Output pre"

    def execute(self):

        file=open(str(self.Readpath_sum), "r")

        lines = file.readlines()
        count=0

        for i in lines:
            count+=1
            if "FEM Eigenvalues [Hz]." in i: #Get FEM eigenvalues
                FEMEigval_start=count+1
                FEMEigvalMX=i.split(' ')
                NFEMEigval=int(FEMEigvalMX[len(FEMEigvalMX)-1])
                self.FEMEigval=[]
                FEMEigval1=[1.]
                FEMEigval2=[2.]
                while FEMEigval_start<=(count+1)<(FEMEigval_start+NFEMEigval):
                    self.FEMEigval.append(lines[count].split(' '))
                    count+=1
                for j in self.FEMEigval:
                    FEMEigval1=np.vstack((FEMEigval1, int(j[(len(j)-4)])))
                    FEMEigval2=np.vstack((FEMEigval2, float(j[(len(j)-1)])))
                    self.FEMEigval=np.hstack((FEMEigval1,FEMEigval2))
                self.FEMEigval=np.delete(self.FEMEigval,(0),axis=0)
            elif "SubDyn's Total Mass (structural and non-structural)=" in i:    #Get Total Mass
                TotMassMX=i.split(' ')
                self.TotMass=float(TotMassMX[len(TotMassMX)-1])
            elif "SubDyn's Total Mass CM coordinates (Xcm,Ycm,Zcm)   =" in i:  #Get coordinates of COM
                CMTotMassMX=i.split(' ')
                self.CMTotMass=[]
                if not self.TEST:
                    if CMTotMassMX[len(CMTotMassMX)-6]=='':
                        if CMTotMassMX[len(CMTotMassMX)-7] != '':
                            CMTotMass1=float(CMTotMassMX[len(CMTotMassMX)-7])
                        else:
                            CMTotMass1=float(CMTotMassMX[len(CMTotMassMX)-5])
                    else:
                        CMTotMass1=float(CMTotMassMX[len(CMTotMassMX)-6])
                    CMTotMass2=float(CMTotMassMX[len(CMTotMassMX)-4])
                    CMTotMass3=float(CMTotMassMX[len(CMTotMassMX)-1])
                else:
                    CMTotMass1=float(CMTotMassMX[len(CMTotMassMX)-5])
                    CMTotMass2=float(CMTotMassMX[len(CMTotMassMX)-3])
                    CMTotMass3=float(CMTotMassMX[len(CMTotMassMX)-1])
                self.CMTotMass.append(CMTotMass1)
                self.CMTotMass.append(CMTotMass2)
                self.CMTotMass.append(CMTotMass3)

        file.close()

        self.SubDynAOuts.FEMeigenvalues=self.FEMEigval
        self.SubDynAOuts.StructureMass=self.TotMass-self.RNAinputs.mass #TODO This removes the RNA lumped mass from the value pySubDyn returns to JacketSE, but the total mass is still written into SubDyn's output file
        #self.SubDynAOuts.StructureMass=self.TotMass
        self.SubDynAOuts.StructureCM=self.CMTotMass

        print self.SubDynAOuts.FEMeigenvalues[0]
        print self.SubDynAOuts.FEMeigenvalues[1]
        print self.SubDynAOuts.StructureMass
        print self.SubDynAOuts.StructureCM

        print "Read Output function"

        if self.Delete_file==True:
            os.remove(self.Driver_path)
            os.remove(self.InputFile_path)
            os.remove(self.Readpath_out)
            os.remove(self.Readpath_sum)


class pySubDynA(Assembly):

    gravity =Float(  units='m/s**2',iotype='in',desc='Gravity Acceleration (ABSOLUTE VALUE!)')
    wdepth =Float(  units='m',iotype='in',desc='Water depth')
    dck_botz=Float(units='m', iotype='in', desc='Distance from bottom of TP to MSL')
    nodes=Array(iotype='in', desc='Nodes from JacketSE')
    mems=Array(iotype='in',desc='Member connections between')
    Reacts=Array(iotype='in', desc='Base joint IDs and fixity')
    SDPropSet=Array(iotype='in', desc='Properties to be used in SubDyn')
    RNAinputs=VarTree(RNAprops(), iotype='in', desc='Basic Inertial Properties of RNA')

    TEST=Bool(iotype='in')

    #INPUT FILE INPUTS----------------------------------------------------------

    InputFile_name=Str(iotype='in')

    #Simulation Control
    Echo=Array(iotype='in')
    SDdeltaT=Array(iotype='in')
    IntMethod=Array(iotype='in')
    SttcSolve=Array(iotype='in')

    #FEA and CRAIG-BAMPTON PARAMETERS
    FEMmod=Array(iotype='in')
    NDiv=Array(iotype='in')
    CBMod=Array(iotype='in')
    Nmodes=Array(iotype='in')
    JDampings=Array(iotype='in')

    #Structure Joints
    SDjointsHeader=Array(iotype='in')
    SDjoints=Array(iotype='in', desc='Leg joints from jacket.py. In SD joints are used to build the structure, not nodes.')

    #Base Reaction Joints
    BaseRxnJointsHeader=Array(iotype='in')
    BaseRxnJoints=Array(iotype='in', desc='Base joints of the structure.')

    #Interface Joints
    InterfaceRxnJointsHeader=Array(iotype='in')
    InterfaceJoints=Array(iotype='in', desc='Interface joints.')

    #Members
    MembersHeader=Array(iotype='in')
    Members=Array(iotype='in', desc='Build members from joints, not nodes.')

    #MEMBER X-SECTION PROPERTY data 1/2
    NPropSets=Array(iotype='in')
    PropSet1Header=Array(iotype='in')
    PropSet1=Array(iotype='in')

    #MEMBER X-SECTION PROPERTY data 2/2
    PropSet2Header=Array(iotype='in')
    PropSet2=Array(iotype='in')

    #MEMBER COSINE MATRICES COSM(i,j)
    COSMHeader=Array(iotype='in')
    COSMs=Array(iotype='in')

    #JOINT ADDITIONAL CONCENTRATED MASSES
    CmassHeader=Array(iotype='in')
    Cmass=Array(iotype='in')

    #OUTPUT: SUMMARY & OUTFILE
    SSSum=Array(iotype='in')
    OutCOSM=Array(iotype='in')
    OutAll=Array(iotype='in')
    OutSwtch=Array(iotype='in')
    TabDelim=Array(iotype='in')
    OutDec=Array(iotype='in')
    OutFmt=Array(iotype='in')
    OutSFmt=Array(iotype='in')

    #MEMBER OUTPUT LIST
    MemOutListHeader=Array(iotype='in')
    MemOutList=Array(iotype='in', desc='Members whose loads and dynamics will be output.')

    #SSOutline
    SSOutlist=Array(iotype='in')

    #DRIVER INPUTS--------------------------------------------------------------

    InputandDriverpath=Str(iotype='in')
    Driver_path=Str(iotype='in')

    EchoD=Array(iotype='in')

    #Environmental Conditions
    Gravity=Array(iotype='in')
    WtrDpth=Array(iotype='in')

    #SubDyn
    SDInputFile=Array(iotype='in')
    OutRootName=Array(iotype='in')
    NSteps=Array(iotype='in')
    TimeInterval=Array(iotype='in')
    TP_RefPoint=Array(iotype='in')
    SubRotateZ=Array(iotype='in')

    #INPUTS
    InputsMod=Array(iotype='in')
    InputsFile=Array(iotype='in')

    #STEADY INPUTS
    uTPInSteady=Array(iotype='in')
    uDotTPInSteady=Array(iotype='in')
    uDotDotTPInSteady=Array(iotype='in')

    #INPUTS TO RUN SUBDYN-------------------------------------------------------

    SDpath=Str(iotype='in')

    #INPUTS TO READ OUTPUT------------------------------------------------------

    Readpath_out=Str(iotype='in')
    Readpath_sum=Str(iotype='in')
    Delete_file=Bool(iotype='in')
    InputFile_path=Str(iotype='in')
    Driver_path=Str(iotype='in')

    #OUTPUTS--------------------------------------------------------------------
    SubDynAOuts=VarTree(SubDynAOutputs(), iotype='out', desc='Selected outputs from SubDyn') #test

    print "pySubDynA pre"

    def configure(self):
        self.add('WIF', WriteInputFile())
        self.add('WD', WriteDriver())
        self.add('RSD', RunSubDyn())
        self.add('RO', ReadOutput())

        self.driver.workflow.add(['WIF', 'WD', 'RSD', 'RO'])

        self.connect('TEST',['WIF.TEST','WD.TEST','RO.TEST']) #test flag

        #INPUT FILE CONNECTIONS-------------------------------------------------

        self.connect('InputFile_name', 'WIF.InputFile_name')
        self.connect('SDInputFile', 'WIF.SDInputFile')

        #Simulation Control
        self.connect('Echo', 'WIF.Echo')
        self.connect('SDdeltaT', 'WIF.SDdeltaT')
        self.connect('IntMethod', 'WIF.IntMethod')
        self.connect('SttcSolve', 'WIF.SttcSolve')

        #FEA and CRAIG-BAMPTON PARAMETERS
        self.connect('FEMmod', 'WIF.FEMmod')
        self.connect('NDiv', 'WIF.NDiv')
        self.connect('CBMod', 'WIF.CBMod')
        self.connect('Nmodes', 'WIF.Nmodes')
        self.connect('JDampings', 'WIF.JDampings')

        #Structure Joints
        self.connect('SDjointsHeader', 'WIF.SDjointsHeader')
        self.connect('SDjoints', 'WIF.SDjoints')

        #Base Reaction Joints
        self.connect('BaseRxnJointsHeader', 'WIF.BaseRxnJointsHeader')
        self.connect('BaseRxnJoints', 'WIF.BaseRxnJoints')

        #Interface Joints
        self.connect('InterfaceRxnJointsHeader', 'WIF.InterfaceRxnJointsHeader')
        self.connect('InterfaceJoints', 'WIF.InterfaceJoints')

        #Members
        self.connect('MembersHeader', 'WIF.MembersHeader')
        self.connect('Members', 'WIF.Members')

        #MEMBER X-SECTION PROPERTY data 1/2
        self.connect('NPropSets', 'WIF.NPropSets')
        self.connect('PropSet1Header', 'WIF.PropSet1Header')
        self.connect('PropSet1', 'WIF.PropSet1')

        #MEMBER X-SECTION PROPERTY data 2/2
        self.connect('PropSet2','WIF.PropSet2')
        self.connect('PropSet2Header', 'WIF.PropSet2Header')

        #MEMBER COSINE MATRICES COSM(i,j)
        self.connect('COSMHeader', 'WIF.COSMHeader')
        self.connect('COSMs', 'WIF.COSMs')

        #JOINT ADDITIONAL CONCENTRATED MASSES
        self.connect('CmassHeader', 'WIF.CmassHeader')
        self.connect('Cmass', 'WIF.Cmass')

        #OUTPUT: SUMMARY & OUTFILE
        self.connect('SSSum', 'WIF.SSSum')
        self.connect('OutCOSM', 'WIF.OutCOSM')
        self.connect('OutAll', 'WIF.OutAll')
        self.connect('OutSwtch', 'WIF.OutSwtch')
        self.connect('TabDelim', 'WIF.TabDelim')
        self.connect('OutDec', 'WIF.OutDec')
        self.connect('OutFmt', 'WIF.OutFmt')
        self.connect('OutSFmt', 'WIF.OutSFmt')

        #MEMBER OUTPUT LIST
        self.connect('MemOutListHeader', 'WIF.MemOutListHeader')
        self.connect('MemOutList', 'WIF.MemOutList')

        #SSOutline
        self.connect('SSOutlist', 'WIF.SSOutlist')

        #DRIVER CONNECTIONS-----------------------------------------------------

        self.connect('InputandDriverpath', 'WD.InputandDriverpath')
        self.connect('Driver_path', 'WD.Driver_path')

        self.connect('EchoD', 'WD.EchoD')

        #Environmental Conditions
        self.connect('Gravity', 'WD.Gravity')
        self.connect('WtrDpth', 'WD.WtrDpth')

        #SubDyn
        self.connect('SDInputFile', 'WD.SDInputFile')
        self.connect('OutRootName', 'WD.OutRootName')
        self.connect('NSteps', 'WD.NSteps')
        self.connect('TimeInterval', 'WD.TimeInterval')
        self.connect('TP_RefPoint', 'WD.TP_RefPoint')
        self.connect('SubRotateZ', 'WD.SubRotateZ')

        #INPUTS
        self.connect('InputsMod', 'WD.InputsMod')
        self.connect('InputsFile', 'WD.InputsFile')

        #STEADY INPUTS
        self.connect('uTPInSteady', 'WD.uTPInSteady')
        self.connect('uDotTPInSteady', 'WD.uDotTPInSteady')
        self.connect('uDotDotTPInSteady', 'WD.uDotDotTPInSteady')

        #INPUTS TO RUN SUBDYN---------------------------------------------------

        self.connect('SDpath', 'RSD.SDpath')

        #INPUTS TO READ OUTPUT--------------------------------------------------
        self.connect('Readpath_out', 'RO.Readpath_out')
        self.connect('Readpath_sum', 'RO.Readpath_sum')
        self.connect('Delete_file', 'RO.Delete_file')
        self.connect('InputFile_path','RO.InputFile_path')
        self.connect('Driver_path','RO.Driver_path')

        #OUTPUT CONNECTIONS-----------------------------------------------------
        self.connect('RO.SubDynAOuts','SubDynAOuts')

        #CONNECTIONS FROM JACKET------------------------------------------------
        self.connect('gravity', 'WIF.gravity') #CJB+
        self.connect('gravity', 'WD.gravity') #CJB+
        self.connect('wdepth', 'WD.wdepth') #CJB+
        self.connect('dck_botz', 'WD.dck_botz') #CJB+
        self.connect('nodes', ['WIF.nodes', 'WD.nodes']) #CJB+
        self.connect('mems', 'WIF.mems') #CJB+
        self.connect('Reacts', 'WIF.Reacts') #CJB+
        self.connect('SDPropSet', 'WIF.SDPropSet') #CJB+
        self.connect('RNAinputs', ['WIF.RNAinputs', 'RO.RNAinputs']) #CJB+

        print "pySubDynA function"


if __name__ == '__main__':

    #Example inputs correspond to Test03 under CertTest, an OC4 'Jacket' SIMPLIFIED, 1-bay

    test = pySubDynA()

    test.TEST=True #Test to see if the code is being run locally

    #INPUTS TO RUN SUBDYN-------------------------------------------------------
    #(PATH INFORMATION FOR INPUT FILE AND DRIVER)

    test.Base_name="pySubDynTest" #Input name here

    #INPUT FILE PATH
    test.InputFile_name=str(test.Base_name)+".txt"
    test.InputandDriverpath="C:\wisdem\plugins\JacketSE\src\jacketse\SubDyn\CertTest"
    test.InputFile_path=str(test.InputandDriverpath)+"\\"+str(test.InputFile_name)

    #DRIVER PATH
    test.Driver_name=str(test.Base_name)+"D"+".txt"
    test.Driver_path=str(test.InputandDriverpath)+"\\"+str(test.Driver_name)

    #PATH TO RUN SUBDYN
    test.SDEXEpath="C:\wisdem\plugins\JacketSE\src\jacketse\SubDyn"+"\\"+"bin\SubDyn_Win32.exe"
    test.SDpath=str(test.SDEXEpath)+' '+str(test.Driver_path)
    #test.SDpath='r'+"'''"+test.SDEXEpath+' '+test.SDDriverpath+"'''"

    #PATH TO READ OUTPUT (INPUTS TO READ OUTPUT)
    test.Readpath_out=str(test.InputandDriverpath)+"\\"+str(test.Base_name)+".SD.out"
    test.Readpath_sum=str(test.InputandDriverpath)+"\\"+str(test.Base_name)+".SD.sum"
    test.Delete_file=False #Deletes driver, input, and output files. Does not delete Echo file.

    #INPUT FILE INPUTS----------------------------------------------------------

    #Simulation Control
    test.Echo=np.array([True, "Echo", "- Echo input data to ""<rootname>.SD.ech"" (flag)"])
    test.SDdeltaT=np.array(["DEFAULT", "SDdeltaT", "- Local Integration Step. If ""default"", the glue-code integration step will be used."])
    test.IntMethod=np.array([4, "IntMethod", "- Integration Method [1/2/3/4 = RK4/AB4/ABM4/AM2]."])
    test.SttcSolve=np.array([False, "SttcSolve", "- Solve dynamics about static equilibrium point"])

    #FEA and CRAIG-BAMPTON PARAMETERS
    test.FEMmod=np.array([3, "- FEM switch: element model in the FEM. [1= Euler-Bernoulli(E-B);  2=Tapered E-B (unavailable);  3= 2-node Timoshenko;  4= 2-node tapered Timoshenko (unavailable)]"])
    test.NDiv=np.array([1, "NDiv", "- Number of sub-elements per member"]) #CJB "HARDWIRED" INTO THE CODE AS 1 TO ALLOW JACKETSE'S NODES TO BE USED AS SUBDYN'S JOINTS
    test.CBMod=np.array([False, "- [T/F] If True perform C-B reduction, else full FEM dofs will be retained. If True, select Nmodes to retain in C-B reduced system."])
    test.Nmodes=np.array([75, "- Number of internal modes to retain (ignored if CBMod=False). If Nmodes=0 --> Guyan Reduction."])
    test.JDampings=np.array([2, "JDampings", "- Damping Ratios for each retained mode (% of critical) If Nmodes>0, list Nmodes structural damping ratios for each retained mode (% of critical), or a single damping ratio to be applied to all retained modes. (last entered value will be used for all remaining modes)."])

    #Structure Joints
    test.SDjointsHeader=np.array([['JointID','JointXss','JointYss','JointZss','[Coordinates of Member joints in SS-Coordinate System]'], ['(-)','(m)','(m)','(m)']])

    #Base Reaction Joints
    test.BaseRxnJointsHeader=np.array([['RJointID','RctTDXss','RctTDYss','RctTDZss','RctRDXss','RctRDYss','RctRDZss','[Global Coordinate System]'], ['(-)',('flag'),('flag'),('flag'),('flag'),('flag'),('flag')]])

    #Interface Joints
    test.InterfaceRxnJointsHeader=np.array([['IJointID','ItfTDXss','ItfTDYss','ItfTDZss','ItfRDXss','ItfRDYss','ItfRDZss','[Global Coordinate System]'], ['(-)',('flag'),('flag'),('flag'),('flag'),('flag'),('flag')]])

    #Members
    test.MembersHeader=np.array([['MemberID','MJointID1','MJointID2','MPropSetID1','MPropSetID2','COSMID'], ['(-)','(-)','(-)','(-)','(-)','(-)']])

    #MEMBER X-SECTION PROPERTY data 1/2
    test.NPropSets=np.array([6, 'NPropSets', '- # of structurally unique x-sections (i.e. # of X-sectional props utilized throughout all of the members)'])
    test.PropSet1Header=np.array([['PropSetID', 'YoungE', 'ShearG', 'MatDens', 'XsecD', 'XsecT'],['(-)', '(N/m2)', '(N/m2)', '(kg/m3)', '(m)', '(m)']])
    test.PropSet1=np.array([[1, 2.10000e+11, 8.07690e+10, 7850.00, 0.800000, 0.020000],\
                   [2, 2.10000e+11, 8.07690e+10, 7850.00, 1.200000, 0.050000],\
                   [3, 2.10000e+11, 8.07690e+10, 7850.00, 1.200000, 0.035000],\
                   [4, 2.10000e+11, 8.07690e+10, 7850.00, 1.200000, 0.040000],\
                   [5, 2.10000e+11, 8.07690e+10, 3339.12, 2.082000, 0.491000],\
                   [6, 2.10000e+11, 8.07690e+10, 7850.00, 2.082000, 0.060000]])

    #MEMBER X-SECTION PROPERTY data 2/2
    test.PropSet2=np.array([])
    test.PropSet2Header=np.array([["PropSetID","YoungE","ShearG","MatDens","XsecA","XsecAsx","XsecAsy","XsecJxx","XsecJyy","XsecJ0"],["(-)","(N/m2)","(N/m2)","(kg/m3)","(m2)","(m2)","(m2)","(m4)","(m4)","(m4)"]])

    #MEMBER COSINE MATRICES COSM(i,j)
    test.COSMHeader=np.array([["COSMID","COSM11","COSMID12","COSMID13","COSMID21","COSMID22","COSMID23","COSMID31","COSMID32","COSMID33"],["(-)","(-)","(-)","(-)","(-)","(-)","(-)","(-)","(-)","(-)"]])
    test.COSMs=np.array([])

    #JOINT ADDITIONAL CONCENTRATED MASSES
    test.CmassHeader=np.array([["CMJointID","JMass","JMXX","JMYY","JMZZ"],["(-)","(kg)" ,"(kg*m^2)","(kg*m^2)","(kg*m^2)"]])
    test.Cmass=np.array([])

    #OUTPUT: SUMMARY & OUTFILE
    test.SSSum=np.array([True, "SSSum", "- Output a Summary File (flag).It contains: matrices K,M  and C-B reduced M_BB, M-BM, K_BB, K_MM(OMG^2), PHI_R, PHI_L. It can also contain COSMs if requested."])
    test.OutCOSM=np.array([True, "OutCOSM", "- Output cosine matrices with the selected output member forces (flag)"])
    test.OutAll=np.array([True, "OutAll", "- [T/F] Output all members' end forces "])
    test.OutSwtch=np.array([1, "OutSwtch", "- [1/2/3] Output requested channels to: 1=<rootname>.SD.out;  2=<rootname>.out (generated by FAST);  3=both files."])
    test.TabDelim=np.array([True, "TabDelim", "- Generate a tab-delimited output in the <rootname>.SD.out file"])
    test.OutDec=np.array([1, "OutDec", "- Decimation of output in the <rootname>.SD.out file"])
    test.OutFmt=np.array(["Es11.4e2", "OutFmt", "- Output format for numerical results in the <rootname>.SD.out file"])
    test.OutSFmt=np.array(["A11", "OutFmt", "- Output format for header strings in the <rootname>.SD.out file"])

    #MEMBER OUTPUT LIST
    test.MemOutListHeader=np.array([['MemberID','NoutCnt','NodeCnt','[NOutCnt=how many nodes to get output for [< 10]; NodeCnt are local ordinal numbers from the start of the member, and must be >=1 and <= NDiv+1] If NMOutputs=0 leave blank as well.]'],\
                                           ['(-)','(-)','(-)']])

    #SSOutline
    test.SSOutlist=np.array([["ReactFXss, ReactFYss, ReactFZss, ReactMXss, ReactMYss, ReactMZss",'-Base reactions (forces onto SS structure)'],\
                    ["IntfFXss,  IntfFYss,  IntfFZss,  IntfMXss, IntfMYss, IntfMZss",'-Interface reactions (forces from SS structure)'],\
                    ["IntfTDXss,  IntfTDYss,  IntfTDZss,  IntfRDXss, IntfRDYss, IntfRDZss",'-Interface deflections '],\
                    ["IntfTAXss,  IntfTAYss,  IntfTAZss,  IntfRAXss, IntfRAYss, IntfRAZss",'Interface accelerations']])

    test.SDjoints=np.array([[1, -5.93900, -5.93900, -43.12700],\
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

    test.BaseRxnJoints=np.array([[1,1,1,1,1,1,1],\
                                 [2,1,1,1,1,1,1],\
                                 [3,1,1,1,1,1,1],\
                                 [4,1,1,1,1,1,1]])

    test.InterfaceJoints=np.array([[5,1,1,1,1,1,1],\
                                   [6,1,1,1,1,1,1],\
                                   [7,1,1,1,1,1,1],\
                                   [8,1,1,1,1,1,1]])

    test.Members=np.array([[1, 1, 5, 2, 2],\
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

    test.MemOutList=np.array([[1,2,1,2]])

    #DRIVER INPUTS--------------------------------------------------------------

    test.EchoD=np.array([True, "Echo", "- Echo the input file data (flag)"])

    #Environmental Conditions
    test.Gravity=np.array([9.81, "Gravity", "- Gravity (m/s^2)"])
    test.WtrDpth=np.array([43.127, "WtrDpth", "- Water Depth (m) positive value"])

    #SubDyn
    test.SDInputFile=np.array([test.InputFile_path, "SDInputFile"])
    test.OutRootName=np.array([str(test.InputandDriverpath)+"\\"+str(test.Base_name), "OutRootName"])
    test.NSteps=np.array([600, "NSteps", "- Number of time steps in the simulations (-)"])
    test.TimeInterval=np.array([0.005, "TimeInterval", "- TimeInterval for the simulation (sec)"])
    test.TP_RefPoint=np.array([0.0, 0.0, 18.15, "TP_RefPoint", "- Location of the TP reference point in global coordinates (m)"])
    test.SubRotateZ=np.array([0.0, "SubRotateZ", "- Rotation angle of the structure geometry in degrees about the global Z axis."])

    #INPUTS
    test.InputsMod=np.array([1, "InputsMod", "- Inputs model {0: all inputs are zero for every timestep, 1: steadystate inputs, 2: read inputs from a file (InputsFile)} (switch)"])
    test.InputsFile=np.array(['""', "InputsFile", "- Name of the inputs file if InputsMod = 2"])

    #STEADY INPUTS
    test.uTPInSteady=np.array([0.1, 0.0, 0.0, 0.0, 0.0, 0.0, "uTPInSteady", "- input displacements and rotations ( m, rads )"])
    test.uDotTPInSteady=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "uDotTPInSteady", "- input translational and rotational velocities ( m/s, rads/s)"])
    test.uDotDotTPInSteady=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "uDotDotTPInSteady", "- input translational and rotational accelerations ( m/s^2, rads/s^2)"])

    test.run()
