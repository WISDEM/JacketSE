#-------------------------------------------------------------------------------
# Name:        Utilization.py
# Purpose: Assembly for the calculation of utilization ratios of Tower And Jacket
#          !!!!!Still to do pile lateral capacity!!!!!
# Author:      rdamiani/sning
#
# Copyright:   (c) rdamiani 2014
# Licence:     <Apache>
#
# History:          17/02/2014 Created
#                    10/2014-found that in the tower utilization the self-weight was missing, added. This must have happened when getting Andrew's utilization from frame3dd which accounts for it.
#                          -utilization is now calculated as such in andrew's code, before tower supplement was returning (utilization-1), this had to be corrected here where the +1 factors have been removed
#-------------------------------------------------------------------------------
import numpy as np
from openmdao.main.api import VariableTree, Component, Assembly
from openmdao.main.datatypes.api import Int, Float, Array, VarTree, Bool

# from commonse.Tube import Tube
# from commonse.RigidMember import RigidMember
import ApiCodeChecks as API
from towerse.tower import AeroLoads
from commonse.rna import RotorLoads
from towerse.towerSupplement import hoopStressEurocode, shellBucklingEurocode, bucklingGL, vonMisesStressUtilization
from VarTrees import JcktGeoOutputs, TwrGeoOutputs
from commonse.Frustum import frustum
#______________________________________________________________________________#


class TowerUtilOutputs(VariableTree):
    """Jacket Utilization Basic Outputs"""
    StressUtil = Array(dtype=np.float, desc='Tower Utilization at each given section (stress)')
    GLUtil = Array(dtype=np.float, desc='Tower Utilization at each given section (global buckling)')
    EUshUtil = Array(dtype=np.float, desc='Tower Utilization at each given section (shell buckling)')


class JacketUtilOutputs(VariableTree):
    """Jacket Utilization Basic Outputs"""

    t_util = Array(dtype=np.float, desc='Member Utilization (tensile strength)')
    cb_util = Array(dtype=np.float, desc='Member Utilization (compression-bending buckling)')
    XjntUtil = Array(dtype=np.float, desc='X-Joint Utilization')
    KjntUtil = Array(dtype=np.float, desc='K-Joint Utilization')


#______________________________________________________________________________#

class JcktUtilization(Component):
    """This Component Calculates ULS Maximum Member and Joint Utilization per API Checks, \n
    additionally calculates ULS Maximum Utilization for Tower following GL and Eurocode.\n
    Finally, it calculates the embedment pile length and mass.
    \n
    INPUT
    JcktGeoOut       - Jacket geometry outpyt fully built through BuildGeometry. \n
    """

    #inputs
    Pilemems    =Array(dtype=int,     iotype='in',desc='Connectivity Array for all members of the Pile 1 portion only.') # TO DO just need the size of these 4 arrays
    Legmems     =Array(dtype=int,     iotype='in',desc='Connectivity Array for all members of the Leg 1 portion only.')
    Twrmems     =Array(dtype=int,     iotype='in',desc='Connectivity Array for all members of the tower portion only (including rigid member at the top if requested).')
    TPmems      =Array(dtype=int,     iotype='in',desc='Connectivity Array for all members of the TP portion only.')


    MbrFrcs   =   Array(dtype=np.float,    iotype='in', desc='Forces and Moments at Element Nodes from Frame3DD')

    JcktGeoOut = VarTree(JcktGeoOutputs(), iotype='in', desc='Geometry of the Jacket -Node Coordinates and Member Connectivity')
    nlegs =Int(             units=None,    iotype='in', desc='Number of Legs')

    XjntIDs     =Array(dtype=int,     iotype='in',desc='X-Joint IDs, as for Frame3DD, used to do checks later')
    KjntIDs     =Array(dtype=int,     iotype='in',desc='K-Joint IDs, as for Frame3DD, used to do checks later')


    #outputs
    Utilouts=   VarTree(JacketUtilOutputs(), iotype='out', desc='Jacket Utilization Basic Outputs')

    def execute(self):
        #Simplify nomenclature
        nlegs=self.nlegs
        allmems=self.JcktGeoOut.mems
        Pilemems=self.Pilemems
        Legmems=self.Legmems
        Twrmems=self.Twrmems
        TPmems=self.TPmems

        #self.JcktGeoOut.mems=np.vstack((self.Pilemems,self.Legmems,LLURmems,ULLRmems,Mudmems,Hmems,Stmpmems,Stemmems,Girdmems,Brcmems,Strtmems,self.Twrmems))
        nfrcs=self.MbrFrcs.shape[-1]
        MbrFrcs=np.empty([nfrcs,9])
        MbrFrcs[:,:-1]=self.MbrFrcs.squeeze().T  #This is now [nElements*2,8] where the cols are [Elem No., Node No, N, Vx,Vy, Mxx,Myy,Mzz]
        #augment MbrFrcs with the ToC flag
        ToC=np.sign(MbrFrcs[:,2])*np.tile([-1.,1.],nfrcs/2) #ToC:  1 or -1 depending on tension or compression.
        MbrFrcs[:,8]=ToC

        #(f1,f2,mass,MbrFrcs,Reacts)=BuildJckt.ReadOutFrame3DD(OutFrame3DD)
        #Nhead=np.max(Reacts[:,3])  #Fz in global coordinates

        #Member checks
        #Number of elements to check (here mems is really elements)
        n_allmems=allmems.shape[0]  #number of mems across the entire OWT
        n_pilemems=Pilemems.shape[0] #number of mems in piles
        n_twrmems=Twrmems.shape[0] #number of mems in tower
        n_TPmems=TPmems.shape[0] #number of mems in TP
        n_jcktmems=n_allmems-n_pilemems-n_twrmems# -n_TPmems #Number of members belonging to jacket alone, no piles no tower no TP



        mbr_idx=n_pilemems+np.arange(0,n_jcktmems) #indices of members for which we will do API member-checks
        #prepare other parameters for the checks
        mbr_strct=self.JcktGeoOut.TubeObjs
        mbr_strct.XnsfM=self.JcktGeoOut.XnsfM  #Augment structure with this, used for joint checks

        t_chk,t_util, cb_chk, cb_util=API.ApiMbrChk(MbrFrcs,mbr_idx,mbr_strct)
        self.Utilouts.t_util=t_util.astype(float)
        #self.Utilouts.t_util[self.Utilouts.t_util<0]=np.NaN  #Remove this nan for optimizer's sake
        max_tutil=np.nanmax(t_util)
        self.Utilouts.cb_util=cb_util.astype(float)
        #self.Utilouts.cb_util[self.Utilouts.cb_util<0]=np.NaN #Remove this nan for optimizer's sake
        max_cbutil=np.nanmax(cb_util)
        #__________________________________________#
        #Joint checks
        #We need to identify node IDs for XJoints and KJoints

        #Take care of only certain members (used to remove TP stuff) via this new list of indices
        junk=mbr_idx*2
        mbr_idx2=np.array([junk,junk+1]).flatten()
        mbr_idx2.sort()

        XjntChk,self.Utilouts.XjntUtil=API.XjntChk(self.XjntIDs,self.JcktGeoOut.mems,MbrFrcs[mbr_idx2,:],mbr_strct)
        max_XjntUtil=np.nanmax(self.Utilouts.XjntUtil)
        #Lets us check the K-joints now
        KjntChk,self.Utilouts.KjntUtil=API.KjntChk(self.KjntIDs,self.JcktGeoOut.mems,MbrFrcs[mbr_idx2,:],mbr_strct)
        max_KjntUtil=np.nanmax(self.Utilouts.KjntUtil)


class IEC_PSFS(VariableTree):
    """IEC Basic PSFs"""
    gamma_f = Float(1.35, desc='Safety factor on loads')
    gamma_m = Float(1.1,  desc='Safety factor on materials')
    gamma_n = Float(1.0,  desc='Safety factor on consequence of failure')
    gamma_b = Float(1.1,  desc='Safety factor for buckling')
    gamma_g = Float(1.1,  desc='Safety factor for gravity loads')


class TwrUtilization(Component):
    "THIS UNTILIZATION WORKS WITH WIND AND MAIN THRUST LOADS ALONG GLOBAL X ONLY"
    #inputs
    towerWindLoads = VarTree(AeroLoads(), iotype='in', desc='Aero loads in inertial coordinate system.')
    towerWaveLoads = VarTree(AeroLoads(), iotype='in', desc='Hydro loads in inertial coordinate system.')
    top_F = Array(iotype='in')
    top_M = Array(iotype='in')
    Dt = Float(iotype='in', units='m', desc='TowerTop OD from Tower.')
    tt = Float(iotype='in', units='m', desc='TowerTop wall thickness from Tower.')
    Twr_data = VarTree(TwrGeoOutputs(), iotype='in', desc='Tower Node data.')  # From Tower
    L_reinforced = Float(30.0, iotype='in', desc='reinforcement length.')
    IECpsfIns= VarTree(IEC_PSFS(), iotype='in', desc='Basic IEC psfs.')
    g = Float(9.81, iotype='in', units='m/s**2', desc='Gravity Acceleration (ABSOLUTE VALUE!).')

    #outputs
    utilization = VarTree(TowerUtilOutputs(), iotype='out')
    # stress = Array(iotype='out', units='N/m**2', desc='von Mises stress along tower on downwind side (yaw-aligned +x).  normalized by yield stress.  includes safety factors.')
    # shell_buckling = Array(iotype='out', desc='shell buckling utilization.  should be <= 1 for feasibility.  includes safety factors')
    # tower_buckling = Array(iotype='out', desc='tower buckling constraint.  should be <= 1 for feasibility.  includes safety factors')

    def execute(self):

        #simplify nomenclature
        gamma_f = self.IECpsfIns.gamma_f
        gamma_m = self.IECpsfIns.gamma_m
        gamma_n = self.IECpsfIns.gamma_n
        gamma_b = self.IECpsfIns.gamma_b
        gamma_g = self.IECpsfIns.gamma_g

        z = self.towerWindLoads.z  # this should be the same for wind and wave
        Px = self.towerWindLoads.Px + self.towerWaveLoads.Px
        Py = self.towerWindLoads.Py + self.towerWaveLoads.Py
        Pz = self.towerWindLoads.Pz + self.towerWaveLoads.Pz

        #Make it all along x
        Px=np.sqrt(Px**2+Py**2)
        Py *=0.

        Twr = self.Twr_data.TwrObj
        # #Augment with z etc as needed by JTwrUtilization
        # Twr.nodes=self.Twrouts.nodes
        # Twr.RNAmass=self.Twrouts.TopMass[0]
        # Twr.CMzoff=self.Twrouts.TopMass[-1]
        # Twr.RNA_Thrust=np.sqrt((self.RNA_F[0:2]**2).sum())
        # Twr.RNA_M=self.RNA_F[3:]
        # #GLUtil,EUshUtil=TwrUtil.JTwrUtilization(Twr,self.wind_dict,self.water_dict,wind_load=self.twr_wndload,water_load=self.twr_wtrload,ploton=False,gravity=self.gravity)#
        # #max_GLUtil=np.nanmax(GLUtil)
        # #max_EUUtil=np.nanmax(EUshUtil)

        twr_z = self.Twr_data.nodes[2, :]
        d_top = self.Dt
        t_top = self.tt
        twr_D = np.hstack((Twr.D, d_top))  # add top diameter
        twr_t = np.hstack((Twr.t, t_top))  # add top thickness
        twr_A = np.hstack((Twr.Area, np.pi/4.*(d_top**2-(d_top-2.*t_top)**2)))  # add top diameter station
        # twr_Amid = np.pi/4.*(twr_D-twr_t)**2  # Mid surface inscribed area (for torsion)
        twr_Jyy = np.hstack((Twr.Jyy, np.pi/64.*(d_top**4-(d_top-2.*t_top)**4)))  # add top diameter station

        # if twr_z.size > twr_D.size:  # this may include a fictitious z for the attachment to RNA (CMzoff)
        #     twr_z = twr_z[:-1]  # pop the last element #Need to remove RNA joint, not part of tower

        # n = twr_z.shape[0]  # number of stations along the tower (removing top one in case there is a rigid link on top)

        #Take care of the fact we get materials spearately for each tower element
        rho = np.array([mat.rho for mat in Twr.mat])  # I wonder whether these 3 could be condensed into 1 loop instead of 3
        E = np.array([mat.E for mat in Twr.mat])
        fy = np.array([mat.fy for mat in Twr.mat])
        #replicate the last item since mat was for elements not nodes
        rho = np.hstack((rho, rho[-1]))
        E = np.hstack((E, E[-1]))
        fy = np.hstack((fy, fy[-1]))



        #Calculate internal loads
        # MTTg = Twr.RNAmass*self.gravity  # [N]  Tower top weight

        n = len(z)

        Vx = np.zeros(n)
        Vy = np.zeros(n)
        Fz = np.zeros(n)
        Mx = np.zeros(n)
        My = np.zeros(n)
        Tz = np.zeros(n)
        Vx[-1] = self.top_F[0]
        Vy[-1] = self.top_F[1]
        Fz[-1] = self.top_F[2]
        Mx[-1] = self.top_M[0]
        My[-1] = self.top_M[1]
        Tz[-1] = self.top_M[2]

        for i in reversed(range(n-1)):
            delta_z=z[i+1]-z[i]
            vol=frustum(twr_D[i],twr_D[i+1],delta_z)[0]-frustum(twr_D[i]-2.*twr_t[i], twr_D[i+1]-2.*twr_t[i+1], delta_z)[0]

            Vx[i] = Vx[i+1] + 0.5*(Px[i] + Px[i+1])*delta_z
            Vy[i] = Vy[i+1] + 0.5*(Py[i] + Py[i+1])*delta_z
            Fz[i] = Fz[i+1] + 0.5*(Pz[i] + Pz[i+1])*delta_z - 0.5*(rho[i]+rho[i+1])*self.g*vol  #This was missing self weight of tower shell

            Mx[i] = Mx[i+1] + Vy[i+1]*(z[i+1]-z[i]) + (Py[i]/6.0 + Py[i+1]/3.0)*(z[i+1]-z[i])**2
            My[i] = My[i+1] + Vx[i+1]*(z[i+1]-z[i]) + (Px[i]/6.0 + Px[i+1]/3.0)*(z[i+1]-z[i])**2
            Tz[i] = Tz[i+1]



        L_reinforced = self.L_reinforced*np.ones_like(twr_D)

        # axial and shear stress (all stress evaluated on +x yaw side)
        axial_stress = Fz/twr_A - np.sqrt(Mx**2+My**2)/twr_Jyy*twr_D/2.0
        shear_stress = 2 * np.sqrt(Vx**2+Vy**2) / twr_A
        hoop_stress = hoopStressEurocode(self.towerWindLoads, self.towerWaveLoads,
            z, twr_D, twr_t, L_reinforced)

        # von mises stress
        VMutil = vonMisesStressUtilization(axial_stress, hoop_stress, shear_stress,
            gamma_f*gamma_m*gamma_n, fy)

        self.utilization.StressUtil = VMutil   #  utilization

        # shell buckling
        shell_buckling = shellBucklingEurocode(twr_D, twr_t, axial_stress, hoop_stress,
            shear_stress, L_reinforced, E, fy, gamma_f, gamma_b)

        self.utilization.EUshUtil = shell_buckling  #utilization

        # global buckling
        tower_height = z[-1] - z[0]
        tower_buckling = bucklingGL(twr_D, twr_t, Fz, My, tower_height, E, fy, gamma_f, gamma_b, gamma_g)

        self.utilization.GLUtil = tower_buckling  #utilization


        #import matplotlib.pyplot as plt
        #plt.plot(self.utilization.EUshUtil, z)
        #plt.plot(self.utilization.GLUtil, z)
        #plt.show()


#____________________________________________________#

class PrepArray(Component):
    "This stupid component to rearrange array as I could not slice them in jacketassembly"
    RNA_FM = Array(np.zeros([3,2]),iotype='in', desc='rotor aerodynamic forces & moments in hub-aligned,yawed coordinate system')
    RNA_F = Array(np.zeros([3]),iotype='out', desc='rotor aerodynamic forces in hub-aligned,yawed coordinate system')
    RNA_M = Array(np.zeros([3]),iotype='out', desc='rotor aerodynamic moments in hub-aligned,yawed coordinate system')
    def execute(self):
        self.RNA_F=self.RNA_FM[:,0]
        self.RNA_M=self.RNA_FM[:,1]

class UtilAssembly(Assembly):

    # inputs
    #RNA_FM = Array(np.zeros([3,2]),iotype='in', desc='rotor aerodynamic forces & moments in hub-aligned,yawed coordinate system')
    RNA_F = Array(np.zeros([3]),iotype='in', desc='rotor aerodynamic forces in hub-aligned,yawed coordinate system')
    RNA_M = Array(np.zeros([3]),iotype='in', desc='rotor aerodynamic moments in hub-aligned,yawed coordinate system')
    rna_weightM = Bool(True, iotype='in', desc='flag to consider or not the RNA weight effect on Moment')

    r_hub = Array(np.zeros([3]),iotype='in', desc='position of rotor hub relative to tower top in yaw-aligned c.s.')
    r_cm  = Array(np.zeros([3]),iotype='in', desc='position of RNA CM relative to tower top in yaw-aligned c.s.')

    tilt = Float(iotype='in', units='deg')
    g = Float(9.81, iotype='in', units='m/s**2', desc='Gravity Acceleration (ABSOLUTE VALUE!)')

    towerWindLoads = VarTree(AeroLoads(), iotype='in', desc='Aero loads in inertial coordinate system')
    towerWaveLoads = VarTree(AeroLoads(), iotype='in', desc='Hydro loads in inertial coordinate system')
    Dt = Float(iotype='in', units='m', desc='TowerTop OD from Tower')
    tt = Float(iotype='in', units='m', desc='TowerTop wall thickness from Tower')
    Twr_data = VarTree(TwrGeoOutputs(), iotype='in', desc='Tower Node data')  # From Tower
    L_reinforced = Float(30.0, iotype='in', desc='reinforcement length')
    IECpsfIns = VarTree(IEC_PSFS(), iotype='in', desc='Basic IEC psfs')


    # outputs
    tower_utilization = VarTree(TowerUtilOutputs(), iotype='out', desc='Tower Utilization Basic Outputs')
    jacket_utilization = VarTree(JacketUtilOutputs(), iotype='out', desc='Jacket Utilization Basic Outputs')


    def configure(self):

        self.add('PrepArray',PrepArray())
        self.add('rotorLoads', RotorLoads())
        self.add('TwrUtil', TwrUtilization())
        self.add('JcktUtil', JcktUtilization())

        #self.driver.workflow.add(['PrepArray', 'rotorLoads', 'JcktUtil', 'TwrUtil'])
        self.driver.workflow.add(['rotorLoads', 'JcktUtil', 'TwrUtil'])

        #connections to PrepArray
        #self.connect('RNA_FM', 'PrepArray.RNA_FM')

        # connections to rotorLoads
        self.connect('RNA_F', 'rotorLoads.F')
        self.connect('RNA_M', 'rotorLoads.M')
        self.connect('rna_weightM','rotorLoads.rna_weightM')

        self.connect('r_hub', 'rotorLoads.r_hub')
        self.connect('r_cm', 'rotorLoads.rna_cm')
        self.connect('Twr_data.TopMass[0]', 'rotorLoads.m_RNA')
        self.connect('tilt', 'rotorLoads.tilt')
        self.connect('g', 'rotorLoads.g')

        # connections to TwrUtil
        self.connect('rotorLoads.top_F', 'TwrUtil.top_F')
        self.connect('rotorLoads.top_M', 'TwrUtil.top_M')
        self.connect('towerWindLoads', 'TwrUtil.towerWindLoads')
        self.connect('towerWaveLoads', 'TwrUtil.towerWaveLoads')
        self.connect('Dt', 'TwrUtil.Dt')
        self.connect('tt', 'TwrUtil.tt')
        self.connect('Twr_data', 'TwrUtil.Twr_data')
        self.connect('L_reinforced', 'TwrUtil.L_reinforced')
        self.connect('IECpsfIns', 'TwrUtil.IECpsfIns')
        self.connect('g', 'TwrUtil.g')


        #Connections to inputs of JcktUtilization
        self.create_passthrough('JcktUtil.Legmems')
        self.create_passthrough('JcktUtil.Pilemems')
        self.create_passthrough('JcktUtil.Twrmems')
        self.create_passthrough('JcktUtil.TPmems')
        self.create_passthrough('JcktUtil.MbrFrcs')
        self.create_passthrough('JcktUtil.nlegs')
        self.create_passthrough('JcktUtil.XjntIDs')
        self.create_passthrough('JcktUtil.KjntIDs')
        self.create_passthrough('JcktUtil.JcktGeoOut')



        #Connections to outputs
        self.connect('TwrUtil.utilization', 'tower_utilization')
        self.connect('JcktUtil.Utilouts', 'jacket_utilization')


