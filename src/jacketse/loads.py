import numpy as np
import math
from openmdao.main.api import VariableTree, Component, Assembly
from openmdao.main.datatypes.api import Int, Float, Array, Bool, VarTree
import sys

from commonse.WindWaveDrag import FluidLoads, TowerWindDrag, TowerWaveDrag
from commonse.environment import PowerWind, LinearWaves
from commonse.utilities import sind, cosd

from VarTrees import RNAprops

#NOTE: ON 5/9/2014 the input of wave height changed to be the full wave height (peak-to-peak) as Andrew's load routine uses it, if i use ww_loads.py, then this would have to be 0.5(peak-to-peak)

# TODO: share a common set of VariableTrees

class WaterInputs(VariableTree):

    """Basic Wave Inputs needed to calculate Jacket Loads"""
    # inputs
    wdepth = Float(units='m', desc='Water Depth')
    z_floor = Float(units='m', desc='reference location of sea floor')
    wlevel = Float(units='m', desc='Water Level: Distance from bottom of structure (including buried pile) to water surface')
    T = Float(units='s', desc='50-yr wave period (single wave here, not peak spectral period)')
    HW = Float(units='m', desc='50-yr wave height PEAK-to-PEAK!!!  (it used to be half-amplitude of sinusoid, changed since)')
    psi = Float(45., units='deg', desc='Wave direction angle relative to global CS.')
    rho = Float(1027., units='kg*m**-3', desc='Water Density')
    Cm = Float(2.0, units=None, desc='Added mass coefficient')
    Cd =Float(iotype='in', units=None, desc='User input drag coefficient, if left blank it will be automatically calculated based on Cylinder Re')
    mu = Float(1.3351e-3, units='kg/(m*s)', desc='dynamic viscosity of water')
    Uc = Float(units='m/s', desc='mean current speed')


class WindInputs(VariableTree):

    """Basic Wind Inputs needed to calculate Jacket Loads"""
    # inputs
    al_shear = Float(0.14, units=None, desc='power law exponent wind from IEC 61400-3')
    psi = Float(45.,units='deg', desc='Wind angle relative global CS. Positive according to RHR with positive z direction.')
    rho = Float(1.225, units='kg*m**-3', desc='Air Density. Set to 0 if no tower wind load is included (thrust from rotor still on)')
    U50HH = Float(70., units='m/s', desc='50-yr return 3-s gust [m/s] :From IEC Class I-III or whatever gust to be associated with DLC under consideration (e.g. 30 m/s for DLC 1.6 under max thrust)')
    HH = Float(units='m', desc='Hub-height above MSL')
    mu = Float(1.7934e-5, units='kg/(m*s)', desc='dynamic viscosity of air')
    Cdj =Float(iotype='in', units=None, desc='User input drag coefficient for jacket members, if left blank it will be automatically calculated based on Cylinder Re')
    Cdt =Float(iotype='in', units=None, desc='User input drag coefficient for tower, if left blank it will be automatically calculated based on Cylinder Re')

class LoadOutputs(VariableTree):

    """Jacket Loading Basic Outputs"""
    wDrag1_3 = Array(dtype=np.float, units='N', desc='Drag forces on leg 1 and 3')
    wDrag2_4 = Array(dtype=np.float, units='N', desc='Drag forces on leg 2 and 4')
    #wl_idx=   Array(dtype=int,      units=None, desc='Indices of the leg/pile nodes where water loading is applied')
    twr_dragf = Array(dtype=np.float, units='N', desc='Drag forces on tower nodes')
    nds_frc = Array(dtype=np.float,
                    desc='Concentrated Forces and Moments: (NODE ID, Fx,Fy,Fz,Mxx,Myy,Mzz)')
    RNAload =Array(np.zeros([3,2]),dtype=np.float, desc='RNA aero-loads rotated into yawed coordinate system, [3,2], 1st col. forces, 2nd col. moments')
# ______________________________________________________________________________#


class JcktLoadPre(Component):

    """This function sets the wind and wave loading onto the jacket/tower structure \n
    for FRAME3DD to work properly. It comes from JcktLoad under PYTHON (nonGIT).\n
    NOTE: 1. It assumes a DLC 1.6 case load at 45 degree, along 4-legged jacket diagonal. \n
          2. Still to do other load conditions and 3-legged case.\n
          3. Buoyancy effect is considered negligible (flooded legs) as well as hydrostatic pressure. \n
    \n
    INPUT & OUTPUT described below in Component Core.
    """
    # inputs
    nlegs = Int(iotype='in')  # comes from JcktGeoIn.nlegs
    nodes = Array(iotype='in')  # comes from JcktGeoOut.nodes

    pillegDs = Array(iotype='in', units='m', desc='ODs of pile and leg #1.')
    twrDs = Array(iotype='in', units='m', desc='ODs of Tower.')
    #twrTs = Array(iotype='in', units='m', desc='shell thickness of Tower.')

    Twrmems = Array(iotype='in', dtype=int, desc='Connectivity Array for all members of the tower portion only (including rigid member at the top if requested)')
    Legmems = Array(iotype='in', dtype=int, desc='Connectivity Array for all members of the Leg 1 portion only')
    Pilemems = Array(iotype='in', dtype=int, desc='Connectivity Array for all members of the Pile 1 portion only ')

    # outputs
    pillegZs_out = Array(iotype='out', units='m', desc='pile and leg z node locations')
    pillegDs_out = Array(iotype='out', units='m', desc='corresponding diameters')
    twrZs_out = Array(iotype='out', units='m', desc='tower z node locations')
    twrDs_out = Array(iotype='out', units='m', desc='corresponding tower diameters')
    #twrTs_out = Array(iotype='out', units='m', desc='corresponding tower shell thickness')

    legndIDs= Array(iotype='out', units=None, desc='Node IDs for the legs')
    twrndIDs= Array(iotype='out', units=None, desc='Node IDs for the tower')
    pilendIDs= Array(iotype='out', units=None, desc='Node IDs for the piles')

    def execute(self):
        # Simplify nomenclature
        nlegs = self.nlegs
        nodes = self.nodes
        Pilemems = self.Pilemems
        Legmems = self.Legmems
        Twrmems = self.Twrmems

        # indices of leg 1
        legidx = np.arange(Legmems[0, 0], Legmems[-1, 1] / nlegs + 1) - 1

        #Indices in JcktGeoOut.nodes
        pileidx = np.array([])  # Initialize in case empty piles
        pileZs = np.array([])
        pileDs = np.array([])
        self.pilendIDs=np.array([]).reshape([-1,1])
        if Pilemems.size:
            # indices of pile 1 (not counting leg joint)
            pileidx = np.arange(Pilemems[0, 0], Pilemems[Pilemems.shape[0]/nlegs-1, 1] ) - 1
            self.pilendIDs=np.arange(Pilemems[0,0],Pilemems[-1,1]+1,legidx.size+1).reshape([-1,1])  #actual node IDs for all pile nodes to be used later
            pileZs = nodes[pileidx, 2]
            pileDs = self.pillegDs[pileidx]

        # legs and Tower must exist

        #actual node IDs of all legs
        self.legndIDs=np.arange(Legmems[0,0],Legmems[-1,1]+1).reshape([-1,1])

        # Indices of tower nodes
        twridx = np.arange(Twrmems[0, 0], Twrmems[-1, 1] + 1) - 1

        self.twrndIDs=(twridx+1).reshape([-1,1]) #Actual IDs of tower nodes, possibly including RNA node in RigidTop=True;np.arange(Twrmems[0,0],Twrmems[-1,1]+1).reshape([-1,1])

        # Zs and OD of pile, leg, tower
        legZs = nodes[legidx, 2]
        legDs = np.hstack((self.pillegDs[legidx[0]:], self.pillegDs[-1]))
                          # ODs are associated with members not nodes, so I
                          # need to replicate for the last (top) node
        self.twrZs_out = nodes[twridx, 2]
        self.twrDs_out = np.hstack(
            (self.twrDs, self.twrDs[-1]))  # ODs are associated with flexible members not nodes, so I need to replicate for the last (top) node;
                                                     # this should be
                                                     # conservative as Dt may
                                                     # be a little less than
                                                     # D[-1]
        #self.twrTs_out = np.hstack((self.twrTs, self.twrTs[-1]))

        # this may include a fictitious z for the attachment to RNA (CMzoff)
        if self.twrZs_out.size > self.twrDs_out.size:
            self.twrZs_out = self.twrZs_out[:-1]  # pop the last element

        if nlegs == 4:  # Diagonal loading condition ONLY for the time being

            self.pillegZs_out = np.hstack((pileZs, legZs))
            self.pillegDs_out = np.hstack((pileDs, legDs))

        else:  # 3legged jacket in progress
            # TODO: THIS CASE
            sys.exit("nlegs <>4 not implemented yet for loading")



# from scipy.integrate import cumtrapz


# class DistributedToPointLoads(Component):

#     FluidLoads = VarTree(FluidLoads(), iotype='in', desc='aero/hydro loads in inertial coordinate system')

#     Tx = Array(iotype='out')
#     Ty = Array(iotype='out')
#     Mxx = Array(iotype='out')
#     Myy = Array(iotype='out')
#     N = Array(iotype='out')


#     def execute(self):

#         loads = self.windLoads

#         adrag = np.sqrt(loads.Px**2 + loads.Py**2)

# z_from_top = loads.z[-1] - loads.z[::-1]  # distance from top, in order
# from top to bottom

# Tx_adrag = cumtrapz(adrag[::-1], z_from_top)[::-1]  # Shear due to wind pressure, from bottom to top; cumtrapz cuts dimension
# Myy_adrag = cumtrapz(adrag[::-1]*z_from_top, z_from_top)[::-1]  #
# Bending due to wind pressure, bottom to top

#         self.Tx = Tx_adrag*np.cos(loads.beta)
#         self.Ty = Tx_adrag*np.sin(loads.beta)
#         self.Myy = Myy_adrag*np.cos(loads.beta)
#         self.Mxx = Myy_adrag*np.sin(loads.beta)
#         self.N = Tx_adrag


class JcktLoadPost(Component):

    # inputs
    towerWindLoads   = VarTree(FluidLoads(), iotype='in', desc='aero loads in inertial coordinate system')
    towerWaveLoads   = VarTree(FluidLoads(), iotype='in', desc='aero loads in inertial coordinate system')
    pileLegWindLoads = VarTree(FluidLoads(), iotype='in', desc='aero loads in inertial coordinate system')
    pileLegWaveLoads = VarTree(FluidLoads(), iotype='in', desc='aero loads in inertial coordinate system')
    nlegs            = Int(units=None, iotype='in', desc='Number of Jacket Legs')  # comes from JcktGeoIn.nlegs
    al_bat3D         = Float(   units='rad',iotype='in', desc='Batter Angle in 3D.')
    VPFlag=      Bool(False,units=None,    iotype='in',desc='Vertical Pile Flag [Y/N]: If True the Mudbrace is put at the expected joint with Piles, i.e. at the bottom of leg.')

    RNA_F       =Array(dtype=np.float,  iotype='in',desc='Unfactored Rotor Forces and Moments excluding weight, aligned with rotor CS. Array(6)')

    RNAinputs=VarTree(RNAprops(), iotype='in', desc='Basic Inertial Properties of RNA')

    TwrRigidTop =Bool( units=None,iotype='in',desc='Rigid Member used in tower top or not')
    #CMzoff=Float(      units='m', iotype='in',desc='RNA CM z offset from Tower Top Flange as It comes out of Tower Processing (trumped by RigidTop incase)')       # RNA CMzoff [m]

    legndIDs= Array(iotype='in', units=None, desc='Node IDs for the legs')
    twrndIDs= Array(iotype='in', units=None, desc='Node IDs for the tower')
    pilendIDs= Array(iotype='in', units=None, desc='Node IDs for the piles')

    wdepth   =Float(iotype='in', units='m', desc='Water Depth, needed to refine deltaz at the node just below surface')
    # outputs
    Loadouts = VarTree(LoadOutputs(), iotype='out', desc='Jacket Loading Basic Outputs')

    def execute(self):
        #simplify nomenclature
        psi_wi=self.towerWindLoads.beta[0]  #psi is stored in beta already, however for legs and piles it needs to be adjusted for psi
        psi_wa=self.towerWaveLoads.beta[0]
        twrZs=self.towerWindLoads.z
        pillegZs=self.pileLegWaveLoads.z
        pillegDs=self.pileLegWaveLoads.d
        nlegs=self.nlegs
        al_bat3D=self.al_bat3D
        VPFlag=self.VPFlag

        TwrRigidTop=(self.RNAinputs.CMoff[2] != 0.) and self.TwrRigidTop  #This makes sure we do not have TwrRigidTop=True with CMzoff=0., no length segment that is.

        pilendIDs=self.pilendIDs
        legndIDs=self.legndIDs
        twrndIDs=self.twrndIDs

        wdepth=self.wdepth
       #COSINE MATRIX FROM LOCAL TO GLOBAL COORDINATE SYSTEMS
        DIRCOSwind=np.array([ [cosd(psi_wi) , -sind(psi_wi) ,0.],
                          [    sind(psi_wi) , cosd(psi_wi),0.],
                          [    0.,                         0.,                      1.]])
        DIRCOSwave=np.array([ [cosd(psi_wa) , -sind(psi_wa) ,0.],
                          [   sind(psi_wa) ,   cosd(psi_wa),0.],
                          [     0.,                         0.,                      1.]])

        #I need to get the right node IDs to be able to assign forces for Frame3DD

        #Get the values of the loads at the nodes of pile and leg 1
        #for the other legs, the wdrag is going to be the same

        if nlegs== 4: #Diagonal loading condition ONLY for the time being
            pileidx=np.array([])  #Initialize in case empty piles
            if pilendIDs.size:
                pileidx=pilendIDs[0:pilendIDs.size/nlegs] #np.arange(Pilemems[0,0],Pilemems[-1,1]/nlegs+1)-1   #indices of pile 1

            legidx=legndIDs[0:legndIDs.size/nlegs]    #legidx=np.arange(Legmems[0,0],Legmems[-1,1]/nlegs+1)-1 #indices of leg 1

            deltaz=((np.roll(pillegZs,-1)-pillegZs)/2)[0:-1]  #these are the appropriate 1/2 deltazs for the entire pile-leg assembly
            deltaz=np.hstack((deltaz[0],(np.roll(deltaz,-1)+deltaz)[0:-1],deltaz[-1]))  #This is the actual DeltaZ to be assigned at each node starting from the 2nd and ending at the one before last

            #Correct the delta z for the node just below water to avoid jumps in loading - use refined mesh however
            idx_bw=np.nonzero(pillegZs<= wdepth)[0][-1]
            #Also Attempt at using load at z=0
            ###Px0=self.pileLegWaveLoads.Px_i0overd2*pillegDs[idx_bw]**2+self.pileLegWaveLoads.Px_d0overd*pillegDs[idx_bw]
            ###Py0=self.pileLegWaveLoads.Py_i0overd2*pillegDs[idx_bw]**2+self.pileLegWaveLoads.Py_d0overd*pillegDs[idx_bw]
            Px0=self.pileLegWaveLoads.Px0
            Py0=self.pileLegWaveLoads.Py0
            deltaz0=(wdepth-pillegZs[idx_bw])

            if (wdepth-pillegZs[idx_bw])> deltaz[idx_bw]/2.:  #point with its deltaz entriely below surface
                deltaz0 -= deltaz[idx_bw]/2.  #note deltaz0 before deltaz as deltaz gets modified
                #deltaz[idx_bw]=deltaz[idx_bw]/2. -(pillegZs[idx_bw]-wdepth)
            else:
                deltaz[idx_bw]=deltaz[idx_bw]/2.

            waDrag0=np.sqrt(Px0**2.+Py0**2.)*deltaz0  #In case add this one to the original wDrag, or alternatively do awhat I do below, increasing Deltaz

            waDrag=np.sqrt(self.pileLegWaveLoads.Px**2.+self.pileLegWaveLoads.Py**2.)*deltaz #[N] forces from waves, considered normal to the sloping 1 leg.
            waDrag[idx_bw] += waDrag0

            wiDrag=np.sqrt(self.pileLegWindLoads.Px**2.+self.pileLegWindLoads.Py**2.)*deltaz * (np.cos(al_bat3D))**2 #[N] forces from wind normal to the sloping 1 leg. Wind is horizontal, that's why cos^2

            #Forces on legs 1 (and 3 though sign of Fz reversed)
            junk=np.zeros(waDrag.size)
            waDrag1_3=  (np.dot(DIRCOSwave , np.vstack((waDrag*np.cos(al_bat3D),junk,waDrag*np.sin(al_bat3D))))).T   #[n,3] [N] This is an approx for air drag, as it is not normal to the leg to begin with, but it makes it easier
            wiDrag1_3=  (np.dot(DIRCOSwind , np.vstack((wiDrag*np.cos(al_bat3D),junk,wiDrag*np.sin(al_bat3D))))).T   #[n,3] [N]

            #Forces on legs 2 (and 4), it is as if they were vertical
            waDrag2_4=  (np.dot(DIRCOSwave ,np.vstack((waDrag,junk, waDrag*0.)))).T #[n,3]  [N]
            wiDrag2_4=  (np.dot(DIRCOSwind ,np.vstack((wiDrag,junk, wiDrag*0.)))).T #[n,3]  [N]

            #Add them together
            wDrag1_3=waDrag1_3+wiDrag1_3 #[N]
            wDrag2_4=waDrag2_4+wiDrag2_4 #[N]
            wDrag=waDrag+wiDrag #[N]

            wl_idx=np.nonzero(wDrag)[0] #indices of the pile-leg nodes where loading is applied--

            if VPFlag and pileidx:  #Need to account for the cos
               #Forces on legs 1 (and 3 though sign of Fz reversed)
                wDrag1_3[pileidx]=wDrag2_4[pileidx]

            n_pillegLds=wl_idx.size
            n_loads=n_pillegLds*nlegs #Total number of wave load nodes

            pilleg_ndsfrc=np.zeros([n_loads,7])
            nNodespile=pileidx.size  #number of nodes in 1 pile
            nNodesleg=legidx.size  #number of nodes in 1 leg

            for ii in range(0,nlegs): #Cycle through legs
                idx0=nNodespile*ii  #start index for pile IDs
                idx1=idx0+nNodespile #stop index for pile IDs
                idx2=nNodesleg*ii   #start index for leg IDs
                idx3=idx2+nNodesleg  #stop index for leg IDs
                nd_ids=np.vstack((self.pilendIDs[idx0:idx1],self.legndIDs[idx2:idx3]))[wl_idx]
                if np.mod(ii,2)==0:
                    lds=wDrag2_4[wl_idx,:]
                elif ii==1:
                    lds=np.hstack((wDrag1_3[wl_idx,0:2],-wDrag1_3[wl_idx,2].reshape(-1,1)))
                else:
                    lds=wDrag1_3[wl_idx,:]

                idx0=n_pillegLds*ii #start index in pilleg_ndsfrc
                idx1=idx0+n_pillegLds  #stop index in pilleg_ndsfrc
                pilleg_ndsfrc[idx0:idx1,:-3]=np.hstack((nd_ids,lds))

                #Store wave loads
                self.Loadouts.wdrag1_3=wDrag1_3
                self.Loadouts.wdrag2_4=wDrag2_4

        else: #3legged jacket in progress
            #TO DO THIS CASE
            sys.exit("nlegs <>4 not implemented yet for loading")

        #-----Need to add hydrostatic loading -BUOYANCY FORCE
        #IT is a distributed force along the non-flooded members
        #Still TO DO

        #________________________________TOWER__________________________________#
        #Now get the values of the loads at the nodes of Tower

        junk=(np.roll(twrZs,-1)-np.roll(twrZs,1))[1:-1]/2.
        deltaz=np.hstack(((twrZs[1]-twrZs[0])/2.,junk,(twrZs[-1]-twrZs[-2])/2.))#these are the appropriate deltazs

        #just wind loading for tower but wtrload perhaps in the future for tripods
        junk=(np.sqrt(self.towerWindLoads.Px**2+self.towerWindLoads.Py**2)*deltaz).reshape([-1,1]) #[N] forces

        #decompose along local x,y and add z component
        junk=np.hstack((junk,np.zeros([junk.size,2]))).T
        self.Loadouts.twr_dragf=np.dot(DIRCOSwind , junk)    #*np.sin(wind_dict['psi']), junk*0.])

        #Tower Loads for Frame3DD
        n_twrlds=(twrndIDs.shape[0]-TwrRigidTop)*(self.Loadouts.twr_dragf.any()) #distributed loads along tower+RNA node load
        twr_ndsfrc=np.zeros([n_twrlds,7])
        twr_ndsfrc[0:n_twrlds,:-3]=np.hstack((twrndIDs[0:len(twrndIDs)-TwrRigidTop].reshape([-1,1]),self.Loadouts.twr_dragf.T))

        #____________ADD CONCENTRATED RNA FORCES________________________________#
        #First account for moment transfer in global coordinate system as if yaw=0
        RNAload=np.copy(self.RNA_F.reshape([2,3]).T)#RNAload: 1columns forces, 2nd column moments

        #Store yaw-rotated forces and moments at the hub still  [3,2]
        self.Loadouts.RNAload=np.dot(DIRCOSwind,RNAload)

        Deltax=self.RNAinputs.Thoff[0]-self.RNAinputs.CMoff[0]
        Deltay=self.RNAinputs.Thoff[1]-self.RNAinputs.CMoff[1]
        Deltaz=self.RNAinputs.Thoff[2]-self.RNAinputs.CMoff[2]

        if  not(TwrRigidTop):
            Deltax=self.RNAinputs.Thoff[0]
            Deltay=self.RNAinputs.Thoff[1]
            Deltaz=self.RNAinputs.Thoff[2]

        #Rotor Loads - no weight yet
        RNAload[0,1] += -RNAload[1,0]*Deltaz + RNAload[2,0]*Deltay #Mxx=-Fy*Deltaz +Fz*deltay
        RNAload[1,1] +=  RNAload[0,0]*Deltaz - RNAload[2,0]*Deltax  #Myy= Fx*deltaz-Fz*deltax
        RNAload[2,1] += -RNAload[0,0]*Deltay + RNAload[1,0]*Deltax #Mzz=-Fx*Deltay +Fy*deltax

        #Then rotate them in the global coordinate system for PYFRAME to account for yaw
        RNAload=np.dot(DIRCOSwind,RNAload) #Gravity is already accounted for by Frame3DD for concentrated masses: Jacket.Twr.RNAmass*aux['g_z'],np.zeros(3)))  self.RNA_F.reshape([2,3]).T

        if not(TwrRigidTop):
            twr_ndsfrc[-1,1:] += RNAload.T.flatten()
        else:
            twr_ndsfrc = np.vstack(( twr_ndsfrc,np.hstack((twrndIDs[-1],RNAload.T.flatten())) ))

        # If we apply this load at tower-top and not tower-top+CMzoff we need to account for different DeltaZ
##        if  not(TwrRigidTop):
##            RNAload[0,1] -= RNAload[1,0]*ThOffset_z  #*CMzoff #Mxx
##            RNAload[1,1] += RNAload[0,0]*CMzoff #Myy
##            twr_ndsfrc[-1,1:] += RNAload.T.flatten()
##        else: #Put forces at the RNA.CM location
##            RNAload[0,1] -= RNAload[1,0]*ThOffset_z  #*CMzoff #Mxx
##            RNAload[1,1] += RNAload[0,0]*CMzoff #Myy
##
##            twr_ndsfrc[-1,:] += np.hstack((twrndIDs[-1],RNAload.T.flatten()))

        #STACK ALLOF THEM TOGETHER
        self.Loadouts.nds_frc=np.vstack((pilleg_ndsfrc, twr_ndsfrc))  #Concentrated loads



        # import matplotlib.pyplot as plt
        # plt.plot(self.towerWindLoads.Px, self.towerWindLoads.z)
        # plt.plot(self.towerWaveLoads.Px, self.towerWaveLoads.z)
        # plt.plot(self.pileLegWindLoads.Px, self.pileLegWindLoads.z)
        # plt.plot(self.pileLegWaveLoads.Px, self.pileLegWaveLoads.z)
        # plt.show()


class JcktLoad(Assembly):

    # normally these connections come from other components
    nlegs = Int(iotype='in')  # comes from JcktGeoIn.nlegs
    nodes = Array(iotype='in')  # comes from JcktGeoOut.nodes
    pillegDs = Array(dtype=np.float, units='m', iotype='in', desc='ODs of pile and leg #1.')
    twrDs = Array(units='m', iotype='in', desc='ODs of Tower.')
    #twrTs = Array(iotype='in')
    Twrmems = Array(dtype=int, iotype='in', desc='Connectivity Array for all members of the tower portion only (including rigid member at the top if requested)')
    Legmems = Array(dtype=int, iotype='in', desc='Connectivity Array for all members of the Leg 1 portion only')
    Pilemems = Array(dtype=int, iotype='in', desc='Connectivity Array for all members of the Pile 1 portion only ')

    VPFlag  = Bool(False,units=None, iotype='in',desc='Vertical Pile Flag [Y/N]: If True the Mudbrace is put at the expected joint with Piles, i.e. at the bottom of leg.')
    al_bat3D=Float(      units='rad',iotype='in',desc='Batter Angle in 3D.')
    RNA_F       =Array(dtype=np.float,  iotype='in',desc='Unfactored Rotor Forces and Moments excluding weight. Array(6).')
    TwrRigidTop =Bool( units=None,iotype='in',desc='Rigid Member used in tower top or not.')

    RNAinputs=VarTree(RNAprops(), iotype='in', desc='Basic Inertial Properties of RNA')
    #CMzoff  =Float(      units='m', iotype='in',desc='RNA CM z offset from Tower Top Flange as It comes out of Tower Processing (trumped by RigidTop incase)')       # RNA CMzoff [m]
    gravity =Float(  units='m/s**2',iotype='in',desc='Gravity Acceleration (ABSOLUTE VALUE!)')

    windIns = VarTree(WindInputs(), iotype='in')
    waterIns = VarTree(WaterInputs(), iotype='in')

    # #Following inputs are for utilization, I would keep this separate
    # topF = Array(iotype='in')
    # topM = Array(iotype='in')
    # n_reinforced = Int(3, iotype='in')
    # E = Float(210e9, iotype='in', units='N/m**2', desc='material modulus of elasticity')
    # sigma_y = Float(450.0e6, iotype='in', units='N/m**2', desc='yield stress')
    # gamma_f = Float(1.35, iotype='in', desc='safety factor on loads')
    # gamma_m = Float(1.1, iotype='in', desc='safety factor on materials')
    # gamma_n = Float(1.0, iotype='in', desc='safety factor on consequence of failure')


    #inputs
#    JcktGeoIn=   VarTree(JcktGeoInputs(),  iotype='in', desc='Jacket Geometry Basic Inputs')
#    JcktGeoOut = VarTree(JcktGeoOutputs(), iotype='in', desc='Geometry of the Jacket -Node Coordinates and Member Connectivity')

#    gravity     =Float(                units='m/s**2',iotype='in',desc='Acceleration Gravity.')

    #outputs
    Loadouts=   VarTree(LoadOutputs(), iotype='out', desc='Jacket Loading Basic Outputs')
    twrWindLoads = VarTree(FluidLoads(), iotype='out', desc='Tower wind loads in inertial coordinate system')
    twrWaveLoads = VarTree(FluidLoads(), iotype='out', desc='Tower wave loads in inertial coordinate system')
    #stress = Array(iotype='out', units='N/m**2', desc='von Mises stress along tower on downwind side (yaw-aligned +x).  normalized by yield stress.  includes safety factors.')
    #z_buckling = Array(iotype='out', units='m', desc='z-locations along tower where shell buckling is evaluted')
    #buckling = Array(iotype='out', desc='a shell buckling constraint.  should be <= 0 for feasibility.  includes safety factors')

    def configure(self):

        self.add('pre', JcktLoadPre())
        self.add('windj', PowerWind())
        self.add('windt', PowerWind())
        self.add('wavej', LinearWaves())
        self.add('wavet', LinearWaves())
        self.add('windLoadsj', TowerWindDrag())
        self.add('windLoadst', TowerWindDrag())
        self.add('waveLoadsj', TowerWaveDrag())
        self.add('waveLoadst', TowerWaveDrag())
        self.add('post', JcktLoadPost())


        self.driver.workflow.add(['pre', 'windj', 'windt', 'wavej', 'wavet',
            'windLoadsj', 'windLoadst', 'waveLoadsj', 'waveLoadst', 'post'])

        # connections to pre
        self.connect('nlegs', 'pre.nlegs')
        self.connect('nodes', 'pre.nodes')
        self.connect('pillegDs', 'pre.pillegDs')
        self.connect('twrDs', 'pre.twrDs')
        #self.connect('twrTs', 'pre.twrTs')
        self.connect('Twrmems', 'pre.Twrmems')
        self.connect('Legmems', 'pre.Legmems')
        self.connect('Pilemems', 'pre.Pilemems')

        # connections to windj/t
        self.connect('windIns.U50HH', ['windj.Uref', 'windt.Uref'])
        self.connect('windIns.HH + waterIns.wdepth + waterIns.z_floor', ['windj.zref', 'windt.zref'])
        self.connect('pre.pillegZs_out', 'windj.z')
        self.connect('pre.twrZs_out', 'windt.z')
        self.connect('waterIns.wdepth + waterIns.z_floor', ['windj.z0', 'windt.z0'])
        self.connect('windIns.psi', ['windj.betaWind', 'windt.betaWind'])
        self.connect('windIns.al_shear', ['windj.shearExp', 'windt.shearExp'])


        # connections to wavej/t
        self.connect('waterIns.Uc', ['wavej.Uc', 'wavet.Uc'])
        self.connect('waterIns.wdepth + waterIns.z_floor', ['wavej.z_surface', 'wavet.z_surface'])
        self.connect('waterIns.HW', ['wavej.hmax', 'wavet.hmax'])
        # self.connect('waterIns.T', ['wavej.T', 'wavet.T'])
        self.connect('waterIns.T', 'wavej.T')
        self.connect('waterIns.T', 'wavet.T')
        self.connect('waterIns.z_floor', ['wavej.z_floor', 'wavet.z_floor'])
        self.connect('gravity', ['wavej.g', 'wavet.g'])
        self.connect('waterIns.psi', ['wavej.betaWave', 'wavet.betaWave'])
        self.connect('pre.pillegZs_out', 'wavej.z')
        self.connect('pre.twrZs_out', 'wavet.z')

        # connections to windLoadsj
        self.connect('windj.U', 'windLoadsj.U')
        self.connect('windj.beta', 'windLoadsj.beta')
        self.connect('windIns.rho', 'windLoadsj.rho')
        self.connect('windIns.mu', 'windLoadsj.mu')
        self.connect('windIns.Cdj', 'windLoadsj.cd_usr')
        self.connect('pre.pillegZs_out', 'windLoadsj.z')
        self.connect('pre.pillegDs_out', 'windLoadsj.d')

        # connections to windLoadst
        self.connect('windt.U', 'windLoadst.U')
        self.connect('windt.beta', 'windLoadst.beta')
        self.connect('windIns.rho', 'windLoadst.rho')
        self.connect('windIns.mu', 'windLoadst.mu')
        self.connect('windIns.Cdt', 'windLoadst.cd_usr')
        self.connect('pre.twrZs_out', 'windLoadst.z')
        self.connect('pre.twrDs_out', 'windLoadst.d')

        # connections to waveLoadsj
        self.connect('wavej.U', 'waveLoadsj.U')
        self.connect('wavej.A', 'waveLoadsj.A')
        self.connect('wavej.beta', 'waveLoadsj.beta')
        self.connect('wavej.U0', 'waveLoadsj.U0')
        self.connect('wavej.A0', 'waveLoadsj.A0')
        self.connect('wavej.beta0', 'waveLoadsj.beta0')

        self.connect('waterIns.rho', 'waveLoadsj.rho')
        self.connect('waterIns.mu', 'waveLoadsj.mu')
        self.connect('waterIns.Cm', 'waveLoadsj.cm')
        self.connect('waterIns.Cd', 'waveLoadsj.cd_usr')
        self.connect('waterIns.wlevel', 'waveLoadsj.wlevel')
        self.connect('pre.pillegZs_out', 'waveLoadsj.z')
        self.connect('pre.pillegDs_out', 'waveLoadsj.d')

        # connections to waveLoadst
        self.connect('wavet.U', 'waveLoadst.U')
        self.connect('wavet.A', 'waveLoadst.A')
        self.connect('wavet.beta', 'waveLoadst.beta')
        self.connect('wavet.U0', 'waveLoadst.U0')
        self.connect('wavet.A0', 'waveLoadst.A0')
        self.connect('wavet.beta0', 'waveLoadst.beta0')
        self.connect('waterIns.rho', 'waveLoadst.rho')
        self.connect('waterIns.mu', 'waveLoadst.mu')
        self.connect('waterIns.Cm', 'waveLoadst.cm')
        self.connect('waterIns.Cd', 'waveLoadst.cd_usr')
        self.connect('waterIns.wlevel', 'waveLoadst.wlevel')
        self.connect('pre.twrZs_out', 'waveLoadst.z')
        self.connect('pre.twrDs_out', 'waveLoadst.d')

        # connections to post
        self.connect('windLoadst.windLoads', 'post.towerWindLoads')
        self.connect('waveLoadst.waveLoads', 'post.towerWaveLoads')
        self.connect('windLoadsj.windLoads', 'post.pileLegWindLoads')
        self.connect('waveLoadsj.waveLoads', 'post.pileLegWaveLoads')
        self.connect('pre.pilendIDs','post.pilendIDs')
        self.connect('pre.legndIDs','post.legndIDs')
        self.connect('pre.twrndIDs','post.twrndIDs')
        self.connect('nlegs','post.nlegs')
        self.connect('TwrRigidTop','post.TwrRigidTop')
        self.connect('RNAinputs','post.RNAinputs')
        self.connect('RNA_F','post.RNA_F')
        self.connect('al_bat3D','post.al_bat3D')
        self.connect('VPFlag','post.VPFlag')
        self.connect('waterIns.wdepth','post.wdepth')

        # connections to outputs
        self.connect('post.Loadouts','Loadouts')
        #self.connect('windLoadst.windLoads','twrWindLoads')
        #self.connect('waveLoadst.waveLoads','twrWaveLoads')
        self.create_passthrough('windLoadst.windLoads')
        self.create_passthrough('waveLoadst.waveLoads')


if __name__ == '__main__':

    test = JcktLoad()

    # parameters from other components
    test.nlegs = 4
    test.nodes = np.array([[-7.7, -7.7, 0.0], [-7.6, -7.6, 1.5], [-6.71723019101, -6.71723019101, 14.7415471348], [-5.93610731855, -5.93610731855, 26.4583902218], [-5.24492719924, -5.24492719924, 36.8260920114], [-4.58, -4.58, 46.8], [7.7, -7.7, 0.0], [7.6, -7.6, 1.5], [6.71723019101, -6.71723019101, 14.7415471348], [5.93610731855, -5.93610731855, 26.4583902218], [5.24492719924, -5.24492719924, 36.8260920114], [4.58, -4.58, 46.8], [7.7, 7.7, 0.0], [7.6, 7.6, 1.5], [6.71723019101, 6.71723019101, 14.7415471348], [5.93610731855, 5.93610731855, 26.4583902218], [5.24492719924, 5.24492719924, 36.8260920114], [4.58, 4.58, 46.8], [-7.7, 7.7, 0.0], [-7.6, 7.6, 1.5], [-6.71723019101, 6.71723019101, 14.7415471348], [-5.93610731855, 5.93610731855, 26.4583902218], [-5.24492719924, 5.24492719924, 36.8260920114], [-4.58, 4.58, 46.8], [8.881784197e-16, -7.1314002458, 8.52899631296], [0.0, -6.30255839886, 20.9616240171], [0.0, -5.569153853, 31.9626922049], [-1.7763568394e-15, -4.88996327105, 42.1505509342], [7.1314002458, 8.881784197e-16, 8.52899631296], [6.30255839886, -1.7763568394e-15, 20.9616240171], [5.569153853, -8.881784197e-16, 31.9626922049], [4.88996327105, -8.881784197e-16, 42.1505509342], [-8.881784197e-16, 7.1314002458, 8.52899631296], [0.0, 6.30255839886, 20.9616240171], [0.0, 5.569153853, 31.9626922049], [1.7763568394e-15, 4.88996327105, 42.1505509342], [-7.1314002458, -8.881784197e-16, 8.52899631296], [-6.30255839886, 1.7763568394e-15, 20.9616240171], [-5.569153853, 8.881784197e-16, 31.9626922049], [-4.88996327105, 8.881784197e-16, 42.1505509342], [0.0, -7.6, 1.5], [7.6, 0.0, 1.5], [0.0, 7.6, 1.5], [-7.6, 0.0, 1.5], [0.0, 0.0, 46.8], [0.0, 0.0, 50.8], [0.0, 0.0, 53.8], [0.0, 0.0, 54.8], [0.0, 0.0, 57.6333333333], [0.0, 0.0, 60.4666666667], [0.0, 0.0, 63.3], [0.0, 0.0, 66.1333333333], [0.0, 0.0, 68.9666666667], [0.0, 0.0, 71.8], [0.0, 0.0, 81.11], [0.0, 0.0, 90.42], [0.0, 0.0, 99.73], [0.0, 0.0, 109.04], [0.0, 0.0, 118.35], [0.0, 0.0, 127.66], [0.0, 0.0, 130.0]])
    test.pillegDs = np.array([2.0, 1.8, 1.6, 1.6, 1.6])
    test.twrDs = np.array([5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.9, 5.53083333333, 5.16166666667, 4.7925, 4.42333333333, 4.05416666667])
    test.Twrmems = np.array([[48, 49], [49, 50], [50, 51], [51, 52], [52, 53], [53, 54], [54, 55], [55, 56], [56, 57], [57, 58], [58, 59], [59, 60], [60, 61]])
    test.Legmems = np.array([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12], [13, 14], [14, 15], [15, 16], [16, 17], [17, 18], [19, 20], [20, 21], [21, 22], [22, 23], [23, 24]])
    test.Pilemems = []
    test.VPFlag  = False
    test.al_bat3D=.78
    test.RNA_F       =[1000.,0.,0.,0.,0.,0.]
    test.TwrRigidTop =True
    test.RNAinputs.CMoff[2]=1.34
    test.gravity=9.81


    # inputs for utilization
    test.topF = [1.3295e6, 0.0, -1e6]
    test.topM = [6.2829e6, 0.0, 0.0]
    test.twrTs = test.twrDs/200.0
    test.n_reinforced = 3
    test.E = 210e9
    test.sigma_y = 450.0e6
    test.gamma_f = 1.35
    test.gamma_m = 1.3
    test.gamma_n = 1.0


    # wind inputs
    test.windIns.al_shear = 0.14
    test.windIns.psi = 0.0
    test.windIns.rho = 1.225
    test.windIns.U50HH = 70.0
    test.windIns.HH = 100.0
    test.windIns.mu = 1.7934e-5

    # wave inputs
    test.waterIns.wdepth = 30.0
    test.waterIns.z_floor = 0.0
    test.waterIns.wlevel = 30.0
    test.waterIns.T = 12.0
    test.waterIns.HW = 10.0
    test.waterIns.psi = 0.0
    test.waterIns.rho = 1027.0
    test.waterIns.Cm = 2.0
    test.waterIns.mu = 1.3351e-3
    test.waterIns.Uc = 2.0

    test.run()

    #print test.stress
    #print test.z_buckling
    #print test.buckling