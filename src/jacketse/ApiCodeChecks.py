#-------------------------------------------------------------------------------
# Name:        ApiCodeChecks.py
# Purpose:     It contains main equations for API code Checks
#
# Author:      rdamiani
#
# Created:     11/09/2012;
# Modified:    12/13:     modified mbr_strct, mbr_idx to accept JacketSE.py in openmdao
#               4/14: Andrew expanded the tutil and modified cbutil so to remove jumps
# Copyright:   (c) rdamiani 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np

from commonse.utilities import cubic_spline_eval, smooth_min

#from jacket import Mip_Mop, Mjnt
#Here are some main constants to be assigned at the beginning, and that potnetially could be changed by the user
FSchord=1.2 #Safety Factor for Cmpilehord/Joint varserification API
FSjAPI=1.6 # Safety factor for joint checks
delta_all=4./3  #Delta allowable for API, to multiply allowables

def main():
    """Sample call: for testing (see also SweepJacket.py) \n"""
    import BuildJckt                                                  #\n
    InFrame3DD='C:\RRD\PYTHON\TURBINEDESIGN\JACKET_CST\jckt_UHPMDD.3DD' #\n
    JcktInfo='C:\RRD\PYTHON\TURBINEDESIGN\JACKET_CST\Jckt_Twr_Inp_PMDDbattered.py' #\n
    f1,f2,mass,MbrFrcs,Jckt=BuildJckt.main(JcktInfo,InFrame3DD)           #\n
    mbr_idx=np.arange(1,Jckt.n_mems+1) #indices of members for which we will do API checks
    print('Jacket Mass (only member mass) is {f5.3}'.format(Jckt.mass()))
    import ApiCodeChecks
    #Member Check
    t_chk,t_util, cb_chk, cb_util=ApiCodeChecks.ApiMbrChk(MbrFrcs,mbr_idx,Jckt.mbr_strct) #\n
    # For joint checks, check out SweepJacket.py as well
    KjntChk,KjntUtil=ApiCodeChecks.KjntChk(Jckt,MbrFrcs)
    XjntChk,XjntUtil=ApiCodeChecks.XjntChk(Jckt,MbrFrcs)

    """"""
    pass

if __name__ == '__main__':
    main()
#______________________________________________________________________________#
def Mip_Mop(MMvec,Mloc,Mjoint):
    """This function computes the In-Plane and Out-Of-Plane portion of MMvec, \n
    where MMvec is a moment (vector) along the local frame element CS. \n
    INPUT \n
    MMvec  -float(3,n), Mxx,Myy,Mzz following Frame3DD reference frame [Nm] \n
    Mloc   -float(3,3,n), direction cosine matrix identifying the local element triad, with \n
            x along axis, y and z along principal axis  \n
    Mjoint -float(3,3,n), direction cosine matrix identifying the joint triad, with \n
            x,y in the joint plane and z normal to it  \n
    OUTPUT    \n
    MIP,|MOoP| -2*float(n), norm of MOoP and signed values of the In-Plane bending moments (Kjnt component) [Nm] \n
                       """
                       #To revise dot produce and vectorization HERE
    from numpy import cos,sin,arctan2,arctan
    out=np.zeros([2,np.shape(MMvec)[1]])
    M_OoP=np.array([MMvec[0,:]*np.dot(Mloc[0,:],Mjoint[0,:])+\
                    MMvec[1,:]*np.dot(Mloc[1,:],Mjoint[0,:])+\
                    MMvec[2,:]*np.dot(Mloc[2,:],Mjoint[0,:])    ,\
                    MMvec[0,:]*np.dot(Mloc[0,:],Mjoint[1,:])+\
                    MMvec[1,:]*np.dot(Mloc[1,:],Mjoint[1,:])+\
                    MMvec[2,:]*np.dot(Mloc[2,:],Mjoint[1,:])    ,\
                    np.zeros_like(MMvec[2,:])])  #[Mxx*iloc.i_jnt+Myy*jloc.i_jnt+Mzz*Kloc.i_jnt ; Mxx*iloc.j_jnt+Myy*jloc.j_jnt+Mzz*Kloc.j_jnt; 0 0 0 ...0 0 0 0] [3,n]
    #Note M_IP will have the right sign for API convention, + if fibers compress brace
    M_IP=MMvec[0,:]*np.dot(Mloc[0,:],Mjoint[2,:])+\
         MMvec[1,:]*np.dot(Mloc[1,:],Mjoint[2,:])+\
         MMvec[2,:]*np.dot(Mloc[2,:],Mjoint[2,:])

    #Now return the magnitude for M_OoP, and the signed M_IP


    return M_IP,np.sqrt(np.sum(M_OoP**2,0))
#______________________________________________________________________________#
def Mjnt(ich,ibrc):
    """This function finds the triad xl,yl,zl for a local joint. \n
    INPUT \n
    ich   -float(3,n), direction cosines (components of unit vector) for chord x-axis unit vector \n
    ibrc  -float(3,n), direction cosines (components of unit vector) for brace x-axis unit vector \n
    OUTPUT    \n
    out   -float(3,3,n), direction cosine matrices of coordinate system having
             x-axis along chord, z-axis normal to the plane of joint, y-axis following RHR \n"""
    ich=ich.reshape(3,-1) #just to the next commands will work
    ibrc=ibrc.reshape(3,-1) #just to the next commands will work
    out=np.zeros([3,3,ich.shape[1]])
    out[2,:,:]=np.cross(ich,ibrc,0,0,0) #ich X ibrc=Kjnt
    out[1,:,:]=np.cross(out[2,:,:],ich,0,0,0) #Kjnt X ich =Jjnt
    out[0,:,:]=ich #Ijnt=ich
    return out

#______________________________________________________________________________#
#API allowables
def ApiFt(fy):
   """Tension Allowable for API"""
   return 0.6*fy

def ApiFv(fy):
   """Shear Allowable for API"""
   return 0.4*fy

def ApiFxe(E,D,t):
    """Elastic local buckling allowable stress"""
    return 2*0.3*E*t/D

def ApiFxc(E,fy,D,t):
    """Inelastic local buckling allowable stress"""
    #make sure we can treat floats as well as arrays as input
    fyy=np.array([fy])
    EE=np.array([E])
    DD=np.array([D])
    tt=np.array([t])

    Fxc = fyy*(1.64-0.23*(DD/tt)**0.25)

    return Fxc

def ApiCmPile(fa,Fe):
    """Reduction factor for buckling allowables for pile and leg"""
    Cm, _, _ = smooth_min(1.-0.4*fa/Fe, 0.85, pct_offset=0.01)
    return Cm

def ApiFa(E,fy,D,t,Klr):
    """Allowable axial stress in compression"""
    #make inputs into arrays just incase they were simple floats
    fyy=np.array([fy])
    EE=np.array([E])
    DD=np.array([D])
    tt=np.array([t])
    KlrK=np.array([Klr])

    n = len(fyy.squeeze())
    fyy2 = np.zeros(n)
    for i in range(n):
        fyy2[i], _ , _ = smooth_min(ApiFxc(EE.squeeze()[i], fyy.squeeze()[i], DD.squeeze()[i], tt.squeeze()[i]), fyy.squeeze()[i])

    Cc=np.sqrt(2.*np.pi**2*EE/fyy2)
    out=(1-KlrK**2/(2.*Cc**2))*fyy2 / (5./3. + 3./8.*KlrK/Cc - KlrK**3/(8.*Cc**3))
    idx= (Klr>=Cc)
    out[idx]=12.*np.pi**2 * EE[idx]/ (23. * KlrK[idx]**2)
    return out.squeeze()

def ApiFb(E,fy,D,t):
    """Allowable axial stress in bending"""
    #Make sure we can process single floats besides arrays
    fyy=np.array([fy])
    EE=np.array([E])
    DD=np.array([D])
    tt=np.array([t])

    n = len(fyy)
    Fb = np.zeros(n)

    for i in range(n):
        dt = (DD/tt).squeeze()[i]
        fy2 = fyy.squeeze()[i]
        E = EE.squeeze()[i]
        delta = 0.2*10340.0/fy2
        x1 = 10340.0/fy2
        x2 = 20680.0/fy2

        if dt < x1-delta:
            Fb[i] = 0.75*fy2
        elif dt < x1+delta:
            Fb[i] = cubic_spline_eval(x1-delta, x1+delta, 0.75*fy2, (0.84 - 1.74*fy2*(x1+delta)/E)*fy2,
                                      0.0, -1.74*fy2/E*fy2, dt)
        elif dt < x2-delta:
            Fb[i] = (0.84 - 1.74*fy2*dt/E)*fy2
        elif dt < x2+delta:
            Fb[i] = cubic_spline_eval(x2-delta, x2+delta, (0.84 - 1.74*fy2*(x2-delta)/E)*fy2,
                                      (0.72 - 0.58*fy2*(x2+delta)/E)*fy2, -1.74*fy2/E*fy2, -0.58*fy2/E*fy2, dt)
        else:
            Fb[i] = (0.72 - 0.58*fy2*dt/E)*fy2

##        Fb=(0.72 - 0.58*fy2*dt/E)*fy2  #initialize
##        idx=(DD/tt) < (x1+delta)
##        Fb[idx]=(cubic_spline_eval(x1-delta, x1+delta, 0.75*fy2, (0.84 - 1.74*fy2*(x1+delta)/E)*fy2,
##                                      0.0, -1.74*fy2/E*fy2, dt))[idx]
##        idx=(DD/tt) < (x1-delta)
##        Fb[idx]=0.75*fy2[idx]
##
##        idx=(DD/tt) < (x2+delta)


    return Fb


#Code Checks for various loading conditions
def chk_tens(f_zz,fy,delta_all):
    """This function performs check on tensioned member"""
    chk=f_zz / (ApiFt(fy)*delta_all)
    return (np.abs(chk)<=1), chk

def chk_shr(tau,fy,delta_all):
    """This function performs check on sheared member"""
    chk= tau / (ApiFv(fy)*delta_all)

    return (chk<=1), chk

def chk_cb(fa,fbx,fby,delta_all,E,fy,D,t,Klr):

    Fb = ApiFb(E, fy, D, t)
    Fa = ApiFa(E, fy, D, t, Klr)
    Fep = 12./23. * np.pi**2 * E * delta_all/Klr**2  #12/23 is a PSF that AISC 335-89 used in teh past, I have reestablished that on 8/18/2015
    Cm = ApiCmPile(fa, Fep)

    chk1 = fa/(Fa*delta_all) + Cm*np.sqrt(fbx**2+fby**2)/ \
        ((1.- fa/Fep)*Fb*delta_all)
    chk2 = fa/(0.6*fy*delta_all) + np.sqrt(fbx**2+fby**2)/(Fb*delta_all)

    chk= (chk1 <= 1) & (chk2 <=1)
    return chk, chk1.squeeze(), chk2

#Member checks
def ApiMbrChk(MbrFrcs,mbr_idx,mbr_strct):
    """Verify tension and compression+bending for the various members.
    For the time being we assume we can just read the main output file, i.e. we assume
    that the largest stresses are at the nodes, which may not be the case,
    but in the future we will have more element per member and so that would be ok.
    INPUT:  MbrFrcs - an 2nM x 9 array, as read from the output file of Frame3DD, with [nM,nN,Nx,Vx,Vy,Mzz,Mxx,Myy,ToC],
                        nM being the member number of the element (repeated for node 1 and 2)
                        nN being the node of the element (node i-th and j-th)
                        Nx,Vx,Vy: axial and shear forces
                        Mxx,Myy:  bending moments (note Frame 3DD calls them Myy and Mzz, as it calls x what I call z)
                        Mzz:      torque
                        ToC:    1=tension, -1=compression element force at that node
            mbr_strct - array(n) of tube objects with all the structural properties needed
            mbr_idx - integer(m), indices of members to be checked (they may be non-contigous) first member =0 for Python, not as in Frame3DD
    \n
    Sample Call: mbr_idx=np.arange(0,Jckt.n_mems) \n#
                 mbr_strct=Jckt.mbr_strct         \n#
                 ApiMbrChk(MbrFrcs,mbr_idx,mbr_strct) \n# """

    #Scan through all of the members (elements) (vectorization is possible)
    #need to account for indices repeated twice per member, given the 2 nodes
    junk=mbr_idx*2
    mbr_idx2=np.array([junk,junk+1]).flatten()
    mbr_idx2.sort()

    #mbr_idx2=mbr_idx.repeat(2)
    #mbr_idx3=mbr_idx2-1 #this is to operate within python's indices
    #mbr_idx=mbr_idx-1 #this is to operate within python's indices

    N=MbrFrcs[mbr_idx2,2] #axial force
    Vx=MbrFrcs[mbr_idx2,3]
    Vy=MbrFrcs[mbr_idx2,4]
    Mzz=MbrFrcs[mbr_idx2,5] #torsion
    Mxx=MbrFrcs[mbr_idx2,6]
    Myy=MbrFrcs[mbr_idx2,7]

    ToC=MbrFrcs[mbr_idx2,8]  #tension or compression

    #Properties
    D=mbr_strct.D[mbr_idx].repeat(2,0) #Diameter
    t=mbr_strct.t[mbr_idx].repeat(2,0) #thickness
    Area=mbr_strct.Area[mbr_idx].repeat(2,0) #x-sec areas
    BdgMxx=mbr_strct.BdgMxx[mbr_idx].repeat(2,0) #Bending modulus for bending about xx
    BdgMyy=mbr_strct.BdgMyy[mbr_idx].repeat(2,0) #Bending modulus for bending about yy
    #dpth=mbr_strct.z[mbr_idx]

    Asx=mbr_strct.Asx[mbr_idx].repeat(2,0) #shear areas along x
    Asy=mbr_strct.Asy[mbr_idx] .repeat(2,0)#shear areas along x
    Jxx=mbr_strct.Jxx[mbr_idx].repeat(2,0) #Jxx
    Jyy=mbr_strct.Klr[mbr_idx].repeat(2,0) #Jyy
    Jzz=mbr_strct.Klr[mbr_idx].repeat(2,0) #Jzz=Jo
    Klr=mbr_strct.Klr[mbr_idx].repeat(2,0) #Jzz=Jo

    E=np.array([oinst.E for oinst in mbr_strct.mat[mbr_idx]]).repeat(2,0).astype(float) #Young's Modulus
    G=np.array([oinst.G for oinst in mbr_strct.mat[mbr_idx]]).repeat(2,0).astype(float) #Shear Modulus
    rho=np.array([oinst.rho for oinst in mbr_strct.mat[mbr_idx]]).repeat(2,0).astype(float) #density
    fy=np.array([oinst.fy for oinst in mbr_strct.mat[mbr_idx]]).repeat(2,0).astype(float) #density

    #Calculate stresses as fa, fbx, fby (fby and fbz for frame3DD ref sys)
    f_N= N/Area  #axial stress due to N alone  # TODO: no abs (also below, multipler times)
    fbx=Mxx/BdgMxx        #axial stress due to Mxx alone
    fby=Myy/BdgMyy        #axial stress due to Myy alone
    #Combine axial stresses
    # f_tens=f_N+(abs(fbx)+abs(fby))*ToC #We are conservatively putting max stresses in the same point though it is not (x vs y)
    f_stress1 = f_N + fbx
    f_stress2 = f_N - fbx
    f_stress3 = f_N + fby
    f_stress4 = f_N - fby

    #tens util accounts for both compressed and tensile, assuming fy applies for both, so use abs in teh end
    # t_chk,t_util=chk_tens(f_tens,fy,delta_all) # this contains results for compressed members too, which is not what we want so we can remove
    t_chk1, t_util1 = chk_tens(f_stress1, fy, delta_all)
    t_chk2, t_util2 = chk_tens(f_stress2, fy, delta_all)
    t_chk3, t_util3 = chk_tens(f_stress3, fy, delta_all)
    t_chk4, t_util4 = chk_tens(f_stress4, fy, delta_all)
    t_chk = np.concatenate([t_chk1, t_chk2, t_chk3, t_chk4])
    t_util = np.concatenate([t_util1, t_util2, t_util3, t_util4])
    # t_util=np.abs(t_util.astype(float)) #remove abs for optimizer's sake
    #t_util=t_util.astype(float)
    #t_util[ToC<1]=np.NaN  #decided to leave t_util to consider both compressend and tensile members against fy strength
    #max_tutil=np.nanmax(t_util)

    #now isolate only compressed members, return nans elsewhere
    # cb_util=np.zeros(np.size(mbr_idx2)) #*np.NaN #remove nans for optimizer's sake
    # cb_chk=np.zeros(np.size(mbr_idx2)) #*np.NaN
    ####debug
    ###a,b,c=chk_cb(119.e6,111.275e6,107.219e6,4./3,2.1e11,3.447e8,1.713,0.034,29.)
    # idx_c=np.nonzero(ToC<1)[0]
    # cb_chk[idx_c],cb_util1,cb_util2=chk_cb(abs(f_N[idx_c]),fbx[idx_c],fby[idx_c],delta_all,E[idx_c],fy[idx_c],D[idx_c],t[idx_c],Klr[idx_c]) #[2n], note we do not care about signs of fbx,fby as tehy are squared in this check
    cb_chk, cb_util1, cb_util2=chk_cb(-f_N,fbx,fby,delta_all,E,fy,D,t,Klr)  # compressive members become positive
    cb_util = np.concatenate([cb_util1, cb_util2])


    # cb_util[idx_c]=np.nanmax(np.vstack((cb_util1,cb_util2)),0)

    return t_chk, t_util, cb_chk, cb_util


#______________________________________________________________________#
#ALLOWABLES FOR JOINTS
def ApiPa(Qus,Qfs,fy_chord,tht,tchord,FSjnt=1.6):
    """Brace-Axial-load allowable capacity of joint.
    INPUT:
        Qus      -float(n), Strength factors
        Qfs      -float(n), Chord load factors
        fy_chord -float(n), Chord yield stress at the joint or 0.8 of tensile strength if less [Pa]
        FSjnt    -float, joint safety factor from API, default=1.6
        tht      -float(n), angle between chord and brace [rad](<=pi/2)
        tchord   -float(n), thickness of the chord [m]"""
    return Qus*Qfs*fy_chord*tchord**2 / (FSjnt*np.sin(tht))

def ApiMa(Qus,Qfs,fy_chord,tht,Dbrc,tchord,FSjnt=1.6):
    """Brace-Axial-load allowable capacity of joint.
    INPUT:
        Qus      -float(n), Strength factors
        Qfs      -float(n), Chord load factors
        fy_chord -float(n), Chord yield stress at the joint or 0.8 of tensile strength if less [Pa]
        FSjnt    -float, joint safety factor from API, default=1.6
        tht      -float(n), angle between chord and brace [rad](<=pi/2)
        Dbrc     -float(n), brace diameter [m]
        tchord   -float(n), thickness of the chord [m]"""
    return Qus*Qfs*fy_chord*tchord**2 *Dbrc / (FSjnt*np.sin(tht))

#JOINT FUNCTIONS
def gam_chord(Dchord,tchord):
    """This function returns the gamma factor from API.
    INPUT
    Dchord -float(n), chord OD [m]
    tchord -float(n), chord thickness [m]"""
    return Dchord/(2.*tchord)

def Qg(gap,Dchord,tchord,tbrc,fy_chord,fy_brc):
    """This function returns the gap factor.
    INPUT
    gap    -float(n), gap size [m]
    Dchord -float(n), chord OD [m]
    tchord -float(n), chord thickness [m]
    tbrc   -float(n), brace thickness [m]
    fy_chord -float(n), chord yield stress at the joint or 0.8 of tensile strength if less [Pa]
    fy_brace -float(n), brace yield stress at the joint or 0.8 of tensile strength if less [Pa]"""
    g2D=gap/Dchord
    Qg1=np.max([1.,1.+0.2*(1.-2.8*g2D)**3])
    phi=tbrc*fy_brc / (tchord*fy_chord)
    Qg2=0.13+0.65*phi*gam_chord(Dchord,tchord)**0.5
    if g2D>=0.05:
        out=Qg1
    elif g2D<=-0.05:
        out=Qg2
    else:
        out=np.interp(g2D,[-0.05,0.05],[Qg2,Qg1])
    return out

def Qbeta(Dchord,Dbrc):
    """This function returns the Geometric Factor
    INPUT
    Dchord -float(n), chord OD [m]
    tchord -float(n), chord thickness [m]"""
    bet=Dbrc/Dchord*np.ones(np.size(Dchord))
    out=1.*np.ones(np.size(Dchord))
    out[bet>0.6]=0.3/ (bet[bet>0.6]*(1.-0.833*bet[bet>0.6]))
    return out

def Qu(gap,Dchord,Dbrc,tchord,tbrc,fy_chord,fy_brc,Jtype,BrcLoad):
    """This function returns the Strength Factor
    INPUT
    gap    -float(n), gap size [m] (significant only if Jtype='K'
    Dchord -float(n), chord OD [m]
    Dbrc   -float(n), brace OD [m]
    tchord -float(n), chord thickness [m]
    tbrc   -float(n), brace thickness [m]
    fy_chord -float(n), chord yield stress at the joint or 0.8 of tensile strength if less [Pa]
    fy_brace -float(n), brace yield stress at the joint or 0.8 of tensile strength if less [Pa]
    Jtype     -string(n), type of joint, "X","K","T/Y"
    BrcLoad   -string(n), brace loading condition, "tension","compression","IP_bending","OoP_bending"""
    Qgs=Qg(gap,Dchord,tchord,tbrc,fy_chord,fy_brc)*np.ones(np.size(tbrc))
    Qbetas=Qbeta(Dchord,Dbrc)*np.ones(np.size(tbrc))
    gams=gam_chord(Dchord,tchord)*np.ones(np.size(tbrc))
    betas=Dbrc/Dchord*np.ones(np.size(tbrc))

    if BrcLoad=='IP_Bending':  #In-plane bending
        out=(5.+0.7*gams)*betas**1.2
    elif BrcLoad=='OoP_Bending':  #Out-of-Plane Bending
        out=(4.5+0.2*gams)*betas**2.6 + 2.5

    elif BrcLoad=='Tension':  #Tension
        if Jtype=='K':
            out=np.min([(16.+1.2*gams)*betas**1.2*Qgs,40.*betas**1.2*Qgs])
        elif Jtype=='T/Y':
            out=30.*betas
        elif Jtype=='X':
            out=20.7+(betas-0.9)*(17.*gams-220.)
            out[betas<=0.9]=23.*betas[betas<=0.9]

    elif BrcLoad=='Compression':  #Compression
        if Jtype=='K':
            out=np.min([(16.+1.2*gams)*betas**1.2*Qgs,40.*betas**1.2*Qgs])
        elif Jtype=='T/Y':
            out=np.min([2.8+36*betas**1.6,2.8+(10.+0.8*gams)*betas**1.6])
        elif Jtype=='X':
            out=(2.8+(12.+0.1*gams)*betas)*Qbetas

    return out

def ApiGap(Dchord,Dbrc1,Dbrc2,a1,a2):
    """This function calculates the Gap Size for a jointwith up to 2 braces in the plane.
    INPUT
    Dchord -float(n), chord OD [m]
    Dbrc1  -float(n), brace1 OD [m]
    Dbrc2  -float(n), brace2 OD [m]
    a1     -float(n), angle between brace 1 and chord [rad]
    a2     -float(n), angle between brace 2 and chord [rad]"""
    return Dchord/2.*(np.tan(np.pi/2.-a1)+np.tan(np.pi/2.-a2)) \
        -Dbrc2/np.sin(a2)/2.-Dbrc1/np.sin(a1)/2.

def Cchord(Dchord,Dbrc,Jtype,BrcLoad):
    """C1,2,3: Chord parameters for load factor.
    INPUT
    Dchord -float(n), chord OD [m]
    Dbrc   -float(n), brace OD [m]
    Jtype     -string(n), type of joint, "X","K","T/Y"
    BrcLoad   -string(n), brace loading condition, "tension","compression","IP_bending","OoP_bending"""
    betas=Dbrc/Dchord*np.ones(np.size(Dbrc))
    if BrcLoad=='IP_Bending' or BrcLoad=='OoP_Bending':
        out=np.array([0.2,0.,0.4]).reshape(3,1).repeat(np.size(Dbrc),1)
    elif Jtype=='K':
        out=np.array([0.2,0.2,0.3]).reshape(3,1).repeat(np.size(Dbrc),1)
    elif Jtype=='T/Y':
        out=np.array([0.3,0.,0.8]).reshape(3,1).repeat(np.size(Dbrc),1)
    elif Jtype=='X':
        out=np.array([0.2,0.,0.5]).reshape(3,1).repeat(np.size(Dbrc),1)
        out1=np.array([-0.2,0.,0.2]).reshape(3,1).repeat(np.size(Dbrc),1)
        idxs=(betas==1)
        out[:,idxs]=out1[:,idxs]
        idxs=(betas>0.9) & (betas<1)
        #Manual interpolation as I could not vetorize otherwise
        out[:,idxs]=out[:,idxs]+(out1[:,idxs]-out[:,idxs])/(1-.9)*(betas[idxs]-0.9)

    return out

def Qf(Dchord,Dbrc,tchord,fy_chord,Pc,M_IPc,Jtype,BrcLoad,FSc=1.2):
    """Chord Load Factor.
    INPUT
    Dchord    -float(n), chord diameter [m]
    Dbrc      -float(n), brace diameter [m]
    tchord    -float(n), chord thickness [m]
    fy_chord  -float(n), chord yield stress [Pa]
    Pc        -float(n), nominal axial load [N]
    M_IPc     -float(n), nominal In_bending load for the chord; note: Total bending = sqrt(M_IPc**2+M_OoPc**2) [Nm],
                         positive if it compresses at the foot of the brace
    Jtype     -string(n), type of joint, "X","K","T/Y"
    BrcLoad   -string(n), brace loading condition, "tension","compression","IP_bending","OoP_bending"
    FSc       -float(n), safety factor for chord load factor from API """
    from commonse.Tube import Tube
    #make sure  FSc is a vector
    FSc=FSc*np.ones(np.size(Dchord)) #to revise, as it may not work when Dchord is a scalar

    chord=Tube(Dchord,tchord)
    Py=fy_chord*chord.Area  #yield capacity of chord
    Mp=chord.Jxx/((Dchord-2.*tchord)/2) * fy_chord #plastic moment capacity of chord
    C=Cchord(Dchord,Dbrc,Jtype,BrcLoad) #joint type coefficients
    A=np.sqrt( (FSc*Pc/Py)**2+(FSc*M_IPc/Mp)**2 ) #parameter
    return 1.+C[0]*FSc*Pc/Py - C[1]*FSc*M_IPc/Mp  - C[2]*A**2

#Joint Checks
def ApiJntChk(Pc,M_IPc,P,M_IP,M_OoP,fy_chord,fy_brc,tht,gap,Dchord,Dbrc,tchord,tbrc,Jtype,DeltaAll=4./3.,FSjnt=1.6,FSc=1.2):
    """The actual check routine on the joint.
    INPUT:
        Pc       -float(n), axial load on chord [N] (average between either side of the joint)
        M_IPc    -float(n), In-Plane bending load on chord [NM] (average between either side of the joint)
        P        -float(n), axial load on brace [N]
        Pa       -float(n), API calculated axial-load allowable for brace [N]
        M_IP     -float(n), In-plane bending moment for brace [Nm]
        M_OoP    -float(n), Out-of-plane bending moment for brace [Nm]
        Ma       -float(n), allowable capacity for bending moment brace [Nm]
        DeltaAll -float, allowable increase per API, (default =4./3)
    OUTPUT:
        chk     -bool(4,n), utilization check flags (1,0) pass(1) or not(0)
        JntUtil -float(4,n), Utilization for the 4 load conditions Tension/Compression/IPbending/OoPBending

    """

    #Explore all load possibilities for brace checks
    braceloads=('Tension','Compression','IP_Bending','OoP_Bending')
    chk=np.ones(np.size(Dbrc)).astype(int)
    JntUtil=np.zeros(np.size(Dbrc))*np.NaN

    #Now get the right allowables for every P and MIP,MOoP
    #First calculate them all
    Pa=np.zeros([2,np.size(Dbrc)])*np.NaN
    Ma=Pa.copy()
    for ild,BrcLoad in enumerate(braceloads):
        Qus=Qu(gap,Dchord,Dbrc,tchord,tbrc,fy_chord,fy_brc,Jtype,BrcLoad)
        Qfs=Qf(Dchord,Dbrc,tchord,fy_chord,Pc,M_IPc,Jtype,BrcLoad,FSc)
        if ild<2:
            Pa[ild,:]=ApiPa(Qus,Qfs,fy_chord,tht,tchord,FSjnt)
        else:
            Ma[ild-2,:]=ApiMa(Qus,Qfs,fy_chord,tht,Dbrc,tchord,FSjnt)
    #take care of allowable in compression/tension when P is compression/tension
    Pall=Pa[0,:] #Tension Allowable first
    idx=np.nonzero(P<0.)[0]
    Pall[idx]=Pa[1,idx] #Now it works for both compression and tension

    #Here is the check
    JntUtil= np.abs(P/(Pall*DeltaAll)) +(M_IP/(Ma[0,:]*DeltaAll))**2 + np.abs(M_OoP/(Ma[1,:]*DeltaAll))
    chk=JntUtil<=1.

    return chk, JntUtil

#___________________________________________________________________________#

def KjntChk(KjntIDs,mems,MbrFrcs,mbr_strct):
    """This function sweeps through the K-joint nodes of a jacket, and performs K-joint-checks.
    INPUT: \n
        KjntIDs  - array(nKjnts), K-Joint IDs, as from Frame3DD
        mems     - array(m,2), connectivity array, for each member ith, mems[i,0:2] are the 1st and 2nd node
        MbrFrcs  -   MbrFrcs - an 2nM x 9 array, as read from the output file of Frame3DD, with [nM,nN,Nx,Vx,Vy,Mzz,Mxx,Myy,Toc],
                        nM being the member number of the element (repeated for node 1 and 2)
                        nN being the node of the element (node i-th and j-th)
                        Nx,Vx,Vy: axial and shear forces
                        Mxx,Myy:  bending moments (note Frame 3DD calls them Myy and Mzz, as it calls x what I call z)
                        Mzz:      torque
                        ToC:    1=tension, -1=compression element force at that node
                        mbr_idx   - integer(m), indices of members to be checked (they may be non-contigous) first member =0 for Python, not as in Frame3DD
        mbr_strct - array(n) of tube objects with all the structural properties needed
    OUTPUT
        KjntUtil -float(nKjnts), utilization for nKjnts k-joint nodes
        KjntChk  -int(nKjnts), flag (1/0), K-joint passes(1) or fails(0)

        """

    #Get all the nodes for K-joints and look in MbrFrcs for the right indices
    #However, since we will have multiple planes, we need to repeat calculations for all involved planes (2 at a time)
    KjntUtil=np.zeros(KjntIDs.size)*np.NaN
    KjntChk=np.ones(KjntIDs.size).astype(int)
    Jtype='K'
    for jj,nj in enumerate(KjntIDs):
        #consider the first 2 hits in MbrFrcs as the chord,
        #the second 2's as the braces in one plane, the third 2's as braces in the other plane; this is for the normal Xbrace joint
        #however, the bottom joint may have the mudbrace attached, in which case we may get mix-bag of hits, first 2 are always chord
        jntfrcs=MbrFrcs[MbrFrcs[:,1]==KjntIDs[jj],:] #this gives me the T/F indices to select forces for K brace checks

        #check if we are at the top of the leg and we may be including TP braces, girders and struts
        #normally we would get 6 hits, but at the top of the leg we may get 7 if TP girders or if Hbraces are present

        mudbraces=False#Initialize

        if np.abs(jntfrcs[1,0]-jntfrcs[0,0]) !=1 : #This means we are at the top of the leg, or possibly at the bottom at the junction with vertical pile
            legtop=True
            brc_mbr=jntfrcs[1:,0] #brace member element numbers, skipping leg element
            nd2=mems[brc_mbr.astype(int)-1,0] #1st nodes of those 2 brace elements (-1 for python)
            idx=np.nonzero(nd2 == KjntIDs[jj])[0] #index of elements having local joint as start node
            idx1=np.nonzero(nd2 != KjntIDs[jj])[0]#index of elements having local joint as end node
            #brc_mbr1=np.array([brc_mbr[idx],brc_mbr[idx1]]) #This may put them in pairs of 2

            jntfrcs[1:,:]=jntfrcs[np.hstack((1+idx,1+idx1)),:]
            ch_idx=0   #chord index for later

            #Now Chord Loads-from the only chord element
            Pc=np.abs(jntfrcs[0,2])*jntfrcs[0,8] #this takes care of the complex signs at the element ends
            Mzzc=-jntfrcs[0,7] #This is to account for the right signs and directions (- instead of + sign, for avg)
            Myyc=-jntfrcs[0,6] #This is to account for the right signs and directions (- instead of + sign, for avg)
            Mxxc=-jntfrcs[0,5] #This is to account for the right signs and directions (- instead of + sign, for avg)

        elif (jntfrcs[:,0].size==3):   # at the base with just mudbraces, thus the plane to investigate is made up of the two braces in the two different jacket planes
            #Use same Pc,Mzzc,Myyc,Mxxc as at the top
            legtop=False
            mudbraces=True

        else: #inner leg joint
            legtop=False
            brc_mbr=jntfrcs[2:,0] #brace member element numbers, skipping first 2 (leg elements)
            nd2=mems[brc_mbr.astype(int)-1,0] #1st nodes of those 2 brace elements (-1 for python)
            idx=np.nonzero(nd2 == KjntIDs[jj])[0] #index of elements having local joint as start node
            idx1=np.nonzero(nd2 != KjntIDs[jj])[0]#index of elements having local joint as end node
            #brc_mbr1=np.array([brc_mbr[idx],brc_mbr[idx1]]) #This may put them in pairs of 2

            jntfrcs[2:,:]=jntfrcs[np.hstack((2+idx,2+idx1)),:] #
            ch_idx=1
            #Now Chord Loads:averaging
            junk=np.array([0,1])
            Pc=np.mean(np.abs(jntfrcs[junk,2])*jntfrcs[junk,8]) #this takes care of the complex signs at the element ends
            Mzzc=0.5*(jntfrcs[junk[1],7]-jntfrcs[junk[0],7]) #This is to account for the right signs and directions (- instead of + sign, for avg)
            Myyc=0.5*(jntfrcs[junk[1],6]-jntfrcs[junk[0],6]) #This is to account for the right signs and directions (- instead of + sign, for avg)
            Mxxc=0.5*(jntfrcs[junk[1],5]-jntfrcs[junk[0],5]) #This is to account for the right signs and directions (- instead of + sign, for avg)

            ##The 2nd brace may be a horizontal brace, which would show up after the x-brc from another plane
            ##brc_mbr=jntfrcs[np.array([3,4]),0] #4th and 5th hit
            ##nd2=mems[brc_mbr.astype(int)-1,0] #1st nodes of those 2 brace elements (-1 for python)
            ##If the 1st node of any of those elements is not the local joint, then that element belongs to a different plane
            ##if np.any(nd2 != np.array([KjntIDs[jj],KjntIDs[jj]])):
            ##jntfrcs[np.array([3,4]),:]=jntfrcs[np.array([4,3]),:]
            ##else:
            ##print 'Problem with brace identification in Kjnt check, KjntIDs= {:d}'.format(KjntIDs[jj])
        mbr_nos=jntfrcs[:,0].astype(int)  #all member(element) hits ordered appropriately by pairs in planes

        if legtop and np.mod(mbr_nos.size,2)==0:
            #This means we are at the top of a leg, we may get 4 hits(leg,Xbrc, no Hbrc but TP stmp), 6 (leg,Xbrc,Hbrc,TPstmp), or 7 (TP sitting there, thus leg,Xbrc,TPgirder,TPstruts,TPbrcs)
            #If 7 hits, no need to change anything, the above logic takes care of 3 pairs: Girder-Xbrc, TPbrc-TPstrut, Girder-Xbrc
            #In the other cases the 4th hit is the stump, but the treatment is sketchy, leave it alonw
            mbr_nos=np.hstack((mbr_nos[0:2],mbr_nos[3:]))  #we are popping the 3rd hit from the ordered list, which should be the stump all the time:

        nhits=mbr_nos.size
        nplanes=np.max([2,(nhits-ch_idx)/2])  #number of planes to investigate, minimum2 xcept for mudbraces in case we have 3 hits only
        nbrcppl=(nhits-(ch_idx+1))/nplanes  #number of braces per plane

        if mudbraces:
            nplanes=1
            nbrcppl=2

        #Establish Chord Element local transformation matrix, direction cosines, that will be used later
        Mloc_c=mbr_strct.XnsfM[:,:,jntfrcs[ch_idx,0]-1]#chord local system [3,3], for the 2nd element (upper chord if possible) as a convention (also -1 for python)
        Dchord=mbr_strct.D[mbr_nos[ch_idx]-1]          #D for the 2nd element (upper chord if possible) as a convention (also -1 for python)
        tchord=mbr_strct.t[mbr_nos[ch_idx]-1]          #t for the 2nd element (upper chord if possible) as a convention (also -1 for python)
        fy_chord=mbr_strct.mat[mbr_nos[ch_idx]-1].fy   #fy for the 2nd element (upper chord if possible) as a convention (also -1 for python)

        for kk in range(0,nplanes): #cycle on planes, or pairs of braces

            brc_idx=np.arange(0,nbrcppl)+kk*nbrcppl+(ch_idx+1)  #brace indices
            #brc_idx=2*kk+np.arange(1,3)+ch_idx #This works because we usually get 2 leg mem, but at the top, we only get one. We need to revise.

            #Now I need to take care of dimensions for the members that are involved
            #the 2 braces in the plane of interest
            Dbrc=mbr_strct.D[mbr_nos[brc_idx]-1]  #normally [2] however it could be [1] if mudbrace joint
            if any(np.isnan(Dbrc)): continue  #this skips rigid members
            tbrc=mbr_strct.t[mbr_nos[brc_idx]-1]  #[2]
            fy_brc=np.array([obj.fy for obj in mbr_strct.mat[mbr_nos[brc_idx]-1]]) #[2]

            #Here is the local transformation matrix for the 2 braces
            Mloc_b=mbr_strct.XnsfM[:,:,jntfrcs[brc_idx,0].astype(int)-1]#brace local systems [2,3,3] (-1 for python)
            #Brace-to-chord included angle
            tht=np.mod(np.arccos(np.dot(Mloc_c[0,:],Mloc_b[0,:])),np.pi)
            tht=np.min([np.pi-tht,tht],axis=0) # NOT SURE  I NEED THIS #HERE!!!!!!!!!!!!!!!!
            Kgap=ApiGap(Dchord,Dbrc[0],Dbrc[-1],tht[0],tht[-1])

            #Since we want to check both braces in case there are,  do another loop
            for ll in range(0,brc_idx.size):
            #we need to figure out which brace needs to have the Mjnt triad flipped about i_chord to get the right sign for M_IPc
                brc_mbr=jntfrcs[brc_idx[ll],0] #
                nd2=mems[brc_mbr.astype(int)-1,0]#1st node of brace member under consideration
                #if this node number does not coincinde with leg-brc intersection node, Mjnt needs to be flipped
                Mjnt_sign= np.int(nd2==nj)-np.int(nd2!=nj) #This gives me 1 (true) or -1(False)
                #Now brace loads
                P=np.abs(jntfrcs[brc_idx[ll],2])*jntfrcs[brc_idx[ll],8] #Axial load with correct sign +=tension
                Mzzb=jntfrcs[brc_idx[ll],7]
                Myyb=jntfrcs[brc_idx[ll],6]
                Mxxb=jntfrcs[brc_idx[ll],5]

                #IP and OoP bending moments of brace (note the sign for ibrc does not matter for the brace component, we keep it at whatever we get
                M_IPb,M_OoPb=Mip_Mop(np.array([[Mxxb],[Myyb],[Mzzb]]),\
                             Mloc_b[:,:,ll],Mjnt(Mloc_c[0,:],Mloc_b[0,:,ll]) )

                #IP and OoP bending moments of chord (note the sign for ibrc matters for the chord loading)
                M_IPc,M_OPc=Mip_Mop(np.array([[Mxxc],[Myyc],[Mzzc]]),\
                       Mloc_c,Mjnt(Mloc_c[0,:],Mjnt_sign*Mloc_b[0,:,ll]) )

                jnt_chk, jnt_util=ApiJntChk(Pc,M_IPc,P,M_IPb,M_OoPb,fy_chord,fy_brc[ll],tht[ll],Kgap,Dchord,Dbrc[ll],tchord,tbrc[ll],Jtype,DeltaAll=4./3.,FSjnt=1.6,FSc=1.2)
                #Here is the check
                KjntUtil[jj]=np.nanmax(np.vstack((jnt_util,KjntUtil[jj]))) #These get updated for the various chord/brc combinations even fot the given node jj
                KjntChk[jj]&=np.prod(jnt_chk)

    return KjntChk,KjntUtil

#___________________________________________________________________________#

def XjntChk(XjntIDs,mems,MbrFrcs,mbr_strct):
    """This function sweeps through the X-joint nodes of a jacket, and performs X-joint-checks.
    INPUT: \n
        XjntIDs  - array(nXjnts), X-Joint IDs, as from Frame3DD
        mems     - array(m,2), connectivity array, for each member ith, mems[i,0:2] are the 1st and 2nd node
        MbrFrcs  -   MbrFrcs - an 2nM x 9 array, as read from the output file of Frame3DD, with [nM,nN,Nx,Vx,Vy,Mzz,Mxx,Myy,ToC],
                        nM being the member number of the element (repeated for node 1 and 2)
                        nN being the node of the element (node i-th and j-th)
                        Nx,Vx,Vy: axial and shear forces
                        Mxx,Myy:  bending moments (note Frame 3DD calls them Myy and Mzz, as it calls x what I call z)
                        Mzz:      torque
                        ToC:    1=tension, -1=compression element force at that node
        mbr_strct - array(n) of tube objects with all the structural properties needed

    OUTPUT
        XjntUtil -float(nXjnts), utilization for nXjnts X-joint nodes
        XjntChk  -int(nXjnts), flag (1/0), X-joint passes(1) or fails(0)

        """

    #Initialize output
    XjntUtil=np.zeros(XjntIDs.size)*np.NaN
    XjntChk=np.ones(XjntIDs.size).astype(int)

    Jtype='X'
    gap=0.  #Need to put something for gap but it is ignored if X-brc
    for jj,nj in enumerate(XjntIDs):
        #consider the first 2 hits in MbrFrcs as the chord, the second 2 as the brace then reverse
        jntfrcs=MbrFrcs[MbrFrcs[:,1]==XjntIDs[jj],:]#this gives me the T/F indices to select forces for X brace checks

        mbr_nos=jntfrcs[:,0]
        for kk in np.array([0,1]):
            junk=2*kk+np.arange(0,2)
            junk1=np.arange(2,4)-2*kk
            #Now I need to take care of dimensions for the members that are involved
            #Since the members involved in the X brace are all of the same dimension, I need to get parameters for 1 member belonging to the joint only
            Dchord=mbr_strct.D[mbr_nos[junk[1]]-1]  #use upper chord, and then -1 for python, though here the chord is for sure uniform, since it is a virtual chord, being the brace acting as chord
            tchord=mbr_strct.t[mbr_nos[junk[1]]-1]
            fy_chord=mbr_strct.mat[mbr_nos[junk[1]]-1].fy
            Dbrc=mbr_strct.D[mbr_nos[junk1[0]]-1]
            tbrc=mbr_strct.t[mbr_nos[junk1[0]]-1]
            fy_brc=mbr_strct.mat[mbr_nos[junk1[0]]-1].fy
            #Establish Chord Element local transformation matrix, direction cosines, that will be used later
            Mloc_c=mbr_strct.XnsfM[:,:,jntfrcs[junk[1],0]-1] #chord local system [3,3], for the 2nd element (upper chord) as a convention (also -1 for python)
            #Since we want to check both braces, on either sides of chord, do another loop
            for ll in np.array([0,1]):
            #we need to figure out which brace needs to have the Mjnt triad flipped about i_chord to get the right sign for M_IPc
                brc_mbr=jntfrcs[junk1[ll],0] #-1 for python
                nd2=mems[brc_mbr+1,0] #1st node of brace member under consideration mbrNiNj[Jckt.mbrNiNj[:,0]==brc_mbr,1]
                #if this node number does not coincinde with X-brc intersection node, Mjnt needs to be flipped
                Mjnt_sign= np.int(nd2==nj)-np.int(nd2!=nj) #This gives me 1 (true) or -1(False)
                #Now brace loads
                P=np.abs(jntfrcs[junk1[ll],2])*jntfrcs[junk1[ll],8] #Axial load with correct sign +=tension
                Mzzb=jntfrcs[junk1[ll],7]
                Myyb=jntfrcs[junk1[ll],6]
                Mxxb=jntfrcs[junk1[ll],5]
                #Here is the local transformation matrix for the brace
                Mloc_b=mbr_strct.XnsfM[:,:,jntfrcs[junk1[ll],0]-1]#brace local system [3,3] (-1 for python)

                #IP and OoP bending moments of brace (note the sign for ibrc does not matter for the brace component, we keep it at whatever we get
                M_IPb,M_OoPb=Mip_Mop(np.array([[Mxxb],[Myyb],[Mzzb]]),\
                             Mloc_b,Mjnt(Mloc_c[0,:],Mloc_b[0,:]) )
                #Now Chord Loads
                Pc=np.mean(np.abs(jntfrcs[junk,2])*jntfrcs[junk,8]) #this takes care of the complex signs at the element ends
                Mzzc=0.5*(jntfrcs[junk[1],7]-jntfrcs[junk[0],7]) #This is to account for the right signs and directions (- instead of + sign, for avg)
                Myyc=0.5*(jntfrcs[junk[1],6]-jntfrcs[junk[0],6]) #This is to account for the right signs and directions (- instead of + sign, for avg)
                Mxxc=np.mean(jntfrcs[junk,5]) #Here we use the usual mean, as torsion is aligned

               #IP and OoP bending moments of brace (note the sign for ibrc matters for the chord loading)
                M_IPc,M_OPc=Mip_Mop(np.array([[Mxxc],[Myyc],[Mzzc]]),\
                       Mloc_c,Mjnt(Mloc_c[0,:],Mjnt_sign*Mloc_b[0,:]) )

                #Brace-to-chord included angle
                tht=np.mod(np.arccos(np.dot(Mloc_c[0,:],Mloc_b[0,:])),np.pi)
                tht=np.min([np.pi-tht,tht])
                jnt_chk, jnt_util=ApiJntChk(Pc,M_IPc,P,M_IPb,M_OoPb,fy_chord,fy_brc,tht,gap,Dchord,Dbrc,tchord,tbrc,Jtype,DeltaAll=4./3.,FSjnt=1.6,FSc=1.2)
                #Here is the check
                XjntUtil[jj]=np.nanmax(np.vstack((jnt_util,XjntUtil[jj]))) #These get updated for the various chord/brc combinations even fot the given node jj
                XjntChk[jj]&=np.prod(jnt_chk)

    return XjntChk,XjntUtil