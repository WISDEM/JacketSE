#-------------------------------------------------------------------------------
# Name:        PlotJacket.Py
# Purpose:
#
# Author:      rdamiani
#
# Created:     14/05/2014
# Copyright:   (c) rdamiani 2014
# Licence:     <your licence>

#All lines marked CJB+, CJB-, or CJBe have been added, "removed" (commented out), or edited respectively by Casey Broslawski. Summer 2016
#-------------------------------------------------------------------------------
import matplotlib.pyplot as plt #first for peregrine's sake
import numpy as np
from commonse.axisEqual3D import axisEqual3D

from mpl_toolkits.mplot3d import Axes3D

def main(myjckt,util=False,savefileroot=[],toscreen=True):

    #CJB+ Prepare utilization data to be plotted
    #TODO The arrays should be formatted in a dedicated module in jacketSE
    t_util_list=np.array([[1,1]])  #CJB+
    count_t=1  #CJB+
    for i in myjckt.LoadFrameOuts.jacket_utilization.t_util:  #CJB+
        if i>=1:  #CJB+
            step=np.array([count_t,'r'])  #CJB+
        elif .8<=i<1:  #CJB+
            step=np.array([count_t,'m'])  #CJB+
        elif .6<=i<.8:  #CJB+
            step=np.array([count_t,'y'])  #CJB+
        elif .4<=i<.6:  #CJB+
            step=np.array([count_t,'g'])  #CJB+
        elif .2<=i<.4:  #CJB+
            step=np.array([count_t,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_t,'k'])  #CJB+
        count_t+=1  #CJB+
        t_util_list=np.vstack((t_util_list, step))  #CJB+
    t_util_list=np.delete(t_util_list,0,0)  #CJB+
    #print t_util_list  #CJB+

    cb_util_list=np.array([[1,1]])  #CJB+
    count_cb=1  #CJB+
    for j in myjckt.LoadFrameOuts.jacket_utilization.cb_util:  #CJB+
        if j>=1:  #CJB+
            step=np.array([count_cb,'r'])  #CJB+
        elif .8<=j<1:  #CJB+
            step=np.array([count_cb,'m'])  #CJB+
        elif .6<=j<8:  #CJB+
            step=np.array([count_cb,'y'])  #CJB+
        elif .4<=j<6:  #CJB+
            step=np.array([count_cb,'g'])  #CJB+
        elif .2<=j<.4:  #CJB+
            step=np.array([count_cb,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_cb,'k'])  #CJB+
        count_cb+=1  #CJB+
        cb_util_list=np.vstack((cb_util_list, step))  #CJB+
    cb_util_list=np.delete(cb_util_list,0,0)  #CJB+
    #print cb_util_list  #CJB+

    tower_Stressutil_list=np.array([[1,1]])  #CJB+
    count_tower=1  #CJB+
    for k in myjckt.LoadFrameOuts.tower_utilization.StressUtil:  #CJB+
        if k>=1:  #CJB+
            step=np.array([count_tower,'r'])  #CJB+
        elif .8<=k<1:  #CJB+
            step=np.array([count_tower,'m'])  #CJB+
        elif .6<=k<.8:  #CJB+
            step=np.array([count_tower,'y'])  #CJB+
        elif .4<=k<.6:  #CJB+
            step=np.array([count_tower,'g'])  #CJB+
        elif .2<=k<.4:  #CJB+
            step=np.array([count_tower,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_tower,'k'])  #CJB+
        count_tower+=1  #CJB+
        tower_Stressutil_list=np.vstack((tower_Stressutil_list, step))  #CJB+
    tower_Stressutil_list=np.delete(tower_Stressutil_list,0,0)  #CJB+
    #print tower_Stressutil_list  #CJB+

    tower_GLutil_list=np.array([[1,1]])  #CJB+
    count_tower=1  #CJB+
    for k in myjckt.LoadFrameOuts.tower_utilization.GLUtil:  #CJB+
        if k>=1:  #CJB+
            step=np.array([count_tower,'r'])  #CJB+
        elif .8<=k<1:  #CJB+
            step=np.array([count_tower,'m'])  #CJB+
        elif .6<=k<.8:  #CJB+
            step=np.array([count_tower,'y'])  #CJB+
        elif .4<=k<.6:  #CJB+
            step=np.array([count_tower,'g'])  #CJB+
        elif .2<=k<.4:  #CJB+
            step=np.array([count_tower,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_tower,'k'])  #CJB+
        count_tower+=1  #CJB+
        tower_GLutil_list=np.vstack((tower_GLutil_list, step))  #CJB+
    tower_GLutil_list=np.delete(tower_GLutil_list,0,0)  #CJB+
    #print tower_GLutil_list  #CJB+

    tower_EUshutil_list=np.array([[1,1]])  #CJB+
    count_tower=1  #CJB+
    for k in myjckt.LoadFrameOuts.tower_utilization.EUshUtil:  #CJB+
        if k>=1:  #CJB+
            step=np.array([count_tower,'r'])  #CJB+
        elif .8<=k<1:  #CJB+
            step=np.array([count_tower,'m'])  #CJB+
        elif .6<=k<.8:  #CJB+
            step=np.array([count_tower,'y'])  #CJB+
        elif .4<=k<.6:  #CJB+
            step=np.array([count_tower,'g'])  #CJB+
        elif .2<=k<.4:  #CJB+
            step=np.array([count_tower,'b'])  #CJB+
        else:  #CJB+
            step=np.array([count_tower,'k'])  #CJB+
        count_tower+=1  #CJB+
        tower_EUshutil_list=np.vstack((tower_EUshutil_list, step))  #CJB+
    tower_EUshutil_list=np.delete(tower_EUshutil_list,0,0)  #CJB+
    #print tower_EUshutil_list  #CJB+

    t_condensed=np.array([[1,1]])  #CJB+
    for m in range(len(t_util_list)):  #CJB+
        counter=(m/8+1)  #CJB+
        if t_util_list[m][1]=='r':  #CJB+
            step=np.array([counter, 'r'])  #CJB+
        elif t_util_list[m][1]=='m':  #CJB+
            step=np.array([counter, 'm'])  #CJB+
        elif t_util_list[m][1]=='y':  #CJB+
            step=np.array([counter, 'y'])  #CJB+
        elif t_util_list[m][1]=='g':  #CJB+
            step=np.array([counter, 'g'])  #CJB+
        elif t_util_list[m][1]=='b':  #CJB+
            step=np.array([counter, 'b'])  #CJB+
        else:  #CJB+
            step=np.array([counter, 'k'])  #CJB+
        t_condensed=np.vstack((t_condensed,step))  #CJB+
    t_condensed=np.delete(t_condensed,0,0)  #CJB+
    t_condensed=t_condensed[0:-1:8]  #CJB+
    #print t_condensed  #CJB+

    cb_condensed=np.array([[1,1]])  #CJB+
    for m in range(len(cb_util_list)):  #CJB+
        counter=(m/4+1)  #CJB+
        if cb_util_list[m][1]=='r':  #CJB+
            step=np.array([counter, 'r'])  #CJB+
        elif cb_util_list[m][1]=='m':  #CJB+
            step=np.array([counter, 'm'])  #CJB+
        elif cb_util_list[m][1]=='y':  #CJB+
            step=np.array([counter, 'y'])  #CJB+
        elif cb_util_list[m][1]=='g':  #CJB+
            step=np.array([counter, 'g'])  #CJB+
        elif cb_util_list[m][1]=='b':  #CJB+
            step=np.array([counter, 'b'])  #CJB+
        else:  #CJB+
            step=np.array([counter, 'k'])  #CJB+
        cb_condensed=np.vstack((cb_condensed,step))  #CJB+
    cb_condensed=np.delete(cb_condensed,0,0)  #CJB+
    cb_condensed=cb_condensed[0:-1:4]  #CJB+
    #print cb_condensed  #CJB+

    jacket_check=np.array([[1,1]])  #CJB+
    for n in range(len(t_condensed)):  #CJB+
        counter=n+1  #CJB+
        if t_condensed[n][1]=='r' or cb_condensed[n][1]=='r':  #CJB+
            step=np.array([counter, 'r'])  #CJB+
        elif t_condensed[n][1]=='m' or cb_condensed[n][1]=='m':  #CJB+
            step=np.array([counter, 'm'])  #CJB+
        elif t_condensed[n][1]=='y' or cb_condensed[n][1]=='y':  #CJB+
            step=np.array([counter, 'y'])  #CJB+
        elif t_condensed[n][1]=='g' or cb_condensed[n][1]=='g':  #CJB+
            step=np.array([counter, 'g'])  #CJB+
        elif t_condensed[n][1]=='b' or cb_condensed[n][1]=='b':  #CJB+
            step=np.array([counter, 'b'])  #CJB+
        else:  #CJB+
            step=np.array([counter, 'k'])  #CJB+
        jacket_check=np.vstack((jacket_check,step))  #CJB+
    jacket_check=np.delete(jacket_check,0,0)  #CJB+
    #print jacket_check  #CJB+

    tower_check=np.array([[1,1]])  #CJB+
    for p in range(len(tower_Stressutil_list)):  #CJB+
        counter=p+len(jacket_check)+1  #CJB+
        if tower_Stressutil_list[p][1]=='r' or tower_GLutil_list[p][1]=='r' or tower_EUshutil_list[p][1]=='r':  #CJB+
            step=np.array([counter, 'r'])  #CJB+
        elif tower_Stressutil_list[p][1]=='m' or tower_GLutil_list[p][1]=='m' or tower_EUshutil_list[p][1]=='m':  #CJB+
            step=np.array([counter, 'm'])  #CJB+
        elif tower_Stressutil_list[p][1]=='y' or tower_GLutil_list[p][1]=='y' or tower_EUshutil_list[p][1]=='y':  #CJB+
            step=np.array([counter, 'y'])  #CJB+
        elif tower_Stressutil_list[p][1]=='g' or tower_GLutil_list[p][1]=='g' or tower_EUshutil_list[p][1]=='g':  #CJB+
            step=np.array([counter, 'g'])  #CJB
        elif tower_Stressutil_list[p][1]=='b' or tower_GLutil_list[p][1]=='b' or tower_EUshutil_list[p][1]=='b':  #CJB+
            step=np.array([counter, 'b'])  #CJB+
        else:  #CJB+
            step=np.array([counter, 'k'])  #CJB+
        tower_check=np.vstack((tower_check,step))  #CJB+
    tower_check=np.delete(tower_check,0,0)  #CJB+
    #print tower_check  #CJB+

    member_check=np.vstack((jacket_check,tower_check))  #CJB+
    #print member_check  #CJB+

    Kjnts_check=np.array([[1,1,1]])  #CJB+
    for i in myjckt.Build.KjntIDs:  #CJB+
        Kjnts_check=np.vstack((Kjnts_check,myjckt.JcktGeoOut.nodes[i-1]))  #CJB+
    Kjnts_check=np.delete(Kjnts_check,0,0)  #CJB+
    #print Kjnts_check  #CJB+
    #Kjnts_check=np.c_[Kjnts_check,myjckt.LoadFrameOuts.jacket_utilization.KjntUtil]  #CJB+
    #print Kjnts_check  #CJB+
    K_count=1  #CJB+
    K_joints_color=np.array([[1,1]])  #CJB+
    for i in myjckt.LoadFrameOuts.jacket_utilization.KjntUtil:  #CJB+
        if i>=1:  #CJB+
            step=np.array([K_count,'r^'])  #CJB+
        elif .8<=i<1:  #CJB+
            step=np.array([K_count,'mo'])  #CJB+
        elif .6<=i<.8:  #CJB+
            step=np.array([K_count,'yo'])  #CJB+
        elif .4<=i<.6:  #CJB+
            step=np.array([K_count,'go'])  #CJB+
        elif .2<=i<.4:  #CJB+
            step=np.array([K_count,'bo'])  #CJB+
        else:  #CJB+
            step=np.array([K_count,'ko'])  #CJB+
        K_count+=1  #CJB+
        K_joints_color=np.vstack((K_joints_color,step))  #CJB+
    K_joints_color=np.delete(K_joints_color,0,0)  #CJB+
    #print K_joints_color  #CJB+

    Xjnts_check=np.array([[1,1,1]])  #CJB+
    for i in myjckt.Build.XjntIDs:  #CJB+
        Xjnts_check=np.vstack((Xjnts_check,myjckt.JcktGeoOut.nodes[i-1]))  #CJB+
    Xjnts_check=np.delete(Xjnts_check,0,0)  #CJB+
    #print Xjnts_check  #CJB+

    X_count=1  #CJB+
    X_joints_color=np.array([[1,1]])  #CJB+
    for i in myjckt.LoadFrameOuts.jacket_utilization.XjntUtil:  #CJB+
        if i>=1:  #CJB+
            step=np.array([X_count,'r^'])  #CJB+
        elif .8<=i<1:  #CJB+
            step=np.array([X_count,'mo'])  #CJB+
        elif .6<=i<.8:  #CJB+
            step=np.array([X_count,'yo'])  #CJB+
        elif .4<=i<.6:  #CJB+
            step=np.array([X_count,'go'])  #CJB+
        elif .2<=i<.4:  #CJB+
            step=np.array([X_count,'bo'])  #CJB+
        else:  #CJB+
            step=np.array([X_count,'ko'])  #CJB+
        X_count+=1  #CJB+
        X_joints_color=np.vstack((X_joints_color,step))  #CJB+
    X_joints_color=np.delete(X_joints_color,0,0)  #CJB+
    #print X_joints_color  #CJB+

    KXcoords=np.vstack((Kjnts_check,Xjnts_check))  #CJB+
    KXcolors=np.vstack((K_joints_color,X_joints_color))  #CJB+
    #print KXcoords  #CJB+
    #print KXcolors  #CJB+

    #toscreen=False allows peregrine to work
    fig = plt.figure(1)
    ax = Axes3D(fig)
    nodes=myjckt.JcktGeoOut.nodes  #coordinates
    mems=myjckt.JcktGeoOut.mems
    XYZ1=nodes[mems[:,0]-1,:]
    XYZ2=nodes[mems[:,1]-1,:]
    #ax.scatter(nodes[:,0],nodes[:,1],nodes[:,2])
    x,y,z=zip(*nodes)
    #ax.scatter(x,y,z)
    Xs=np.vstack((XYZ1[:,0],XYZ2[:,0])).T
    Ys=np.vstack((XYZ1[:,1],XYZ2[:,1])).T
    Zs=np.vstack((XYZ1[:,2],XYZ2[:,2])).T
    #ax.plot([XYZ1[1:5,0],XYZ2[1:5,0]],[XYZ1[1:5,1],XYZ2[1:5,1]],[XYZ1[1:5,2],XYZ2[1:5,2]])
    for i in range(0,mems.shape[0]): #mems.shape[0])
        ax.plot([XYZ1[i,0],XYZ2[i,0]],[XYZ1[i,1],XYZ2[i,1]],[XYZ1[i,2],XYZ2[i,2]],member_check[i][1]) #CJBe Member check gives memebrs utilization colors
        #CJB TODO If RigidTop=False, then the code does not make the top member. However, the utilization plot (Figure 2) will still show the utiliazation of the top member. (There is a different number of tower members being plotted)
    axisEqual3D(ax)
#    ax.set_aspect('equal')
#    ax.auto_scale_xyz([min(XYZ1[:,0]),max(XYZ1[:,0])],[min(XYZ1[:,1]),max(XYZ1[:,1])],[min(XYZ1[:,2]),max(XYZ1[:,2])])

    if savefileroot:
        plt.savefig(savefileroot+'_config.png',format='png')
    elif toscreen:
        plt.show()

    #Plot utilization of Tower if requested
    if util:
        RigidTop=(myjckt.RNAinputs.CMoff[2] != 0.) and myjckt.Tower.RigidTop
        fig2=plt.figure(2);
        ax1 = fig2.add_subplot(111)
        twr_zs=myjckt.Tower.Twrouts.nodes[2,0:myjckt.Tower.Twrouts.nNodes-RigidTop]-myjckt.Tower.Twrouts.nodes[2,0]
        #Add one diameter at the very top since the tube object does not have the final D
        twr_D=np.hstack((myjckt.Tower.Twrouts.TwrObj.D,myjckt.Tower.Twrins.Dt))

        yrange=(np.min(twr_zs),np.max(twr_zs));
        #hh=plt.Axes(ox','off',\
        ##    'Ylim',yrange,'Nextplot','Add','Visible','On');
        ax1.set_xlabel('VonMises, GL, and EU Utilization Ratios');
        ax1.set_ylabel('z from base of tower (incl. AF) [m]');
        #ax1.set_ylabel('z from base of structure (incl. AF) [m]');
        ax1.set_title('Utilization Ratio');
        ax1.set_xlim([0,2]);
        ax1.set_ylim(yrange);
        #set(plt.gca(),'ylim',yrange)
        #ax1.plot([0, 2], [water_dict['wlevel'],water_dict['wlevel']],'.-r')
        #ax1.plot([0, 2], np.array([water_dict['wlevel']-water_dict['wdpth']]).repeat(2),':k')

        ax1.plot([],[],'k',label='Tower Profile') #CJB+
        ax1.plot(myjckt.LoadFrameOuts.tower_utilization.StressUtil,twr_zs , label='VonMises Util1')
        ax1.plot(myjckt.LoadFrameOuts.tower_utilization.GLUtil,twr_zs, label='GL Util1')
        ax1.plot(myjckt.LoadFrameOuts.tower_utilization.EUshUtil, twr_zs, label='EUsh Util1')
        if myjckt.LoadFrameOuts2.tower_utilization.StressUtil[0] != -9999.:
            ax1.plot([],[],'k',label='Tower Profile') #CJB+
            ax1.plot(myjckt.LoadFrameOuts2.tower_utilization.StressUtil,twr_zs , label='VonMises Util2')
            ax1.plot(myjckt.LoadFrameOuts2.tower_utilization.GLUtil,twr_zs, label='GL Util2')
            ax1.plot(myjckt.LoadFrameOuts2.tower_utilization.EUshUtil, twr_zs, label='EUsh Util2')

        ax1.grid()
        #ax1.legend(['0 m MSL','SeaBed','GL Global Buckling', 'EU Shell Buckling'])
        #ax1.legend(['GL Global Buckling', 'EU Shell Buckling'])
        ax1.legend(bbox_to_anchor=(0.6, 1.0), loc=2)

        ax2=plt.twiny(ax1)# .axes(ax1.get_position(),frameon=False)
        #hh2=plt.axes('Position', pos,'NextPlot','Add','XtickLabel','','Xtick',[],frameon=False);
        ax2.plot(twr_D/2,twr_zs,'k') #CJBe Add color
        ax2.plot(-twr_D/2,twr_zs, 'k') #CJBe Add color
        ax2.set_aspect('equal')
        ax2.set_frame_on(False)
        ax2.set_xticklabels('')
        ax2.set_xticks([])
        ax1.set_xlim([0,2]);
        #hh2.set_aspect('equal')
        ##ax2.axis('equal')

        if savefileroot:
            plt.savefig(savefileroot+'_util.png',format='png')
        elif toscreen:
            plt.show()


    fig3= plt.figure(3) #CJB+
    ax = Axes3D(fig3) #CJB+

    for i in range(myjckt.JcktGeoOut.nmems-len(myjckt.LoadFrameOuts.TPmems)-len(myjckt.LoadFrameOuts.Twrmems)): #CJB+
        #if count<=myjckt.JcktGeoOut.nmems:
        plt.plot([XYZ1[i,0],XYZ2[i,0]],[XYZ1[i,1],XYZ2[i,1]],[XYZ1[i,2],XYZ2[i,2]],'.75') #CJBe Member check gives memebrs utilization colors
    for i in range(0,KXcoords.shape[0]): #CJB+
        plt.plot([KXcoords[i,0]],[KXcoords[i,1]],[KXcoords[i,2]], KXcolors[i][1]) #CJB+
    axisEqual3D(ax) #CJB+
    plt.show() #CJB+
#        plt.close(all)



if __name__ == '__main__':
                             #batter,Dpile,tpile,  Lp   ,Dleg,      tleg,  Dbrc, tbrc,  Dbrc_mud,tbrc_mud   ,Db,  DTRb, Dt, DTRt,   Htwr2frac
    myjckt=SetJacketInputs.main(10.,2.0, 0.035,    40.,   2.2,  0.04,     1.2,  0.03,   1.2,     0.03,       6.5,130.,   4.,130.,0.2)
    main(myjacket)
