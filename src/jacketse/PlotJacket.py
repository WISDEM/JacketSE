#-------------------------------------------------------------------------------
# Name:        PlotJacket.Py
# Purpose:
#
# Author:      rdamiani
#
# Created:     14/05/2014
# Copyright:   (c) rdamiani 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np
from commonse.axisEqual3D import axisEqual3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main(myjckt,util=False,savefileroot=[]):

    fig = plt.figure()
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
        ax.plot([XYZ1[i,0],XYZ2[i,0]],[XYZ1[i,1],XYZ2[i,1]],[XYZ1[i,2],XYZ2[i,2]])
    axisEqual3D(ax)
#    ax.set_aspect('equal')
#    ax.auto_scale_xyz([min(XYZ1[:,0]),max(XYZ1[:,0])],[min(XYZ1[:,1]),max(XYZ1[:,1])],[min(XYZ1[:,2]),max(XYZ1[:,2])])

    if savefileroot:
        plt.savefig(savefileroot+'_config.png',format='png')
    else:
        plt.show()

    #Plot utilization of Tower if requested
    if util:
        RigidTop=(myjckt.RNAinputs.CMoff[2] != 0.) and myjckt.Tower.RigidTop
        fig2=plt.figure();
        ax1 = fig2.add_subplot(111)
        twr_zs=myjckt.Tower.Twrouts.nodes[2,0:myjckt.Tower.Twrouts.nNodes-RigidTop]-myjckt.Tower.Twrouts.nodes[2,0]
        #Add one diameter at the very top since the tube object does not have the final D
        twr_D=np.hstack((myjckt.Tower.Twrouts.TwrObj.D,myjckt.Tower.Twrins.Dt))

        yrange=(np.min(twr_zs),np.max(twr_zs));
        #hh=plt.Axes(ox','off',\
        ##    'Ylim',yrange,'Nextplot','Add','Visible','On');
        ax1.set_xlabel('GL and EU Utilization Ratios');
        ax1.set_ylabel('z from base of tower (incl. AF) [m]');
        #ax1.set_ylabel('z from base of structure (incl. AF) [m]');
        ax1.set_title('Utilization Ratio');
        ax1.set_xlim([0,2]);
        ax1.set_ylim(yrange);
        #set(plt.gca(),'ylim',yrange)
        #ax1.plot([0, 2], [water_dict['wlevel'],water_dict['wlevel']],'.-r')
        #ax1.plot([0, 2], np.array([water_dict['wlevel']-water_dict['wdpth']]).repeat(2),':k')

        ax1.plot(myjckt.LoadFrameOuts.tower_utilization.StressUtil,twr_zs , label='VonMises Util')
        ax1.plot(myjckt.LoadFrameOuts.tower_utilization.GLUtil,twr_zs, label='GL Util')
        ax1.plot(myjckt.LoadFrameOuts.tower_utilization.EUshUtil, twr_zs, label='EUsh Util')
        if myjckt.LoadFrameOuts2.tower_utilization.StressUtil[0] != -9999.:
            ax1.plot(myjckt.LoadFrameOuts2.tower_utilization.StressUtil,twr_zs , label='VonMises Util2')
            ax1.plot(myjckt.LoadFrameOuts2.tower_utilization.GLUtil,twr_zs, label='GL Util2')
            ax1.plot(myjckt.LoadFrameOuts2.tower_utilization.EUshUtil, twr_zs, label='EUsh Util2')

        ax1.grid()
        #ax1.legend(['0 m MSL','SeaBed','GL Global Buckling', 'EU Shell Buckling'])
        #ax1.legend(['GL Global Buckling', 'EU Shell Buckling'])
        ax1.legend(bbox_to_anchor=(1.05, 1.0), loc=2)

        ax2=plt.twiny(ax1)# .axes(ax1.get_position(),frameon=False)
        #hh2=plt.axes('Position', pos,'NextPlot','Add','XtickLabel','','Xtick',[],frameon=False);
        ax2.plot(twr_D/2,twr_zs);
        ax2.plot(-twr_D/2,twr_zs);
        ax2.set_aspect('equal')
        ax2.set_frame_on(False)
        ax2.set_xticklabels('')
        ax2.set_xticks([])
        ax1.set_xlim([0,2]);
        #hh2.set_aspect('equal')
        ##ax2.axis('equal')

        if savefileroot:
            plt.savefig(savefileroot+'_util.png',format='png')
        else:
            plt.show()


if __name__ == '__main__':
                             #batter,Dpile,tpile,  Lp   ,Dleg,      tleg,  Dbrc, tbrc,  Dbrc_mud,tbrc_mud   ,Db,  DTRb, Dt, DTRt,   Htwr2frac
    myjckt=SetJacketInputs.main(10.,2.0, 0.035,    40.,   2.2,  0.04,     1.2,  0.03,   1.2,     0.03,       6.5,130.,   4.,130.,0.2)
    main(myjacket)
