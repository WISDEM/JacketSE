#-------------------------------------------------------------------------------
# Name:        printJacketres.py
# Purpose:     Prints data (often results) for a jacket
#
# Author:      rdamiani
#
# Created:     19/01/2015
# Copyright:   (c) rdamiani 2015
# Licence:     Apache 2015
#-------------------------------------------------------------------------------
import numpy as np

def main(myjckt):
    print "Dpile=%f, tpile=%f  Lp=%f " % (myjckt.Piles.Pileinputs.Dpile,myjckt.Piles.Pileinputs.tpile,myjckt.Piles.Pileinputs.Lp)
    print "Dbrc=%f, tbrc=%f, beta2D=%f  " % (myjckt.Xbraces.Xbrcouts.LLURObj.D[0],myjckt.Xbraces.Xbrcouts.LLURObj.t[0],myjckt.beta2D)
    print "Dbrcmud=%f, tbrcmud=%f  " % (myjckt.Mudbraces.Mbrcouts.brcObj.D[0],myjckt.Mudbraces.Mbrcouts.brcObj.t[0])
    print "batter=%f, dckwidth=%f, footprint=%f, Dleg=%f, tleg=%f,  " % (myjckt.JcktGeoIn.batter,myjckt.dck_width, myjckt.PreBuild.wbase, myjckt.leginputs.Dleg[0],myjckt.leginputs.tleg[0])
    print "Dgir=%f, tgir=%f " % (myjckt.TPinputs.Dgir,myjckt.TPinputs.tgir)
    print "Db=%f DTRb=%f Dt=%f DTRt=%f H2frac=%f " % (myjckt.Tower.Twrins.Db,myjckt.Tower.Twrins.DTRb,myjckt.Tower.Twrins.Dt,myjckt.Tower.Twrins.DTRt,myjckt.Tower.Twrins.Htwr2frac)
    print "\n"
    # print component masses
    print('jacket+TP(structural+lumped) mass (no tower, no piles) [kg] = {:6.0f}'.format(myjckt.LoadFrameOuts.Frameouts.mass[0]+myjckt.TP.TPlumpinputs.mass-myjckt.Tower.Twrouts.mass))
    print('tower mass [kg] = {:6.0f}'.format(myjckt.Tower.Twrouts.mass))
    print('TP mass structural + lumped mass [kg] = {:6.0f}'.format(myjckt.TP.TPouts.mass+myjckt.TP.TPlumpinputs.mass))
    print('piles (all) mass for assigned (not optimum, unless optimization is run) Lp [kg] = {:6.0f}'.format(myjckt.LoadFrameOuts.Mpiles))
    print('frame3dd model mass (structural + TP lumped) [kg] = {:6.0f}'.format(myjckt.LoadFrameOuts.Frameouts.mass[0]+myjckt.TP.TPlumpinputs.mass))
    print('frame3dd model mass (structural + TP lumped) + Pile Mass [kg] = {:6.0f}'.format(myjckt.LoadFrameOuts.Frameouts.mass[0]+myjckt.TP.TPlumpinputs.mass+myjckt.Mpiles))
    print "\n"
    # modal analysis
    print('First two Freqs.= {:5.4f} and {:5.4f} Hz \n'.format(*myjckt.LoadFrameOuts.Frameouts.Freqs))
    #print tower top displacement
    print('Tower Top Displacement in Global Coordinate System DLC1.6 [m] ={:5.4f}'.format(myjckt.LoadFrameOuts.Frameouts.top_deflection[0]))
    print('Tower Top Displacement in Global Coordinate System DLC6.1 [m] ={:5.4f}'.format(myjckt.LoadFrameOuts2.Frameouts.top_deflection[0]))

    print "\n"
    #print Pile Lp0rat
    print ('Lp0 ={:5.4f} and Lp0rat={:5.4f} DLC 1.6').format(myjckt.Piles.Pileinputs.Lp/(1.+myjckt.LoadFrameOuts.Lp0rat),myjckt.LoadFrameOuts.Lp0rat)
    print ('Lp0 ={:5.4f} and Lp0rat={:5.4f} DLC 6.1').format(myjckt.Piles.Pileinputs.Lp/(1.+myjckt.LoadFrameOuts2.Lp0rat),myjckt.LoadFrameOuts2.Lp0rat)

    print "\n"
    #print GL EU utilizations
    print "GLutil=%f EUutil=%f DLC1.6"  % (np.nanmax(myjckt.LoadFrameOuts.tower_utilization.GLUtil),np.nanmax(myjckt.LoadFrameOuts.tower_utilization.EUshUtil))
    print "GLutil=%f EUutil=%f DLC6.1"  % (np.nanmax(myjckt.LoadFrameOuts2.tower_utilization.GLUtil),np.nanmax(myjckt.LoadFrameOuts2.tower_utilization.EUshUtil))
    print "Mudline Footprint=%f"  % (myjckt.PreBuild.wbase)
    #print max API code checks
    print('MAX member compression-bending utilization at joints DLC1.6 = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts.jacket_utilization.cb_util)))
    print('MAX member compression-bending utilization at joints DLC6.1 = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts2.jacket_utilization.cb_util)))
    print('MAX member tension utilization at joints DLC1.6 = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts.jacket_utilization.t_util)))
    print('MAX member tension utilization at joints DLC6.1 = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts2.jacket_utilization.t_util)))
    print('MAX X-joint  utilization at joints DLC1.6 = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts.jacket_utilization.XjntUtil)))
    print('MAX X-joint  utilization at joints DLC6.1 = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts2.jacket_utilization.XjntUtil)))
    print('MAX K-joint  utilization at joints DLC1.6 = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts.jacket_utilization.KjntUtil)))
    print('MAX K-joint  utilization at joints DLC6.1 = {:5.4f}'.format(np.max(myjckt.LoadFrameOuts2.jacket_utilization.KjntUtil)))

if __name__ == '__main__':
    main()



