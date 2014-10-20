# INPUT PARAMETERS FOR Jacket/TOWER JckTwr Class (Attributes for Python's Class)
# UNITS ARE SI - strictly
#JACKET PROPERTIES  (This is exactly as for MAtlab, but with '#' comment character)

# The "(Python)" indicates that the definition is used by Python's Program as well, but not necessarily by Matlab
##Python needs the following
##Jacket prms={'Db':6.7, 'nlegs':4, 'nbays':4, 'batter':7,\
       ##'dck_botz':16, 'wdpth':50, 'pil_TOPz':2, 'leg_BOTz':2, 'weld2D':.2,\
##       'TPmass':666e3, 'TPlth':4, 'gussets':True}
##Twr  twrprms={'Db':6.7, 'Dt':3.4, 'Htwr':80, 'Htwr2':2/25*80, 'DTR':120,\
##             'RNAmass':573800, 'RNA_Ixx':86.579E+6, 'RNA_Iyy':53.530E+6,\
##             'RNA_Izz':58.112E+6 ,'CMzoff':2.34, 'n_elems':10}
#This is for 5MW turbine  VERTICAL PILES

#TOWER PROPERTIES- FIXED TRUNCATED TOWER
JckTwr.Db=6.00 #         [m]   base outer diameter (Python)
JckTwr.Dt=3.87;        # [m] top outer diameter

#!!!!!!!!!!!!!!!!!!!!!!!!!!TOWER HEIGHT LOOK BELOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#JckTwr.DTR=120;                  # D to thickness ratio
JckTwr.DTRb=180;#JckTwr.Db/0.025;      # D to thickness ratio 6t base
JckTwr.DTRt=180;#JckTwr.Dt/0.019;      # D to thickness ratio at top
JckTwr.Htwr2=0.2; # [m] height of untapered section of tower

JckTwr.Dleg=1.8;           #[m] Leg OD
JckTwr.tleg=0.0254;        #[m] Leg thickness
JckTwr.batter=10;          #[m] apparent 2D batter (Python)
JckTwr.Dpile=JckTwr.Dleg;       #[m] pile diameter [m]
JckTwr.tpile=0.037;       #[m] pile thickness [m]
JckTwr.Dbrc=0.457;         #[m] X brace diameter [m]
JckTwr.tbrc=0.0121        #[m] X brace thickness [m]
JckTwr.Dbrch=0.457;       #[m] 1st horizontal brace diameter [m] (Python)
JckTwr.tbrch=0.0133;       #[m] 1st horizontal brace thickness [m] (Python)

JckTwr.nbays=4;  #number of bays (Python)
JckTwr.nlegs=4;  #number of legs (Python)
JckTwr.weld2D=.2;  #weldment allowance, i.e. total clearance=(1+weld2D)*D (Python)
JckTwr.gussets=True; #whetehr or not TP is represented with gussets or box (Python)


#*** !!! The following property may need to be modified to get the right !!! ***
#         frequencies in BMODES below the seabed level                   !!! ***
JckTwr.r_gyr=27/2.0;  #[m] gyration radius for the jacket  below seabed; use r_gyr=wbase/2 to start
#                               this value will be used to calculate an effective Area
#                               of the cross section, to be varied to match frequency as desired;
#*** !!!                                                                 !!! ***
#               TRANSITION PIECE PROPERTIES
JckTwr.Dstt=JckTwr.Dleg;                   # Strut OD (TP) [m]
JckTwr.tstt=0.0254;                        # Strut t (TP) [m]
JckTwr.al_stt=90.-48.403;             # Strut Angle w.r.t. vertical (TP) [deg] 90-beta_stt in mathcad
JckTwr.stthgt=4.392*tan(JckTwr.al_stt*pi/180.);# height above base of TP to where strutts terminate [m] =l.stt.h * tan(90-beta.stt) in mathcad
JckTwr.TPthk=0.25;                    # Thickness of deck [m], from bottom of stringers to top of plate
JckTwr.Dtp= JckTwr.Db;               # TP cylinder OD, usually set equal to tower base D [m]
JckTwr.ttp=JckTwr.Dtp/JckTwr.DTRb;   # TP cylider thickness, cylidrical part (usually set equal to 1.5*tower base t) [m] for 5 MW set =to tb
JckTwr.TPwidth=2*JckTwr.Dtp;         # Width of platform deck [m] (square), usually  2*Db=2*Dtp
JckTwr.ptfm_ow=-1.06+2.22+2.13+1.5+11.36;     #From GL, and K-13 deep water site                 # Height above MSL of the platform bottom, i.e. tower base is at ptfm_ow+TPlth [m]
JckTwr.dck_botz=JckTwr.ptfm_ow       # other name for the same variable. Do not modify. (Python)
JckTwr.TPlth=6.25;#6.3                    # overall Length of transition piece [m] (Python) (from bottom of stringer to top of cylinder)
JckTwr.TPmass=105.e3#190.e3;#


#RNA PROPERTIES (PYTHON ONLY)
JckTwr.Twr.RNAmass=3.5E+05;    # RNA mass [kg]
JckTwr.Twr.RNA_Ixx=4.37e7;   # RNA Ixx, [kg*m^2]
JckTwr.Twr.RNA_Iyy=2.353E7;   # RNA Iyy, [kg*m^2]
JckTwr.Twr.RNA_Izz=2.542E7;   # RNA Izz, [kg*m^2]
JckTwr.Twr.CMzoff=1.97;          # RNA CMzoff [m] to get to hub-Height=90m, even though the CM is at 1.97m above tower-top, close enough
JckTwr.Twr.RNA_Thrust=1185.E+3;   # UNFactored Rotor max thrust [N] (Python); it will be factored in utilization calc (Python)
JckTwr.Twr.RNA_Myy=0.;           # Rotor Myy moment associated with max thrust [N] (Python)
JckTwr.Twr.RNA_Mxx=0.;           # Rotor Myy moment associated with max thrust [N] (Python)
JckTwr.Twr.n_elems=40;           # number of elements for tower

JckTwr.Htwr=90.-(JckTwr.TPlth+JckTwr.dck_botz)-JckTwr.Twr.CMzoff;   #chopped 5 MW tower, we are assuming CMzoff is at hub-height, which is not quite, but small error  # [m] Height of the tower (PYTHON)

# IMPORTANT!!!!! :Auxiliary data (PYTHON)
JckTwr.pile_bat=0.;    # integer,[0/1] flag to state whether piles are vertical (0) or battered (1) (PYTHON)
JckTwr.dck_width=JckTwr.TPwidth; #[m] width of platform (Python): in Matlab it may be the width at joints with legs (usually platform width -2*(weld2D+1)*Dleg)
JckTwr.wdpth=30.;      #    water depth needed by buildMPtower only, loads are given through windwaterinp
JckTwr.pileZtop=1.;    #    where the pile head terminates above sea-bed [m]
JckTwr.legZbot=1.;     #    where the leg foot terminates above sea-bed [m]
