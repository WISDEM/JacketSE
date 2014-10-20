#INPUT PARAMETERS FOR Jacket/TOWER JckTwr Class (Attributes for Python's Class)
#  JACKET PROPERTIES  (This is exactly as for MAtlab, but with '#' comment character)
# The "(Python)" indicates that the definition is used by Python's Program as well, but not necessarily by Matlab
##Python needs the following
##Jacket prms={'Db':6.7, 'nlegs':4, 'nbays':4, 'batter':7,\
       ##'dck_botz':16, 'wdpth':50, 'pil_TOPz':2, 'leg_BOTz':2, 'weld2D':.2,\
##       'TPmass':666e3, 'TPlth':4, 'gussets':True}
##Twr  twrprms={'Db':6.7, 'Dt':3.4, 'Htwr':80, 'Htwr2':2/25*80, 'DTR':120,\
##             'RNAmass':573800, 'RNA_Ixx':86.579E+6, 'RNA_Iyy':53.530E+6,\
##             'RNA_Izz':58.112E+6 ,'CMzoff':2.34, 'n_elems':10}



JckTwr.Dleg=1.278;      # [m] Leg OD
JckTwr.tleg=0.0254*2;   # [m] Leg thickness
JckTwr.batter=7.;        #[m] apparent 2D batter (Python)
JckTwr.Dpile=JckTwr.Dleg; #[m] pile diameter [m]
JckTwr.tpile=0.014;      #[m] pile thickness [m]
JckTwr.Dbrc=0.66;        #[m] X brace diameter [m]
JckTwr.tbrc=0.014;      #[m] X brace thickness [m]

JckTwr.nbays=4;  #number of bays (Python)
JckTwr.nlegs=4;  #number of legs (Python)
JckTwr.weld2D=.2;  #weldment allowance, i.e. total clearance=(1+weld2D)*D (Python)
JckTwr.gussets=True; #whetehr or not TP is represented with gussets or box (Python)


JckTwr.pile_bat=0.;   #integer,[0/1] flag to state whether piles are vertical (0) or battered (1)
JckTwr.dck_width=13.; #[m] width of platform at joints with legs (usually platform width -2*(weld2D+1)*Dleg) (Python)
JckTwr.wdpth=50.;  #water depth needed by buildMPtower only, loads are given through windwaterinp
JckTwr.pileZtop=1.;#    where the pile head terminates above sea-bed [m]
JckTwr.legZbot=1.; #    where the leg foot terminates above sea-bed [m]

#*** !!! The following property may need to be modified to get the right !!! ***
#         frequencies in BMODES below the seabed level                   !!! ***
JckTwr.r_gyr=27/2.0;  #[m] gyration radius for the jacket  below seabed; use r_gyr=wbase/2 to start
#                               this value will be used to calculate an effective Area
#                               of the cross section, to be varied to match frequency as desired;
#*** !!!                                                                 !!! ***
#               TRANSITION PIECE PROPERTIES
JckTwr.Dstt=1.75;             # strut OD (TP) [m]
JckTwr.tstt=1.375*0.0254;       # strut t (TP) [m]
JckTwr.al_stt=90-46.057;      #strut angle w.r.t. vertical (TP) [deg]
JckTwr.stthgt=5.025*tan(JckTwr.al_stt*pi/180);# height above base of TP to where strutts terminate [m]
JckTwr.TPwidth=2*6.7;         # width of platform deck [m] (square)
JckTwr.TPthk=0.3;             # Thickness of deck [m]
JckTwr.Dtp= 6.7;              # TP cylinder OD, usually set equal to tower base D [m]
JckTwr.ttp=1.5*6.7/120;       # TP cylider thickness, cylidrical part (usually set equal to 1.5*tower base D) [m]
JckTwr.ptfm_ow=16.;            # height above MSL of the platform bottom, i.e. tower base is at ptfm_ow+TPlth [m]
JckTwr.dck_botz=JckTwr.ptfm_ow #Other name for the same variable. Do not modify. (Python)
JckTwr.TPlth=7.3;             # overall length of transition piece [m] (Python)
JckTwr.TPmass=300e3;          # TP mass [kg] (Python)

#TOWER PROPERTIES
JckTwr.Db=JckTwr.Dtp; #[m]   base outer diameter (Python)
JckTwr.Dt=0.55*JckTwr.Db;  # [m] top outer diameter
JckTwr.Htwr=90.4-2;  # Height of the tower
JckTwr.DTR=120;#  # D to thickness ratio
JckTwr.Htwr2=2./25*JckTwr.Htwr; # [m] height of untapered section of tower

#RNA PROPERTIES (PYTHON ONLY)
JckTwr.Twr.RNAmass=1.072E+06;    #RNA mass [kg]
JckTwr.Twr.RNA_Ixx=236.582E+6; # RNA Ixx, [kg*m^2]
JckTwr.Twr.RNA_Iyy=131.783E+6; # RNA Iyy, [kg*m^2]
JckTwr.Twr.RNA_Izz=139.101E+6; # RNA Izz, [kg*m^2]
JckTwr.Twr.CMzoff=2.94;       # RNA CMzoff [m]
JckTwr.Twr.n_elems=20;        # number of elements for tower