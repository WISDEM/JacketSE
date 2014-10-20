.. _tutorial-label: 5 MW Jacket Tower in 30 m Waters

.. currentmodule:: rotorse.rotoraero

Tutorial
--------

You can refer to classes (:class:`AeroBase`), or methods (:meth:`CCAirfoil.initFromAerodynFile`) assuming you have set the module with ..currentmodule or ..module.  You can also refer to modules :mod:`rotorse.rotoraero`.

You might want to include a figure (:num:`Figure #somelabel-fig`)

.. _somelabel-fig:

.. figure:: /images/figurename.*
    :height: 4in
    :align: center

    Caption goes here

You can also include code from an example.  This is preferable to writing the actual code here in this page, because then you can run the example and test it.

.. literalinclude:: ../src/rotorse/rotoraerodefaults.py
    :start-after: # --- rotor geometry
    :end-before: # ---

Print out some results from the code

>>> CP = [0.48329808]
>>> CT = [0.7772276]
>>> CQ = [0.06401299]

Tutorial
--------

Two examples are included in this tutorial section: simulation of a land-based tower, and optimization of a land-based tower.

5-MW Jacket Tower in 30 m Waters
================================

.. currentmodule:: jacketse.jacket

This example demonstrates how to setup and run *analysis* for a jacket-tower support structure.  

First, we set up the jacket geometry paramters and instantiate the main jacket assembly:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: # --- Set Jacket Parameters ---#
    :end-before: #Soil inputs
 
Then we set up various component inputs. Start with the soil stratigraphy:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Soil inputs
    :end-before: #Water and wind inputs
    
Then assign water and wind environmental parameters:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #ater and wind inputs
    :end-before: #RNA loads

Then assign the RNA aerodynamic loads:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #RNA loads
    :end-before: #Legs data

Now it is time to assign properties for the various member classes of the structure. 

We will follow a bottom up approach, where the definitions are offered from the bottom members towards the top ones.


Start with the piles:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Pile data
    :end-before: #Legs data

Then continue with the main legs:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Legs data
    :end-before: #Xbrc data

Then move on to the X-braces:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Xbrc data
    :end-before: #Mbrc data

In this example we are considering the Mud-braces. They normally mitigate stress concentrations at the pile-leg connections and provide some torsional stiffeness:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Mbrc data
    :end-before: #Hbrc data

In this case we are also considering a top horizontal brace below the main perimeter girder of the TP:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Mbrc data
    :end-before: #TP data

The next step consists of declaring the transition piece properties.

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #TP data
    :end-before: #RNA data


Then the tower geometry is assigned:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Tower data
    :end-before: #Frame3DD 
    
Finally auxiliary parameters for the Frame3DD solver may be assigned:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Frame3DD
    :end-before: #--- RUN Jacket


`JacketSE` was designed to be modular, and so the specific modules you wish to use in the analysis must be specified. 
There are `OpenMDAO Slots <http://openmdao.org/docs/basics/variables.html#slot-variables>`_ for wind1, wind2, wave1, wave2, soil, tower1, and tower2.  The first five make use of the commonse.environment module.  We use multiple wind/wave/tower modules because we are considering two separate loading conditions in this simulation.  The slots wind1 and wind2 can be any component that inherits from :class:`WindBase`, wave1 and wave2 require components that inherit from :class:`WaveBase`, and soil uses components that inherit from :class:`SoilBase`.

,which defines a power-profile for the wind distribution.  A component, :class:`LogWind` for a logarithmic profile is also available.  We are simulating a land-based turbine so we do not load any wave modules.  The default :class:`WaveBase` module is for no wave loading.  A :class:`LinearWaves` component is available which uses linear wave theory.  A simple textbook-based soil model is provided at :class:`TowerSoil`.  For all slots, users may define any custom component as desired.

.. currentmodule:: jacketse.jacket

Start by importing selecting the main parameters for the jacket geometry, In this case we are using :class:`JcktGeoInputs()` .

Tor tower1 and tower2, a component that implements :class:`TowerBase` must be used.  Two default implementations are provided.  The first, :class:`TowerWithpBEAM`, uses the beam finite element code `pBEAM <https://github.com/WISDEM/pBEAM>`_.  The second, :class:`TowerWithFrame3DD` uses the frame finite element code `Frame3DD <http://frame3dd.sourceforge.net/>`_.

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: # --- Set Jacket Parameters ---#
    :end-before: #--- RUN JACKET ---#

With the jacket configuration setup, we define some of the geometric parameters.  Some of these parameters are specific to the module we have chosen to load, :class:`TowerWithpBEAM`, and may differ for other modules.  Some of the geometric parameters are seen in :num:`Figure #tower-fig`.  The tower is not restricted to 3 sections, any number of sections can be defined.  The array `z` is given in coordinates nondimensionalized by the tower height.  The array `n`, should of length len(tower.z)-1 and represents the number of finite elements to be used in each tower can.  The float `L_reinforced` is a reinforcement length used in the buckling calculations.  Yaw and tilt are needed to handle to mass/load transfer.  For offshore applications, monopile geometry can also be defined (see :class:`TowerSE`).

.. _tower-fig:

.. figure:: /images/tower.*
    :height: 4in
    :align: center

    Example of tower geometric parameterization.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- geometry
    :end-before: # ---


We now define the mass properties for the rotor-nacelle-assembly.  The center of mass locations are defined relative to the tower top in the yaw-aligned coordinate system.  Blade and hub moments of inertia should be defined about the hub center, nacelle moments of inertia are defined about the center of mass of the nacelle.


Environmental properties are defined below, note that the parameters that need to be defined depend on which modules were loaded.  For the power-law wind profile, the only parameter needed is the shear exponent.  For the soil, shear and modulus properties for the soil can be defined, but in this example we assume that all directions are rigid (3 translation and 3 rotation).  In addition, some geometric parameters for the wind profile's extend must be defined, the base (or no-slip location) at `z0`, and the height at which a reference velocity will be defined.



As mentioned earlier, we are allowing for two separate loading cases.  The wind speed, and rotor force/moments for those two cases are now defined.  The wind speed location corresponds to the reference height defined previously as `wind_zref`.  In this simple case, we include only thrust and torque, but in general all 3 components of force and moments can be defined in the hub-aligned coordinate system.  The assembly automatically handles translating the forces and moments defined at the rotor to the tower top.


Safety factors for loading, material, consequence of failure, and buckling are defined


A simplified fatigue analysis is available for the tower.  This requires running an aeroelastic code, like FAST, before hand and inputing the damage equivalent moments.  The locations of the moments are given on a nondimensional tower.  A safety factor, lifetime (in years), and slope of the S-N curve can be defined.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- fatigue
    :end-before: # ---

Finally, some additional parameters used for constraints can be defined.  These include the minimum allowable taper ratio of the tower (from base to top), and the minimum diameter-to-thickness ratio allowed at any section.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- constraints
    :end-before: # ---

In the folllowing specification, we used the default values for wind density and viscosity, material properties (for steel), and acceleration of gravity.  By examining, :class:`TowerSE` the user can see all possible parameters and their defaults.  We can now run the assembly and display some of the outputs.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- run
    :end-before: # ---

The resulting output is shown below.  The outputs f1 and f2 and the first two natural frequencies.  
The outputs top_deflection are the deflection of the tower top in the yaw-aligned +x direction for the two loading cases.  Weldability is an array of outputs that should be <0 for feasibility.  They represent margin against the minimum diameter-to-thickness ratio.  Manufacturability represents margin against the minimum taper ratio and should be <0 for feasibility.


>>> 7 nodes in the constant-OD segment of the tower
('>>>>>>>>>>  needed embedment Lp0=', 25.941872262448666)
>>> First two Freqs.= 0.2298 and 0.2309 Hz
>>> jacket+TP(structural+lumped) mass (no tower, no piles) [kg] = 910177
>>> tower mass [kg] = 342737
>>> TP mass structural + lumped mass [kg] = 428872
>>> piles (all) mass (for assigned (not optimum) Lp [kg] =  48844
>>> frame3dd model mass (structural + TP lumped) [kg] = 1252914


>>> MAX member compression-bending utilization at joints = 0.8486
>>> MAX member tension utilization at joints = 0.9463
>>> MAX X-joint  utilization at joints = 0.2618
>>> MAX K-joint  utilization at joints = 0.8625

>>> mass (kg) = 349486.79362
>>> f1 (Hz) = 0.331531844509
>>> f2 (Hz) = 0.334804545737
>>> top_deflection1 (m) = 0.691606748192
>>> top_deflection2 (m) = 0.708610880714
>>> weldability = [-0.42450142 -0.37541806 -0.30566802]
>>> manufactuability = -0.245

The stress, buckling, and damage loads are shown in :num:`Figure #utilization-fig`.  Each is a utilization and so should be <1 for feasibility.

.. _utilization-fig:

.. figure:: /images/utilization.*
    :width: 6in
    :align: center

    Utilization along tower height for: Von-Mises/yield; shell buckling; global buckling.


Land-Based Tower Optimization
=============================

We begin with the same setup as the previous section, but now import additional modules for optimization.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- optimizer imports
    :end-before: # ---

The optimizer must first be selected and configured, in this example I use SNOPT.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- Setup Pptimizer
    :end-before: # ---

We now set the objective, and in this example it is normalized to be of order 1 for better convergence behavior.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- Objective
    :end-before: # ---

The tower diameters, thickness, and waist location are added as design variables.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- Design Variables
    :end-before: # ---

A recorder is added to display each iteration to the screen.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- recorder
    :end-before: # ---

Finally, constraints are added.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- Constraints
    :end-before: # ---

Now the optimization can be run.

.. literalinclude:: ../src/towerse/tower.py
    :language: python
    :start-after: # --- run opt
    :end-before: # ---








