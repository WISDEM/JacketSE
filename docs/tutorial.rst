.. _tutorial-label: 5 MW Jacket Tower in 30 m Waters

Tutorial
--------

Two examples are included in this tutorial section: simulation of a jacket-tower, and optimization of a jacket-tower support structure.

5-MW Jacket Tower in 30 m Waters
================================

.. currentmodule:: jacketse.jacket

This example demonstrates how to setup and run *analysis* for a jacket-tower support structure.
The structure looks as shown in :num:`Figure jacketTut-fig` :

.. _jacketTut-fig:

.. figure:: /images/jacket_tower_tut.*
    :width: 6in
    :align: center

    Jacket-tower structure for the tutorial example. 

First, we set up the jacket geometry parameters:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #--- Set Jacket Input Parameters ---#
    :end-before: #Soil inputs

Then we set up various component inputs. Start with the soil stratigraphy:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Soil inputs
    :end-before: #Water and wind inputs

Then assign water and wind environmental parameters:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Water and wind inputs
    :end-before: #Pile data

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

In this example we are considering the Mud-braces. They normally mitigate stress concentrations at the pile-leg connections and provide some torsional stiffness:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Mbrc data
    :end-before: #Hbrc data

In this case we are not considering a top horizontal brace below the main perimeter girder of the TP, so we will set its *ndiv* to 0:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Hbrc data
    :end-before: #TP data

The next step consists of declaring the transition piece properties.

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #TP data
    :end-before: #Tower data

Then the tower geometry is assigned:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Tower data
    :end-before: #RNA data

Then the RNA mass properties are specified:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #RNA loads
    :end-before: #Frame3DD

Then assign the RNA aerodynamic loads:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #RNA loads
    :end-before: #Frame3DD

Finally auxiliary parameters for the Frame3DD solver may be assigned:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Frame3DD
    :end-before: #-----Launch the assembly-----#

It is then time to launch the assembly and pass all the inputs to it; note that the assembly is called with parameters depending on the selected inputs:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #-----Launch the assembly-----#
    :end-before: #--- RUN OPTIMIZATION ---#

Then run the assembly:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #--- RUN JACKET ---#
    :end-before: # ---------------

You may print some of the results of the analysis:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #Now show results of modal analysis
    :end-before: #Plot geometry

They should look as shown here:

::

>>>7 nodes in the constant-OD segment of the tower
>>>('>>>>>>>>>>  needed embedment Lp0=', 29.499564957098116)
>>>First two Freqs.= 0.1983 and 0.1994 Hz
>>>jacket+TP(structural+lumped) mass (no tower, no piles) [kg] = 1154582
>>>tower mass [kg] = 335018
>>>TP mass structural + lumped mass [kg] = 335354
>>>piles (all) mass (for assigned (not optimum) Lp [kg] =  48844
>>>frame3dd model mass (structural + TP lumped) [kg] = 1489600
>>>Tower Top Displacement in Global Coordinate System [m] =0.4843
>>>MAX member compression-bending utilization at joints = 0.5608
>>>MAX member tension utilization at joints = 0.6413
>>>MAX X-joint  utilization at joints = 0.2600
>>>MAX K-joint  utilization at joints = 0.3667

If you plot the utilization of the tower, you should get something as in :num:`Figure #utilization-fig`, where
the Von Mises stress, global and shell buckling utilizations are shown along the tower span.  Each curve represents material utilization and so should be <1 for feasibility.

.. _utilization-fig:

.. figure:: /images/util_tut.*
    :width: 6in
    :align: center

    Utilization along tower height for: Von-Mises/yield; shell buckling; global buckling.

.. <!-- ________________________________________ !>

Jacket-Tower Optimization - OpenMDAO Internal Optimization Driver
=================================================================

We begin with the same setup as the previous section, but now we need to set up the optimizer and thus also offer bounds for the design variables:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after:  #Set Optimization Bounds for the various variables:
    :end-before: # --- optimizer imports ---


Also import additional modules for optimization:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: # --- optimizer imports
    :end-before: # ---

The optimizer must first be selected and configured, in this example use SNOPT.

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: # --- Setup Optimizer ---
    :end-before: if SNOPTflag:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: if SNOPTflag:
    :end-before: else:


We now set the objective, and in this example it is normalized to be of order 1 for better convergence behavior.

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: # --- Objective ---
    :end-before: # ----------------------

The batter,pile OD and thickness, Embedment Length, Leg OD and thickness, X-brace OD and thickness,Mud-brace OD and thickness,Tower base OD and DTR, tower-top OD and DTR, and height of constant cross-section segment are added as design variables.

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: # --- Design Variables
    :end-before: #--- Constraints ---#

Constraints are then added; note that we are after a target first eigenfrequeccy of 0.22 Hz:

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #--- Constraints ---#
    :end-before: # ---

A recorder is added to display each iteration to the screen.

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: # --- recorder ---
    :end-before: # ---


Now the optimization can be run.

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: #--- RUN JACKET ---#
    :end-before: # ---------------
    

THe final configuration skeleton is shown in :num:`Figure #jacket_tower_tutOpt`.
If you plot the utilization of the tower, you should get something as in :num:`Figure #utilopt-fig`, where
the Von Mises stress, global and shell buckling utilizations are shown along the tower span.  Each is a utilization and so should be <1 for feasibility.

::

>>>7 nodes in the constant-OD segment of the tower
>>>('>>>>>>>>>>  needed embedment Lp0=', 41.227393227538592)
>>>First two Freqs.= 0.2200 and 0.2212 Hz
>>>jacket+TP(structural+lumped) mass (no tower, no piles) [kg] = 1170273
>>>tower mass [kg] = 410096
>>>TP mass structural + lumped mass [kg] = 349641
>>>piles (all) mass (for assigned (not optimum) Lp [kg] =  98476
>>>frame3dd model mass (structural + TP lumped) [kg] = 1580369
>>>Tower Top Displacement in Global Coordinate System [m] =0.4196
>>>MAX member compression-bending utilization at joints = 0.5483
>>>MAX member tension utilization at joints = 0.6600
>>>MAX X-joint  utilization at joints = 0.2581
>>>MAX K-joint  utilization at joints = 0.3815

.. _jacket_tower_tutOpt-fig:

.. figure:: ./images/jacket_tower_tutOpt.*
    :width: 6in
    :align: center

    Jacket-tower structure for the tutorial example after OpenMDAO driven optimization via SNOPT. 

    
.. _utilopt-fig:

.. figure:: /images/util_tutOpt.*
    :width: 6in
    :align: center

    Utilization along tower height for: Von-Mises/yield; shell buckling; global buckling. Optimization obtained via OpenMDAO pyOPT driver with SNOPT.
        
        
.. <!-- ________________________________________ !>

Jacket-Tower Optimization - External Optimization via PyOPT
===========================================================

.. currentmodule:: jacketse.JacketOpt_PyOPT

For this tutorial, we use an auxiliary module *JacketOpt_PyOPT.py* and the auxiliary input file *MyJacketInputs.py*. 
The new module is just a wrapper for jacketSE and contains: calls to the pyOPT optimization package to perform the optimization, 
objective function and constraints. The input file contains the same jacket input information as in the previous example.

We simply run a SNOPT optimization case by issuing: ::

    >>>python JacketOpt_PyOPT.py MyJacketInputs.py True

The run terminates with the following results:

    
    .. _utilextopt-fig:
    
    .. figure:: /images/util_ExtPySnopt_tut.*
        :width: 6in
        :align: center
    
        Utilization along tower height for: Von-Mises/yield; shell buckling; global buckling. Optimization obtained via pyOPT SNOPT.

.. <!-- ________________________________________ !>

Jacket-Tower Optimization - External Optimization via Python Cobyla
===================================================================

.. currentmodule:: jacketse.JacketOpt_ExtCobyla

For this tutorial, we use an auxiliary module *JacketOpt_ExtCobyla.py* and the same auxiliary input file *MyJacketInputs.py* as above. 
The new module is just a wrapper for jacketSE and contains: calls to the python function *scipy.optimize.fmin_cobyla* to perform the optimization, 
objective function and constraints. The input file contains the same jacket input information as in the previous example.

We start by simply running: ::

    >>>python JacketOpt_ExtCobyla.py MyJacketInputs.py 

The run terminates with the following results:

::

>>>Minimum mass Mjacket, MPiles, TPmass = 1251631.941889 132121.306097 142330.589696
>>>Minimum mass Tower, Jacket(no tower no piles) = 364691.517594 886940.424295
>>>Minimum found at Dpile=1.000000, tpile=0.025400  Lp=48.589424
>>>Minimum found at Dbrc=1.000000, tbrc=0.025400
>>>Minimum found at Dbrcmud=1.000000, tbrcmud=0.025400
>>>Minimum found at batter=14.668193, dckwidth=14.000000, Dleg=1.735573, tleg=0.025400,
>>>Minimum found at Dgir=1.094296, tgir=0.025400
>>>Minimum found at Db=7.000000 DTRb=200.000000 Dt=3.941684 DTRt=200.000000 H2frac=0.250000
>>>Minimum found at Freq 0.230000
>>>Minimum found at GLutil=0.483653 EUutil=0.276818
>>>Minimum found at Mudline Footprint=17.668715 beta3D=54.489871
>>>Elapsed time:  684.977999926 seconds
>>>Execution count:  547
>>>7 nodes in the constant-OD segment of the tower
>>>('>>>>>>>>>>  needed embedment Lp0=', 48.589424211400853)

The overall mass amounts to 1383753.248 kg, which is less than the optimum found in the previous optimizations.

Cobyla seems to perform better than the other optimization options.

 .. _JacketOpt_CobyOptConfig-fig:
 
    .. figure:: ./images/JacketOpt_CobyOptConfig.*
     :width: 6in
     :align: center
 
     Jacket-tower structure for the tutorial example after OpenMDAO driven optimization via SNOPT. 

 
 .. _util_tutCobyOpt-fig:
    
    .. figure:: /images/util_tutCobyOpt.*
     :width: 6in
     :align: center
 
     Utilization along tower height for: Von-Mises/yield; shell buckling; global buckling. Optimization obtained via fmin_cobyla.
