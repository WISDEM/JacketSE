# DEPRECATED
------------

**THIS REPOSITORY IS DEPRECATED AND WILL BE ARCHIVED (READ-ONLY) IN NOVEMBER 2019.**

WISDEM has moved to a single, integrated repository at https://github.com/wisdem/wisdem

---------------
# JacketSE

JacketSE is a systems engineering model for 3 and 4-legged lattice structures (jackets) supporting tubular steel towers for offshore wind turbines.  The analysis uses beam finite element theory, cylinder drag data, linear wave theory, shell and global buckling methods, and API RP-2A standard code checks.  The module is developed as an OpenMDAO assembly.

Author: [NREL WISDEM Team](mailto:systems.engineering@nrel.gov) 

## Documentation

See local documentation in the `docs`-directory or access the online version at <http://wisdem.github.io/JacketSE/>

## Installation

For detailed installation instructions of WISDEM modules see <https://github.com/WISDEM/WISDEM> or to install JacketSE by itself do:

    $ python setup.py install


## Run Unit Tests

To check if installation was successful try to import the package:

    $ python
    > import jacketse.jacket

You may also run the unit tests.

    $ python src/jacketse/test/test_jacketse.py

For software issues please use <https://github.com/WISDEM/JacketSE/issues>.  For functionality and theory related questions and comments please use the NWTC forum for [Systems Engineering Software Questions](https://wind.nrel.gov/forum/wind/viewtopic.php?f=34&t=1002).



