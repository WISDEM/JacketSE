JacketSE
========

JacketSE is a systems engineering model for 3 and 4-legged lattice structures (jackets) supporting tubular steel towers for offshore wind turbines.  The analysis uses beam finite element theory, cylinder drag data, linear wave theory, 
shell and global buckling methods, and API RP-2A standard code checks.  The module is developed as an OpenMDAO assembly.

Author: [RRD](mailto:nrel.wisdem+jacketse@gmail.com)

## Version

This software is a beta version 0.1.0.

## Detailed Documentation

For detailed documentation see <http://wisdem.github.io/JacketSE/>

## Prerequisites

General: NumPy, SciPy, PyFrame3DD, OpenMDAO

## Dependencies:

Wind Plant Framework: [FUSED-Wind](http://fusedwind.org) (Framework for Unified Systems Engineering and Design of Wind Plants)

Sub-Models: CommonSE, PyFrame3DD

Supporting python packages: Pandas, Algopy, Zope.interface, Sphinx, Xlrd, PyOpt, py2exe, Pyzmq, Sphinxcontrib-bibtex, Sphinxcontrib-zopeext, Numpydoc, Ipython

## Installation

First, clone the [repository](https://github.com/WISDEM/JacketSE)
or download the releases and uncompress/unpack (JacketSE.py-|release|.tar.gz or JacketSE.py-|release|.zip) from the website link at the bottom the [WISDEM site](http://nwtc.nrel.gov/WISDEM).

Install JacketSE with the following command:

    $ python setup.py install

or if in an activated OpenMDAO environment:

    $ plugin install


## Run Unit Tests

To check if installation was successful try to import the module from within an activated OpenMDAO environment:

    $ python
    > import jacketse.jacket

You may also run the unit tests.

    $ python src/jacketse/test/test_jacketse.py

For software issues please use <https://github.com/WISDEM/JacketSE/issues>.  For functionality and theory related questions and comments please use the NWTC forum for [Systems Engineering Software Questions](https://wind.nrel.gov/forum/wind/viewtopic.php?f=34&t=1002).



