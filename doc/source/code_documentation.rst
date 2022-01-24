###########
Python code
###########

The following section shows the documentation of the python API for using the
FLT. It is heavily using the OOP aspect of python in hope to be intuitive and
reusable.

**********
l2g module
**********

The top level of the module setups the logging formatting.

.. automodule:: l2g
   :members:
   :undoc-members:
   :private-members:

****************
l2g.utils module
****************

This module contains utility functions:

 - Read/Write data from MED files
 - Helper functions when running l2g with JSON files.

.. automodule:: l2g.utils
   :members:
   :undoc-members:
   :private-members:

***************
l2g.plot module
***************

Contains functions used extensively in plotting equilibriums.

.. automodule:: l2g.plot
   :members:
   :undoc-members:
   :private-members:

****************
l2g.equil module
****************

This module contains functions for handling equilibrium data:

 - Functions for retrieving and writing equilibrium data in EQDSK G and IMAS
   format
 - EQ class for performing evaluations on a equilibrium
 - Equilibrium class, the main class for storing relevant equilibrium data
   (for a fixed reference when handling different input data types). If new
   storage formats are introduced for loading equilibriums, then a function has
   to be prepared that copies the data to an instance of this class.
 - EquilibriumIterator, a class that can hold multiple equilibriums and acts
   as an iterable object.

.. note::

   The classes are accessible directly from :py:mod:`l2g.equil`, but to avoid
   cluttered and massive single python files, the code has been chopped and
   "hidden" in underscore-prefixed python submodules.

.. automodule:: l2g.equil
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.equil._eq
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.equil._eqdskg
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.equil._equilibrium
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.equil._iter

**************
l2g.hlm module
**************

Contains functions for applying the heat load mapping of different phases
of a scenario (e.g., Ramp-Down, Steady-State, ...)

.. automodule:: l2g.hlm.general
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.hlm.ramp_down
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.hlm.steady_state
   :members:
   :undoc-members:
   :private-members:

***************
l2g.comp module
***************

The heart or component or compute module. It holds the main classes that
enables the run of FLT.

.. note::

   Most of these classes are "hidden" in python files that have underscore "_" at
   the start of the name. All of them are mainly accessible, intentionally via
   through the :py:mod:`l2g.comp` namespace.

l2g.comp._data
==============

This module contains data classes. They only hold quantities based on what
phase or role and maybe contain functions for generating file-like objects
for visualization or storage (e.g. generating a VTK object for storage or
visualization in ParaView).

.. automodule:: l2g.comp._data
   :members:
   :undoc-members:
   :private-members:

l2g.comp._field_line_tracer
===========================

The main class that is used for running FLT. It is equipped with functions for
loading input data, setting parameters, running FLT and post-compute
application of heat load maps.

.. automodule:: l2g.comp._field_line_tracer
   :members:
   :undoc-members:
   :private-members:

l2g.comp._io
============

This module contains functions for input output handling of FLT result data.

.. automodule:: l2g.comp._io
   :members:
   :undoc-members:
   :private-members:

l2g.comp._settings
==================

This module holds all settings that are and will be introduced to control the
flow of a FLT study.

.. automodule:: l2g.comp._settings
   :members:
   :undoc-members:
   :private-members:


********************
l2g.comp.core module
********************

This module holds the "lowest" level functions that operate with external
libraries (L2G_cpp and Embree) for running FLT. It is written in Cython to
enable high performance of computation and OpenMP parallelization.

.. automodule:: l2g.comp.core._core
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.comp.core.ext_wrapper
   :members:
   :undoc-members:
   :private-members: