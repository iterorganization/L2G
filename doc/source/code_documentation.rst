###########
Python code
###########

The following section shows the documentation of the python API for using the
FLT. It is heavily using the OOP aspect of python in hope to be intuitive and
reusable.

********************
Top module functions
********************

The top level of the module setups the logging formatting.

.. autofunction:: l2g.enableLogging

.. autofunction:: l2g.enableDebugging

.. autofunction:: l2g.addStreamHandler


******************
Plotting submodule
******************

Contains functions used in creating 2d graphs and movies of a plasma
equilibrium.

.. automodule:: l2g.plot
   :members:
   :undoc-members:
   :private-members:

*************************************
Plasma equilibrium handling submodule
*************************************

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

.. autoclass:: l2g.equil.EQ
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: l2g.equil.Equilibrium
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: l2g.equil.EQDSKIO
   :members:
   :undoc-members:
   :private-members:


.. autofunction:: l2g.equil.correct_equilibrium_helicity

***************************
Head load mapping submodule
***************************

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

*********************
Computation submodule
*********************

The heart or component or compute module. It holds the main classes that
enables the run of FLT.

.. note::

   Most of these classes are "hidden" in python files that have underscore "_" at
   the start of the name. All of them are mainly accessible, intentionally via
   through the :py:mod:`l2g.comp` namespace.

Data storage classes
====================

This module contains data classes. They only hold quantities based on what
phase or role and maybe contain functions for generating file-like objects
for visualization or storage (e.g. generating a VTK object for storage or
visualization in ParaView).

.. autoclass:: l2g.comp.L2GResults
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: l2g.comp.L2GPointResults
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: l2g.comp.L2GFLs
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: l2g.comp.L2GResultsHLM
   :members:
   :undoc-members:
   :private-members:

FieldLineTracer class
=====================

The main class that is used for running FLT. It is equipped with functions for
loading input data, setting parameters, running FLT and post-compute
application of heat load maps.

.. autoclass:: l2g.comp.FieldLineTracer
   :members:
   :undoc-members:
   :private-members:

MEDCoupling I/O class
=====================

This module contains functions for input output handling of FLT result data.

.. autoclass:: l2g.comp.MEDMeshIO
   :members:
   :undoc-members:
   :private-members:

.. autofunction:: l2g.comp.dump_flt_mesh_results_to_med

.. autofunction:: l2g.comp.load_flt_mesh_results_from_med

.. autofunction:: l2g.comp.save_fls_to_vtk

.. autofunction:: l2g.comp.save_mesh_to_vtk

.. autofunction:: l2g.comp.save_results_to_vtk

Settings classes
================

This module holds all settings that are and will be introduced to control the
flow of a FLT study.

.. autoclass:: l2g.settings.Parameters
   :members:
   :undoc-members:
   :private-members:

.. autoclass:: l2g.settings.Options
   :members:
   :undoc-members:
   :private-members:


****************************************************
Core computation wrapper code of external C++ kernel
****************************************************

This module holds the "lowest" level functions that operate with external
libraries (L2G_cpp and Embree) for running FLT. It is written in Cython to
enable high performance of computation and OpenMP parallelization.

Not all functions, especially the ones that are only called with Cython context
are shown in this section. Best way is to delve into the files directly to see
the structure.

There you will see that in some cases one function can have up to three
different signatures, based on the use case: usable in Python, GIL-free for
use in Cython context and as previous case + overloading.

.. automodule:: l2g.comp.core._core
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.comp.core.ext_wrapper
   :members:
   :undoc-members:
   :private-members:
