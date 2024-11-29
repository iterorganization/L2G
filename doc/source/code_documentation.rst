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

.. automodule:: l2g.equil
   :members:
   :undoc-members:
   :private-members:

***************************
Head load mapping submodule
***************************

Contains functions for applying the heat load mapping of different phases
of a scenario (e.g., Ramp-Down, Steady-State, ...).

.. automodule:: l2g.hlm.general
   :members:
   :undoc-members:
   :private-members:


.. automodule:: l2g.hlm.ramp_down
   :members:
   :undoc-members:
   :private-members:


.. automodule:: l2g.hlm.elm_plm
   :members:
   :undoc-members:
   :private-members:


*********************
Computation submodule
*********************

The heart or component or compute module. It holds the main classes that
enables the run of FLT.

Data storage classes
====================

This module contains data classes. They only hold quantities based on what
phase or role and maybe contain functions for generating file-like objects
for visualization or storage (e.g. generating a VTK object for storage or
visualization in ParaView).

The :py:class:`l2g.comp.L2GResults` ontains the arrays for storing all
quantities that are related to FLT.

.. automodule:: l2g.comp.L2GResults
   :members:
   :undoc-members:
   :private-members:

The :py:class:`l2g.comp.L2GFLs` contains the arrays for storing the points of
specific fieldline we want to plot.

.. automodule:: l2g.comp.L2GFLs
   :members:
   :undoc-members:
   :private-members:


The :py:class:`l2g.comp.L2GResultsHLM` is the class for containing data
related to the heat load mapping.
.. automodule:: l2g.comp.L2GResultsHLM
   :members:
   :undoc-members:
   :private-members:

FieldLineTracer class
=====================

The main class that is used for running FLT. It is equipped with functions for
loading input mesh data, equilibrium data, setting parameters, running FLT and
post-compute application of heat load maps.

.. autoclass:: l2g.comp.FieldLineTracer
   :members:
   :undoc-members:
   :private-members:


Settings classes
================

This module holds all settings and parameters the user can set.


Parameters are settings that affect the numerical aspect.

.. autoclass:: l2g.settings.Parameters
   :members:
   :undoc-members:
   :private-members:

Options are settings that affect the global behaviour (i.e., activating RT or
not).

.. autoclass:: l2g.settings.Options
   :members:
   :undoc-members:
   :private-members:


************************************************
Cython implementations and wrapped external code
************************************************

This module holds Cython implementation of algorithms and the wrappers for
external code.

.. note::

   The classes :py:class:`l2g.external.bfgs_2d` and
   :py:class:`l2g.external.equilibrium_analysis.EQA` are implemented with speed
   in mind. They perform the same as the python equivalents (find_minimum of
   scipy and the l2g.equil.EQ class)


.. automodule:: l2g.external.bfgs_2d
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.external.equilibrium_analysis
   :members:
   :undoc-members:
   :private-members:

The following submodules are wrappers to the external L2G c++ kernel code.
They are paramount, without them there is no FLT and you can see the usage
in the :py:class:`l2g.comp.FieldLineTracer` code as example.

.. automodule:: l2g.external.bicubic
   :members:
   :undoc-members:
   :private-members:


.. automodule:: l2g.external.embree
   :members:
   :undoc-members:
   :private-members:

.. automodule:: l2g.external.flt
   :members:
   :undoc-members:
   :private-members:


And finally, the RKF45 in this case is an implementation of rht :term:`RKF45`
algorithm but in this case specifically for solving the fieldline equations in
the plasma axisymmetric case.

.. automodule:: l2g.external.rkf45
   :members:
   :undoc-members:
   :private-members:

******************
Mesh file handling
******************

This module provides a minimal functionality for reading/writing unstructured
grid mesh files with fields. There is a single class that based on the file
extension utilizes different backends.

An unstructured grid mesh here is considered a mesh built from triangles. Other
cell types are ignored.

Mesh class
==========

This class is the main interface used for reading/writing data from/to a
MEDCOUPLING med file.

.. code-block:: python

   import l2g.mesh
   file_path = '/path/to/file.med'

   mesh = l2g.mesh.Mesh(file_path)

   # Get mesh data
   mesh.getMeshData() # Returns vertices, triangles

   # Get field data
   # Set the index if there are multiple time series
   mesh.setIndex(2)
   # Get a field
   mesh.getField("Field name")

   # Now the opposite. Let's write the mesh to a new file.
   new_mesh = mesh.writeMeshTo('/new/file/path.ext')

   number_of_triangles = len(new_mesh.triangles)

   # Let's add some time steps
   import numpy as np

   new_mesh.setNumberOfTimeSteps(5)

   for i in range(5):
      # Set the index
      new_mesh.setIndex(i)
      # Provisionally set the time
      new_mesh.setTime(i * 0.25)
      array = np.random.random(number_of_triangle)
      new_mesh.addField("Field name", array)

   # Now write the fields
   new_mesh.writeFields()

.. autoclass:: l2g.mesh.Mesh
   :members:
   :undoc-members:
   :private-members:

#######
Scripts
#######

In the bin/ directory contains various scripts. Here is a short description of
each:

 - flat
     Script that takes YAML configuration files and run FLT cases.
 - get_disruption_profile_from_imas
     Script that obtains the disruption heat load profile from the IMAS
     Disruption IDS.
 - imas2eqdsk
     From an Equilibirum IDS time slice create a G-EQDSK file.
 - mat2med
     Read matlab fields and store it to a MED file (that contains a mesh of same
     size as fields in matlab)
 - med2mat
     Convert the fields of MED file to a matlab file.
 - mkEqdskMovie
     Create a GIF out of a set of EQDSK files. It plots the boundary flux
     surfaces
 - mkImasMovie
     Creates a GIF out of an IMAS Equilibrium IDS by taking the boundary flux
     values stored there.
 - mkImasMovieFromPsi
     Creates a GIF out of an IMAS Equilibrium IDS by analyzing the equilibrium
     data and obtaining the boundary flux surfaces. In case if the Equilibrium IDS
     do not have stored the boundary flux surfaces.
 - plotIP
     Creates a plot of the plasma current and power from the IMAS Summary IDS.
 - plotMFlux
     Plots the poloidal magnetic flux from either an IMAS or G-EQDSK source.
 - submitFLAT
     Same as flat but it runs the FLT case on the cluster (via SLURM).

In order to see how to use them simply run each script with ``-h`` or read
the code itself.