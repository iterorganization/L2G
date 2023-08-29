########
Examples
########

***********
First steps
***********

This chapter will show the minimum number of steps required to run a
:term:`FLT` case.

Mesh data
=========

Essentially when it comes to :term:`FLT` all we require is the mesh data, a
target and shadowing mesh and an equilibrium source. If we use scripting there
is no strict format on the mesh files, what is important is only the data type
of the arrays and their shape. For instance the format of the mesh data is
supposed to be as following:

.. code-block:: python

   # We obtain mesh data from a file. We only handle vertices and triangles.
   # The triangles are simple, plane type and *not* higher orders (quadratic)

   import numpy as np

   # It is essential the vertices are supplied as a 1D numpy array with the
   # data type np.float32.

   # Example dummy data
   vertices = np.zeros(9, dtype=np.float32)

   # Next are triangles. Again it is a 1D numpy array of unsigned integers.
   triangles = np.array([0,1,2,3,4,5,6,7,8], dtype=np.uint32)

   # In this case the first triangle is the first trio of vertices ID's,
   # index starts with 0, second is the second trio, etc...

   # Again we create here dummy variable for shadowing mesh
   shadow_vertices_1 = np.zeros(9, dtype=np.float32)
   shadow_triangles_1 = np.array([0,1,2,3,4,5,6,7,8], dtype=np.uint32)

   # We can have multiple shadowing meshes.
   shadow_vertices_2 = np.zeros(9, dtype=np.float32)
   shadow_triangles_2 = np.array([0,1,2,3,4,5,6,7,8], dtype=np.uint32)


.. note::

   It is important that the units of the meshes (vertices) is always in meters!

Equilibrium data
================

Now that we have our mesh data, in this case we will have variables only for
one dummy mesh, we require a source of equilibrium. In L2G the class for
handling equilibrium data is the `l2g.equil.Equilibrium`. But the source of
the data can be either the EQDSK-G or from IMAS format:

EQDSK-G format
--------------

.. code-block:: python

   import l2g.equil

   # For instance if we have a single EQDSK file we may use
   eqdsk_file = "/path/to/EQDSK/file.eqdsk"

   # There is a helper function that can gather the equilibrium data from a
   # EQDSK format file

   equilibrium = l2g.equil.getEquilibriumFromEQDSKFile(eqdsk_file)
   # Check the function on what are the other options available, for instance
   # for ITER there is a fixed helicity, by default the signs of the magnetic
   # poloidal flux function is changed to mimic the correct helicity.

   # The equilibrium object is now ready to be used  for :term:`FLT`

IMAS format
-----------

.. code-block:: python

   import l2g.equil

   # For instance we wish to gather a single time slice equilibrium data from
   # the IMAS database.

   equilibrium = l2g.equil.getEquilibriumFromIMASSlice(run=135011, shot=7,
      user='public', machine='iter', time=400.0)

   # The equilibrium object is now ready to be used  for :term:`FLT`

Now we wish to setup the :term:`FLT` class in order to start it.

.. code-block:: python

   import l2g.comp
   import l2g.comp.core

   flt_obj = l2g.comp.FieldLineTracer()

   # First we set the equilibrium and mesh data
   flt_obj.setTargetData(vertices, triangles)

   # Now we load the Shadowing mesh data
   embree_obj = l2g.comp.core.PyEmbreeAccell()
   # We can also include the target mesh as shadow, in case of self shadowing
   # effects. Otherwise if the shadow meshes already includes the target, we
   # can skip this step
   embree_obj.commitMesh(vertices, triangles)

   embree_obj.commitMesh(shadowing_vertices_1, shadowing_triangles_1)
   embree_obj.commitMesh(shadowing_vertices_2, shadowing_triangles_2)

   # Tie the embree_obj to the flt_obj
   flt_obj.setEmbreeObj(embree_obj)

   # Now we can set the case parameters. Check the code PUT_LINK_HERE to see
   # the documentation of available options
   ...

   # Now we have to call two function that propagate the settings and
   # equilibrium data to the C++ kernel. This has to be done at least once
   # after an equilibrium or settings are changed.
   flt_obj.applyParameters()
   flt_obj.loadEq()

   # Now we can run the :term:`FLT` on the target mesh.
   flt_obj.runFltOnMesh()

   # After the computation is finished, the results are stored in the
   # flt_obj.mesh_results. We can save it by using the following function

   l2g.comp.save_results_to_vtk(flt_obj.mesh_results, "/path/to/save.vtu")

   # If we wish to plot fieldlines, we can do it as following:

   # Set the list of the FieldlineTracer object `fl_ids`. This list contains
   # the ID's of the triangles, with the index starting from 0.
   flt_obj.fl_ids = [0, 1, 2]

   # If applyParameters or loadEq weren't called, do it again.
   flt_obj.applyParameters()
   flt_obj.loadEq()

   # Now get the FLs
   flt_obj.getFL()

   # The FL results are stored in the fl_results attribute of the
   # FieldlineTracer object.
   l2g.comp.save_results_to_vtk(flt_obj.fl_results, "/path/to/save.vtu")

This is in essence the functions required to run a FLT case. To perform
additional steps, such as applying misalignments or changing parameters and/or
options, you can read the following examples:

.. literalinclude:: ../../examples/example_1.py
   :language: python

A "synthetic" case, where an imaginary equilibrium and meshes are generated
all in one script:

.. literalinclude:: ../../examples/synthetic_case.py
   :language: python

**************************
Text based (YAML) workflow
**************************

In the :py:mod:`l2g.workflow` is the code for a `YAML <https://yaml.org>`_ text
based workflow. The idea is to use the YAML file as a storage for information
on where the program can find the meshes, equilibrium data and/or which heat
load mapping to use.

We can have multiple YAML files or a single **master** file, as the YAML format
allows for multiple YAML documents inside a single file.

Then with the help of the executable scripts ``flat`` and ``submitFLAT``, which
take the YAML file as the primary input argument, we can run a case which is
the combination of the data described in the provided YAML file.

YAML file content
=================

A YAML file should contain the following type of a document block:

 - geometry
 - equilibrium
 - hlm

.. note::

   This is only a convention, there is nothing special between these different
   YAML documents, except the information they hold.

What is common for these three documents, is that they all have the following
common entries:

 - name: Name or label of the document. It is used when writing results to files.
 - type: What kind of information it holds. Whether information about the
     meshes, equilibriums or heat load mapping profiles.


To obtain the list of all fields that could be filled in these yaml blocks:

.. code-block:: python

   import l2g.settings

   # For obtaining Parameters fields, to be filled in the Geometry block under
   # the parameters field
   print(l2g.settings.Parameters().dump())

   # For obtaining Options fields, to be filled in the Geometry block under
   # the parameters field
   print(l2g.settings.Options().dump())

   # And to obtain the YAML fields
   import l2g.workflow
   print(l2g.workflow.GEOMETRY())
   print(l2g.workflow.EQUILIBRIUM())
   print(l2g.workflow.HLM())

   # For more information check the code directly


Geometry
--------

The geometry block in a YAML file should contain the name of the mesh files
to be used as target mesh, shadowing mesh or even the artificial fieldline
catcher mesh.

Additionally it contains FLT options and parameters. Such as to what length
fieldlines are followed, the resolution of the fieldlines, etc...

Finally it can also contain entries for applying rotational or longwave
misalignment. These are special cases.

An example YAML document for geometry:


.. code-block:: yaml

   ---
   type: geometry
   name: geometry_label
   target_mesh:
     file: /path/to/the/MEDCoupling/mesh/file.med
     include_target_in_shadow: false/true
   shadow_meshes:
     files:
       - /1/path/to/a/shadowing/MEDCoupling/mesh/file.med
       - /2/path/to/a/shadowing/MEDCoupling/mesh/file.med
       - /3/path/to/a/shadowing/MEDCoupling/mesh/file.med
       - ...
     exclude_pattern: # Optional
       - pattern_1
       - pattern_2
       - ...
     artificial_fieldline_catcher_files: # Optional
       - /1/path/to/afl/catcher/MEDCoupling/mesh/file.med
       - /2/path/to/afl/catcher/MEDCoupling/mesh/file.med
       - /3/path/to/afl/catcher/MEDCoupling/mesh/file.med
       - ...
   cutoff_conlen: 4.3 # In meters.
   rotational_misalignment:
     vectors: []
     angles: []
   longwave_misalignment:
     vectors: []
     lengths: []
   parameters:
     ...
   options:
     ...
   ...

The entries ``type``, ``name``, ``target_mesh`` and ``shadow_meshes`` are
self-explanatory.

Target mesh
^^^^^^^^^^^

This entry contains two entries ``file`` and ``include_target_in_shadow``.
The first one is self-explanatory while the other one *means* that the target
mesh should also be loaded along the shadowing meshes for the intersection
tests. The reason why this option is here is to avoid loading the same meshes
twice.

.. note::

   The structure is that there can only be one mesh assigned as target, but
   for shadow there can be more than one.



For the parameters and options entries see :py:mod:`l2g.settings.Parameters`
and :py:mod:`l2g.settings.Options`.


.. note::

   One of the more important parameters to set is the scaling factor that
   transforms the dimension unit of the meshes to meters. By default the
   program expects meshes to be in millimeters:

     - parameter
       * target_to_m = 1e-3 # Default from mm to m
       * shadow_to_m = 1e-3 # Default from mm to m

   Of course it is expected when providing multiple shadowing meshes that they
   all have the same dimension unit. Otherwise you have to properly transform
   the units.

Shadow meshes
^^^^^^^^^^^^^

This entry contains 3 entries.

Equilibrium
-----------

The equilibrium block specifies the input data of the equilibrium data. Usually
this only means either a list of EQDSK files or parameters to the IMAS
database. The equilibrium_type field specifies in what form is the equilibrium
data stored.

Example of YAML equilibrium:

.. code-block:: yaml

   ---
   type: equilibrium
   equilibrium_type: eqdsk_files
   name: collection_of_equilibriums
   eqdsk_files:
    - /path/to/eqdsk/file1.eqdsk
    - /path/to/eqdsk/file2.eqdsk
    - ...
   ...
   ---
   type: equilibrium
   equilibrium_type: imas
   name: equilibrium_from_imas
   imas:
     user: public
     shot: 135013
     run: 2
     device: iter
     version: '3'
     time_start: 630
     time_end: 673
   ...

HLM
---

The hlm block specifies what heat load mapping to apply to the study:

.. code-block:: yaml

   ---
   type: hlm
   name: ramp-down
   hlm_type: ramp-down
   ip_transition: 10.e+6 # Default
   ...

   ---
   type: hlm
   name: flat-top
   hlm_type: elm
   shadow_meshes:
     - /home/ITER/simicg/MESH_DIRECTORY/FULL_BLANKET_MESH/FullTokamak.med
     - /home/ITER/simicg/MESH_DIRECTORY/DIVERTOR/Divertor.med
   parameters:
     max_connection_length: 1000 # default
     time_step: 0.01
     abs_error: 0.0001
     rel_error: 0.0001
   r_break: 0.025 # default
   ...

   ---
   type: hlm
   name: single_exponential
   hlm_type: single
   p_sol: 7.5E+6
   lambda_q: 0.012
   ...

   ---
   type: hlm
   name: double_exponential
   hlm_type: double
   p_sol: 7.5E+6
   lambda_q_near: 0.005
   lambda_q_main: 0.17
   ratio: 4
   ...

   ---
   type: hlm
   hlm_type: custom
   name: list_of_custom_profiles
   profile_files:
     - /path/to/profile/file1.txt
     - /path/to/profile/file2.txt
     - /path/to/profile/file3.txt
     - /path/to/profile/file4.txt
     - /path/to/profile/file5.txt
     - /path/to/profile/file6.txt
     - /path/to/profile/file7.txt
     - /path/to/profile/file8.txt
     - /path/to/profile/file9.txt
   ...
