########
Examples
########

***********
First steps
***********

This chapter will show the minimum number of steps required to run a
:term:`FLT` case using the python API.

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
   tlas_obj = l2g.comp.core.PyTLAS()
   # We can also include the target mesh as shadow, in case of self shadowing
   # effects. Otherwise if the shadow meshes already includes the target, we
   # can skip this step
   tlas_obj.commitMesh(vertices, triangles)

   tlas_obj.commitMesh(shadowing_vertices_1, shadowing_triangles_1)
   tlas_obj.commitMesh(shadowing_vertices_2, shadowing_triangles_2)

   # Tie the tlas_obj to the flt_obj
   flt_obj.setTLASObj(tlas_obj)

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
