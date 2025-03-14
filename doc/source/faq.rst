############################
(Frequently) Asked Questions
############################

What is the difference between cutoff connection length (old naming of max_fieldline_length) and maximum connection length?
===========================================================================================================================

There are two parameters called cutoff_conlen and maximum_connection_length
that a user can set for a FLT case.

The maximum connection lengths, means to what length are the fieldlines
followed before stopping. With this you reduce computation time as some
fieldlines can be kilometers long in certain cases. This also causes confusion,
if we folow a fieldline and it spans to a certain maximum reachable length, set
by us, then it cannot really be treated as a connection fieldline as it is not
really connecting to any physical geometry.

The cutoff connection length is used to determine, conservatively, the wetted
areas on your analyzed geometry. In order for us to determine whether an area
is wetted, we take a fieldine from that area, follow it and see if the length
is larger than this cutoff connection length value. If it is, the area is
considered wetted.

I changed my target meshes in YAML case file but the results do not change
==========================================================================

The program checks if a resulting file exists. It does not check the contents
entirely as that is too complex to perform. Therefore when in doubt that the
resulting file does not match your expectation, delete the resulting MED file
so that it is created from scratch based on the YAML case file.

What is a .med file and why is it used
======================================

This is a binary HDF5 based file type for storing meshes and quantities related
to meshes. The python interface is fairly simple and strong with a numpy
coupling. While the library can use whatever file format as the it is not tied
directly to specific file formats, the .med format has shown great success for
reading/writing, editing, manipulating with great performance as well as
showing in visualization software like ParaView.

I am using align_lcfs but the limiter plasma is unnaturally pushed towards the geometry
=======================================================================================

With the align_lcfs option the limiter plasma is actually seated on the
geometry we want to analyze, or as we say target geometry. But be careful, if
in the target geometry you do not have the actual geometry where the plasma is
supposed to touch, the software will quite simply push it to the next closest
cell in the mesh. So it is extremely ``important``, when we are analyzing
limiter plasma and we are using align_lcfs options, to have ENOUGH geometry,
including the part where the plasma is really supposed to touch. Otherwise the
result will be garbage.

Usually in the axisymmetric sources of plasma equilibrium the silhouette of the
wall represents the heated wall during operation and the models that we analyze
are usually cold geometry. Because of the thermal expansion there is a slight
disconnection between the geometry we want to analyze and the plasma
equilibrium source.

Additionally, to understand what align_lcfs really does: In the actual
computation of the FLT it introduces an additional offset that zeros the
limiter plasma directly to the geometry. You can still manually shift the
plasma vertically and radially, what this option does is then finally clips the
plasma to the closest part of the geometry.

.. noted::

   With the advanced parameter ``lcfs_max_align_dist`` you tell the program
   that the LCFS can be aligned only if the alignment distance is smaller than
   the value of the parameter.
