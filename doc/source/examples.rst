########
Examples
########

This section shows examples of how to use the code. They are not complete
examples.

To find complete examples, check the synthetic case, where the whole case is
done inside the python file:

.. literalinclude:: ../../examples/synthetic_case.py
   :language: python

************************************
YAML structure for describing a case
************************************

The commands ``flat`` and ``submitFLAT`` uses a YAML type file as input, where
the input data is described in different YAML documents with the following type
fields:

 - geometry
 - equilibrium
 - hlm

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
========

The geometry document describes what input meshes to be used as the target,
shadowing and the meshes for filtering out fieldlines that escape the tokamak
chambre (also called artificial fieldlines). Additionally parameters and
options are set in this block for setting the tracing tool.

Example YAML geometry block:

.. code-block:: yaml

   ---
   type: geometry
   name: name_for_target
   target_mesh: /path/to/target/mesh.med
   shadow_meshes:
    - /path/to/shadowing/mesh1.med
    - /path/to/shadowing/mesh2.med
    - ...
   afl_catchers_meshes:
    - /path/to/afl/catching/mesh1.med
    - /path/to/afl/catching/mesh2.med
    - ...
   cutoff_conlen: 4.3
   parameters:
       time_step: 0.01
       time_end: 31.41 # 5 full revolutions
       max_connection_length: 10000 # in meters
       self_intersection_avoidance_length: 0.001 # in meters
       side: owl
   ...

Equilibrium
===========

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

Hlm
===

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

The command ``flat`` and ``submitFLAT`` uses a JSON file, which describes a
FLT case to run. In this JSON format file we describe the case, the input
files, parameters, ... The commands parse the file, runs the case and if all
goes well produces a result file along the auxiliary data files and images.

An example JSON file can be seen at the end of this section
:ref:`here <example json>`.

Parameters
==========

The following list of parameters used for running a case. The names of the
parameters are shown verbatim as they are used either in a JSON file or as
attributes in python.

.. note::

   Unless specified, every parameter is required to be set!.

name
----

Name of the case. This affects the name displayed in
the SLURM job list and the output name of files.


output_directory
----------------

Specify the path for dumping result files.

eq_type
-------

Specify in what format is the input equilibrium
data.

 * "eqdsk"
 * "imas"

eqdsk_files
-----------

**REQUIRED** if ``eqdsk_type==eqdsk``. Specify a list of EQDSK files (paths).
The order of the EQDSK files corresponds to the order of FLT results in the
resulting MED file.

imas
----

**REQUIRED** if ``eqdsk_type==imas``. Specify a set of IMAS parameters from
which to obtain equilibriums data:

 * user
 * shot
 * run
 * device
 * version
 * time_series - True if more than 1 time slice
 * time_start - From which time to take equilibrium data.
 * time_end - To which time to take equilibrium data.


.. _json-global-parameters:

parameters
----------

**OPTIONAL**. This block globally sets the FLT settings for all FLT objects and
act as the default settings. Check the :py:class:`l2g.comp._settings.Settings`
for all available parameters to set and what are their defaults.

.. _json-global-options:

options
-------

**OPTIONAL**. This block contains options on how FLT is performed. For instance
if the fieldlines are followed until a set maximum length or toroidal angle or
if intersection  tests should be activated during tracing, etc... See
:py:class:`l2g.comp._settings.Options` for more information.



.. _main-flt-object:

flt
---

Main FLT object. This block contains information on what input mesh is used as
target geometry, shadowing geometry. Essentially it contains the following
settings:

 * target_mesh: Path to a mesh used as target mesh or shadowed mesh
 * shadow_meshs: List of paths, pointing to meshes which are used in the
   intersection tests.
 * afl_catcher_meshes: List of paths, pointing to meshes which are also used in
   the intersection tests, but if a FL penetrate them they are marked finally
   as shadowed.
 * parameters: FLT parameters just for this FLT object (see
   :ref:`json-global-parameters`)
 * options: FLT options just for this FLT object
   (see :ref:`json-global-options`)


elm
---

**OPTIONAL**. This blocks applies the flat-top Steady-State plasma profile to
the :ref:`main FLT object<main-flt-object>`. Additionally, it also runs FLT,
but in order to obtain OWL connection length graph.

Additionally it requires the following settings:

 * shadow_meshes: List of meshes to be used for the shadowing. In this case the
   whole tokamak + divertor is recommended.
 * parameters: FLT parameters, but again just for this FLT object (see
   :ref:`json-global-parameters`)
 * options: FLT options just for this FLT object
   (see :ref:`json-global-options`)
 * "r_break": Breakpoint location in meters.


ramp-down
---------

**OPTIONAL**. Presence of this block applies the ramp-down plasma profile to
the :ref:`main FLT object<main-flt-object>`. Used only when we have ``IMAS``
equilibrium input as from the IMAS we can gather more data necessary for
calculating the heat load.

Additionally it should have the ``Ip transition`` parameter which tells us when
does the H to L mode transition happens. Every parameter necessary for
calculating the heat load comes from the IMAS dabatase.

custom-hlm-*
------------

**OPTIONAL**. Presence of this block (notice the **asterisk**) applies a plasma
heat load profile on the :ref:`main FLT object<main-flt-object>`. The
**asterisk** is replaced by a user specified name (for instance:
custom-hlm-plasma-1, where the name of the heat load array will be "plasma-1").

The following custom profiles are avilable:

 * single-exp
 * double-exp

See :py:func:`l2g.hlm.general.single_exponential_qpar`,
:py:func:`l2g.hlm.general.single_exponential_psol`,
:py:func:`l2g.hlm.general.double_exponential_psol` for a list of parameters
you can set. In order to change the parameters you simply set in your custom
block

.. note::

   Drsep, Bt, Bpm, Rb are automatically supplied by the program and are not to
   be defined. All the other parameters must be written in the custom plasma
   profile block.

See the following example on the custom :term:`hlm` blocks.


.. code-block:: yaml

   ---
   type: geometry
   name: name_for_target
   target_mesh: /path/to/target/mesh.med
   shadow_meshes:
    - /path/to/shadowing/mesh1.med
    - /path/to/shadowing/mesh2.med
    - ...
   afl_catchers_meshes:
    - /path/to/afl/catching/mesh1.med
    - /path/to/afl/catching/mesh2.med
    - ...
   cutoff_conlen: 4.3
   parameters:
       time_step: 0.01
       time_end: 31.41 # 5 full revolutions
       max_connection_length: 10000 # in meters
       self_intersection_avoidance_length: 0.001 # in meters
       side: owl
   ...
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
   ---
   type: hlm
   name: name_for_heat_load_mapping
   hlm_type: single
   p_sol: 7.5E+6
   lambda_q: 0.012
   ...


wall_limiter
------------

Optional ``wall_limiter``. Specify a custom set of points, describing the wall
silhouette. In meters.

The wall limiter variable block contains two sub-arrays called
``r`` and ``z``. These r, z points describe the profile of the
blanket+divertor. These points are used for analysis in the
:py:class:`l2g.equil._eq.EQ` class. Within the analysis class, equilibrium data
is analyzed used these points as reference. In rare cases the equilibrium data
and the actual geometry used in the FLT analysis might not be aligned, either
by a fixed major radius offset or something else. In this case the user can
specify inside the JSON file the points of the wall silhouette which will be
used then for equilibrium analysis.

