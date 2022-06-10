########
Examples
########

This section houses a myriad of examples on how to run l2g. From the simple
"hello world" approach to examples showing the available features and
functionalities.

*********
Hello FLT
*********


************************************
JSON structure for describing a case
************************************

The command ``runL2G`` and ``submitL2G`` uses a JSON file, which describes a
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

.. code-block:: json

   // Example, Custom profile
   "custom-hlm-example-1": // Start blocks with "custom-hlm-" and add your own name for the profile
   {
       "type": "single-exp", // or "double-exp", ...
       "P_sol": 10e6, // In watts. One can specify P_sol or directly
       "lambda_q": 0.050 // In meters.
   },
   "custom-hlm-example-2": // Start blocks with "custom-hlm-" and add your own name for the profile
   {
       "type": "single-exp", // or "double-exp", ...
       "q_parallel": 10e6, // In W/m^2.
       "lambda_q": 0.050 // In meters.
   },
   "custom-hlm-example-3": // Start blocks with "custom-hlm-" and add your own name for the profile
   {
       "type": "double-exp", // or "double-exp", ...
       "q_parallel": 10e6, // In W/m^2.
       "lambda_q_main": 0.050 // In meters.
       "lambda_q_near": 0.050 // In meters.
       "Rq": 4 // No units.
   },
   "custom-hlm-example-4": // Start blocks with "custom-hlm-" and add your own name for the profile
   {
       "type": "double-exp", // or "double-exp", ...
       "q_parallel": 10e6, // In W/m^2.
       "lambda_q_main": 0.050 // In meters.
       "lambda_q_near": 0.050 // In meters.
       "Rq": 4 // No units.
   },




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


Example JSON
------------

.. _example json:

.. code-block:: json
   :caption: "JSON example"

   {
       "name": "case_name",
       "output_directory": "/path/for/output/files",

       // If we have eqdsk files
       "eq_type": "eqdsk",
       "eqdsk_files":
       [
           "/path/to/eqdsk1",
           "/path/to/eqdsk2",
       ],

       // If we have IMAS
       "eq_type": "imas",
       "imas":
       {
           "user": "public",
           "shot": 999999,
           "run": 1,
           "device": "iter",
           "version": "3",
           "time_series": true,
           "time_start": 630,
           "time_end": 673
       },

       "parameters": {}, // Change parameters defined at :class:`Settings`
       "options": {}, // Change options defined at :class:`Options`

       // Specify HLM profiles
       // Example, flat-top steady state
       "elm":
       {
            "shadow_meshes": // List of meshes used for OWL connection length graph
            [
                "/path/to/mesh1",
                "/path/to/mesh2",
                ...
            ]
            "parameters": {}, // Change parameter options, similar to FLT parameters
            "r_break": 0.025 // Specify the breakpoint position in meters.
       },

       // Example, Ramp-Down profile
       "ramp-down":
       {
            "Ip transition": 10e6
       },

       // Example, Start-Up
       "start-up":
       {
            "lambda_q_main": 0.17,
            "lambda_q_near": 0.005,
            "Rq": 4,
            "p_sol": 6e6
       },

       // Example, Custom profile
       "custom-hlm-example-1": // Start blocks with "custom-hlm-" and add your own name for the profile
       {
           "type": "single-exp", // or "double-exp", ...
           "P_sol": 10e6, // In watts. One can specify P_sol or directly
           "lambda_q": 0.050 // In meters.
       },
       "custom-hlm-example-2": // Start blocks with "custom-hlm-" and add your own name for the profile
       {
           "type": "single-exp", // or "double-exp", ...
           "q_parallel": 10e6, // In W/m^2.
           "lambda_q": 0.050 // In meters.
       },
       "custom-hlm-example-3": // Start blocks with "custom-hlm-" and add your own name for the profile
       {
           "type": "double-exp", // or "double-exp", ...
           "q_parallel": 10e6, // In W/m^2.
           "lambda_q_main": 0.050 // In meters.
           "lambda_q_near": 0.050 // In meters.
           "Rq": 4 // No units.
       },
       "custom-hlm-example-4": // Start blocks with "custom-hlm-" and add your own name for the profile
       {
           "type": "double-exp", // or "double-exp", ...
           "q_parallel": 10e6, // In W/m^2.
           "lambda_q_main": 0.050 // In meters.
           "lambda_q_near": 0.050 // In meters.
           "Rq": 4 // No units.
       },
   }