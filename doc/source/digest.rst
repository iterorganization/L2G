
#############
Digesting L2G
#############

********
Workflow
********

The way L2G is structured is to have a singular python interface, providing all
the tools a user needs or might need in order to perform their FLT studies.
From the, obviously, running a :term:`FLT` case to plotting utilities for
intermediate data and diagnostics of their input equilibrium data.

Python is used for interface due to it's availability and broad specter of
libraries, which brings functionality and relatively high productivity.

For high performance and compatibility with python, C++ is used for the kernel,
e.g., the heart of the FLT computation and Cython is used to glue together
externally used C++ code with python. Parallelization is performed on
Python/Cython side with the use of OpenMP. Reason for using OpenMP is discussed
in :ref:`why_openmp`.

.. graphviz::
   :name: l2g workflow
   :caption: L2G scheme
   :align: center

    digraph "l2g-workflow" {
        size="8,6";
        rankdir="TD";
        splines=true;
        compound=true;

        graph [fontname="Verdana", fontsize="9"];
        node [fontname="Verdana", fontsize="9"];
        edge [fontname="Sans", fontsize="6"];

        l2g_py [label="l2g Py", shape="component", target="_blank"];
        l2g_cpp [label="l2g Cpp", shape="component", target="_blank"];

        output_graphics [label="graphics", shape="note",];
        output_data [label=".MED data", shape="box3d"]

        input_equil [label="Equilibrium data", shape="note"]
        input_mesh [label="Mesh data", shape="box3d"]

        ext_embree [label="Embree", shape="component", href="https://www.embree.org/", target="_blank"]
        ext_alglib [label="Alglib++", shape="component", href="https://www.alglib.net/", target="_blank"]
        ext_rkf45 [label="RKF45", shape="component"]

        py_plot [label="Plotting", shape="component"]
        py_cython [label="Cython+OpenMP", shape="component", href="https://cython.org/"]
        py_mesh [label="Mesh I/O (med, vtk)", shape="component"]

        subgraph cluster_0 {
            color=lightgrey;

            input_equil;
            input_mesh;
            label="Input data";
        }

        subgraph cluster_1 {
            rank=UD;

            l2g_py -> l2g_cpp [label="Call FLT", constraint=false];
            l2g_cpp -> l2g_py [label="Return results", constraint=false];
            label="L2G code";
        }
        subgraph cluster_2 {
            style=filled;
            color=yellow;

            "output_graphics";
            "output_data";
            label="Output results";
        }
        subgraph cluster_3 {
            color=lightgreen;

            {ext_embree ext_alglib ext_rkf45};
            label="External libraries";
        }

        subgraph cluster_4 {
            {rank = same; py_plot py_cython py_mesh};
            label="Python modules"
        }

        input_equil -> l2g_py [ltail=cluster_0]
        l2g_py -> output_graphics [lhead=cluster_2, label="Write"]
        ext_alglib -> l2g_cpp [ltail=cluster_3, style=dashed, arrowhead=none]
        py_cython -> l2g_py [ltail=cluster_4, style=dashed, arrowhead=none]
    }


******************
External libraries
******************

There are many external libraries used by the python module, from the
:term:`FLT` kernel to support for mesh formats (mainly MED).

MEDCOUPLING
===========

The main data format for meshes or geometries is the MED format. The MED
format, powered by HDF5, offers an intuitive interface and performance for
reading/writing data. This means that through Python we can write or read a
massive amount of data with the use of numpy arrays and we leave to
MEDCOUPLING to write/read numpy arrays. From experience MEDCOUPLING is more
intuitive to use with better performance than VTK.

FLT Kernel
==========

The kernel of the code (where field line tracing, :term:`FLT`, is performed) is
written in C++. Parallelization is performed in the Cython wrapping. In order
to achieve thread safety in the kernel, each thread has it's local data and
objects stored in vectors or containers. Each thread uses it's OpenMP ID as
address for accessing and storing data.


Embree
------

For :term:`FLT` we require Finite Ray-Tracing, since we do not have infinite
rays, but segmented rays (field-lines) for which we would like to see if during
the tracing it hits any of the shadow geometry. Embree is a Ray-Tracing
library, with an impressive performance and simple API to use in code.

.. todo::

   Benchmark Embree performance in order to justify the word impressive
   performance.

Alglib
------

For tracing a field-line interpolation methods are used on the input data,
mainly the poloidal magnetic flux map. The 2D bicubic interpolation method of
Alglib is used in this regard.

RKF45
-----

An implementation of the method RKF45 (Runge-Kutta-Fehlberg 45) is used for
solving the field-line equations. The reason bpehind is its performance and
accuracy and the feature of it's step-adaptivity. Since in :term:`FLT` curves
are being traced, depending on the input data, user should have the power to
set the resolution, e.g., the distance between each point being traced on a
field-line. Using RKF45, we can specify at which parametric time steps (in this
case the toroidal angle) we wish to obtain the next point on the field-line
trajectory. How many steps the method might actually need to go from the
current parametric time to the next one is handled by the algorithm, but in the
results we will obtain consistent field-line points.

.. _why_openmp:

Why OpenMP
==========

OpenMP is used for parallelization of the C++ code. The calls for OpenMP is
performed on Cython side and not in the C++ code. With this the C++ code can be
a simple, yet smart enough implementation that can be called from OpenMP
threads. With identifying which data is required locally by each thread, we can
create vectors in which each thread uses it's own designated thread ID as
location for reading and writing data.

As would be in C++ the way to activate parallel blocks in the code is simple in
C++. This is also one of the reasons for using OpenMP instead of relying on
OpenMPI. Even if with OpenMP we sacrifice the option of having multiple compute
nodes running one case, since FLT is an embarrassingly parallel problem, this
can be easily solved by partitioning the input target geometry and running a
case as a multi-part case on separate compute nodes. Of course this comes with
the drawbacks of multiple loading of the same geometry, again this can be
mitigated with either having a background service waiting for order or
different implementation of the workflow.