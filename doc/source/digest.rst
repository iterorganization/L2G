
##################
Digesting the code
##################

********
Workflow
********

The way the code is structured is to have a singular python interface,
providing all the tools a user needs or might need in order to perform their
FLT studies. From the, obviously, running a :term:`FLT` case to plotting
utilities for intermediate data and diagnostics of their input equilibrium data.

Python is used for interface due to it's availability and broad specter of
libraries, which brings functionality and relatively high productivity.

For high performance and compatibility with python, C++ is used for the kernel,
e.g., the heart of the FLT computation and Cython is used to glue together
externally used C++ code with python. Parallelization is performed on
C++ side with the use of OpenMP. Reason for using OpenMP is discussed
in :ref:`why_openmp`.

.. graphviz::
   :name: l2g workflow
   :caption: L2G scheme
   :align: center

    digraph "l2g-workflow" {
        size="8,6";
        rankdir="TD";
        splines=false;
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

        ext_tinybvh [label="TinyBvh", shape="component", href="https://github.com/jbikker/tinybvh", target="_blank"]
        ext_rkf45 [label="RKF4(5)", shape="component"]

        py_plot [label="Plotting", shape="component"]
        py_cython [label="Cython", shape="component", href="https://cython.org/"]
        py_mesh [label="Mesh I/O (med, vtk)", shape="component"]
        py_analysis [label="Plasma equilibrium analysis", shape="component"]

        subgraph cluster_1 {
            color=lightgreen;

            {ext_tinybvh ext_rkf45};
            label="External libraries";
        }

        subgraph cluster_2 {
            {py_plot py_cython py_mesh py_analysis};
            label="Python modules"
        }

        subgraph cluster_3 {
            /*rank=UD;*/

            l2g_py -> l2g_cpp [label="Call FLT", constraint=false];
            l2g_cpp -> l2g_py [label="Return results", constraint=false];
            label="L2G code";
        }

        subgraph cluster_4 {
            color=lightgrey;

            input_equil;
            input_mesh;
            label="Input data";
        }

        subgraph cluster_5 {
            style=filled;
            color=yellow;

            "output_graphics";
            "output_data";
            label="Output results";
        }


        /* Connections */
        input_equil -> l2g_py [ltail=cluster_0]
        /*l2g_py -> input_equil [dir=back, lhead=cluster_0]*/
        l2g_py -> output_graphics [lhead=cluster_2, label="Write"]
        py_cython -> l2g_py [ltail=cluster_4, style=dashed, arrowhead=none]
        ext_tinybvh -> l2g_cpp [ltail=cluster_3, style=dashed, arrowhead=none]

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
----------

The kernel of the code (where field line tracing, :term:`FLT`, is performed) is
written in C++. Parallelization is achieved in the C++ code using OpenMP.

tinybvh
-------

The tinybvh project is used to test for intersection with 2D surface triangular
meshes during :term:`FLT`.


RKF45
-----

An implementation of the method RKF45 (Runge-Kutta-Fehlberg 45) is used for
solving the field-line equations. It is a robust and stable algorithm for
solving non-stiff problems and it is shown to be perfect.

.. _why_openmp:

Why OpenMP
==========

OpenMP is used for parallelization of the C++ code. The calls for OpenMP is
performed on C++ code and via Cython we can control the number of threads to
run on.

The implementation is done as efficient as possible to achieve, at least,
decent scalability on a single compute node.