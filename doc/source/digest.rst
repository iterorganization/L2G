
####################
Digesting the module
####################

.. graphviz::
   :name: l2g workflow
   :caption: Workflow of the l2g code
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
