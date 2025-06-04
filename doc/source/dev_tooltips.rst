####################
Development tooltips
####################

This section provides tooltips when developing this package.

***********************
Python module structure
***********************

Currently the python module structure of the L2G package is as following:

 * l2g

    - comp

       + core

    - equil
    - hlm
    - plot
    - utils
    - settings

In the directories we can see python files that have an ``_`` prefix or not.
Files with ``_`` prefix contains classes, which are then imported or exposed
via the ``__init__.py`` in the same place.

The reason for that is to split large classes into separate files, but to avoid
imports like ``from module.class_a import class_a`` the ``__init__.py`` in the
same directory exposes the classes so that it becomes
``from module import class_a``.

We can also see that there are python files, without the prefix ``_`` which
are in the submodules. These python files, usually a collection of functions
(see for example ``l2g.hlm.general``) act as their own submodule. This means
there is no need to expose the functions from the ``__init__.py`` in the same
directory. Of course this is not necessary, we can play with namespace
designing as much as we want.

The name of all the modules is in lowercase and the submodules have short
names, but long enough for the user to understand the role:

 - comp

    + The component or computation submodule is so called the heart of l2g,
      containing calls to the external C++ kernel and python objects for
      constructing cases.
 - comp.core

    + Contains the Cython implementation, most of the numerical work is
      performed here, either by implementing functions in Cython or calling the
      external C++ code.
 - equil

    + Contains the classes for handling equilibrium data and its diagnostic.
 - hlm

    + Contains the functions to be used for heat load mapping. The reason for
      a separate submodule is to have well documented functions for applying
      the plasma profiles and to separate them of course from other parts of
      the code
 - plot

    + Plotting functions used for creating movies or display of the equilibrium
      data
 - utils

    + Utility functions, helpers for setting up cases, handling lower levels of
      i/o and the handling of files which describe cases.
 - settings

    + Contains the switches, options and parameters used when running a FLT
      case.

***********
Adding code
***********

If we wish to add a class, and it is fairly big, it should be added into it's
appropriate submodule with the "_" prefix and then exposed in the submodule
namespace by importing it in the submodule ``__init__.py``.

If we wish to add a set of functions, according to what these functions are,
for example if they complement a submodule, we can add it into the submodule
``__init__.py`` or in python files that act as submodules (python files
without the ``_`` prefix). Otherwise, we can create a new submodule, by simply
creating the correct namespace or putting the python file in its place in the
code structure.
