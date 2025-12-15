############
Introduction
############

This module is used for performing FLT studies, primarily for ITER cases. It
includes:

 - :term:`FLT`, following :term:`FL` on a triangular surface mesh or points.
 - Obtaining characteristic graphs, such as outer midplane (:term:`OMP`)
   connection  length graphs.
 - Calculating relevant quantities for use in :term:`FLT` study.
 - Applying heat loads for many different phases (Start-Up, Flat-Top,
   Ramp-Down, custom, ...)

The core of the library is held within an external C++ library, with the same
name (L2G-CPP). The external library performs the numerical steps (solving
:term:`FL` equations, checking for intersection with the shadow geometry, etc)
while this library provides a high-level interface for rapid run of different
cases, via the helper binary scripts or by integrating it's API to different
workflows. The parallelization is performed inside this module with the use
of OpenMP inside Cython.

The following technologies are used for :term:`FLT` inside the external CPP
library:

  - RKF45 (Runge-Kutta-Fehlberg 45 adaptive step) for :term:`FL` equations
  - `TinyBVH <https://github.com/jbikker/tinybvh/>`_, used for intersection tests
  - `Cython <https://cython.org/>`_, used for wrapping and enabling the use of
    external C++ libraries, enabling the use of OpenMP multithreading.
  - Implementaion of `bicubic interpolation <https://en.wikipedia.org/wiki/Bicubic_interpolation/>`_,
    used for 2D spline interpolation of axisymmetric equilibrium data.