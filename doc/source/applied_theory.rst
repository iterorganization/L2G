
##############
Applied theory
##############

This section describes the core numerical equations used for :term:`FLT`.

In the following sub-sections a :term:`FLT` case is to be considered as having
a target geometry, i.e., a set of points or a mesh from which we wish to
analyze.

Additionally to this we also have a shadowing geometry, a set of meshes or it
can be nothing as well, depending what we need. So in other words a geometry
which will act as obstruction to the field lines, which has an effect to the
wetted area.

The field lines are traced backwards, meaning, we actually start from the end
of the trajectory and follow backwards to observe how long the field lines is
before it is intercepted.

.. note::

   In case of triangular surface meshes, the origin or actually end of a
   field line starts at the barycenter of the triangle. There is no sampling
   used.

.. note::

   For these sections consider the following list of abbreviations:

    - target element: One of the triangles in mesh.
    - target point: One of the specified target points.

*******************
Tracing field lines
*******************

In order to trace magnetic field lines in a tokamak the following quantities
are used from standard axisymmetric equilibrium data:

 - Poloidal magnetic flux :math:`\Psi`. This quantity is presented with single
   matrix of size N x M where N and M are positive integers. The definition
   domain is a cylindrical R x Z domain, where R represents the major radius,
   usually defined with N points and where Z represents vertical position,
   usually defined with M points. The units of the poloidal flux map is in
   Webb/rad and the positions are in meters.
 - Poloidal current function F. Defined with the following relation
   :math:`F=B_t R`, where F is the poloidal current function, :math:`B_t` is
   the toroidal component of the magnetic field and R is the major radius
   position. In the vaccum F is considered constant, meaning that in order to
   obtain the toroidal component of the magnetic field, we simply need the
   radial position of the point where we wish to evaluate the value.

.. note::

   For ITER in order to have correct helicity, the toroidal component of the
   magnetic field points Counter-Clockwise, while the gradient of the poloidal
   magnetic flux should point outside (global minimum at O-point).

For poloidal magnetic flux we use interpolation methods, usually 2D bicubic
interpolation method with the ability to evaluate derivatives inside the
computation domain. The computation domain is defined with the R and Z points.

The magnetic field lines are traced in a Cylindrical coordinate system with the
following differential equations:

.. math::

   \frac{\partial s_r}{\partial \zeta} & = - \frac{\partial \Psi(R, Z)}{\partial Z} \frac{R}{F(R)} \\\\
   \frac{\partial s_z}{\partial \zeta} & =   \frac{\partial \Psi(R, Z)}{\partial R} \frac{R}{F(R)}

:math:`s_r` and :math:`s_z` denotes the R and Z position of the trajectory,
while :math:`\zeta` is the toroidal angle and in this case treated as
parametric time. Following the magnetic surface or a contour of the poloidal
magnetic flux, means going perpendicular in the direction of the gradient at
any (R, Z) point. The term :math:`\frac{R}{F(R)}` or in other words
:math:`\frac{1}{B_t}` comes from following the field lines or streamlines in a
Cylindrical coordinate system.

As mentioned, the :math:`\frac{\partial \Psi(R, Z)}{\partial Z}` and
:math:`\frac{\partial \Psi(R, Z)}{\partial R}` terms are obtained by obtaining
derivatives via an interpolation method on the :math:`\Psi(R, Z)`. The set of
differential equations are solved using RKF45.

This covers just following the field lines, so the next step, obviously, is
checking if there are any collisions with the shadowing geometry. Embree is
used, where we trace a finite ray to see if we obtain hits with the
shadowing geometry.

So in order to follow a field line and check for intersections with the
shadowing geometry, the following steps are performed:

 #. Start from initial position :math:`P_0=(R_0, Z_0, \zeta_0)`.
 #. Obtain the next point :math:`P=(R, Z, \zeta)` via solving the ODE.
 #. Convert the points to Cartesian coordinate system to obtain
    :math:`P_0=(X_0, Y_0, Z_0)` and :math:`P=(X, Y, Z)`
 #. Check for intersection with the shadowing geometry using the segment
    :math:`\vec{P_0 P}` as a part of the field line. If there
    is a collision, stop tracing.
 #. If no collision, go to step 2.

Of course, we also need to check if the field line is longer than the specified
maximum connection length (to cut time on computation, as sometimes we only
need to follow a few meters) or if a field line leaves the computation (R, Z)
domain. This can happen when we have an input geometry, with a magnetic surface
going outside the tokamak with no shadowing geometry to intercept.

*****************************************
Determining if a field line wets the area
*****************************************

Now that we have a field line traced we have to determine if it wets the area.

One way is to see if a field line passes the midplane. So along the trajectory
we see if it passes the midplane, opening a path of particles to land to our
target. That involves additional check during tracing and it is open to
discussion whether this is a final solution.

So the easy way and the conservative way is that the wetted area is defined
by having a filter based on the connection lengths. This means, we specify a
minimum connection length, i.e., a few meters and *IF* a field line is longer
than this value we consider the target wetted.