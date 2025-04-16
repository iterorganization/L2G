"""Using BFGS to find the minimum of the magnetic poloidal flux.

This example shows (and tests) the BFGS by using different starting position
and showing how the method progresses.

Run, for exapmles, as python3 example_bfgs.py
"""

import l2g.equil

eqdsk_file_path = "../shot_135011_run_7_t_399.927598s.eqdsk"

# Get equilibrium
equilibrium = l2g.equil.getEquilibriumFromEQFile(eqdsk_file_path)

from l2g.external.bicubic import PyBicubic

bicubic = PyBicubic(equilibrium.grid_r, equilibrium.grid_z, equilibrium.psi)

from l2g.external.bfgs_2d import PyBfgs2d

min_obj = PyBfgs2d()
min_obj.setInterpolator(bicubic)

# Get bounds
min_r = min(equilibrium.grid_r)
max_r = max(equilibrium.grid_r)
min_z = min(equilibrium.grid_z)
max_z = max(equilibrium.grid_z)
min_obj.setBounds(min_r, max_r, min_z, max_z)

half_r = (min_r + max_r) * 0.5
half_z = (min_z + max_z) * 0.5

low_half_r = (min_r + half_r) * 0.5
low_half_z = (min_z + half_z) * 0.5

upp_half_r = (half_r + max_r) * 0.5
upp_half_z = (half_z + max_z) * 0.5

points = [
    (half_r, half_z), (low_half_r, low_half_z), (upp_half_r, upp_half_z)
]

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt


f = plt.figure()
ax = f.add_subplot()

cnt = ax.contourf(equilibrium.grid_r, equilibrium.grid_z, equilibrium.psi, 30, cmap="jet")
plt.colorbar(cnt, orientation='vertical')

minimum_points = []

# Plot the solution of the algorithm and all the intermediate points.
for (r, z) in points:
    out_of_bounds, *p = min_obj.findMinimum(r, z, 100)
    if out_of_bounds:
        raise Exception()
    minimum_points.append(p)
    points = min_obj.getPoints()
    points_x = [_[0] for _ in points]
    points_y = [_[1] for _ in points]
    ax.plot(points_x[0], points_y[0], 'ro', ms=10)
    ax.plot(points_x, points_y, 'o-')

print(minimum_points)

ax.set_xlabel("R [m]")
ax.set_ylabel("R [m]")
ax.axis("equal")
plt.show()

