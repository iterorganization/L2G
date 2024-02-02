import l2g.equil

equilibrium = l2g.equil.getEquilibriumFromEQFile("../shot_135011_run_7_t_399.927598s.eqdsk")

from l2g.external.bicubic import PyBicubic

bicubic = PyBicubic(equilibrium.grid_r, equilibrium.grid_z, equilibrium.psi)

print(bicubic.getValues(4.0, 4.0))
from l2g.external.bfgs_2d import PyBgfs2d

min_obj = PyBgfs2d()
min_obj.setInterpolator(bicubic)

# Get bounds
min_r = min(equilibrium.grid_r)
max_r = max(equilibrium.grid_r)
min_z = min(equilibrium.grid_z)
max_z = max(equilibrium.grid_z)
min_obj.setBounds(min_r, max_r, min_z, max_z)

x = []
# x.append(min_obj.findMinimum(6.3,0.5, 2))
# print("asd")
x.append(min_obj.findMinimum(5.0,-3.4, 10))
x.append(min_obj.findMinimum(equilibrium.mag_axis_r,equilibrium.mag_axis_z, 10))

minR = min(equilibrium.grid_r)
maxR = max(equilibrium.grid_r)
minZ = min(equilibrium.grid_z)
maxZ = max(equilibrium.grid_z)
lower_z = minZ + 0.15 * (maxZ - minZ)
print((minR+maxR)*0.5, minZ + 0.15 * (maxZ - minZ))
x.append(min_obj.findMinimum((minR+maxR)*0.5, minZ + 0.15 * (maxZ - minZ), 2))

print(x)
for _ in x:
    print(min_obj.secondDerivativeTest(_[0], _[1]))
    print(f"{_=}")