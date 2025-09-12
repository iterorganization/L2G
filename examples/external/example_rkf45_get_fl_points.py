from math import cos
import l2g.equil

equilibrium = l2g.equil.getEquilibriumFromEQFile("../shot_135011_run_7_t_399.927598s.eqdsk")
from l2g.external.bicubic import PyBicubic

bicubic = PyBicubic(equilibrium.grid_r, equilibrium.grid_z, equilibrium.psi)

from l2g.external.rkf45 import PyRKF45FLT
import numpy as np
rkf45_obj = PyRKF45FLT()
rkf45_obj.setInterpolator(bicubic)
# rkf45_obj.set_vacuum_fpol(equilibrium.fpol_vacuum)
rkf45_obj.set_vacuum_fpol(equilibrium.fpol_vacuum)

R = []
Z = []
Th = []

sr = 8.18
sz = 0.5
sth = 0.0

turns = 0
pcareR = []
pcareZ = []
length = 0.0
while sth <= 100*2*np.pi:

    nr, nz, nth = rkf45_obj.run_step(sr, sz, sth, 0.01)
    length += np.sqrt(nr*nr + sr*sr - 2*nr*sr*cos(nth-sth) + (nz - sz)*(nz - sz))
    R.append(nr)
    Z.append(nz)
    Th.append(nth)

    if abs(nth - turns * 2 * np.pi) < 1e-1:
        turns += 1
        pcareR.append(nr)
        pcareZ.append(nz)
    sr = nr
    sz = nz
    sth = nth
print(length)

import matplotlib

import matplotlib.pyplot as plt
plt.plot(R, Z)
plt.axis("equal")
plt.show()

print(len(pcareR), turns)
plt.plot(pcareR, pcareZ, "o")
plt.axis("equal")
plt.show()

R = np.array(R)
Z = np.array(Z)
Th = np.array(Th)

X = np.cos(Th) * R
Y = np.sin(Th) * R

dX = X[1:] - X[:-1]
dY = Y[1:] - Y[:-1]
dZ = Z[1:] - Z[:-1]

dl = np.sum(np.sqrt(dX*dX + dY*dY + dZ*dZ))
print(f"Length from cartesian={dl}")

f = plt.figure()
ax = f.add_subplot(projection="3d")
ax.plot(X, Y, Z)
ax.plot(R, np.zeros(R.shape), Z)
ax.plot(X[0], Y[0], Z[0], "ro")
ax.plot(X[-1], Y[-1], Z[-1], "yo")
ax.set_box_aspect([1,1,1])
plt.show()
