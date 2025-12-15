if __name__ != "__main__":
    import sys
    sys.exit(0)

import matplotlib
matplotlib.use('QtAgg')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay

#Construct the 'poloidal flux map'
N = 100
x = np.linspace(0, 4, N)
y = np.linspace(0, 4, N)
xx, yy = np.meshgrid(x, y)

# Simple square function positioned in point (2, 2)

z = 5*(xx - 2)**2 + 3.5*(yy - 2)**2

# Set the strength of the FPOL value. (F=R Bt)
FPOL_vaccuum = 5.3
# Usually we only need one value, since the central solenoid magnetic field
# corresponds with 1 / r.

# Dummy value filled just for completeness. Otherwise it contains the
# information of the poloidal current function inside the plasma core.

fpolArray = np.array([FPOL_vaccuum for _ in range(N)], dtype=np.float64)

# Plot how the poloidal flux map looks
plt.contour(x, y, z, 30)
plt.xlabel("X[m]")
plt.ylabel("Z[m]")
plt.title("Synthetic poloidal flux map")
plt.show()

def createTarget():
    """Creates a target, spanning 20 degrees toroidally (-10, to 10). Angle
    resolution is 0.1 degree.

    Vertical height is 1 unit, from 1.5 to 2.5. Vertical resolution is 0.1.
    Radial position is 1 unit.

    In order to take advantage of the scipy.spatial.Delaunay, we create a 2D
    surface (Z, Phi), so that we get triangles. Otherwise we will get
    tetrahedrons if we create (x, y, z) points
    """
    degree_start = -10
    degree_step = 0.1
    degree_end = 10
    degree_num = int((degree_end - degree_start) / degree_step) + 1
    degrees = np.linspace(degree_start, degree_end, degree_num)

    vertical_start = 1.5
    vertical_step = 0.01
    vertical_end = 2.5
    vertical_num = int((vertical_end - vertical_start) / vertical_step) + 1
    verticals = np.linspace(vertical_start, vertical_end, vertical_num)

    points_z_phi = np.empty((len(degrees) * len(verticals), 2), dtype=np.float32)
    N_degrees = len(degrees)
    for i in range(len(verticals)):
        for j in range(N_degrees):
            points_z_phi[i * N_degrees + j, 0] = verticals[i]
            points_z_phi[i * N_degrees + j, 1] = degrees[j]

    tri = Delaunay(points_z_phi)

    R = 1
    points_x_y_z = np.empty((len(points_z_phi), 3), dtype=np.float32)
    for i in range(len(points_z_phi)):
        points_x_y_z[i, 0] = R * np.cos(np.deg2rad(points_z_phi[i, 1]))
        points_x_y_z[i, 1] = R * np.sin(np.deg2rad(points_z_phi[i, 1]))
        points_x_y_z[i, 2] = points_z_phi[i, 0]

    return points_x_y_z, np.asarray(tri.simplices, np.uint32)

def createShadow():
    """Creates a target, spanning 60 degrees toroidally (-30, to 30). Angle
    resolution is 0.1 degree. It does not span over target area!

    Vertical height is 3 units, from 0.5 to 3.5. Vertical resolution is 0.1.

    Radial position is 1 unit.

    In order to take advantage of the scipy.spatial.Delaunay, we create a 2D
    surface (Z, Phi), so that we get triangles. Otherwise we will get
    tetrahedrons if we create (x, y, z) points
    """
    degree_start = -30
    degree_step = 0.1
    degree_end = 30
    degree_num = int((degree_end - degree_start) / degree_step) + 1
    degrees = np.linspace(degree_start, degree_end, degree_num)

    vertical_start = 0.5
    vertical_step = 0.01
    vertical_end = 3.5
    vertical_num = int((vertical_end - vertical_start) / vertical_step) + 1
    verticals = np.linspace(vertical_start, vertical_end, vertical_num)

    points_z_phi = []
    N_degrees = len(degrees)
    for i in range(len(verticals)):
        for j in range(N_degrees):
            if 1.4 <= verticals[i] <= 2.6 and -9.9 <= degrees[j] <= 10.1:
                continue
            points_z_phi.append([verticals[i], degrees[j]])

    points_z_phi = np.asarray(points_z_phi, dtype=np.float32)

    tri = Delaunay(points_z_phi)
    R = 0.9
    points_x_y_z = np.empty((len(points_z_phi), 3), dtype=np.float32)
    for i in range(len(points_z_phi)):
        points_x_y_z[i, 0] = R * np.cos(np.deg2rad(points_z_phi[i, 1]))
        points_x_y_z[i, 1] = R * np.sin(np.deg2rad(points_z_phi[i, 1]))
        points_x_y_z[i, 2] = points_z_phi[i, 0]

    return points_x_y_z, np.asarray(tri.simplices, np.uint32)
tg_points, tg_tri = createTarget()

tr = tg_tri[0]
sh_points, sh_tri = createShadow()

# f = plt.figure()
# ax = f.add_subplot(111, projection="3d")
# ax.scatter(tg_points[:, 0], tg_points[:, 1], tg_points[:, 2]) #tg_tri.simplices.reshape(len(tg_tri.simplices) * 3))
# ax.scatter(sh_points[:, 0], sh_points[:, 1], sh_points[:, 2])

# max_range = np.array([sh_points[:, 0].max()-sh_points[:, 0].min(), sh_points[:, 1].max()-sh_points[:, 1].min(), sh_points[:, 2].max()-sh_points[:, 2].min()]).max() / 2.0

# mid_x = (sh_points[:, 0].max()+sh_points[:, 0].min()) * 0.5
# mid_y = (sh_points[:, 1].max()+sh_points[:, 1].min()) * 0.5
# mid_z = (sh_points[:, 2].max()+sh_points[:, 2].min()) * 0.5
# ax.set_xlim(mid_x - max_range, mid_x + max_range)
# ax.set_ylim(mid_y - max_range, mid_y + max_range)
# ax.set_zlim(mid_z - max_range, mid_z + max_range)
# plt.show()

import l2g
from l2g.external.tlas import PyTLAS
import l2g.comp
import l2g.equil
import l2g.mesh
l2g.enableLogging()
l2g.addStreamHandler()

# Construct fake Equilibrium

equilibrium = l2g.equil.Equilibrium()

# Grid and poloidal data
equilibrium.grid_dim_r = 100
equilibrium.grid_dim_z = 100

equilibrium.grid_r = np.linspace(0, 4, 100)
equilibrium.grid_z = np.linspace(0, 4, 100)
equilibrium.psi = z

# Toroidal component
equilibrium.fpol = fpolArray
equilibrium.fpol_flux = fpolArray
equilibrium.fpol_vacuum = FPOL_vaccuum

# Custom wall silhouette.
# Spanning in radial direction from 1 - 5 units and heigh of -2 to 8
equilibrium.wall_contour_r = [1.0, 1.0, 5.0, 5.0]
equilibrium.wall_contour_z = [-2.0, 8.0, 8.0, -2.0]

# Boundary flux value
equilibrium.psi_boundary = 1.0
# Flux at magnetic axis
equilibrium.psi_axis = 0.0
# Position of magnetic axis
equilibrium.mag_axis_r = 2
equilibrium.mag_axis_z = 2
# Direction of toroidal current
equilibrium.psi_sign = 1

case = l2g.comp.FieldLineTracer()
case.name = "synthetic"

# set target data
# Target is already in meters!
case.parameters.target_to_m = 1
case.parameters.shadow_to_m = 1
case.setTargetData(tg_points.flatten(), tg_tri.flatten())

# Set equilibrium data
case.setEquilibrium(equilibrium)

# Construct shadow
tlasObj = PyTLAS()
# Add the shadow and target meshes to TLAS.
tlasObj.commitMesh(sh_points.flatten(), sh_tri.flatten())
tlasObj.commitMesh(tg_points.flatten(), tg_tri.flatten())

# Attach the TLAS object to the case.
case.setTLASObj(tlasObj)

case.parameters.time_step = 0.01 # toroidal angle resolution
case.parameters.max_fieldline_length = 1.0 # in meters
# Shadow dimension multiplier to be used when processing data
case.parameters.num_of_threads = 8 # Number of OMP threads default
# On which side do we take the LCFS parameters. Inner midplane or outer
case.parameters.self_intersection_avoidance_length = 0.0005


# Finally, to commit changes
case.applyParameters() # Propagates parameters

# Prepare the magnetic data
case.processMagneticData()

# Run FLT
case.runFLT()

# Save results to vtk.
case.results.vertices = tg_points.flatten()
case.results.triangles = tg_tri.flatten()
l2g.mesh.save_results_to_vtk(case.results, f"{case.name}.vtu")


# FL ids
fl_ids = [4662,813,5057,2625,8627,21260,14270]

case.fl_ids = fl_ids
case.options.switch_getFL_with_FLT = False # This is the default setting
# Same settings as before, nothing to change. Now let's get those Fls
case.getFL()

l2g.mesh.save_results_to_vtk(case.fl_results, f"{case.name}_fls.vtk")