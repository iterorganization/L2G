import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay

N = 100
x = np.linspace(0, 4, N)
y = np.linspace(0, 4, N)
xx, yy = np.meshgrid(x, y)

z = 5*(xx - 2)**2 + 3.5*(yy - 2)**2

# z = 20 * np.exp(-(xx - 2)**2) * np.exp(-(yy - 2)**2)


FPOL_vaccuum = 5.3
fpolArray = np.array([FPOL_vaccuum for i in range(N)], dtype=np.float64)

plt.contour(x, y, z, 30)
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

from L2G.core import PyEmbreeAccell
from L2G import FieldLineTracer, enableLogging, EQDSKIO
enableLogging()

# Construct fake EQDSK
eqdsk = EQDSKIO()
eqdsk.setNH(100)
eqdsk.setNW(100)

eqdsk.setRLEFT(0)
eqdsk.setRDIM(4)

eqdsk.setZMID(2)
eqdsk.setZDIM(4)

eqdsk.setPSIRZ(z)
eqdsk.setFPOL(fpolArray)

# Necessary values
# Boundary flux value
eqdsk.setSIBRY(1.0)
# Flux at magnetic axis
eqdsk.setSIMAG(0.0)
# Position of magnetic axis
eqdsk.setRMAXIS(2)
eqdsk.setZMAXIS(2)

case = FieldLineTracer()
case.name = "synthetic"

# set target data
case.setTargetData(tg_points.flatten(), tg_tri.flatten())

# Set equilibrium data
case.eqdskio = eqdsk
case.eq.setEqdsk(case.eqdskio)

# Construct shadow
embreeObj = PyEmbreeAccell()
embreeObj.commitMesh(sh_points.flatten(), sh_tri.flatten())
embreeObj.commitMesh(tg_points.flatten(), tg_tri.flatten())

case.setEmbreeObj(embreeObj)

case.parameters.time_step = 0.01 # toroidal angle resolution
case.parameters.time_end = 2.5 # in radians, toroidally
# Shadow dimension multiplier to be used when processing data
case.parameters.target_dim_mul = 1 # default
case.parameters.num_of_threads = 4 # Number of OMP threads default
# On which side do we take the LCFS parameters. Inner midplane or outer
case.parameters.side = "iwl"
case.parameters.P_sol = 7.5e6
case.parameters.F_split = 0.5 # default
case.parameters.q_parallel = None # Deactivated as we have P_sol
case.parameters.lambda_q_main = 0.012 # Main lambda, also default
case.parameters.self_intersection_avoidance_length = 0.0005
case.parameters.time_end = 30


# Custom wall silhouette.
# Spanning in radial direction from 1 - 5 units and heigh of -2 to 8

case.eq.setCustomLimiterSilhouette([1.0, 1.0, 5.0, 5.0], [-2, 8, 8, -2])

# Finally, to commit changes
case.applyParameters() # Propagates parameters
case.loadEq() # Loads the equilibrium data to the FLT kernel


# Run FLT
case.runFltOnMesh()

# Save results
case.saveMeshToVTK(f"{case.name}.vtu")


# FL ids
fl_ids = [4662,813,5057,2625,8627,21260,14270]

case.fl_ids = fl_ids
case.options.switch_getFL_with_FLT = False # This is the default setting
# Same settings as before, nothing to change. Now let's get those Fls
case.getFL()

case.saveFlToVTK(f"{case.name}_fls.vtk")