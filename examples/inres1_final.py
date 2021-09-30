#!/usr/bin/env python3
from L2G import FieldLineTracer
import L2G.meshio_utils
import L2G.core
import numpy as np
import os

targetFile = os.path.expanduser('~/smiter-aux/Data/VTK/inres1.vtk')
shadowFile = os.path.expanduser('~/smiter-aux/Data/VTK/inrshad1.vtk')
eqdskFile = os.path.expanduser('~/smiter-aux/Data/Equilibrium/EQ3.eqdsk')

tg_vertices, tg_cells = L2G.meshio_utils.readVtkMesh(targetFile)
sh_vertices, sh_cells = L2G.meshio_utils.readVtkMesh(shadowFile)

embreeObj = L2G.core.PyEmbreeAccell()
# Convert the units of the mesh from mm to m
embreeObj.commitMesh(sh_vertices * 1e-3, sh_cells)

case = FieldLineTracer()
case.name = "inres1"

case.target_vertices = tg_vertices
case.target_triangles = tg_cells
# or, the following resets any previous saved results. So use the function
# when you switch target mesh.
#case.setTargetData(tg_vertices, tg_cells)

# Set shadowing object
case.setEmbreeObj(embreeObj)

# Set parameters
# Some are default but all are to be present just to show.
# Those that are marked default can be removed and will not cause a change
case.parameters.plasma_r_displ = -0.006
case.parameters.plasma_z_displ = 0.0 # default
case.parameters.time_step = 0.01 # toroidal angle resolution
case.parameters.time_end = 2.5 # in radians, toroidally
# Shadow dimension multiplier to be used when processing data
case.parameters.target_dim_mul = 1e-3 # default
case.parameters.num_of_threads = 1 # Number of OMP threads default
# On which side do we take the LCFS parameters. Inner midplane or outer
case.parameters.side = "iwl"
case.parameters.P_sol = 7.5e6
case.parameters.F_split = 0.5 # default
case.parameters.q_parallel = None # Deactivated as we have P_sol
case.parameters.lambda_q_main = 0.012 # Main lambda, also default
case.parameters.lambda_q_near = 0.005 # Not used in this case (single exponential)
#The following settings are not used
case.parameters.R_q = 1
case.parameters.R_breakpoint = 0.05



case.loadEqdsk(eqdskFile)

# Set custom Limiter silhouette for calculating LCFS
smiterauxWallMeshPoints = np.asarray([
 [4.173, -2.512],
 [4.092,-1.507],
 [4.062, -.492],
 [4.040, .523],
 [4.040, 1.538],
 [4.070, 2.553],
 [4.120, 3.568],
 [4.320, 4.324],
 [4.897, 4.696],
 [5.746, 4.517],
 [6.578, 3.879],
 [7.390, 3.165],
 [7.894, 2.452],
 [8.259, 1.670],
 [8.383, .623],
 [8.295, -.430],
 [7.888, -1.349],
 [7.273, -2.263],
 [6.258, -3.501]])

case.eq.setCustomLimiterSilhouette(smiterauxWallMeshPoints[:,0],
                                   smiterauxWallMeshPoints[:,1])

# Finally, to commit changes
case.applyParameters() # Propagates parameters
case.loadEq() # Loads the equilibrium data to the FLT kernel

# Run FLT
case.runFltOnMesh()

# Evaluate LCFS values
case.evaluateEq()
# Apply single exponential profile
# Add drsep, q and q_par to the results.
case.applySingleExponential()

# Save results
case.saveMeshToVTK(f"{case.name}.vtu")
