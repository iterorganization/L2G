#!/usr/bin/env python3
from L2G import FieldLineTracer, geom_utils, enableLogging
import L2G.core
import os

enableLogging()

targetFile = os.path.expanduser('~/smiter-compute/KFWP4/target/target.vtk')
shadowFile = os.path.expanduser('~/smiter-compute/KFWP4/shadow/shadow.vtk')
eqdskFile = os.path.expanduser('~/smiter-compute/KFWP4/eqdsk/16_97s.eqdsk')

tg_vertices, tg_cells = L2G.geom_utils.readVtkMesh(targetFile)
sh_vertices, sh_cells = L2G.geom_utils.readVtkMesh(shadowFile)

embreeObj = L2G.core.PyEmbreeAccell()
# Convert the units of the mesh from mm to m
embreeObj.commitMesh(sh_vertices * 1e-3, sh_cells)

case = FieldLineTracer()
case.name = "nf55"

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
case.parameters.plasma_r_displ = 0.025
case.parameters.time_step = 0.01 # in radians, toroidally
case.parameters.time_end = 2.5 # toroidal angle resolution
# Shadow dimension multiplier to be used when processing data
case.parameters.num_of_threads = 8 # Number of OMP threads
# On which side do we take the LCFS parameters. Inner midplane or outer
case.parameters.side = "owl"
case.parameters.P_sol = 5e6
case.parameters.lambda_q_main = 0.0379 # Main lambda

case.loadEqdsk(eqdskFile)

# Last two fpol values are way too noisy
fpol = case.eqdskio.getFPOL()
fpol[-1] = fpol[-3]
fpol[-2] = fpol[-3]

case.eqdskio.setFPOL(fpol)

# Set custom Limiter silhouette for calculating LCFS
RLIM = case.eqdskio.getRLIM()
RLIM = [_ + 0.025 for _ in RLIM] # For Midplane parameters additional shift
ZLIM = case.eqdskio.getZLIM()
case.eq.setCustomLimiterSilhouette(RLIM, ZLIM)

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

# Now switch to calculating by some maximum connection length
case.parameters.max_connection_length = 5.0
case.options.switch_runFLT = 1
case.applyParameters()
# No need to load Eq.

# Run FLT
case.runFltOnMesh()
# No need to evaluate Eq
case.applySingleExponential()
# Save results
case.saveMeshToVTK(f"{case.name}_maxconlen.vtu")

# Now lets get some FLs
fl_ids = [92875, 92866, 94701, 82728, 56885, 52760, 49089, 43575]
case.fl_ids = fl_ids
case.options.switch_getFL_with_FLT = True # This is the default setting
# When plotting FLs, apply FLT

# Same settings as before, nothing to change. Now let's get those Fls
case.getFL()

case.saveFlToVTK(f"{case.name}_fls.vtk")