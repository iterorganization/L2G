import l2g.comp.core
import l2g.comp
import l2g.equil


# Path to MED file! Filled by you!
mesh = "/path/to/MEd/file.med"
medmesh = l2g.comp.MEDMeshIO()
medmesh.readMeshFromMedFile(mesh)
vertices, cells = medmesh.getMeshData()

embree_obj = l2g.comp.core.PyEmbreeAccell()
embree_obj.commitMesh(vertices, cells)

equilibrium = l2g.equil.Equilibrium()
# Optional
equilibrium.Area = 1
equilibrium.Ip = 1e6
equilibrium.Psol = 1e8
equilibrium.a = 10
equilibrium.fpol = []
equilibrium.fpol_flux = []
equilibrium.psi_axis = -1
equilibrium.psi_boundary = -1
equilibrium.psi_sign = -1

# Needed, filled by you!!
equilibrium.fpol_vacuum = 15.6
equilibrium.psi = [[]]
equilibrium.grid_dim_r = N_r
equilibrium.grid_r = []
equilibrium.grid_dim_z = N_z
equilibrium.grid_z = []

# Set the field line tracer
flt = l2g.comp.FieldLineTracer()
flt.setEquilibrium(equilibrium)
flt.setTargetData(vertices, cells)
flt.setEmbreeObj(embree_obj)

flt.runFltOnMesh()