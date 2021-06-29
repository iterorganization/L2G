#!/usr/bin/env python3
import L2G.core
from L2G import EQDSKIO, EQ, geom_utils
import numpy as np
from time import time
import os

targetFile = os.path.expanduser('~/smiter-aux/Data/VTK/inres1.vtk')
shadowFile = os.path.expanduser('~/smiter-aux/Data/VTK/inrshad1.vtk')
eqdskFile = os.path.expanduser('~/smiter-aux/Data/Equilibrium/EQ3.eqdsk')

eqdsk = EQDSKIO(eqdskFile)
# Wall used for calculating LCFS
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

eqdsk.setRLIM(smiterauxWallMeshPoints[:,0])
eqdsk.setZLIM(smiterauxWallMeshPoints[:,1])

radialDisplacement = -0.006
# Additional limiter displacement

# RLIM = eqdsk.getRLIM()
# RLIM = [_ + 0.025 for _ in RLIM] # For Midplane parameters additional shift
# eqdsk.setRLIM(RLIM)
eq = EQ(eqdsk, rDispl=radialDisplacement)

tg_vertices, tg_cells = geom_utils.readVtkMesh(targetFile)
N_vertices = len(tg_vertices) // 3
N_cells = len(tg_cells) // 3
print(f'Vertices #: {N_vertices}')
print(f'Cells #: {N_cells}')

sh_vertices, sh_cells = geom_utils.readVtkMesh(shadowFile)

x = L2G.core.PyFLT()
x.setNDIM(eqdsk.getNW(), eqdsk.getNH())
# Poloidal component interpolator data
x.setRARR(eq._R)
x.setZARR(eq._Z)
x.setPSI(eq._psi.flatten())
x.setShift(radialDisplacement, 0.0)
# Toroidal component interpolator data
x.setFARR(eq._psi_fpol)
x.setFPOL(eq._fpol)
x.setVacuumFPOL(eq._fpol[-1])
x.prepare()

x.setAbsError(1e-4)
x.setRelError(1e-4)

embreeObj = L2G.core.PyEmbreeAccell()
# Convert the units of the mesh from mm to m
embreeObj.commitMesh(sh_vertices * 1e-3, sh_cells)
x.applyRT(embreeObj)
# Prepare Threads AFTER applying embree
x.prepareThreadContainers()
print(N_cells)

print('Starting FLT')
start = time()
timeStop = 2.5 # In radians, toroidally
timeStep = 0.01
x.setTimeSpan(timeStop, timeStep)
results = L2G.core.runFLT(x, tg_vertices, tg_cells, 8)

print(f'Finished. It took {time()-start} seconds.')


# Calculate the power deposition.


# The eqdsk is a limiter equilibrium

eq.evaluate()
print(f'Eqdsk type: {eq.type_}.')
print(f'LCFS: {eq.psiLCFS}. Contact point: {eq.contactPoint}')
print(f'O point: {eq.oPoint}')
print(f'Evaluated BTCen: {eqdsk.getFPOL()[0] / eq.oPoint[0]}')
psiLCFS = eq.psiLCFS
# psiLCFS = 11.104981885689330
Rb, Btotal, Bpm = eq.getIWL_midplane(psiLCFS)
print(f'Rb = {Rb}, Btotal = {Btotal}, Bpm = {Bpm}')
drsep = (psiLCFS - results.flux) / (Rb * Bpm)
P_sol = 7.5e6
q = np.zeros(N_cells, dtype=np.float32)
lambda_q = 0.012
K = P_sol * 0.5 / (2 * np.pi * Rb * Bpm * lambda_q)
print(K * Btotal)
q = np.where(results.mask == 0,
             K * results.direction * results.Bdot * np.exp(-drsep / lambda_q),
             0)


# Setting incident angle
results.angle = np.rad2deg(results.angle)
results.angle = np.where(results.angle > 90.0, results.angle - 90.0,
                         90.0 - results.angle)

import vtk
from vtk.util import numpy_support


# Write mask
data = vtk.vtkPolyData()
points = vtk.vtkPoints()
cells = vtk.vtkCellArray()

scalarBaryCent = numpy_support.numpy_to_vtk(results.baryCent.reshape(N_cells, 3))
scalarBaryCent.SetName('Bary center')
scalarNormals = numpy_support.numpy_to_vtk(results.normals.reshape(N_cells, 3))
scalarNormals.SetName('Normals')
scalarFlux = numpy_support.numpy_to_vtk(results.flux)
scalarFlux.SetName('Flux')
scalarBVec = numpy_support.numpy_to_vtk(results.BVec.reshape(N_cells, 3))
scalarBVec.SetName('BVec')
scalarBDot = numpy_support.numpy_to_vtk(results.Bdot)
scalarBDot.SetName('BdotN')
scalarAngle = numpy_support.numpy_to_vtk(results.angle)
scalarAngle.SetName('Angle')
scalarMask = numpy_support.numpy_to_vtk(results.mask)
scalarMask.SetName('Mask')
scalarDirection = numpy_support.numpy_to_vtk(results.direction)
scalarDirection.SetName('Direction')
scalarQ = numpy_support.numpy_to_vtk(q)
scalarQ.SetName('Q')
scalarDrsep = numpy_support.numpy_to_vtk(drsep * 1000)
scalarDrsep.SetName('drsep')
scalarConlen = numpy_support.numpy_to_vtk(results.conlen)
scalarConlen.SetName('conlen')
for i in range(N_vertices):
    points.InsertPoint(i, tg_vertices[i*3:i*3+3])

for i in range(N_cells):
    triangle = vtk.vtkTriangle()
    pointIds = triangle.GetPointIds()

    pointIds.SetId(0, tg_cells[3*i])
    pointIds.SetId(1, tg_cells[3*i+1])
    pointIds.SetId(2, tg_cells[3*i+2])
    cells.InsertNextCell(triangle)


data.SetPoints(points)
data.SetPolys(cells)
cellData = data.GetCellData()
cellData.AddArray(scalarBaryCent)
cellData.AddArray(scalarNormals)
cellData.AddArray(scalarFlux)
cellData.AddArray(scalarBVec)
cellData.AddArray(scalarBDot)
cellData.AddArray(scalarAngle)
cellData.AddArray(scalarMask)
cellData.AddArray(scalarDirection)
cellData.AddArray(scalarQ)
cellData.AddArray(scalarDrsep)
cellData.AddArray(scalarConlen)

writer = vtk.vtkPolyDataWriter()
writer.SetInputData(data)
writer.SetFileName('test.vtk')
writer.Write()

