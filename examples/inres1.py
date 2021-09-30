#!/usr/bin/env python3
import L2G.core
from L2G.eq import EQDSKIO, EQ
import L2G.meshio_utils
import L2G.geom_utils
import numpy as np
from time import time
import os

targetFile = os.path.expanduser('~/smiter-aux/Data/VTK/inres1.vtk')
shadowFile = os.path.expanduser('~/smiter-aux/Data/VTK/inrshad1.vtk')
eqdskFile = os.path.expanduser('~/smiter-aux/Data/Equilibrium/EQ3.eqdsk')

eqdsk = EQDSKIO(eqdskFile)
eq = EQ(eqdsk, rDispl=0.006)

tg_vertices, tg_cells = L2G.meshio_utils.readVtkMesh(targetFile)
N_vertices = len(tg_vertices) // 3
N_cells = len(tg_cells) // 3
print(f'Vertices #: {N_vertices}')
print(f'Cells #: {N_cells}')

sh_vertices, sh_cells = L2G.meshio_utils.readVtkMesh(shadowFile)

x = L2G.core.PyFLT()
x.setNDIM(eqdsk.getNW(), eqdsk.getNH())
# Poloidal component interpolator data
x.setRARR(eq._R)
x.setZARR(eq._Z)
x.setPSI(eq._psi.flatten())
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
mask = np.zeros(N_cells, dtype=np.bool)

print('Starting FLT')
start = time()

polylines = []

for ID in range(N_cells):
    # ID=1683
    trianglePointIds = tg_cells[ID * 3 : ID * 3 + 3]
    trianglePoints = np.empty(9, np.float)
    p1 = 3*tg_cells[ID*3]
    p2 = 3*tg_cells[ID*3 + 1]
    p3 = 3*tg_cells[ID*3 + 2]
    trianglePoints[:3] = tg_vertices[p1: p1 + 3]
    trianglePoints[3:6] = tg_vertices[p2: p2 + 3]
    trianglePoints[6:] = tg_vertices[p3: p3 + 3]
    # Get starting position or barycenter
    elementBary = L2G.geom_utils.getBaryCenter(trianglePoints, coord='cyl', mul=1e-3)
    # Get normal for direction guess
    elementNorm = L2G.geom_utils.getTriangleNormal(trianglePoints)
    elementFlux = eq._psiSpl.ev(elementBary[1], elementBary[0])
    elementBVec = eq.getBCart(*elementBary)
    elementBdot = np.dot(elementNorm, elementBVec)
    # Initial direction for FL from toroidal component
    direction = np.sign(eq._fpolSpl(elementFlux))
    # Direction correction determined by dot(BVec, Norm)
    if elementBdot < 0:
        direction = -direction

    x.setTimeSpan(0.5, 0.01)
    x.setDirection(direction)
    x.setIV(*elementBary)
    points = x.getFL() # [r1, z1, phi1, r2, z2, phi2, ...]
    polylines.append(points)
    x.runFLT()
    mask[ID] = x.isHit()
print(f'Finished. It took {time()-start} seconds.')

import vtk

# Write mask
data = vtk.vtkPolyData()
points = vtk.vtkPoints()
cells = vtk.vtkCellArray()
scalar = vtk.vtkFloatArray()
scalar.SetName('Mask')

for i in range(N_vertices):
    points.InsertPoint(i, tg_vertices[i*3:i*3+3])

for i in range(N_cells):
    triangle = vtk.vtkTriangle()
    pointIds = triangle.GetPointIds()

    pointIds.SetId(0, tg_cells[3*i])
    pointIds.SetId(1, tg_cells[3*i+1])
    pointIds.SetId(2, tg_cells[3*i+2])
    cells.InsertNextCell(triangle)

    scalar.InsertTuple(i, [mask[i]])

data.SetPoints(points)
data.SetPolys(cells)
data.GetCellData().SetScalars(scalar)

writer = vtk.vtkPolyDataWriter()
writer.SetInputData(data)
writer.SetFileName('test.vtk')
writer.Write()

# Write polylines

N_polylines = len(polylines)
N_points = np.sum([len(_) for _ in polylines])

polylineData = vtk.vtkPolyData()
polylinePoints = vtk.vtkPoints()
polylinePoints.SetNumberOfPoints(N_points)
polylineCells = vtk.vtkCellArray()
id_offset = 0
for line in polylines:
    r = np.asarray(line[::3]) * 1000
    z = np.asarray(line[1::3]) * 1000
    phi = np.asarray(line[2::3])

    x = r * np.cos(phi)
    y = r * np.sin(phi)
    polyline = vtk.vtkPolyLine()
    polylineIds = polyline.GetPointIds()
    polylineIds.SetNumberOfIds(len(r))
    for i in range(len(r)):
        polylinePoints.SetPoint(id_offset + i, [x[i], y[i], z[i]])
        polylineIds.SetId(i, id_offset + i)
    polylineCells.InsertNextCell(polyline)
    id_offset += len(r)
polylineData.SetPoints(polylinePoints)
polylineData.SetLines(polylineCells)

writer = vtk.vtkPolyDataWriter()
writer.SetInputData(polylineData)
writer.SetFileName('fieldlines.vtk')
writer.Write()
