from L2G.core import PyEmbreeAccell, PyFLT
import numpy as np

def readVtkMesh(inputFileLoc):
    '''This method reads .vtk file and stores it to Mesh module. It reads
    only points and cells and stores them into Mesh module. Any
    additional fields on triangles are skipped.

    Args:
        inputFileLoc (str): Path to VTK file.
        meshName (str): New name of mesh.

    '''
    vertices = []
    edges = []
    elements = []
    with open(inputFileLoc) as f:
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        points = f.readline()
        if points.split()[0] == "POINTS":
            nr_points = int(points.split()[1])
        else:
            nr_points = int(f.readline().split()[1])
        for point in range(nr_points):
            pt = f.readline().split()
            vertices.append(float(pt[0]))
            vertices.append(float(pt[1]))
            vertices.append(float(pt[2]))
        cells = f.readline()
        if cells.split():
            nr_cells = int(cells.split()[1])
        else:
            nr_cells = int(f.readline().split()[1])
        for cell in range(nr_cells):
            cl = f.readline().split()
            # cl[0] holds the number (string) defining how many
            # points compose the element (vtkCell)
            if cl[0] == '1':  # Vertex
                edges.append([int(cl[1]) + 1])

            if cl[0] == '2':  # Edge
                edges.append([int(cl[1]) + 1,
                              int(cl[2]) + 1])

            elif cl[0] == '3':  # Triangle
                elements.append([int(cl[1]) + 1,
                                 int(cl[2]) + 1,
                                 int(cl[3]) + 1])

            elif cl[0] == '4':  # Quad
                elements.append([int(cl[1]) + 1,
                                 int(cl[2]) + 1,
                                 int(cl[3]) + 1,
                                 int(cl[4]) + 1])

    return vertices, edges, elements

x = PyEmbreeAccell()

vertices = np.array([0,0,0,1,0,0,0,1,0], dtype=np.float64)
triangles = np.array([1,2,3], dtype=np.uint32)

shadow = '/home/ITER/simicg/smiter-aux/Data/VTK/inrshad1.vtk'
vertices, _, triangles = readVtkMesh(shadow)
x.commitMesh_fromList(vertices, triangles)
print('Commited mesh')

y = PyFLT()
y.applyRT(x)
