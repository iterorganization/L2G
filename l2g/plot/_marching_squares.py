from typing import Dict, Tuple, List, Optional
import numpy as np
from scipy.optimize import bisect
from l2g.external.bicubic import PyBicubic

def calculate_length(paths: list, which="all") -> float:
    l = 0.0

    path_types = paths[1]

    for i, path in enumerate(paths[0]):

        # Sometimes we wish to calculate the length of only closed contours.
        if not which=="all" and path_types[i] != which:
            continue

        pathx = path[0]
        pathy = path[1]
        x0 = pathx[0]
        y0 = pathy[0]
        for i in range(1, len(pathx)):
            x = pathx[i]
            y = pathy[i]
            dx = x - x0
            dy = y - y0
            l += np.sqrt(dx*dx + dy*dy)

            x0 = x
            y0 = y
    return l

def calculate_length_longest(paths: list) -> float:
    l = 0.0

    path_types = paths[1]

    for i, path in enumerate(paths[0]):
        pathx = path[0]
        pathy = path[1]
        x0 = pathx[0]
        y0 = pathy[0]
        current_length = 0.0
        for i in range(1, len(pathx)):
            x = pathx[i]
            y = pathy[i]
            dx = x - x0
            dy = y - y0
            current_length += np.sqrt(dx*dx + dy*dy)

            x0 = x
            y0 = y
        if current_length > l:
            l = current_length
    return l

class Segment:
    def __init__(self):
        self.points: Tuple = None, None
        self.next: Optional[Tuple[int, int]] = None
        self.id: Tuple[int, int] = None

        # Whether this is a starting segment
        self.start = False

    def __repr__(self):
        return f"{self.i}, {self.j}\n{self.points}"

class Marching(object):
    def __init__(self):
        pass
        self.c_interpolator: PyBicubic
        self.c_interpolator_set: bool = False
        self.c_data_set: bool = False

        self.c_x: np.ndarray
        self.c_y: np.ndarray
        self.c_f: np.ndarray

    def setData(self, x, y, f):
        self.c_x = x
        self.c_y = y
        self.c_f = f
        self.setInterpolator(PyBicubic(x, y, f))
        self.c_data_set = True

    def setInterpolator(self, i: PyBicubic):
        self.c_interpolator = i
        self.c_interpolator_set = True

    def getContourPath(self, val: float) -> Tuple[list, list]:
        r"""Obtain a collection of all contours on a 2D grid.

        The following orientation is used for obtaining segments from cells.
        Direction of segments point so that the vertex with the higher value is
        always on the right.



                         e1
                     1       2
                     o-------o
          -          |       |
          j       e4 |       | e2
          v          |       |
                     o-------o
                     4       3
                         e3

                        |i>

        For reconstructing segments into a continuous polyline, all cells will
        determine the next neighboring segment the same way. A positive
        rotating arrow in the cell and the last vertice that has the function
        value higher than the contour specify the direction of the next
        neighboring segment.

        For instance. If we take the first case 1, where the vertice 1 is
        full (marked x), then the next sell will be the cell on the left
        (in the sketch):

                        i - 1, j

        On the diagram below.
                             e1
                     1                2
                     X---------------O
          -          |   _______     |
          j       e4 |  /       \    |
          v          | |         |   | e2
                     | v         |   |
                     |           |   |
                     |   _______/    |
                     |               |
                     O---------------O
                     4               3
                         e3

                             |i>

        Let's take the case where the mask index is 14, meaning vertices
        2,3 and 4 have higher function value than the value of the contour.

        Then the next neighboring segment index will be the cell above (in
        the sketch):

                        i, j-1

                                 e1
                     1                2
                     O---------------X
          -          |   _______     |
          j       e4 |  /       \    |
          v          | v         |   | e2
                     |           |   |
                     |           |   |
                     |  \_______/    |
                     |               |
                     X---------------X
                     4               3
                         e3

                             |i>

        Arguments:
            val (float): Value of the contour
        """
        if not self.c_data_set:
            return []

        mask = np.zeros(self.c_f.shape, dtype=np.uint8)
        # Set elements to 1, when the value of the data array is larger than value.
        mask[self.c_f > val] = 1

        # Go over every 2x2 cell

        def fun_hor(x, y):
            return self.c_interpolator(x, y)[0] - val

        def fun_ver(y, x):
            return self.c_interpolator(x, y)[0] - val


        segment_map: Dict[Tuple, Segment] = {}
        saddle_segments: List[Segment] = []

        x = self.c_x
        y = self.c_y

        # Edges indexes
        columns_x = x.size - 1
        rows_y = y.size - 1

        for i in range(x.shape[0] - 1):
            for j in range(y.shape[0] - 1):
                cell = mask[j: j+2, i: i+2]
                ind = 0

                # print(i, j, x.shape[0] -1, y.shape[0] -1)

                # Get the index mask
                ind += cell[0, 0]
                ind += cell[0, 1] << 1
                ind += cell[1, 1] << 2
                ind += cell[1, 0] << 3

                # Obtain segments
                if ind == 0 or ind == 15:
                    continue
                if ind == 5 or ind == 10:
                    # Create Saddle segments
                    sad_obj1 = Segment()
                    sad_obj1.i = i
                    sad_obj1.j = j
                    sad_obj2 = Segment()
                    sad_obj2.i = i
                    sad_obj2.j = j
                    # Do not map it in the segment map

                    if ind == 5:
                        sad_obj1.next = i, j + 1
                        sad_obj2.next = i, j - 1

                        p1 = x[i], bisect(fun_ver, y[j], y[j+1], args=(x[i]))
                        p2 = bisect(fun_hor, x[i], x[i+1], args=(y[j+1])), y[j+1]
                        sad_obj1.points = p1, p2

                        p1 = x[i+1], bisect(fun_ver, y[j], y[j+1], args=(x[i+1]))
                        p2 = bisect(fun_hor, x[i], x[i+1], args=(y[j])), y[j]
                        sad_obj2.points = p1, p2
                    else: # ind == 10
                        sad_obj1.next = i + 1, j
                        sad_obj2.next = i - 1, j
                        p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j+1])), y[j+1]
                        p2 = x[i+1], bisect(fun_ver, y[j], y[j+1], args=(x[i+1]))
                        sad_obj1.points = p1, p2

                        p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j])), y[j]
                        p2 = x[i], bisect(fun_ver, y[j], y[j+1], args=(x[i]))
                        sad_obj2.points = p1, p2
                    saddle_segments.append(sad_obj1)
                    saddle_segments.append(sad_obj2)

                    continue

                obj = Segment()
                obj.id = (i, j)

                segment_map[obj.id] = obj

                if ind == 1:
                    # Segment going from e1 to e4
                    obj.next = i - 1, j
                    p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j])), y[j]
                    p2 = x[i], bisect(fun_ver, y[j], y[j+1], args=(x[i]))
                elif ind == 2:
                    # Segment going from e2 to e1
                    obj.next = i, j - 1
                    p1 = x[i+1], bisect(fun_ver, y[j], y[j+1], args=(x[i+1]))
                    p2 = bisect(fun_hor, x[i], x[i+1], args=(y[j])), y[j]
                elif ind == 3:
                    # Segment going from e2 to r4
                    obj.next = i - 1, j
                    p1 = x[i+1], bisect(fun_ver, y[j], y[j+1], args=(x[i+1]))
                    p2 = x[i], bisect(fun_ver, y[j], y[j+1], args=(x[i]))
                elif ind == 4:
                    # Segment going from e3 to r2
                    obj.next = i + 1, j
                    p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j+1])), y[j+1]
                    p2 = x[i+1], bisect(fun_ver, y[j], y[j+1], args=(x[i+1]))
                elif ind == 6:
                    obj.next = i, j - 1
                    p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j+1])), y[j+1]
                    p2 = bisect(fun_hor, x[i], x[i+1], args=(y[j])), y[j]
                elif ind == 7:
                    obj.next = i - 1, j
                    p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j+1])), y[j+1]
                    p2 = x[i], bisect(fun_ver, y[j], y[j+1], args=(x[i]))
                elif ind == 8:
                    obj.next = i, j + 1
                    p1 = x[i], bisect(fun_ver, y[j], y[j+1], args=(x[i]))
                    p2 = bisect(fun_hor, x[i], x[i+1], args=(y[j+1])), y[j+1]
                elif ind == 9:
                    obj.next = i, j + 1
                    p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j])), y[j]
                    p2 = bisect(fun_hor, x[i], x[i+1], args=(y[j+1])), y[j+1]
                elif ind == 11:
                    obj.next = i, j + 1
                    p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j+1])), y[j+1]
                    p2 = x[i+1], bisect(fun_ver, y[j], y[j+1], args=(x[i+1]))
                elif ind == 12:
                    obj.next = i + 1, j
                    p1 = x[i], bisect(fun_ver, y[j], y[j+1], args=(x[i]))
                    p2 = x[i+1], bisect(fun_ver, y[j], y[j+1], args=(x[i+1]))
                elif ind == 13:
                    obj.next = i + 1, j
                    p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j])), y[j]
                    p2 = x[i+1], bisect(fun_ver, y[j], y[j+1], args=(x[i+1]))
                elif ind == 14:
                    obj.next = i, j - 1
                    p1 = x[i], bisect(fun_ver, y[j], y[j+1], args=(x[i]))
                    p2 = bisect(fun_hor, x[i], x[i+1], args=(y[j])), y[j]
                obj.points = p1, p2
                obj.ind = ind


        # print(segment_map)
        # segments = list(segment_map.keys())
        # print(segments)
        # Start with one

        # Construct paths.

        paths = []
        types = []
        starting_ids_of_paths = []

        # Process the starting segments on edge
        while segment_map:
            segment_id = next(iter(segment_map))
            segment = segment_map.pop(segment_id)
            # print(f'Starting from: {segment_id}')
            points_x = []
            points_y = []
            starting_id = segment_id
            path_type = "segment"
            new_path = True

            while 1:
                points_x.append(segment.points[0][0])
                points_y.append(segment.points[0][1])
                # If everything works correctly then the second point is not
                # actually needed.
                # points_x.append(el.points[1][0])
                # points_y.append(el.points[1][1])
                if segment.next is None:
                    break

                if segment.next not in segment_map:
                    # print(f"Element {segment.next} not in segment map")
                    # Closing
                    if segment.next == starting_id:
                        # print("Closing circle")
                        # Adding the first point. If the element was processed
                        # It most certainly is the first.
                        points_x.append(points_x[0])
                        points_y.append(points_y[0])
                        path_type = "closed"
                    else:
                        # Check other paths if we connect to them.
                        for i, previous_starting_id in enumerate(starting_ids_of_paths):
                            if previous_starting_id == segment.next:
                                new_path = False
                                # Let's append and insert the path there.

                                points_x += paths[i][0]
                                points_y += paths[i][1]

                                starting_ids_of_paths[i] = starting_id
                                paths[i] = (points_x, points_y)
                                break

                        if not new_path:
                            break
                    break
                segment_id = segment.next
                segment = segment_map.pop(segment.next)

            if new_path:
                paths.append((points_x, points_y))
                types.append(path_type)
                starting_ids_of_paths.append(starting_id)

        return paths, types

def plot_paths(ax, paths: Tuple[list, list], *args, **kwargs) -> None:
    label = False
    if 'label' in kwargs:
        label = True
    for path in paths[0]:
        ax.plot(path[0], path[1], *args, **kwargs)

        # ax.plot(path[0][0], path[1][0], "bo", ms=12)
        # ax.plot(path[0][-1], path[1][-1], "yo", ms=12)

        if label:
            label = False
            kwargs.pop('label')
    return None


if __name__ == "__main__":

    r = np.linspace(0, 50, 100)
    z = np.linspace(0, 50, 100)

    # Use the following when overlaying with IMSHOW
    OVERLAY_IMSHOW = True

    if OVERLAY_IMSHOW:
        r = np.array([i for i in range(51)])
        z = np.array([i for i in range(51)])
    rr, zz = np.meshgrid(r, z)
    print(f'rr shape={rr.shape}')

    # data = 200 / np.sqrt(((rr - 25)**2 + (zz-25)**2 + 0.1))

    # Creating a function with a saddle point.

    # First some elipse shapped contours
    data = np.sqrt((0.5*(rr - 35)**2 + 0.1*(zz-25)**2)) * \
           5e-3*((rr - 15)**2 + (zz - 25)**2)

    march = Marching()
    march.setData(r, z, data)
    print(data.shape)


    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt

    f = plt.figure()
    ax = f.add_subplot(111, projection="3d")
    ax.contour3D(rr, zz, data, 50)
    plt.show()

    # contour_paths = getContourPath(r, z, data, val)
    # val = 15
    # contour_path = march.getContourPath(val)

    # for x in contour_path:
    #     print(x)

    # print(f"{len(contour_path)=}")
    # # print(f"{calculate_length(contour_path)=}")
    # # print(contour_path)
    # plot_paths(plt, contour_path, 'r-')

    val = 4.5
    contour_paths = march.getContourPath(val)
    print(f"{len(contour_paths)=}")
    print(f"{calculate_length(contour_paths)=}")
    # print(contour_paths)

    f = plt.figure()
    ax = f.add_subplot()

    mask = np.zeros(data.shape, dtype=np.uint8)
    mask[data > val] = 1
    if OVERLAY_IMSHOW:
        ax.matshow(mask)
    else:
        # plt.contourf(r, z, mask)
        pass
    plot_paths(ax, contour_paths, 'r-')
    ax.grid()
    ax.axis("equal")
    plt.show()
