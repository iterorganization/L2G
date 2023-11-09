from typing import Dict, Tuple
import numpy as np
from scipy.optimize import bisect
from l2g.external.bicubic import PyBicubic

def calculate_length(paths: list) -> float:
    l = 0.0
    for path in paths:
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

class Segment:
    def __init__(self):
        self._id = None, None
        self.p1 = None
        self.p2 = None
        self.points: Tuple = None, None
        self.next = None
        self.i = None
        self.j = None

        # Whether this is a starting segment
        self.start = False

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

    def getContourPath(self, val: float) -> list:
        """Obtain a collection of all contours on a 2D grid.

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

            Arguments:
                x (np.ndarray): X axis data
                y (np.ndarray): Y axis data
                data (np.ndarray): 2D array of data
                val (float): Value of the contour
        """
        if not self.c_data_set:
            return []

        mask = np.zeros(self.c_f.shape, dtype=np.uint8)
        # Set elements to 1, when the value of the data array is larger than value.
        mask[self.c_f > val] = 1

        # Go over every 2x2 cell
        c = 0

        def fun_hor(x, y):
            return self.c_interpolator(x, y)[0] - val

        def fun_ver(y, x):
            return self.c_interpolator(x, y)[0] - val

        segment_map: Dict[Tuple, Segment] = {}
        saddle_segments: List[Segment] = []

        x = self.c_x
        y = self.c_y
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

                c += 1
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
                        p2 = bisect(fun_hor, x[i], x[i+1], args=(y[i])), y[j]
                        sad_obj2.points = p1, p2
                    else: # ind == 10
                        sad_obj1.next = i + 1, j
                        sad_obj2.next = i - 1, j
                        p1 = bisect(fun_hor, x[i], x[i+1], args=(y[j+1])), y[j+1]
                        p2 = x[i+1], bisect(fun_ver, y[j], y[j+1], args=(x[i+1]))
                        sad_obj1.points = p1, p2

                        p1 = bisect(fun_hor, x[i], x[i+1], args=(y[i])), y[j]
                        p2 = x[i], bisect(fun_ver, y[j], y[j+1], args=(x[i]))
                        sad_obj2.points = p1, p2
                    saddle_segments.append(sad_obj1)
                    saddle_segments.append(sad_obj2)

                    continue

                obj = Segment()
                obj.i = i
                obj.j = j
                segment_map[(i, j)] = obj
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
        while segment_map:
            segment_id = next(iter(segment_map))
            el = segment_map.pop(segment_id)
            # print(f'Starting from: {segment_id}')
            x_out = []
            y_out = []
            processed_segments = [segment_id]
            while 1:
                x_out.append(el.points[0][0])
                y_out.append(el.points[0][1])
                # If everything works correctly then the second point is not
                # actually needed.
                # x_out.append(el.points[1][0])
                # y_out.append(el.points[1][1])
                if el.next is None:
                    break
                if el.next not in segment_map:
                    # print(f"Element {el.next} not in segment map")
                    # Closing
                    if el.next in processed_segments:
                        # print("Closing circle")
                        # Adding the first point. If the element was processed
                        # It most certainly is the first.
                        x_out.append(x_out[0])
                        y_out.append(y_out[0])
                    break
                # print(f"Next: {el.next}")
                processed_segments.append(el.next)
                el = segment_map.pop(el.next)

            paths.append((x_out, y_out))

        return paths

def plot_paths(ax, paths, *args, **kwargs) -> None:
    label = False
    if 'label' in kwargs:
        label = True

    for path in paths:
        ax.plot(path[0], path[1], *args, **kwargs)

        ax.plot(path[0][0], path[1][0], "bo", ms=12)
        ax.plot(path[0][-1], path[1][-1], "yo", ms=12)

        if label:
            label = False
            kwargs.pop('label')
        break
    return None


if __name__ == "__main__":

    r = np.linspace(0, 50, 59)
    z = np.linspace(0, 50, 129)

    # Use the following when overlaying with IMSHOW
    OVERLAY_IMSHOW = False

    if OVERLAY_IMSHOW:
        r = np.array([i for i in range(51)])
        z = np.array([i for i in range(51)])
    rr, zz = np.meshgrid(r, z)
    print(f'rr shape={rr.shape}')

    # data = 200 / np.sqrt(((rr - 25)**2 + (zz-25)**2 + 0.1))
    data = np.sqrt(((rr - 25)**2 + (zz-25)**2))

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

    val = 5
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
    plot_paths(ax, contour_paths, 'ro-')
    ax.grid()
    ax.axis("equal")
    plt.show()
