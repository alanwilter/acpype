import array
import math
from . import v3


class SpaceHash(object):
    """
    A geometrical object to sort a set of vertices into disjoint
    boxes to optimize pair sorting by distance.

    The boxes are defined by div, and close pairs generates listing
    of indices of the vertices that are within div apart.
    """

    def __init__(self, vertices, div=5.3, padding=0.05):
        self.vertices = vertices
        self.div = div
        self.inv_div = 1.0 / self.div
        self.padding = padding

        zero3 = lambda: [0.0] * 3
        self.minima = zero3()
        self.maxima = zero3()
        self.spans = zero3()
        self.sizes = zero3()

        for i in range(3):
            self.minima[i] = min([v[i] for v in self.vertices])
            self.maxima[i] = max([v[i] for v in self.vertices])
            self.minima[i] -= self.padding
            self.maxima[i] += self.padding
            self.spans[i] = self.maxima[i] - self.minima[i]
            self.sizes[i] = int(math.ceil(self.spans[i] * self.inv_div))

        self.size1_size2 = self.sizes[1] * self.sizes[2]
        self.size2 = self.sizes[2]

        self.cells = {}
        self.spaces = []
        for i_vertex, vertex in enumerate(self.vertices):
            space = self.vertex_to_space(vertex)
            self.spaces.append(space)
            space_hash = self.space_to_hash(space)
            cell = self.cells.setdefault(space_hash, array.array('L'))
            cell.append(i_vertex)

    def vertex_to_space(self, v):
        return [int((v[i] - self.minima[i]) * self.inv_div) for i in range(3)]

    def space_to_hash(self, s):
        return s[0] * self.size1_size2 + s[1] * self.size2 + s[2]

    def neighbourhood(self, space):
        def neighbourhood_in_dim(space, i_dim):
            i = max(0, space[i_dim] - 1)
            j = min(self.sizes[i_dim], space[i_dim] + 2)
            return list(range(i, j))
        for s0 in neighbourhood_in_dim(space, 0):
            for s1 in neighbourhood_in_dim(space, 1):
                for s2 in neighbourhood_in_dim(space, 2):
                    yield [s0, s1, s2]

    def close_pairs(self):
        n_vertex = len(self.vertices)
        for i_vertex0 in range(n_vertex):
            space0 = self.spaces[i_vertex0]
            for space1 in self.neighbourhood(space0):
                hash1 = self.space_to_hash(space1)
                for i_vertex1 in self.cells.get(hash1, []):
                    if i_vertex0 < i_vertex1:
                        yield i_vertex0, i_vertex1


