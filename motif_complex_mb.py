import itertools as it
import numpy as np
import math
import operator


def main():
    tri_lattice = {"a": 1,
               "b": 1,
               "alpha": 60}

    rect_lattice = {"a": 1,
               "b": 3,
               "alpha": 90}

    motif_1 = {1: (0.0, 0.0)}
    motif_2 = {1: (0.0, 0.0),
               2: (0.5, 0.5)}
    motif_3 = {1: (0.0, 0.0),
               2: (0.5, 0.5),
               3: (0.8, 0.8)}

    xtal = Crystal2d(tri_lattice, motif_1)
    print(xtal.motif_graph_2d(3))


class Crystal2d():
    """Generates 2d crystal as lattice translations and labelled motif points

    Attributes:
        lattice (dict): Dictinary of the lattice attibutes, a, b and alpha
        motif (dict): Dictionary of labelled motif points in fractional coordinates

    """
    def __init__(self, lattice, motif):
        """
        Args: 
            lattice (dict): Dictinary of the lattice attibutes, a, b and alpha
            motif (dict): Dictionary of labelled motif points in fractional coordinates
        """
        self.lattice = lattice
        self.motif = motif
        self.cartesian_motif = self.create_cartesian(lattice, motif)
        self.motif_ord = len(motif)
        self.vertices = list(motif.keys())
        self.unit_cell_volume = lattice["a"] * lattice["b"]
        self.expanded_motif = {}

    '---------------------------------------------------------------------------'
    'Outputs a list of just all point positions in a crystal'
    'Used to test for distance maximum and count translations'
    '---------------------------------------------------------------------------'
    def getpoints(self):
        allpoints = []
        for i in self.vertices:
            if i in self.expanded_motif:
                allpoints.extend(self.expanded_motif[i])
            else: allpoints.extend(self.motif[i])

        return allpoints

    '---------------------------------------------------------------------------'
    'Adds all points with distinct REDUCED translation indices in layers s to e'
    '---------------------------------------------------------------------------'
    def expand(self, e, s=1):

        'Generate a list of reduced translation indices in 2d starting at layer s and ending at layer e'
        translist = self._ti_list(2,e,s)

        cart_vect = self.shear(self.lattice, np.array([(self.lattice["a"], 0), (0, self.lattice["b"])]))
        expanded_motif = {}
        'Calculate lattice translate vectors from indices and add new points to crystal'
        for i in self.vertices:
            for j in translist:
                j_vect = [sum(z) for z in zip([(j[0]*i) for i in cart_vect[0]], 
                                                [j[1]*i for i in cart_vect[1]])]
                move_point = [sum(x) for x in zip(self.motif[i], j_vect)]

                if i in expanded_motif:
                    expanded_motif[i].append(move_point)
                else:
                    expanded_motif[i] = [move_point]

        self.expanded_motif = expanded_motif

        return translist

    '---------------------------------------------------------------------------'
    'Creates motif graph MG(C,s) for a particular 2d crystal and distance s'
    '---------------------------------------------------------------------------'

    def motif_graph_2d(self, s):

        'Initialise vertex and edge data'
        vert_list = {i:[[0,0]] for i in self.motif}
        edge_list = {j:[] for j in it.combinations_with_replacement(self.vertices,2)}
        # print(edge_list)

        vertrans = []

        'Expand graph in layers up to the layer beyond the target distance'
        i = 1
        while self._maxdist(self.getpoints()) < s:
            vertrans.extend(self.expand(i, i-1))
            i+=1

        'Check in-cell motif distances'
        for i in self.vertices:
            for j in self.vertices:
                d = np.linalg.norm(np.array(self.expanded_motif[j][0]) - np.array(self.expanded_motif[i][0]))
                # print(d)
                if d < s and (i,j) in edge_list.keys():
                        edge_list[(i,j)].append([[0,0],d])

        'Check all other motif distances'
        for i in self.vertices:
            for j in self.vertices:
                mot_point = self.expanded_motif[i][0]
                mot_test = self.expanded_motif[j]
                'Check distance between motif point and each translated point'
                for k in mot_test:
                    p = mot_test.index(k)-1
                    d = np.linalg.norm(np.array(k) - np.array(mot_point))
                    'Add vertex and edge data if distance less than threshold'
                    if d != 0 and d < s:
                        'Avoid duplicate reverse translatons'
                        if vertrans[p] not in vert_list[j]:
                            vert_list[j].append(vertrans[p])
                        'Avoid duplication of (i,j) and (j,i) edges (again, from reverse translations)'
                        if (i,j) in edge_list.keys():
                            edge_list[(i,j)].append([vertrans[p],d])


        return [vert_list, edge_list]

    def shear(self, lattice, points):
        transform_mat = np.array([[lattice["a"], lattice["b"] * np.cos(lattice["alpha"] * np.pi / 180)],
                                [0, lattice["a"] * np.sin(lattice["alpha"] * np.pi / 180)]])

        return np.dot(transform_mat, points.T).T

    def create_cartesian(self, lattice, motif):
        """Convert the fractional coordinates and lattice to cartesian coordinates
        """
        labels, points = zip(*[(k, p) for k, p in motif.items()])
        sheared_points = self.shear(lattice, np.array(points))

        return {l: sheared_points[i] for i, l in enumerate(labels)}

        # CURRENTLY UNUSED, MAY BE USEFUL LATER?
        # def tile(motif, lattice, n=3):
        #     """ Return the points of an nxn cell as a list of coordinates for
        #     each labelled point. The most interior cell is marked as a separate
        #     label."""
        #     return_points = []
        #     point_labels = list(motif.keys())

        #     for i in range(n):
        #         for j in range(n):
        #             for k, p in enumerate(motif):
        #                 return_points.append((p[0] + i * lattice["a"], p[1] + j * lattice["b"]))
        #                 label = point_labels[k]
        #                 if i == j == int(n/2): point_labels.append(label)
        #                 else: point_labels.append(label + "_")

        #     sheared_points = shear(lattice, np.array(return_points))
        #     return_motif = {}
            
        #     for i, p in enumerate(sheared_points):
        #         if point_labels[i] not in return_motif:
        #             return_motif[point_labels[i]] = [p]
        #         else:
        #             return_motif[point_labels[i]].append(p)

        #     return return_motif




    '-----------------------------------------------------------'
    'Creates a standardized list of reduced translation indices'
    'for concentric layers in a dim-dimensional crystal'
    'starting at layer s and ending at layer e (s defaults to 1)'
    '-----------------------------------------------------------'

    def _ti_list(self, dim, e, s=1):


        j = [k for k in range (1,e)]

        if s == 1:
            tlist = [[1,0],[0,1],[1,1],[-1,1]]
        else:
            tlist = []


        i = s
        while i <= e:
            for b in j:
                'Check that this is not a reducible index'
                if i != 1 and math.gcd(i,b) == 1:
                    'Append all distinct index permutations'
                    tlist.append([i,b])
                    tlist.append([b,i])
                    tlist.append([-b,i])
                    tlist.append([-i,b])
            i+=1

        return tlist

    
    '-----------------------------------------------------------'
    'Calculate the maximum pairwise distance of a list of points'
    'Used to test whether distance threshold reached in expansion'
    '-----------------------------------------------------------'
    def _maxdist(self, l):

        distlist = []
        if len(l) == 1:
            return 0

        else:
            for dist_pairs in it.combinations(l,2):
                vectors = [np.array(dist_pairs[0]), np.array(dist_pairs[1])]
                distlist.append(np.linalg.norm(vectors[1]-vectors[0]))

        return max(distlist)

if __name__ == "__main__":
    main()