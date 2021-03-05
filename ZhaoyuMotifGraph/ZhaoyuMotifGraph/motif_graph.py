"""
This crystal object class is an application of motif calculation and other related operations, like graph
generation and spectrum calculation.
_author_ = "Matt Bright"
_credits_ = "Matt Bright,Zhaoyu Han"
"""

import itertools as it
import numpy as np
import math
import ase
import ase.io
from ase.visualize import view
import networkx as nx
import matplotlib.pyplot as plt

def main():
    """
    examples for this class
    """
    crystal = Crystal_3d()
    crystal.read_cif(".../data/T2_0.cif")
    crystal.plot()

class crystal_object:
    """Crystal object creator
    This class create a crystal object， on which can perform a series of
    operations using methods from this class. You can choose either to generate
    this object directly from a cif file using file_path, or to add it mannualy by
    fullfiling lattice and fmotif variables.
    
    Attributes:
        file_path: A String that indicates the file_path of the cif file, which is uesd for generating crystal
        lattice: A list that contains 3 list indicating the lattice of one crystal
        fmotif: A dictonary that contains fractional coordinates for points in crystal. Keys should be integers that start from 1 to 2,3,...
    """
    def __init__(self, file_path = None, lattice = None, fmotif=None):
        if file_path is not None:
            crystal_file = self.read_cif(file_path)
            self.lattice = crystal_file[0]
            self.fmotif = crystal_file[1]
        elif lattice is not None and fmotif is not None:
            self.lattice = lattice
            self.fmotif = fmotif
        else:
            print("at least one of the file_path and lattice-fmotif pair should be given")
        
        self.file_path = file_path
        self.motif_ord = len(self.fmotif)
        self.vertices = list(self.fmotif.keys())
        self.lens = [np.linalg.norm(i) for i in self.lattice]
        self.alpha = np.arccos(np.dot((self.lattice[1] / np.linalg.norm(self.lattice[1])),
                                      (self.lattice[2] / np.linalg.norm(self.lattice[2]))))
        self.beta = np.arccos(np.dot((self.lattice[0] / np.linalg.norm(self.lattice[0])),
                                     (self.lattice[2] / np.linalg.norm(self.lattice[2]))))
        self.gamma = np.arccos(np.dot((self.lattice[0] / np.linalg.norm(self.lattice[0])),
                                      (self.lattice[1] / np.linalg.norm(self.lattice[1]))))
    
    def read_cif(self, file_path):
        """CIF file reader
           
        It will read the cif file from file_path and give the required values back. Note this is
        a private method and generally is not intended to be called.
           
        Attributes:
            file_path: A String that indicates the file path of the cif file.
        
        Returns:
            returns a list that contains both lattices and fractional coordinates of each points for
            the crystal.
        """
        crystal = ase.io.read(file_path)

        lattice = crystal.get_cell()
        points = crystal.get_positions()
        fmotif = {}

        for i in range(len(points)):
            fmotif[i + 1] = points[i]

        return [lattice, fmotif]

           
    def cartesianize(self):
        """points cartisianizer
        convert points in crystal from fractional coordinates to cartesian coordinates. Note this is
        a private method and generally is not intended to be called.
        """
        cmotif = dict.fromkeys(self.vertices, [])
        for i in self.fmotif:
            frac = np.array(self.fmotif[i])
            carts = np.multiply(frac, self.lattice[0]) + np.multiply(frac, self.lattice[1]) + np.multiply(frac,
                                                                                                          self.lattice[
                                                                                                              2])
            cmotif[i] = [list(carts)]

        return [self.lattice, cmotif]

    
    def int_3d(self, l, translist=[(0, 0, 0)]):
        """ 3d layer creators
        creates a list of tuples which are all non-parallel integer vectors with positive first
        non-zero entries whose maximal distance is equal to the input layer l. Note this is a
        private method and generally is not intended to be called.
        """
        layerlist = [i for i in range(-l, l + 1)]

        layerlist.sort(reverse=True)

        allvecs = (list(item for item in it.combinations_with_replacement(layerlist, 3)))

        for i in allvecs:
            if translist == [(0, 0, 0)]:
                translist.append(i)
            else:
                for j in translist[1:]:
                    parcheck = np.cross(np.array(i), np.array(j))
                    if i not in translist and np.count_nonzero(parcheck) != 0:
                        translist.append(i)
                    else:
                        break

        for k in translist[1:]:
            nz = np.nonzero(k)
            fnz = nz[0][0]
            if k[fnz] < 0:
                flip = tuple([-1 * i for i in k])
                translist.remove(k)
                translist.append(flip)

        return translist

    def expand3d(self, s, trans=[(0, 0, 0)]):
        """
        Generate a list of reduced translation indices in 2d starting with a given input layer
        and ending at layer s. Note this is a private method and generally is not intended to
        be called.
        """
        tlist = self.int_3d(s, trans)

        c = self.cartesianize()
        c.append(tlist)

        # Calculate lattice translate vectors from indices and add new points to crystal

        for i in c[1]:
            origo = c[1][i][0]
            for j in tlist[1:]:
                j_vect = j[0] * np.array(self.lattice[0]) + j[1] * np.array(self.lattice[1]) + j[2] * np.array(
                    self.lattice[2])
                movepoint = list(np.add(j_vect, origo))
                c[1][i].append(movepoint)

        return c

    def motif_graph_3d(self, cutoff):
        """ Motif graph creators
        
        Make a motif graph up to a specific cutoff distance for the crystal. Now the cutoff is a fixed length value but for afterwards there might be methods to automatically generate an expected cutoff.

        Attributes:
            cutoff: A float that indicates cutoff distance for generating the motif graph
        
        Returns:
            returns a list that contains two values. The first is a list contains all vertices
            for the motif graph. The second is a dictionary contains whether each two vertices
            have edges between them. Keys of dictionary are points pairs and values of dictionary
            are list contains layer and their edge weight(null if there's no edge). For example:
            
            [[1, 2, 3],
            {(1, 1): [], (1, 2): [[(0, 0, 0), 0.6066300355241241]],
             (1, 3): [], (2, 2): [], (2, 3): [], (3, 3): []}]
        """

        distcheck = 0
        layer = 1
        tlist = [(0, 0, 0)]
        vertlist = self.vertices
        edgelist = {j: [] for j in it.combinations_with_replacement(self.vertices, 2)}
        originlist = [self.cartesianize()[1][i][0] for i in self.cartesianize()[1]]
        'Cartesian co-ordinates of points in the cell with the lattice at the origin'

        dset = [0]
        'Check distances within motif'
        first_el = [j for j in it.combinations_with_replacement(self.vertices, 2)]
        for j in first_el:
            d = np.linalg.norm(np.array(originlist[j[1] - 1]) - np.array(originlist[j[0] - 1]))
            dset.append(d)
            if 0 < d <= cutoff:
                edgelist[j].append([(0, 0, 0), d])

        distcheck = max(dset)
        if distcheck > cutoff:
            print('Still within motif')

            return [vertlist, edgelist]

        else:
            'Expand until first point where distances greater than scale appear'
            while distcheck < cutoff:
                startcheck = len(tlist)
                scaffold = self.expand3d(layer, trans=tlist)

                'Check distances in layers out from motif'
                for i in originlist:
                    for j in scaffold[1]:
                        edgekey = (originlist.index(i) + 1, j)
                        if edgekey in list(edgelist.keys()):
                            for k in scaffold[1][j][startcheck:]:
                                place = scaffold[1][j].index(k)
                                d = np.linalg.norm(np.array(k) - np.array(i))
                                dset.append(d)
                                if 0 < d <= cutoff:
                                    edgelist[edgekey].append([scaffold[2][place], d])

                distcheck = max(dset)
                tlist = scaffold[2]
                layer += 1

        motiflist = []
        for j in edgelist:
            if j[0] != j[1] and edgelist[j]:
                motiflist.append(edgelist[j][0][1])

        loopcount = 0
        looplist = []
        for j in edgelist:
            if j[0] == j[1] and edgelist[j]:
                looplist.append(edgelist[j][0][1])

        if motiflist:
            print('Motif distance achieved at ' + repr(max(motiflist)))

        if len(looplist) == len(vertlist):
            print('Loop distance acheived at ' + repr(min(looplist)))

        return [vertlist, edgelist]
    
    def create_graph(self,cutoff):
        """Graph creator

        It will create a graph(V,E) for the crystal with its motif result edges. Weights for edges
        are now remains null and it's under consideration to add weights in following versions.

        Attributes:
            cutoff: A float that indicates cutoff ditances for the motif result.

        Returns:
            A graph G with vertices and edges
        """

        vertices, edge_list = self.motif_graph_3d(cutoff)

        G = nx.Graph()
        G.add_nodes_from(vertices)

        for key in edge_list.keys():
            if edge_list[key]:
                for edge in edge_list[key]:
                    G.add_edge(key[0], key[1], translation_label=edge[0], distance=edge[1])

        return G
    
    def plot(self,cutoff=None,G=None):
        """Plot a graph G
        
        It will plot the graph for the motif result of the crystal in a relatively optimized form
        for the crystal. Also it can be used to plotting for outer graph G.
        
        Attributes:
            G: A graph G with edges and vertices
            cutoff: A float that indicates cutoff distance for the motif result
        
        Returns:
            A plotted graph with circular display
        """
        if G is not None:
            pos = nx.circular_layout(G)
            result_graph = nx.draw_networkx(G, pos, with_labels=True, font_weight='bold')
        elif cutoff is not None:
            G = self.create_graph(cutoff=cutoff)
            pos = nx.circular_layout(G)
            result_graph = nx.draw_networkx(G, pos, with_labels=True, font_weight='bold')
        else:
            print("At least one of the cutoff and G should be given")
        return result_graph

if __name__ == "__main__":
    main()
