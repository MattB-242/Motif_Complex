import itertools as it
import numpy as np
import math
import operator

'-----------------------------------------------------------'
'Calculate the maximum pairwise distance of a list of points'
'Used to test whether distance threshold reached in expansion'
'-----------------------------------------------------------'
def maxdist(l):

    distlist = []
    if len(l) == 1:
        return 0

    else:
        for dist_pairs in it.combinations(l,2):
            vectors = [np.array(dist_pairs[0]), np.array(dist_pairs[1])]
            distlist.append(np.linalg.norm(vectors[1]-vectors[0]))

    return max(distlist)

'---------------------------------------------------------------------------'
'Check the number of edges between each vertex in a Motif Graph'
'---------------------------------------------------------------------------'
    
def edgecount(motif_graph):
        
    edge_count_dict = {}
        
    for (i,j) in list(motif_graph[1].keys()):
        edge_count_dict[(i,j)] = len(motif_graph[1][(i,j)])
        
    
    return edge_count_dict

'-----------------------------------------------------------'
'Creates a standardized list of reduced translation indices'
'for concentric layers in a dim-dimensional crystal'
'starting at layer s and ending at layer e (s defaults to 1)'
'-----------------------------------------------------------'

def ti_list(dim, e, s=1):


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

'----------------------------------------------------'
'Generates 2d crystal as list of lattice translations'
'and dictionary of labelled motif points'
'----------------------------------------------------'
class crystal_2d:

    import numpy as np

    def __init__(self, lattice, motif):

        self.lattice = lattice
        self.motif = motif
        self.motif_ord = len(motif)
        self.vertices = list(motif.keys())
        self.unit_cell_volume = np.linalg.norm(np.cross(lattice[0],lattice[1]))
        self.ab_angle = np.arccos(np.dot(lattice[0],lattice[1]))

    '---------------------------------------------------------------------------'
    'Outputs a list of just all point positions in a crystal'
    'Used to test for distance maximum and count translations'
    '---------------------------------------------------------------------------'
    def getpoints(self):
        allpoints = []
        for i in self.vertices:
                allpoints.extend(self.motif[i])

        return allpoints

    '---------------------------------------------------------------------------'
    'Adds all points with distinct REDUCED translation indices in layers s to e'
    '---------------------------------------------------------------------------'
    def expand(self, e, s=1):

        'Generate a list of reduced translation indices in 2d starting at layer s and ending at layer e'
        translist = ti_list(2,e,s)

        'Calculate lattice translate vectors from indices and add new points to crystal'
        for i in self.vertices:
            for j in translist:
                j_vect = [sum(z) for z in zip([(j[0]*i) for i in self.lattice[0]], [j[1]*i for i in self.lattice[1]])]
                move_point = [sum(x) for x in zip(self.motif[i][0], j_vect)]
                self.motif[i].append(move_point)


        return translist

    '---------------------------------------------------------------------------'
    'Creates motif graph MG(C,s) for a particular 2d crystal and distance s'
    '---------------------------------------------------------------------------'

    def motif_graph_2d(self, s):

        'Initialise vertex and edge data'
        vert_list = {i:[[0,0]] for i in self.motif}
        edge_list = {j:[] for j in it.combinations_with_replacement(self.vertices,2)}
        print(edge_list)

        vertrans = []

        'Expand graph in layers up to the layer beyond the target distance'
        i = 1
        while maxdist(self.getpoints()) < s:
            vertrans.extend(self.expand(i, i-1))
            i+=1

        'Check in-cell motif distances'
        for i in self.vertices:
            for j in self.vertices:
                d = np.linalg.norm(np.array(self.motif[j][0]) - np.array(self.motif[i][0]))
                print(d)
                if d < s and (i,j) in edge_list.keys():
                        edge_list[(i,j)].append([[0,0],d])

        'Check all other motif distances'
        for i in self.vertices:
            for j in self.vertices:
                mot_point = self.motif[i][0]
                mot_test = self.motif[j]
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
    
        '--------------------------------------------------------------------------------------'
        'Run a filtration on Motif Graphs to a given distance s, and list values of s at which'
        'new edges form between each vertex pair i,j, as a dictionary keyed by (i,j)'
        '--------------------------------------------------------------------------------------'
    
        def MG_2d_Filter(self,s_max,increment):
        
            edge_filter_list = {j:[] for j in it.combinations_with_replacement(self.vertices,2)}
            s = 0
            mgraph = self.motif_graph_2d(s)
        
            while s <= s_max:
                try:
                    edge_1 = edgecount(mgraph)
                    s += increment
                    mgraph = self.motif_graph_2d(s)
                    edge_2 = edgecount(mgraph)
                
            
                    for (i,j) in list(edge_filter_list.keys()):
                        if edge_1[(i,j)] != edge_2[(i,j)]:
                            print('New edge formed between vertices '+(i,j)+' at distance '+s +'.')
                            edge_filter_list[(i,j)].append(s)
            
                except:
                    print('Problem at distance {:03f}'.format(s))
                    print('The graph here is')
                    print (mgraph)
                    raise
            
            
        return edge_filter_list
