import itertools as it
import numpy as np
import math

tolerance = 1e-14


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

def int_3d(l, m=0):
    raw_result = set((0,0,0))
    
    if m==0:
        result = [(0,0,0)]
    else:
        result = []
        
    for i in range(3): #the index the first non-zero term will appear
        a_range,b_range,c_range = [(0,1) for _ in range(i)] +[(1,l+1)] + [(-l,l+1) for _ in range(2-i)]
        raw_result.update([(a,b,c) for a in range(*a_range) for b in range(*b_range) for c in range (*c_range) if a==1 or not a==b==c])

    for j in raw_result:
        if np.gcd.reduce(j) == 1 and np.abs(max(j)) >= m:
            result.append(j)
            result.append(tuple([-1*i for i in j]))
    
    result.sort(key = lambda j : max(np.abs(j)))
    
    return result

'---------------------------------------------'
'Return outer translation layer only'
'---------------------------------------------'
def int_3d_outer(l):
    return int_3d(l,m=l)

'----------------------------------------------------'
'Generates 2d crystal as list of lattice translations and dictionary of labelled motif points.'
'Format of lattice entry is an array'
'Format of a motif entry is a dictionary of *integer* vertex labels with co-ordinates expressed as fractions '
'of the lattice basis.'
'----------------------------------------------------'
class Crystal_3d:

    import numpy as np

    def __init__(self, lattice, motif):

        self.lattice = lattice
        self.lattice_array = np.transpose(np.array(lattice))
        self.fmotif = motif
        self.fmotif_ord = len(motif)
        self.vertices = list(motif.keys())
        self.unit_cell_volume = np.abs(np.linalg.det(self.lattice_array))
        self.lens = [np.linalg.norm(i) for i in lattice]
        self.alpha = np.arccos(np.dot(lattice[1]/self.lens[1],lattice[2]/self.lens[2]))
        self.beta = np.arccos(np.dot(lattice[0]/self.lens[0],lattice[2]/self.lens[2]))
        self.gamma = np.arccos(np.dot(lattice[0]/self.lens[0],lattice[1]/self.lens[1]))

    '---------------------------------------------------------------------------'
    'Returns a version of the crystal where all co-ordinates are cartesian'
    '---------------------------------------------------------------------------'    
    def cart_3d(self):
        
        cverts = {i:[np.matmul(self.lattice_array, self.fmotif[i])] for i in self.fmotif}
        
        return(self.lattice_array, cverts)
        

    '---------------------------------------------------------------------------'
    'Creates motif graph MG(C,s) for a particular 2d crystal and distance s'
    '---------------------------------------------------------------------------'

    def motif_graph_3d(self, s):
        
        arr = self.cart_3d()[0]
        pts = self.cart_3d()[1]
        
        'Initialise vertex and edge data'
        vert_list = [i for i in self.fmotif.keys()]
        edge_list = {j:[] for j in it.product(self.vertices,repeat=2)}
        dedge_list = {j:[] for j in it.product(self.vertices,repeat=2)}

        layer = 0
        distcheck = False
            
        while not distcheck:
            distlist = []
            tlist = int_3d_outer(layer)
            for i in tlist:
                neg_i = tuple([-1*k for k in i])
                for j in edge_list:
                    rev_j = (j[1],j[0])
                    if i not in edge_list[j] and neg_i not in edge_list[rev_j]:
                        p_one = np.array(pts[j[0]])
                        p_two = np.array(pts[j[1]] + np.matmul(arr, i))
                        dist = np.linalg.norm(np.subtract(p_two, p_one))
                        distlist.append(dist)
                        
                        'Add edge and distance label'
                        if 0 < dist- tolerance and dist+tolerance < s:
                            edge_list[j].append(i)
                            dedge_list[j].append([i,dist])
            
            
            if min(distlist) > s:
                distcheck = True
            else:
                layer+=1
                
        'sort all edges by distance'
        for j in dedge_list:
            dedge_list[j].sort(key = lambda i: i[1])
            
        return [vert_list, dedge_list, arr]
