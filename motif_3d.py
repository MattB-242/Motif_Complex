import itertools as it
import numpy as np
import math

def int_3d(l, translist = [(0,0,0)]):
    
    layerlist = [i for i in range (-l, l+1)]
    
    layerlist.sort(reverse = True)
    
    allvecs = (list(item for item in it.combinations_with_replacement(layerlist,3)))

    for i in allvecs:
        if translist == [(0,0,0)]:
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
            flip = tuple([-1*i for i in k])
            translist.remove(k)
            translist.append(flip)

    return translist

class Crystal_3d:
    
    def __init__(self, lattice, fmotif):
        self.lattice = lattice
        self.fmotif = fmotif
        self.motif_ord = len(fmotif)
        self.vertices = list(fmotif.keys())
        self.lens = [np.linalg.norm(i) for i in lattice]
        self.alpha = np.arccos(np.dot((self.lattice[1]/np.linalg.norm(self.lattice[1])), (self.lattice[2]/np.linalg.norm(self.lattice[2]))))
        self.beta = np.arccos(np.dot((self.lattice[0]/np.linalg.norm(self.lattice[0])), (self.lattice[2]/np.linalg.norm(self.lattice[2]))))
        self.gamma = np.arccos(np.dot((self.lattice[0]/np.linalg.norm(self.lattice[0])), (self.lattice[1]/np.linalg.norm(self.lattice[1]))))

    def cartesianize(self):
        
        cmotif = dict.fromkeys(self.vertices,[])
        for i in self.fmotif:
            frac = np.array(self.fmotif[i])
            carts = np.multiply(frac,self.lattice[0]) + np.multiply(frac,self.lattice[1]) + np.multiply(frac,self.lattice[2])
            cmotif[i] = [list(carts)]
                
        return [self.lattice, cmotif]
    
    '--------------------------------------------------------------------'
    'EXPAND GRAPH UP TO A GIVEN LAYER'
    '--------------------------------------------------------------------'
    
    def expand3d(self, s, trans = [(0,0,0)]):
        
        'Generate a list of reduced translation indices in 2d starting with a given input layer and ending at layer s'
        tlist = int_3d(s, trans)
        
        c = self.cartesianize()
        c.append(tlist)

        'Calculate lattice translate vectors from indices and add new points to crystal'
        
        for i in c[1]:
            origo = c[1][i][0]
            for j in tlist[1:]:
                j_vect = j[0]*np.array(self.lattice[0]) + j[1]*np.array(self.lattice[1]) + j[2]*np.array(self.lattice[2])
                movepoint = list(np.add(j_vect,origo))
                c[1][i].append(movepoint)
        
        return c
    
    '-------------------------------------'
    'MAKE THE MOTIF GRAPH UP TO DISTANCE s'
    '-------------------------------------'
    def motif_graph_3d(self,s):
        
        distcheck = 0 
        layer = 1 
        tlist = [(0,0,0)]
        vertlist = self.vertices
        edgelist = {j:[] for j in it.combinations_with_replacement(self.vertices,2)}
        originlist = [self.cartesianize()[1][i][0] for i in self.cartesianize()[1]]
        'Cartesian co-ordinates of points in the cell with the lattice at the origin'
            
        dset = [0] 
        'Check distances within motif'
        first_el = [j for j in it.combinations_with_replacement(self.vertices,2)]
        for j in first_el:
            d = np.linalg.norm(np.array(originlist[j[1]-1]) - np.array(originlist[j[0]-1]))
            dset.append(d)
            if 0 < d <= s:
                edgelist[j].append([(0,0,0), d])
        
        
        distcheck = max(dset)
        if distcheck > s:
            print('Still within motif')
        
            return [vertlist, edgelist]
        
        else:
            'Expand until first point where distances greater than scale appear'
            while distcheck < s:
                startcheck = len(tlist)
                scaffold = self.expand3d(layer, trans = tlist)
            
            
                
            

            
                'Check distances in layers out from motif'
                for i in originlist:
                    for j in scaffold[1]:
                        edgekey = (originlist.index(i)+1, j)
                        if edgekey in list(edgelist.keys()):
                            for k in scaffold[1][j][startcheck:]:
                                place = scaffold[1][j].index(k)
                                d = np.linalg.norm(np.array(k) - np.array(i))
                                dset.append(d)
                                if 0 < d <= s:
                                    edgelist[edgekey].append([scaffold[2][place], d])
                                        
                distcheck = max(dset)
                tlist = scaffold[2]
                layer+=1
        
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
                        
        return [vertlist,edgelist]
                          
crystest = Crystal_3d([[1,0,0],[0.1,1,0],[0.3,0.2,1]], {1:[0,0,0], 2:[0.2, 0.3, 0.4], 3:[0.7,0.9,0.8]})

print(crystest.gamma)
print(crystest.motif_graph_3d(1))
