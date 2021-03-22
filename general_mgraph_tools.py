'Create the 3D motif Graph as a Class'
class Motif_Graph:
    
    def __init__ (self, graph):
        
        self.dimension = len(graph[2])
        self.vertices = graph[0]
        self.edges = graph[1]
        self.array = graph[2]
        self.loops = [(key, val) for key, val in self.edges.items() if key[0] == key[1]]
        self.notloops = [(key, val) for key, val in self.edges.items() if key[0] != key[1]]
        self.vertcount = len(graph[0])
        self.notloopcount = sum([len(i[1]) for i in self.notloops])
        self.loopcount = sum([len(i[1]) for i in self.loops])
        self.edgecount = self.notloopcount + self.loopcount
        
    'Extract all distinct translation labels in an input graph'
    def trans_set(self):
        
        translist = []
        for i in self.edges:
            if self.edges[i]:
                for j in self.edges[i]:
                    if np.count_nonzero(j[0]) != 0:
                        translist.append(j[0])
            
        return set(translist)
                
    
    'Get ordered list of distances out of an input graph'
    def distance_list(self):
        
        dlist = []
        
        for i in self.edges:
            if self.edges[i]:
                for j in self.edges[i]:
                    dlist.append(j[1])
                    
        dlist.sort()
        
        return dlist
    
    'Get Molecular Graph'
    
    def molgraph(self):
        mol_graph = nx.Graph()
        mol_graph.add_nodes_from(self.vertices)
        
        for i in self.edges:
            if i[0] != i[1] and self.edges[i]:
                mol_graph.add_edge(*i)
        
        return mol_graph
    
    'Get Molecular Graph Edge Complement'
    def mol_complement(self):
        molc_edges = {j:[] for j in self.edges.keys()}
        
        for i in self.edges:
            if self.edges[i] and len(self.edges[i]) > 1:
                for k in self.edges[i][1:]:
                    molc_edges[i].append(k)
        
        return [self.vertices, molc_edges, self.array]
                
    'Check that motif distance has been reached'
    def check_mdist(self):
        
        for i in self.notloops:
            if not i[1]:
                return False
                break
        
            else:
                return True
    
    'Check that loop distance has been reached'
    def check_ldist(self):
        
        for i in self.loops:
            if not i[1]:
                return False
                break
        
            else:
                return True
    
    def check_bdist(self):
        
        cgraph = Motif_Graph(self.mol_complement())
        
        if nx.is_connected(self.molgraph()) and len(cgraph.trans_set()) >= self.dimension:
            return True

        else:
            return False
