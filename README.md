# Motif_Graph
A topological invariants of periodic point clouds in R3

The function *int_3d(layer, list = (0,0,0)* creates a list of tuples which are all non-parallel integer vectors with positive first non-zero entries whose maximal distance is equal to the input layer l. If a pre-existing list of integer tuples is given, it will add only the vectors not in that tuple (that list defaults to just the zero vector).

The object Crystal_3d takes a list of lattice vectors (in the standard R3 basis) and a dictionary of integer labelled motif points. Motif co-ordinates are given **fractionally**. 

The function *cartesianize()* converts the fractional co-ordinates to cartesian ones. *expand3d (self, layer, list)* expands the graph in layers (again, out from a pre-existing list which is set at the zero vector)

Finally, *motif_graph_3d(self, s)* creates the motif graph up to the distance value s

