
# General_MGraph_Tools
Creates a general  motif graph class, which can be n-periodic, as a list of vertices, a dictionary with edges as keys, each associated with a list of edge voltages and distances. and a lattice translation array.

Functions extracts/exampines a number of useful features:

*trans_set* will output all of the edge voltages in a graph - this is useful in itself but also serves as a utility function for testing whether a graph has achieved the bridging distance

*distance_list* will produce an ordered list of edge distances. Also useful in itself, but serves as a utility function by allowing the actual motif/bridging/loop distance to be extracted from a graph when a filtration on the scale has achieved this.

*molgraph* outputs a networkx simple graph object with vertices equal to the motif graph and a single edge where there is at least one edge in the motif graph. 

*mol_complement* outputs the complement of the molecular graph - that is, it gives the motif graph which has the shortest edge between any two vertices (where such an edge exists) removed.

# Motif_Graph_3d
he object Crystal_3d takes a list of lattice vectors (in the standard R3 basis) and a dictionary of integer labelled motif points. Motif co-ordinates are given **fractionally**. Lattice edge lengths and alpha, beta, gamma angle values are calculated directly from the object

The function *motif_graph_3d(self, s)* creates the motif graph up to the distance value *s* by expanding lattice translates. The output is a list of vertices, a dictionary with edges as keys with each edge associated with a list of edge voltages and their associated distances, sorted by edge distance, and the matrix of lattice vectors. 

It required the following utility functions:

    *cart_3d()* converts the fractional co-ordinates to cartesian ones. 
    *int_3d(l, m=0)* takes integer parameters a list of integer 3-tuples  which are all non-parallel integer vectors with positive first non-zero entries whose       maximal L_infty distance is equal to layer and whose minimal L_infty distance is equal to *m*. 

    *int_3d_outer(l)* takes a single integer parameter and outputs *int_3d(l,l)*, which gives just those 3-tuples containing the maximum value of *l*. Setting *l = 0* gives just the zero tuple. 

Note that in terms of viewing the graph as a voltage graph (Ross qv) only **plus directed** edges are added to the output. 

The functions *get_mdist_3d, get_bdist_3d* and *get_ldist_3d* output the graph at the motif, bridging and loop distance respectively (see note above on bridging distance - proof that this is indeed the minimal distance is TBD). 

The function *get_critical_distances_3d* outputs a list of the motif, bridge and loop distance in that order. 



