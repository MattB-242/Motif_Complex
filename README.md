# Motif_Complex
Calculating topological invariants from periodic crystal structures

Project Aim: 3D calculation of the Abstract Motif Complex from the Motif Graph

Current stage: 2D calculation of Motif Graph

The 2d Motif Graph Calculation has the following utility functions:

- MaxDist(l) takes a list of n-element lists as its input, calculates the pairwise distance between them as co-ordinates in R^n and outputs the maximum distance
- Ti_list (dim, e, s=0) takes a dimension (currently needs to be hardcoded at 2) an end layer. It generates a list of all REDUCED lattice translation vectors among those required to build concentric 2d layers of motif points from layer s (defaults to 0) to layer e. That is, it generates only vectors which are not scalar multiples of already existing vectors in the list.

The rest of the auxiliary functions are all methods on the 2D Crystal Object, the data of which is a pair of lattice translations and a list of motif points, all listed as co-ordinates in R^2:

- getpoints() extracts an unordered list of all currently existing motif points in the crystal to which maxdist(l) can be applied
- expand(e, s=1) adds translates of the motif points to the object in layers s to e (NB the object itself is modified - any subsequent methods called on the crystal will operate on the expanded version)

Finally, motif_graph_2d(s) adds expansion layers using the expand method until the new layer exceeds a threshold distance s. It then searches pairwise through all the distances between all the points in the expanded crystal. 

For each distance d < s between vertices i and j (with possibly i=j) it adds the relevant reduced translation vector to the list of labels whose key is the vertex i, and adds the same translation label plus the distance to the list of labels attached to the edge (i,j). 

The output is therefore a vertex dictionary consisting of entries of type i:[[a_1,a_2]]....] and an edge dictionary consisting of entries ij;[[a_1,a_2],d],...]

*NB There are a couple of kludgy fixes to avoid some duplication and some wierdness with index counting that probably have more elegant solutions...
