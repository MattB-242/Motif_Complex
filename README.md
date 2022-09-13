This motif graph project is an application for generating motif results from crystals and to perform other related functions based on motif results.

The motif calculation part are from Matt Bright's work, which can calculate the motif result from crystal's lattice and fractional coordinates. The original project can be found on "https://github.com/MattB-242/Motif_Complex".


* ZhaoyuMotifGrap.motif_graph module contains following methods:
  * crystal_object(file_path = None, lattice = None, fmotif=None) to create a crystal object. To generate it, you have to either give a cif file path or give the lattice and fractional coordinates for each points in crystal.
    * crystal_object.motif_graph_3d(cutoff) to return motif result of the crystal object.
    * crystal_object.create_graph(cutoff) to return a graph G for the crystal object.
    * crystal_object.plot(cutoff=None,G=None) can return the plotted graph. you need to give it either a cutoff to calculate the graph itself and plot it, or you can give it a outer G to plot.

* ZhaoyuMotifGrap.calculation module contains following methods:
  * compute_spectrum(G,spectrum_method=spectrum.modularity_spectrum) returns both the 42 largest descending-order elements of the vector and unique ones, for the graph G.
  * process_folder(file_directory,cutoff=25,spectrum_method=spectrum.modularity_spectrum) returns both the 42 largest descending-order elements of the vector and unique ones, for each file in a given directory and name

* Also the example folder contains following files:
  * split file contains split(input_file_path,output_file_directory,keyword='END',outfile_name='Crystal') to split a given big file containing multiple crystals. Note that it only works when each elements ends with a general keyword.
  * process_cif can run a sample of calculation module
