"""
This project is an application of motif calculation and other related operations, like graph
generation and spectrum calculation.
_author_ = "Zhaoyu Han"
_credits_ = "Matt Bright,Zhaoyu Han"
"""
import ase.io
import numpy as np
from ase.visualize import view
import ase.lattice.orthorhombic
import networkx as nx
from networkx.linalg import spectrum
import matplotlib.pyplot as plt
import pickle as pk
import ZhaoyuMotifGraph.motif_graph
import os
#def quick_sort(li, start, end):
#    """
#    Takes in a vector and sorts in place , np.sort() ?
#    """
#    # start=end ,assume there's only on value
#    # start>end ,assume there's no value on the right
#    if start >= end:
#        return
#    # define two pointer to point to 0 and the end
#    left = start
#    right = end
#    # regard 0 as mid
#    mid = li[left]
#    while left < right:
#        # let right pointer move to the left, to find value smaller than mid, and put to the left pointer
#        while left < right and li[right] >= mid:
#            right -= 1
#        li[left] = li[right]
#        # let left pointer move to the right, to find value bigger than mid, and put to the right pointer
#        while left < right and li[left] < mid:
#            left += 1
#        li[right] = li[left]
#    # when wile ends, put mid to the middle
#    li[left] = mid
#    # iterate data on the left
#    quick_sort(li, start, left-1)
#    # iterate data on the right
#    quick_sort(li, left+1, end)

def compute_spectrum(G, spectrum_method=spectrum.modularity_spectrum):
    """Computer the spectrum for one single motif graph
    This method calculate the spectrum result for one single motif graph G.
    
    Attributes:
        G: a networkx graph that contains vertices and edges.
        spectrum_method: the spectrum type for this calculation. The default one is modularity. you can change it into other networkx spectrum methods.
    Returns:
        a list that contains two vectors that are array type. The first is the calculated spectrum
        vector which contains the largest 42 elements from the original one with descent order.
        The second is containing the 42 largest unique elements from the original vector.
    """
    try:
        vector = np.nan_to_num(spectrum_method(G))
    except:
        print('Errors happens within spectrum')
    else:
        sorted_vector = np.sort(vector,kind='quick_sort')
        required_vector = sorted_vector[-42:len(sorted_vector)]
        real_required_vector = list(reversed(required_vector))


        required_vector_unique = np.array(list(set(np.round(np.real(vector[-42:len(vector)]), 4))))
        if len(required_vector_unique) != len(required_vector):
            sorted_vector_unique = np.sort(required_vector_unique, kind='quick_sort')
            zeros = [0 for list_length in range(42 - len(sorted_vector_unique))]
            zeros.extend(sorted_vector_unique.tolist())
            vector_unique = list(reversed(zeros))
        else:
            vector_unique = real_required_vector

        return [real_required_vector, vector_unique]


def process_folder(file_directory,cutoff=25,spectrum_method=spectrum.modularity_spectrum):
    """Spectrum calculator for more than one file
    
    Given a folder path we will compute the vector representations of all cif files
    in the folder, and gives sorted vectors with 42-biggest elements for each file.
    
    Attributes:
        file_diectory: a String that indicates the directory where all cif files are stored.
        cutoff: a float that indicates the cutoff distance of the motif cauculation.
        spectrum_method: the spectrum type for this calculation. The default one is modularity. you can change it into other networkx spectrum methods.
    Returns:
        A list than contains three elements. error_vector_file is a list that contains the name
        of files that cannot be processed by the modular method. vectors is a array that contains vectors
        that contains 42 largest elements from their original vector for each elements, vectors_unique is an array
        contains vectors that contains 42 largest unique elements from the original vector for each
        element.
    """
    
    for index, file in enumerate(os.listdir(file_directory)):
        if file[-4:] == ".cif":
            file_crystal = motif_graph.Crystal_3d(file_path=file)
            G = file_crystal.create_graph(cutoff=cutoff)

            error_vector_file = []
            vectors = []
            vectors_unique = []
            
            try:
                # Can output infinity and null values...
                vector = np.nan_to_num(spectrum_method(G))
            except:
                error_vector_file.append(str(file))

                vectors.append(vectors[index - 1])
                print('this one is wrong! ', index, file)
            else:
                sorted_vector = np.sort(vector,kind='quick_sort')
                required_vector = sorted_vector[-42:len(sorted_vector)]
                real_required_vector = list(reversed(required_vector))
                vectors.append(real_required_vector)

                required_vector_unique = np.array(list(set(np.round(np.real(vector[-42:len(vector)]), 4))))
                if len(required_vector_unique) != len(required_vector):
                    sorted_vector_unique = np.sort(required_vector_unique, kind='quick_sort')
                    zeros = [0 for list_length in range(42 - len(sorted_vector_unique))]
                    zeros.extend(sorted_vector_unique,tolist())
                    vectors_unique.append(list(reversed(zeros)))
                else:
                    vectors_unique.append(real_required_vector)

                print(index + ' has been processed')
    
    vectors_array = np.array(vectors)
    vectors_unique_array = np.array(vectors_unique)
    
    return [error_vector_index,vectors_array,vectors_unique_array]

def save_result(vectors,file_path):
    """save the result
    
    Methods to save to result in pk binary format, which makes it faster to proceed.
    
    Attributes:
        vectors: the vecotr you need to save
        file_path: the file_path that you want to save to
    """
    pk.dump(vectors, open(file_path, "wb"))


