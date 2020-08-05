# distutils: language = c++
# distutils: sources = ../src/MassDawg.cpp ../src/utils.cpp ../src/MassDawgNode.cpp

from libcpp.string cimport string 
from libcpp.vector cimport vector

from MassDawg cimport MassDawg

# Create a Cython extension type which holds a C++ instance
# as an attribute and create a bunch of forwarding methods
# Python extension type.
cdef class PyMassDawg:
    cdef MassDawg * m_dawg    # holds the c++ instance that is wrapped

    def __cinit__(self):
        self.m_dawg = new MassDawg()

    def __dealloc__(self):
        del self.m_dawg

    def show(self) -> None:
        '''
        Print the graph to the console as a tree
        '''
        self.m_dawg.show()

    def insert(self, singly_sequence: list, doubly_sequence: list, kmer: str) -> None:
        '''
        Insert a singly and doubly sequence into the graph associated with a kmer

        Inputs:
            singly_sequence:    (list) singly charged masses (floats) to add to the graph
            doubly_sequence:    (list) doubly charged masses (floats) to add to the graph
            kmer:               (str) the kmer associated with these masses
        Outputs:
            None
        '''
        cdef vector[float] singly_vec = singly_sequence
        cdef vector[float] doubly_vec = doubly_sequence

        cdef string input_kmer = str.encode(kmer)

        try:
            self.m_dawg.insert(singly_vec, doubly_vec, input_kmer)
        except:
            print("ERROR: Items must be sorted smallest to largest for insertion")

    def fuzzy_search(self, search_sequence: list, gap_allowance: int, ppm_tol: int) -> list:
        '''
        Search for a sequence in the graph allowing for up to gap_allowance missed masses in the search

        Inputs:
            search_sequence:    (list) the sequence of masses (floats) to search for
            gap_allowance:      (int) the number of gaps allowed in the search
            ppm_tol:            (int) the allowed difference (in parts per million) allowed 
                                      between an observed and theoretical mass to be called a match
        Outputs:
            (list) kmers (strings) found in the recursive search
        '''
        cdef vector[float] search_seq_vec = search_sequence

        cdef vector[string] results = self.m_dawg.fuzzySearch(search_seq_vec, gap_allowance, ppm_tol)

        return list(set(
            [result.decode() for result in results]
        ))

    def search(self, search_sequence: list, ppm_tol: int) -> list:
        '''
        Search for a sequence in the graph with no gaps allowed

        Inputs:
            search_sequence:    (list) the sequence of masses (floats) to search for
            ppm_tol:            (int) the allowed difference (in parts per million) allowed 
                                      between an observed and theoretical mass to be called a match
        Outputs:
            (list) kmers (strings) found in the search
        '''
        cdef vector[float] search_seq_vec = search_sequence
        
        cdef vector[string] results = self.m_dawg.search(search_seq_vec, ppm_tol)
        
        return list(set(
            [result.decode() for result in results]
        ))

    def finish(self):
        '''
        Final compression of any leftover nodes
        '''
        self.m_dawg.finish()