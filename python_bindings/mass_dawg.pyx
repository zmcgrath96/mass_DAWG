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
        self.m_dawg.show()

    def insert(self, singly_sequence: list, doubly_sequence: list, kmer: str) -> None:
        cdef vector[float] singly_vec = singly_sequence
        cdef vector[float] doubly_vec = doubly_sequence

        cdef string input_kmer = str.encode(kmer)

        try:
            self.m_dawg.insert(singly_vec, doubly_vec, input_kmer)
        except:
            print("ERROR: Items must be sorted smallest to largest for insertion")

    def fuzzy_search(self, search_sequence: list, gap_allowance: int, ppm_tol: int) -> list:
        cdef vector[float] search_seq_vec = search_sequence

        cdef vector[string] results = self.m_dawg.fuzzySearch(search_seq_vec, gap_allowance, ppm_tol)

        return list(set(
            [result.decode() for result in results]
        ))

    def finish(self):
        self.m_dawg.finish()