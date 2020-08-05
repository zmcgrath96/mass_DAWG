from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "../src/MassDawg.hpp":
    cdef cppclass MassDawg: 
        MassDawg() except +
        void show()
        void insert(vector[float], vector[float], string) except +
        vector[string] fuzzySearch(vector[float], int, int)
        vector[string] search(vector[float], int)
        void finish()