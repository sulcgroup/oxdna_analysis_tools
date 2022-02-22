import cython
import numpy as np
cimport numpy as numpy
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport strtok, strcpy
from libc.stdlib cimport atoi, atof, malloc, free
from oxDNA_analysis_tools.UTILS.RyeReader import Configuration

@cython.wraparound(False)
@cython.boundscheck(False)

def parse_conf(bytes lines, int nbases):
    cdef int THREE = 3
    cdef int time
    cdef numpy.ndarray[numpy.float32_t, ndim=1] box = np.zeros(THREE, dtype=np.float32, order='C')
    cdef numpy.ndarray[numpy.float32_t, ndim=1] energy = np.zeros(THREE, dtype=np.float32, order='C')
    cdef numpy.ndarray[numpy.float32_t, ndim=2] poses = np.zeros((nbases, THREE), dtype=np.float32, order='C')
    cdef numpy.ndarray[numpy.float32_t, ndim=2] a1s = np.zeros((nbases, THREE), dtype=np.float32, order='C')
    cdef numpy.ndarray[numpy.float32_t, ndim=2] a3s = np.zeros((nbases, THREE), dtype=np.float32, order='C')
    
    cdef int j = 0
    cdef int i = 0

    # Create a pointer to the start of the string containing the configuration
    cdef char * ctext
    ctext = <char *> malloc(len(lines)+1 * sizeof(char))
    if not ctext:
        raise MemoryError("Could not allocate memory for configuration")
    strcpy(ctext, lines)

    # Get the time
    cdef char * sub = strtok(ctext, '= ')
    time = atoi(sub)

    # Get the box and energy
    sub = strtok(NULL, '= ')

    for j in range(THREE):
        box[j] = atof(sub)
        sub = strtok(NULL, ' ')
    sub = strtok(NULL, ' \n')
    
    energy[0] = atof(sub)
    sub = strtok(NULL, ' ')
    energy[1] = atof(sub)
    sub = strtok(NULL, ' \n')
    energy[2] = atof(sub)

    for i in range(nbases):
        for j in range(THREE):
            sub = strtok(NULL, ' ')
            poses[i,j] = atof(sub)
        for j in range(THREE):
            sub = strtok(NULL, ' ')
            a1s[i,j] = atof(sub)
        for j in range(THREE):
            sub = strtok(NULL, ' ')
            a3s[i,j] = atof(sub)
        strtok(NULL, '\n')
        
    cdef out  = Configuration(time, box, energy, poses, a1s, a3s)
    free(ctext)
    return out