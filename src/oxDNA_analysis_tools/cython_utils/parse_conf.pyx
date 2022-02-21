import cython
import numpy as np
cimport numpy as numpy
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport strtok
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
    
    cdef char * text
    text = <char *> PyMem_Malloc(len(lines) * sizeof(char))
    if not text:
        raise MemoryError("Could not allocate memory for configuration")
    cdef char * sub = strtok(&text, '\n')
    print(sub)
    #cdef char * sub2 = strtok(sub, '=')
    #print(sub2)
    #sub2 = strtok(NULL, '=')
    #time = int(sub2)
    split = lines.split(b'\n')
    cdef int i = 0
    cdef int j = 0
    time = int(split[i].split(b'=')[1].strip().decode('UTF-8'))
    i += 1
    for j in range(THREE):
        box[j] = float(split[i].split(b'=')[1].strip().split(b' ')[j].decode('UTF-8'))
    i += 1
    for j in range(THREE):
        energy[j] = float(split[i].split(b'=')[1].strip().split(b' ')[j].decode('UTF-8'))
    i += 1

    cdef bytes line = split[i]
    while i < nbases+THREE:
        for j in range(THREE):
            poses[i-THREE, j] = float(line.split(b' ')[j].decode('UTF-8'))
            a1s[i-THREE, j] = float(line.split(b' ')[THREE+j].decode('UTF-8'))
            a3s[i-THREE, j] = float(line.split(b' ')[THREE+THREE+j].decode('UTF-8'))
            
        i += 1
        line = split[i]

    cdef out  = Configuration(time, box, energy, poses, a1s, a3s)

    return out