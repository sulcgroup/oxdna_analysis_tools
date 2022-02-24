import cython
import numpy as np
cimport numpy as numpy
from cpython.bytes cimport PyBytes_Size 
from libc.stdio cimport fopen, fclose, fread, fseek, FILE
from libc.string cimport strtok, strcpy, memcpy
from libc.stdlib cimport atoi, atof, malloc, free
from oxDNA_analysis_tools.UTILS.RyeReader import Configuration

@cython.wraparound(False)
@cython.boundscheck(False)

def get_confs(list idxs, str traj_path, int start, int nconfs, int nbases):
    # Number of configurations to read
    cdef int conf_count = len(idxs)
    if (start+nconfs >= conf_count):
        nconfs = conf_count - start

    # Configuration start/size markers within the chunk
    cdef numpy.ndarray[numpy.int32_t, ndim=1] sizes = np.zeros(nconfs, dtype=np.int32)
    cdef numpy.ndarray[numpy.int32_t, ndim=1] conf_starts = np.zeros(nconfs, dtype=np.int32)
    cdef int chunk_size = 0
    for i in range(start,start+nconfs):
        sizes[i] = idxs[i].size
        chunk_size += sizes[i]
        conf_starts[i] = idxs[i].offset - idxs[start].offset

    # Convert the path to something C can open
    cdef char *traj_path_c = <char *>malloc(len(traj_path)+1)
    strcpy(traj_path_c, traj_path.encode('utf-8'))
    traj_path_c[len(traj_path)] = b'\0'
    cdef FILE *traj_file = fopen(traj_path_c, "rb")
    if traj_file == NULL:
        print("Could not open trajectory file %s" % traj_path)
        return
    fseek(traj_file, idxs[start].offset, 1)

    # Read in the current chunk
    cdef const char *chunk = <char *>malloc(chunk_size)
    fread(chunk, chunk_size, 1, traj_file)

    # Parse the chunk into Configurations
    cdef list confs = [None]*nconfs
    for i in range(nconfs):
        c = parse_conf(chunk, conf_starts[i], sizes[i], nbases)
        confs[i] = c

    fclose(traj_file)
    free(chunk)
    free(traj_path_c)

    return confs

def parse_conf(char *chunk, int start_byte, int size, int nbases):
    cdef int THREE = 3
    cdef int time
    cdef numpy.ndarray[numpy.float32_t, ndim=1] box = np.zeros(THREE, dtype=np.float32, order='C')
    cdef numpy.ndarray[numpy.float32_t, ndim=1] energy = np.zeros(THREE, dtype=np.float32, order='C')
    cdef numpy.ndarray[numpy.float32_t, ndim=2] poses = np.zeros((nbases, THREE), dtype=np.float32, order='C')
    cdef numpy.ndarray[numpy.float32_t, ndim=2] a1s = np.zeros((nbases, THREE), dtype=np.float32, order='C')
    cdef numpy.ndarray[numpy.float32_t, ndim=2] a3s = np.zeros((nbases, THREE), dtype=np.float32, order='C')
    
    cdef int j = 0
    cdef int i = 0

    # Get a pointer to the start of the configuration
    cdef const char *ptr = chunk + start_byte

    # Get the time
    ptr = strtok(ptr, 't = ')
    time = atoi(ptr)

    # Get the box and energy
    # The energy can't be in a loop because of the format change between it and the conf lines.
    ptr = strtok(NULL, '= ')

    for j in range(THREE):
        box[j] = atof(ptr)
        ptr = strtok(NULL, ' ')
    ptr = strtok(NULL, ' \n')
    
    energy[0] = atof(ptr)
    ptr = strtok(NULL, ' ')
    energy[1] = atof(ptr)
    ptr = strtok(NULL, ' \n')
    energy[2] = atof(ptr)

    # Parse the configuration itself
    for i in range(nbases):
        for j in range(THREE):
            ptr = strtok(NULL, ' ')
            poses[i,j] = atof(ptr)
        for j in range(THREE):
            ptr = strtok(NULL, ' ')
            a1s[i,j] = atof(ptr)
        for j in range(THREE):
            ptr = strtok(NULL, ' ')
            a3s[i,j] = atof(ptr)
        ptr = strtok(NULL, '\n')
        
    cdef out  = Configuration(time, box, energy, poses, a1s, a3s)
    return out