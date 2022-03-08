import cython
import numpy as np
cimport numpy as numpy
from cpython.bytes cimport PyBytes_Size 
from libc.stdio cimport fopen, fclose, fread, fseek, FILE
from libc.string cimport strtok, strcpy
from libc.stdlib cimport atoi, atof, malloc, free
from oxDNA_analysis_tools.UTILS.RyeReader import Configuration

@cython.wraparound(False)
@cython.boundscheck(False)

def get_confs(list idxs, str traj_path, int start, int nconfs, int nbases):
    """
    A string!
    """
    # Number of configurations to read
    cdef int conf_count = len(idxs)
    if (start+nconfs >= conf_count):
        nconfs = conf_count - start

    # Configuration start/size markers within the chunk
    cdef int *sizes = <int *> malloc(nconfs * sizeof(int))
    cdef int *conf_starts = <int *> malloc(nconfs * sizeof(int))
    if not sizes or not conf_starts:
        raise MemoryError("Could not allocate memory for the configuration sizes and starts")
    cdef int chunk_size = 0
    for i in range(nconfs):
        sizes[i] = idxs[start+i].size
        chunk_size += sizes[i]
        conf_starts[i] = idxs[start+i].offset - idxs[start].offset


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
    free(sizes)
    free(conf_starts)

    return confs

cdef parse_conf(char *chunk, int start_byte, int size, int nbases):
    cdef int THREE = 3
    cdef int time
    
    #allocate some memory for our configuration
    cdef double *cbox = <double *> malloc(THREE * sizeof(double))
    cdef double *cenergy = <double *> malloc(THREE * sizeof(double))
    cdef double *cposes = <double *> malloc(nbases * THREE * sizeof(double))
    cdef double *ca1s = <double *> malloc(nbases * THREE * sizeof(double))
    cdef double *ca3s = <double *> malloc(nbases * THREE * sizeof(double))

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
        cbox[j] = atof(ptr)
        ptr = strtok(NULL, ' ')
    ptr = strtok(NULL, ' \n')
    
    cenergy[0] = atof(ptr)
    ptr = strtok(NULL, ' ')
    cenergy[1] = atof(ptr)
    ptr = strtok(NULL, ' \n')
    cenergy[2] = atof(ptr)

    # Parse the configuration itself
    for i in range(nbases):
        for j in range(THREE):
            ptr = strtok(NULL, ' ')
            cposes[i*THREE+j] = atof(ptr)
        for j in range(THREE):
            ptr = strtok(NULL, ' ')
            ca1s[i*THREE+j] = atof(ptr)
        for j in range(THREE):
            ptr = strtok(NULL, ' ')
            ca3s[i*THREE+j] = atof(ptr)
        ptr = strtok(NULL, '\n')

    # Convert the configuration information into numpy arrays and store in a Configuration
    box = np.asarray(<double[:3]>cbox)
    energy = np.asarray(<double[:3]>cenergy)
    poses = np.asarray(<double[:nbases*3]>cposes).reshape(nbases, THREE)
    a1s = np.asarray(<double[:nbases*3]>ca1s).reshape(nbases, THREE)
    a3s = np.asarray(<double[:nbases*3]>ca3s).reshape(nbases, THREE)
        
    cdef out  = Configuration(time, box, energy, poses, a1s, a3s)

    return out