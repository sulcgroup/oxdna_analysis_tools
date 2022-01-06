"""
This parallelizer attaches multiple ErikReaders to a single trajectory.
This is less memory-intensive than multi-file, but is a significant performance decrease on certain systems.
"""

import pathos.multiprocessing as pp
from os import getenv
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2
import numpy as np

#actually unused these days, but just in case...
def get_n_cpu():
    """
    Gets the number of CPUs available from either the slum environment or pathos' check of available resources.

    Returns:
        available_cpus (int): Either the number of slurm tasks assigned to the process or the count of CPUs on the machine.
    """
    try:
        available_cpus = int(getenv('SLURM_NTASKS'))
    except:
        available_cpus = int(pp.cpu_count()/2)

    return available_cpus

#partitions the trajectory file to the number of workers defined by n_cpus
#each worker runs the given function on its section of the trajectory
def fire_multiprocess(traj_file, top_file, function, num_confs, n_cpus, *args, **kwargs):
    """
    Splits a trajectory file into temporary files and attaches a reader to each file.

    Parameters:
        traj_file (str): Name of the trajectory file to split.  
        top_file (str): Name of the topology file associated with the trajectory. 
        num_confs (int): The number of configurations in the trajectory.  
        n_cpus (int): The number of chunks to split the trajectory into.  
        conf_per_processor (int): The number of configurations per chunk (equivalent to floor(num_confs/n_cpus))  

    Returns:
        readers (list of LorenzoReader2s): A list of readers with each one on a unique chunk of the file.
    """
    confs_per_processor = int(np.floor(num_confs/n_cpus))

    reader_pool = []
    processor_pool = pp.Pool(n_cpus)

    #for calculations on symmetric matricies (eRMSD)
    #can't just hand each line to the parallelizer
    if ("matrix", True) in kwargs.items():
        total_calculations = sum([(num_confs-i) for i in range(1, num_confs)])
        calcs_per_cpus = total_calculations/n_cpus
        split_ends = []
        i = 0
        while i < num_confs:
            e = 0
            calcs = 0
            while calcs < calcs_per_cpus:
                calcs += num_confs-i
                e += 1
                i += 1
                if i >=num_confs: break
            split_ends.append(e)
    
    #define sizes of trajectory chunks
    else:
        split_ends = [confs_per_processor for _ in range(n_cpus)]
        split_ends[-1] += num_confs%n_cpus #last chunk gets all the leftovers
    
    #now figure out which configuration each chunk starts on
    split_starts = [0]
    for i in range(n_cpus):
        reader_pool.append(LorenzoReader2(traj_file, top_file))
        #rint(split_starts[i-1], split_ends[i-1])
        if i!= 0:
            split_starts.append(split_starts[i-1]+split_ends[i-1])

    #staple everything together, send it out to the workers, and collect the results as a list
    results = []
    lst = [(r, *args, num_confs, s, e) for r, s, e in zip(reader_pool, split_starts, split_ends)]
    results = processor_pool.starmap_async(function, lst).get()
    processor_pool.close()

    return(results)
