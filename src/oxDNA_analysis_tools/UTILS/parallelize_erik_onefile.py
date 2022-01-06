"""
This parallelizer attaches multiple ErikReaders to a single trajectory.
This is less memory-intensive than multi-file, but is a significant performance decrease on certain systems.
"""

import pathos.multiprocessing as pp
from os import getenv
import numpy as np
from oxDNA_analysis_tools.UTILS.readers import ErikReader

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
def fire_multiprocess(traj_file, function, num_confs, n_cpus, *args, **kwargs):
    """
    Distributes a function over a given number of processes

    Parameters:
        traj_file (str): The name of the trajectory file to analyze.
        function (function): The analysis function to be parallelized.
        num_confs (int): The number of configurations in the trajectory.
        n_cpus (int): The number of processes to launch.
        *args: The arguments for the provided function.
        **kwargs: Used to pass additional parameters to the parallelizer (currently only used by eRMSD.py to correctly allocate tasks to processors)

    Returns:
        results (list): The results from each individual processor's run.

    Note: The manner in which to concatenate the results is function-specific so should be handled in the calling module.
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
        reader_pool.append(ErikReader(traj_file))
        #rint(split_starts[i-1], split_ends[i-1])
        if i!= 0:
            split_starts.append(split_starts[i-1]+split_ends[i-1])

    #staple everything together, send it out to the workers, and collect the results as a list
    results = []
    lst = [(r, *args, num_confs, s, e) for r, s, e in zip(reader_pool, split_starts, split_ends)]
    results = processor_pool.starmap_async(function, lst).get()
    processor_pool.close()

    return(results)
