"""
This parallelizer splits a trajectory into multiple tempfiles and attaches an ErikReader to each one.
This is more memory-intensive than one-file, but is a significant performance increase on certain systems.
"""

import pathos.multiprocessing as pp
from os import getenv, remove
import numpy as np
from tempfile import NamedTemporaryFile
from oxDNA_analysis_tools.UTILS.readers import blocks, ErikReader

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

def split_trajectory(traj_file, num_confs, n_cpus, confs_per_processor):
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
    n_files = 0
    readers = []
    files = []
    rem = num_confs % n_cpus
    
    with open(traj_file, "rb") as f:
        it = blocks(f)
        chunk = next(it) # iterator producing 1 MB chunks of the trajectory
        last_conf_byte = 0

        #create a number of temporary file equal to the number of CPUs
        while n_files < n_cpus:
            out = NamedTemporaryFile(mode='w+b', delete=False)
            conf_count = 0
            
            #If there is a remainder after dividing the number of configurations by the number of CPUs
            #Add one extra configuration to the first rem files
            if n_files < rem:
                a = 1
            else:
                a = 0

            #Find successive configuration start points and write them out to the tempfiles
            while conf_count < confs_per_processor+a:
                next_conf_byte = chunk.find(b"t", last_conf_byte+1)
                if next_conf_byte == -1:
                    out.write(chunk[last_conf_byte:])
                    try:
                        chunk = next(it)
                    except: #next() throws an error if there isn't another chunk
                        break
                    last_conf_byte = 0
                else:  
                    out.write(chunk[last_conf_byte:next_conf_byte])     
                    conf_count += 1
                    last_conf_byte = next_conf_byte 

            #create a reader from the newly created trajectory chunk
            readers.append(ErikReader(out.name))
            files.append(out)
            n_files += 1
        
    return(readers, files)
            

#partitions the trajectory file to the number of workers defined by n_cpus
#each worker runs the given function on its section of the trajectory
def fire_multiprocess(traj_file, function, num_confs, n_cpus, *args):
    """
    Distributes a function over a given number of processes

    Parameters:
        traj_file (str): The name of the trajectory file to analyze.
        function (function): The analysis function to be parallelized.
        num_confs (int): The number of configurations in the trajectory.
        n_cpus (int): The number of processes to launch.
        *args: The arguments for the provided function.

    Returns:
        results (list): The results from each individual processor's run.

    Note: The manner in which to concatenate the results is function-specific so should be handled in the calling module.
    """

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["pathos"])
    confs_per_processor = int(np.floor(num_confs/n_cpus))

    reader_pool = []
    processor_pool = pp.Pool(n_cpus)

    #split_starts and split_ends are around for backwards compatability with the old parallelize algorithm
    reader_pool, tmpfiles = split_trajectory(traj_file, num_confs, n_cpus, confs_per_processor)
    split_starts = [0 for  r in reader_pool]
    split_ends = [confs_per_processor for r in reader_pool]
    rem = num_confs % n_cpus
    for i in range(rem):
        split_ends[i] += 1

    #Staple everything together, send it out to the workers, and collect the results as a list
    #Functions passed to this parallelizer must have the argument order defined by the lst variable (reader, <unique args>, number of configurations total, starting conf id, number of confs for this processor)
    #This args unpacking method was added in Python 3.6, so if you have an older version of Python that's why this isn't working
    results = []
    lst = [(r, *args, num_confs, s, e) for r, s, e in zip(reader_pool, split_starts, split_ends)]
    
    #starmap allows you to have arguments that themselves are iterables
    #async because we don't actually care what order stuff finishes in.
    results = processor_pool.starmap_async(function, lst).get()
    processor_pool.close()
    for f in tmpfiles:
        f.close()
        remove(f.name)

    return(results)
