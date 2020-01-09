import pathos.multiprocessing as pp
from os import getenv
from UTILS.readers import LorenzoReader2
import numpy as np
from tempfile import NamedTemporaryFile
from UTILS.readers import blocks

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

def split_trajectory(traj_file, top_file, num_confs, n_cpus, confs_per_processor):
    """
        Splits a trajectory file into 
    """
    n_files = 0
    readers = []
    files = []
    rem = num_confs % n_cpus
    
    with open(traj_file, "rb") as f:
        it = blocks(f)
        chunk = next(it)
        last_conf_byte = 0
        while n_files < n_cpus:
            out = NamedTemporaryFile(mode='w+b')
            conf_count = 0
            if n_files < rem:
                a = 1
            else:
                a = 0
            while conf_count < confs_per_processor+a:
                next_conf_byte = chunk.find(b"t", last_conf_byte+1)
                if next_conf_byte == -1:
                    out.write(chunk[last_conf_byte:])
                    try:
                        chunk = next(it)
                    except Exception as e:
                        break
                    last_conf_byte = 0
                else:  
                    out.write(chunk[last_conf_byte:next_conf_byte])     
                    conf_count += 1
                    last_conf_byte = next_conf_byte 
            readers.append(LorenzoReader2(out.name, top_file))
            files.append(out)
            n_files += 1
        
    return(readers, files)
            

#partitions the trajectory file to the number of workers defined by n_cpus
#each worker runs the given function on its section of the trajectory
def fire_multiprocess(traj_file, top_file, function, num_confs, n_cpus, *args, **kwargs):
    confs_per_processor = int(np.floor(num_confs/n_cpus))

    reader_pool = []
    processor_pool = pp.Pool(n_cpus)

    reader_pool, tmpfiles = split_trajectory(traj_file, top_file, num_confs, n_cpus, confs_per_processor)
    split_starts = [0 for  r in reader_pool]
    split_ends = [confs_per_processor for r in reader_pool]
    rem = num_confs % n_cpus
    for i in range(rem):
        split_ends[i] += 1

    #staple everything together, send it out to the workers, and collect the results as a list
    results = []
    lst = [(r, *args, num_confs, s, e) for r, s, e in zip(reader_pool, split_starts, split_ends)]
    results = processor_pool.starmap_async(function, lst).get()
    processor_pool.close()

    return(results)
