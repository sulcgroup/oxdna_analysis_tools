from oxDNA_analysis_tools.UTILS.micha_reader import MichaReader
from oxDNA_analysis_tools.UTILS.base_array import base_array
try:
    from Bio.SVDSuperimposer import SVDSuperimposer
except:
    from bio.SVDSuperimposer import SVDSuperimposer
from oxDNA_analysis_tools.UTILS.micha_reader import MichaReader
import numpy as np
from os import getenv, path
from oxDNA_analysis_tools.UTILS.micha_reader import partition, flatten, handle_confs,parse_conf
from multiprocessing import cpu_count, Pool
from functools import partial


import time
start_time = time.time()

#figure out how many chunks we are working on 

def get_n_cpu():
    try:
        available_cpus = int(getenv('SLURM_NTASKS'))
    except:
        available_cpus = cpu_count()
    return available_cpus


def align(centered_ref_coords, cms_ref_cords ,coords):
    # center on centroid
    av1, reference_coords = cms_ref_cords, centered_ref_coords
    av2 = np.mean(coords[:, 0:3], axis=0)
    coords[:, 0:3]= coords[:, 0:3] - av2
    # correlation matrix
    a = np.dot(np.transpose(coords[:, 0:3]), reference_coords)
    u, _, vt = np.linalg.svd(a)
    rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    # check if we have found a reflection
    if np.linalg.det(rot) < 0:
        vt[2] = -vt[2]
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    tran = av1 - np.dot(av2, rot)
    return  (np.dot(coords[:, 0:3], rot) + tran,
             np.dot(coords[:, 3:6], rot),
             np.dot(coords[:, 6:9], rot))


def compute(top, traj,centered_ref_coords, cms_ref_cords, ntopart, ptr):
    reader = MichaReader(top,traj)
    confs = reader._get_confs(ptr*ntopart,ntopart)

    # convert to numpy repr
    aligned_coords = reader._parse_confs(confs)
    aligned_coords = [align(centered_ref_coords,cms_ref_cords, c) for c in aligned_coords]
    sub_mean = np.sum(aligned_coords, axis=0)
    return sub_mean


if __name__ == "__main__":
    top, traj = "/home/erik/Simulations/hinge_all/hinge1/hinge_correct_seq.top", "/home/erik/Simulations/hinge_all/hinge1/aligned.dat"

    reader = MichaReader(top,traj)
    ref_conf = reader.read(0) # say it's always the 1st one
    ref_conf.inbox()

    ncpus = get_n_cpu()
    ntopart = 10

    #prepare to fire multiple processes
    pool = Pool(ncpus)

    n_confs = len(reader.idxs)
    n_chunks = int(n_confs / ntopart +
                         (1 if n_confs % ntopart else 0))

    print(n_chunks)

    reference_coords = ref_conf.positions
    av1 = np.mean(reference_coords, axis=0)
    reference_coords = reference_coords - av1

    compute_func = partial(compute, top, traj,reference_coords, av1, ntopart)

    results = [pool.apply_async(compute_func,(i,))
                            for i in range(n_chunks)]

    #print(
    #    results[0].get()
    #)
    sum_array = np.zeros((3, reader.top_info.bases, 3))


    for r in results:
         sum_array += r.get()
    #pos, a1s, a3s=  np.sum([r.get() for r in results],axis=0)/n_confs
    sum_array = sum_array/n_confs
    pos = sum_array[0]
    a1s = sum_array[1]
    a3s = sum_array[2]
    base_array(
        ref_conf.time,ref_conf.box,ref_conf.energy, pos, a1s,a3s
    ).write_new("/home/erik/Simulations/hinge_all/hinge1/mean.dat")

    print("--- %s seconds ---" % (time.time() - start_time))