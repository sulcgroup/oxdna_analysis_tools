#!/usr/bin/env python3
#Multidimensional_scaling_mean
#Written by: Erik Poppleton
#Date: 3/28/2022
#Computes the RMSD of the contact map for a structure.  The average structure is determined
#by Scikit.learn's MDS algorithm, then subtracts the contact map of each individual structure from the man
#This is used to compute a per-nucleotide deviation in the contact map, which can be visualized with oxView

from matplotlib.pyplot import box
from nbformat import write
import numpy as np
import argparse
from os import path
from sys import exit, stderr
from json import dumps
from collections import namedtuple
from multiprocessing import Pool
from sklearn.manifold import MDS
from oxDNA_analysis_tools.config import check_dependencies
from oxDNA_analysis_tools.rye_contact_map import contact_map
from oxDNA_analysis_tools.distance import vectorized_min_image
from oxDNA_analysis_tools.UTILS.RyeReader import describe, write_conf
from oxDNA_analysis_tools.UTILS.data_structures import Configuration
from oxDNA_analysis_tools.UTILS.get_confs import get_confs

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "ntopart"])

DevsContext = namedtuple("DevsContext",["traj_info",
                                        "top_info",
                                        "ntopart",
                                        "masked_mean_coords"])

#at 2.5 you start to see the hard edges caused by end-loops and see some loop interactions
CUTOFF = 2.5

def make_heatmap(contact_map):
    """
    Convert a matrix of contact distances to a visual contact map.

    Parameters:
        contact_map (numpy.array): An array of all pairwise distances between nucleotides.
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    a = ax.imshow(contact_map, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
       ylabel="nucleotide id",
       xlabel="nucleotide id")
    b = fig.colorbar(a, ax=ax)
    b.set_label("distance", rotation = 270)
    plt.show()

def devs_mds(ctx:DevsContext, chunk_id:int, ):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)
    
    np_poses = np.asarray([c.positions for c in confs])

    devs = np.zeros((ctx.top_info.nbases, ctx.top_info.nbases))

    for c in np_poses:
        c_map = vectorized_min_image(c, c, confs[0].box)
        masked_distances = np.ma.masked_array(c_map, ~(c_map < CUTOFF))

        # Fill the masked values with the cutoff.  Not sure if this is the best practice here.
        masked_distances = np.ma.filled(masked_distances, CUTOFF)
        masked_mean = np.ma.filled(ctx.masked_mean_coords, CUTOFF)

        diff = masked_distances - masked_mean
        diff = np.square(diff)
        devs += diff

    return devs



def main():
    #get commandline arguments
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Calculate molecular contacts, and assembles an average set of contacts based on MDS")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('topology', type=str, nargs=1, help='the topology file associated with the trajectory')
    parser.add_argument('-o', '--output', metavar='output', type=str, nargs=1, help='the name of the .dat file where the mean will be written')
    parser.add_argument('-d', '--dev_file', metavar='dev_file', type=str, nargs=1, help='the name of the .json file where the devs will be written')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    args = parser.parse_args()
    traj = args.trajectory[0]
    top = args.topology[0]
    top_info, traj_info = describe(top, traj)
    example_conf = get_confs(traj_info.idxs, traj_info.path, 1, 1, top_info.nbases)[0]

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    check_dependencies(['python', 'numpy'])

    ntopart = 20
    pool = Pool(ncpus)

    # deduce how many chunks we have to run in parallel
    n_confs  = traj_info.nconfs 
    n_chunks = int(n_confs / ntopart +
                         (1 if n_confs % ntopart else 0))

    ctx = ComputeContext(traj_info, top_info, ntopart)

    # Distribute jobs to the worker processes
    print(f"Starting up {ncpus} processes for {n_chunks} chunks")
    results = [pool.apply_async(contact_map,(ctx,i)) for i in range(n_chunks)]
    print("All spawned")

    # Accumulate the results
    distances = np.zeros((top_info.nbases, top_info.nbases))
    for i, r in enumerate(results):
        print(f"Finished {i+1}/{n_chunks}")
        distances += r.get()
    pool.close()
    pool.join()

    mean_distances = distances / n_confs
    masked_mean = np.ma.masked_array(mean_distances, ~(mean_distances < CUTOFF))

    print("INFO: fitting local distance data", file=stderr)
    mds = MDS(n_components=3, metric=True, max_iter=3000, eps=1e-12, dissimilarity="precomputed", n_jobs=1, n_init=1)
    out_coords = mds.fit_transform(masked_mean, init=example_conf.positions) #without the init you can get a left-handed structure.
    a1s = np.zeros((top_info.nbases, 3))
    a3s = np.zeros((top_info.nbases, 3))

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else:
        outfile = "mean_mds.dat"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    write_conf(outfile,Configuration(0,example_conf.box, np.array([0,0,0]), out_coords, a1s , a3s))
    print("INFO: Wrote mean to {}".format(outfile), file=stderr)

    # Compute the deviations from the mean
    ctx = DevsContext(traj_info, top_info, ntopart, masked_mean)
    pool = Pool(ncpus)

    # Distribute jobs to the worker processes
    print(f"Starting up {ncpus} processes for {n_chunks} chunks")
    results = [pool.apply_async(devs_mds,(ctx,i)) for i in range(n_chunks)]
    print("All spawned")

    # Accumulate the results
    devs = np.zeros((top_info.nbases, top_info.nbases))
    for i, r in enumerate(results):
        print(f"Finished {i+1}/{n_chunks}")
        devs += r.get()
    pool.close()
    pool.join()

    np.ma.masked_array(devs, ~(devs != 0.0))
    devs = devs / n_confs
    devs = np.mean(devs, axis=0)
    devs = np.sqrt(devs)

    #-d names the deviations file
    if args.dev_file:
        devfile = args.dev_file[0].split(".")[0] + ".json"
    else:
        devfile = "devs_mds.json"
        print("INFO: No deviations file name provided, defaulting to \"{}\"".format(devfile), file=stderr)

    with open(devfile, "w") as file:
        file.write(
            dumps({"contact deviation" : list(devs)})
        )
    print("INFO: wrote file {}".format(devfile), file=stderr)

if __name__ == '__main__':
    main()