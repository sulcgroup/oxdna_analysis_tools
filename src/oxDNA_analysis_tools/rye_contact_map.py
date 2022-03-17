import os
import numpy as np
import argparse
import matplotlib.pyplot as plt
from sys import stderr
from multiprocess import Pool
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.RyeReader import describe
from oxDNA_analysis_tools.UTILS.get_confs import get_confs
from oxDNA_analysis_tools.config import check_dependencies
from oxDNA_analysis_tools.distance import vectorized_min_image

from time import time

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "ntopart"])

def contact_map(ctx:ComputeContext,chunk_id:int):
    """
    Computes the average distance between every pair of nucleotides and creates a matrix of these distances.

    Parameters:
        ctx (ComputeContext): A named tuple containing trajectory info, topology info, and the number of configurations to process.
        chunk_id (int): The id of the chunk to process.
    
    Returns:
        distances (numpy.array): A NxN matrix containing pairwise distances between every pair of nucleotides.
    """
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)

    np_poses = np.asarray([c.positions for c in confs])
    distances = np.zeros((ctx.top_info.nbases, ctx.top_info.nbases))
    for c in np_poses:
        distances += vectorized_min_image(c, c, confs[0].box)

    return distances

def main():
    check_dependencies(["python", "numpy", "matplotlib"])

    #get commandline arguments
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Calculate and display the contact map for a structure")
    parser.add_argument('trajectory', type=str, nargs=1, help="The file containing the configurations of which the contact map is needed")
    parser.add_argument('topology', type=str, nargs=1, help="The topology associated with the configuration")
    parser.add_argument('-g', metavar='graph', dest='graph', nargs=1, type=str, help='Filename for the plot')
    parser.add_argument('-d', metavar='data', dest='data', nargs=1, help='The name of the file to save the contact map as a pickle.')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")

    # Get arguments and file metadata
    args = parser.parse_args()
    traj = args.trajectory[0]
    top = args.topology[0]
    top_info, traj_info = describe(top, traj)

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    # how many confs we want to distribute between the processes
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

    # Normalize the distances and convert to nm
    distances /= n_confs
    distances *= 0.8518

    # Plot the contact map
    if args.graph:
        graph_name = args.graph[0]
    else:
        print("INFO: No graph name provided, defaulting to 'contact_map.png'", file=stderr)
        graph_name = "contact_map.png"

    fig, ax = plt.subplots()
    a = ax.imshow(distances, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
    ylabel="nucleotide id",
    xlabel="nucleotide id")
    b = fig.colorbar(a, ax=ax)
    b.set_label("distance (nm)", rotation = 270, labelpad=15)
    plt.tight_layout()
    print("INFO: Saving contact map to '{}'".format(graph_name), file=stderr)
    plt.savefig(graph_name)

    # Save the contact map as a pickle
    if args.data:
        data_name = args.data[0]
    else:
        print("INFO: No data name provided, defaulting to 'contact_map.pkl'", file=stderr)
        data_name = "contact_map.pkl"
    print("INFO: Saving contact map to '{}'".format(data_name), file=stderr)
    np.save(data_name, distances)

if __name__ == "__main__":
    main()