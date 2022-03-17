import argparse
import os
import numpy as np
from sys import stderr
from multiprocessing import Pool
from collections import namedtuple
from json import dumps
from oxDNA_analysis_tools.UTILS.RyeReader import describe, inbox, write_conf
from oxDNA_analysis_tools.UTILS.data_structures import Configuration
from oxDNA_analysis_tools.UTILS.get_confs import get_confs
from oxDNA_analysis_tools.rye_mean import align
import matplotlib.pyplot as plt
import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "mean_coords",
                                              "indexes",
                                              "ntopart"])

def compute(ctx:ComputeContext,chunk_id:int):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)
    confs = (inbox(c, center=True) for c in confs)
    confs = np.asarray([[c.positions, c.a1s, c.a3s] for c in confs])

    SFs = np.empty((ctx.ntopart, ctx.top_info.nbases))
    for i, c in enumerate(confs):
        aligned_conf = align(ctx.mean_coords.positions[ctx.indexes], c, ctx.indexes)[0]
        SFs[i] = np.power(np.linalg.norm(aligned_conf - ctx.mean_coords.positions, axis=1), 2)

    return SFs

def main():
    #handle commandline arguments
    #the positional arguments for this are: 
    # 1. the mean structure from compute_mean.py in json format
    # 2. the trajectory from which to compute the deviations
    # 3. the topology associated with the trajectory
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Compute the RMSD of each nucleotide from the mean structure produced by compute_mean.py")
    parser.add_argument('mean_structure', type=str, nargs=1, help="The mean structure .json file from compute_mean.py")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('topology', type=str, nargs=1, help='the topology file associated with the trajectory')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the RMSF json file to')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Compute mean structure of a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-r', metavar='rmsd_plot', dest='rmsd_plot', nargs=1, help='The name of the file to save the RMSD plot to.')
    parser.add_argument('-d', metavar='rmsd_data', dest='rmsd_data', nargs=1, help='The name of the file to save the RMSD data in json format.')
    args = parser.parse_args()

    #system check
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy", "matplotlib"])

    # Get metadata about input files
    mean = args.mean_structure[0]
    traj = args.trajectory[0]
    top = args.topology[0]
    _, mean_info = describe(top, mean)
    top_info, traj_info = describe(top, traj)

    # -i comes with a list of particles indices representing a subset to compute the mean against.
    # Get the index list which is a space-separated list of particle ids.
    if args.index_file:
        index_file = args.index_file[0]
        with open(index_file, 'r') as f:
            indexes = f.readline().split()
            try:
                indexes = [int(i) for i in indexes]
            except:
                print("ERROR: The index file must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button")
    else:
        indexes = list(range(top_info.nbases))

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

    # get the mean structure from the file path
    mean_conf = get_confs(mean_info.idxs, mean_info.path, 0, 1, top_info.nbases)[0]
    mean_conf = inbox(mean_conf)
    ref_cms = np.mean(mean_conf.positions[indexes], axis=0)
    mean_conf.positions -= ref_cms

    # create a ComputeContext which defines the problem to pass to the worker processes
    ctx = ComputeContext(
        traj_info, top_info, mean_conf, indexes, ntopart
    )

    # Distribute jobs to the worker processes
    print(f"Starting up {ncpus} processes for {n_chunks} chunks")
    results = [pool.apply_async(compute,(ctx,i)) for i in range(n_chunks)]
    print("All spawned")

    # get the results from the workers
    MFs = np.empty((traj_info.nconfs, top_info.nbases))
    for i,r in enumerate(results):
        print(f"finished {i+1}/{n_chunks}",end="\r")
        MFs[i*ctx.ntopart:i*ctx.ntopart+ntopart] = r.get()
    pool.close()
    pool.join()

    # Compute RMSDs and RMSF
    RMSDs = np.sqrt(np.mean(MFs, axis=1)) * 0.8518
    RMSFs = np.sqrt(np.mean(MFs, axis=0)) * 0.8518

    #-o names the output file
    if args.output:
        outfile = args.output[0].strip()
        if not outfile.split(".")[-1] == 'json':
            outfile += ".json"
    else: 
        outfile = "devs.json"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    #-r names the file to print the RMSD plot to
    if args.rmsd_plot:
        plot_name = args.rmsd_plot[0]
    else:
        plot_name = 'rmsd.png'

    # -d names the file to print the RMSD data to
    if args.rmsd_data:
        data_file = args.rmsd_data[0]
    else:
        data_file = 'rmsd_op.json'

    # Save the RMSDs and RMSFs to json files
    print("INFO: writing deviations to {}".format(outfile), file=stderr)
    with open(outfile, 'w') as f:
        f.write(
            dumps({
                "RMSF (nm)" : RMSFs.tolist()
            })
        )

    print("INFO: writing RMSDs to oxView order parameter file, {}".format(data_file))
    with open(data_file, 'w') as f:
        f.write(
            dumps({
                "RMSD (nm)" : RMSDs.tolist()
            })
        )

    print("INFO: writing RMSD plot to {}".format(plot_name))
    plt.plot(RMSDs)
    plt.axhline(np.mean(RMSDs), color='red')
    plt.xlabel('Configuration')
    plt.ylabel('RMSD (nm)')
    plt.savefig(plot_name)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()