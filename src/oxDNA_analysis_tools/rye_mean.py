import argparse
import os
import time
from nbformat import write
import numpy as np
from sys import stderr
from multiprocessing import Pool
from collections import namedtuple
from random import randrange
from oxDNA_analysis_tools.rye_align import align
from oxDNA_analysis_tools.UTILS.RyeReader import describe, inbox, write_conf
from oxDNA_analysis_tools.UTILS.data_structures import Configuration
from oxDNA_analysis_tools.UTILS.get_confs import get_confs
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "centered_ref_coords",
                                              "indexes",
                                              "ntopart"])
def compute(ctx:ComputeContext,chunk_id:int):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)
    confs = (inbox(c, center=True) for c in confs)
    # convert to numpy repr
    np_coords = np.asarray([[c.positions, c.a1s, c.a3s] for c in confs])
    sub_mean = np.zeros(shape=[3,ctx.top_info.nbases,3])
    for c in np_coords:
        sub_mean += align(ctx.centered_ref_coords, c, ctx.indexes)
    
    return sub_mean

def main():
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Computes the mean structure of a trajectory file")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('topology', type=str, nargs=1, help='the topology file associated with the trajectory')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the mean structure to')
    parser.add_argument('-d', '--deviations', metavar='deviation_file', nargs=1, help='Immediatley run compute_deviations.py from the output')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Compute mean structure of a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-a', '--align', metavar='alignment_configuration', nargs=1, help='The id of the configuration to align to, otherwise random')
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    # Get metadata about input files
    traj = args.trajectory[0]
    top = args.topology[0]
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

    # -a specifies the id of the reference configuration, defaults to random
    # Get the reference configuration
    if args.align:
        ref_conf_id = int(args.align[0])
    else:
        ref_conf_id = int(randrange(0, traj_info.nconfs))
    ref_conf = get_confs(traj_info.idxs, traj_info.path, ref_conf_id, 1, top_info.nbases)[0]
    ref_conf = inbox(ref_conf)

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

    # alignment requires the ref to be centered at 0
    reference_coords = ref_conf.positions[indexes]
    ref_cms = np.mean(reference_coords, axis=0) # cms prior to centering
    reference_coords = reference_coords - ref_cms

    # Create a ComputeContext which defines the problem to pass to the worker processes 
    ctx = ComputeContext(
        traj_info, top_info, reference_coords, indexes, ntopart
    )

    # Distribute jobs to the worker processes
    print(f"Starting up {ncpus} processes for {n_chunks} chunks")
    results = [pool.apply_async(compute,(ctx,i)) for i in range(n_chunks)]
    print("All spawned")

    # get the results from the workers
    acc = np.zeros([3, top_info.nbases, 3])
    for i,r in enumerate(results):
        print(f"finished {i+1}/{n_chunks}",end="\r")
        acc += r.get()
    pool.close()
    pool.join()

    print()

    # compute the mean 
    acc /= n_confs
    pos, a1s, a3s = acc

    # renormalize
    a1s = np.array([v/np.linalg.norm(v) for v in a1s])
    a3s = np.array([v/np.linalg.norm(v) for v in a3s])

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else:
        outfile = "mean.dat"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    write_conf(outfile,Configuration(0,ref_conf.box,np.array([0,0,0]), pos, a1s , a3s))
    print("--- %s seconds ---" % (time.time() - start_time))

    # -d runs compute_deviations.py
    if args.deviations:
        from sys import argv
        from oxDNA_analysis_tools import rye_deviations
        dev_file = args.deviations[0]
        print("INFO: Launching compute_deviations")

        #this is probably horrible practice, but to maintain the ability to call things from the command line, I cannot pass arguments between main() calls.
        #so instead we're gonna spoof a global variable to make it look like compute_deviations was called explicitally
        argv.clear()
        argv.extend(['rye_deviations.py', '-o', dev_file, "-r", dev_file.split('.')[0]+"_rmsd.png", "-d", dev_file.split('.')[0]+"_rmsd_data.json"])
        if args.index_file:
            argv.append("-i")
            argv.append(index_file)
        if args.parallel:
            argv.append("-p") 
            argv.append(str(ncpus))
        argv.append(outfile)
        argv.append(traj)
        argv.append(top)

        rye_deviations.main()

if __name__ == '__main__':
    main()