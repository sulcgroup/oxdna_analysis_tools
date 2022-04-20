from sys import stderr
from multiprocessing import Pool
from collections import namedtuple
import numpy as np
import argparse
import os
from oxDNA_analysis_tools.UTILS.RyeReader import describe, inbox, write_conf, write_conf
from oxDNA_analysis_tools.UTILS.get_confs import get_confs
from oxDNA_analysis_tools.UTILS.data_structures import Configuration
from oxDNA_analysis_tools.rye_align import align

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "ref_coords",
                                              "indexes",
                                              "ntopart"])


def compute_centroid(ctx:ComputeContext, chunk_id:int):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)
    confs = [inbox(c) for c in confs]
    np_confs = np.asarray([[c.positions, c.a1s, c.a3s] for c in confs])
    centroid_candidate = np.zeros_like(np_confs[0])
    min_RMSD = np.inf
    centroid_id = -1

    for i, c in enumerate(np_confs):
        c[0] -= np.mean(c[0][ctx.indexes], axis=0) #didn't center earlier because you have to center on the indexed particles
        aligned_conf = align(ctx.ref_coords.positions[ctx.indexes], c, ctx.indexes)[0]
        RMSD = np.sqrt(np.mean(np.linalg.norm(aligned_conf[ctx.indexes] - ctx.ref_coords.positions[ctx.indexes], axis=1)**2))
        if RMSD < min_RMSD:
            min_RMSD = RMSD
            centroid_candidate = c
            centroid_id = i

    t = confs[centroid_id].time

    return (centroid_candidate, min_RMSD, t)

def main():
    #handle commandline arguments
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Find the configuration in a trajectory closest to a provided reference configuration")
    parser.add_argument('reference_structure', type=str, nargs=1, help="The reference structure to search against")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('topology', type=str, nargs=1, help='the topology file associated with the trajectory')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the centroid to')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Alignment and RMSD based on a subset of particles given in a space-separated list in the provided file')
    args = parser.parse_args()

    #system check
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    #Get file paths
    ref = args.reference_structure[0].strip()
    traj = args.trajectory[0].strip()
    top = args.topology[0].strip()
    _, ref_info = describe(top, ref)
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
    ref_conf = get_confs(ref_info.idxs, ref_info.path, 0, 1, top_info.nbases)[0]
    ref_conf = inbox(ref_conf)
    ref_cms = np.mean(ref_conf.positions[indexes], axis=0)
    ref_conf.positions -= ref_cms

    # create a ComputeContext which defines the problem to pass to the worker processes
    ctx = ComputeContext(
        traj_info, top_info, ref_conf, indexes, ntopart
    )

    # Distribute jobs to the worker processes
    print(f"Starting up {ncpus} processes for {n_chunks} chunks")
    results = [pool.apply_async(compute_centroid,(ctx,i)) for i in range(n_chunks)]
    print("All spawned")

    min_RMSD = np.inf
    centroid_candidate = Configuration(0, ref_conf.box, np.zeros(3), np.zeros_like(ref_conf.positions), np.zeros_like(ref_conf.positions), np.zeros_like(ref_conf.positions))
    centroid_time = -1

    # Collect results from the worker processes
    for i, r in enumerate(results):
        print(f"finished {i+1}/{n_chunks}",end="\r")
        centroid, RMSD, t = r.get()
        if RMSD < min_RMSD:
            min_RMSD = RMSD
            centroid_time = t
            centroid_candidate.positions = centroid[0]
            centroid_candidate.a1s = centroid[1]
            centroid_candidate.a3s = centroid[2]
    
    min_RMSD *= 0.8518
    centroid_candidate.time = centroid_time

    #-o names the output file
    if args.output:
        outfile = args.output[0].strip()
    else: 
        outfile = "centroid.dat"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    write_conf(outfile, centroid_candidate)
    print("INFO: Wrote centroid to {}".format(outfile), file=stderr)
    print("INFO: Min RMSD: {} nm".format(min_RMSD), file=stderr)
    print("INFO: Centroid time: {}".format(centroid_time), file=stderr)

if __name__ == '__main__':
    main()