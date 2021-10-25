#!/usr/bin/env python3

from sys import stderr
try:
    from Bio.SVDSuperimposer import SVDSuperimposer
except:
    from bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import argparse
import os
from json import load
from oxDNA_analysis_tools.UTILS import parallelize_erik_onefile
from oxDNA_analysis_tools.UTILS.readers import ErikReader, cal_confs

def compute_centroid(reader, mean_structure, indexes, num_confs, start=None, stop=None):
    """
        Compares each structure to the mean and returns the one with the lowest RMSF

        Parameters:
            reader (readers.LorenzoReader2): An active reader on the trajectory file to analyze.
            mean_structure (numpy.array): The position of each particle in the mean configuration.  A 3xN array.
            num_confs (int): The number of configurations in the reader.  
            <optional> start (int): The starting configuration ID to begin averaging at.  Used if parallel.
            <optional> stop (int): The configuration ID on which to end the averaging.  Used if parallel.

        Returns:
            centroid (numpy.array): The positions corresponding to the structure with the lowest RMSF to the mean.
    """
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)
    confid = 0

    # Use the single-value decomposition method for superimposing configurations 
    sup = SVDSuperimposer()
    lowest_rmsf = 100000 #if you have a larger number than this, we need to talk...
    centroid_candidate = np.zeros_like(mean_structure)
    centroid_a1 = np.zeros_like(mean_structure)
    centroid_a3 = np.zeros_like(mean_structure)

    mysystem = reader.read(n_skip = start)
    
    while mysystem != False and confid < stop:
        mysystem.inbox()
        # calculate alignment transform
        cur_conf = mysystem.positions
        indexed_cur_conf = mysystem.positions[indexes]
        cur_conf_a1 = mysystem.a1s
        cur_conf_a3 = mysystem.a3s
        sup.set(mean_structure, indexed_cur_conf)
        sup.run()
        rot, tran = sup.get_rotran()

        cur_conf = np.einsum('ij, ki -> kj', rot, cur_conf) + tran
        cur_conf_a1 = np.einsum('ij, ki -> kj', rot, cur_conf_a1)
        cur_conf_a3 = np.einsum('ij, ki -> kj', rot, cur_conf_a3)
        RMSF = sup.get_rms()
        print("Frame number:",confid, "RMSF:", RMSF)
        if RMSF < lowest_rmsf:
            centroid_candidate = cur_conf
            centroid_a1 = cur_conf_a1
            centroid_a3 = cur_conf_a3
            lowest_rmsf = RMSF
            centroid_t = mysystem.time
        
        confid += 1
        mysystem = reader.read()

    return centroid_candidate, centroid_a1, centroid_a3, lowest_rmsf, centroid_t

def main():
    #handle commandline arguments
    #the positional arguments for this are: 
    # 1. the mean structure from compute_mean.py in json format
    # 2. the trajectory from which to compute the centroid
    # 3. the name of the file to write out the centroid to.  Should be a .dat because oxView uses file extensions
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Compute the RMSD of each nucleotide from the mean structure produced by compute_mean.py")
    parser.add_argument('mean_structure', type=str, nargs=1, help="The mean structure .json file from compute_mean.py")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the centroid to')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Compute mean structure of a subset of particles from a space-separated list in the provided file')
    args = parser.parse_args()

    #system check
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "Bio", "numpy"])

    #-o names the output file
    if args.output:
        outfile = args.output[0].strip()
    else: 
        outfile = "centroid.dat"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    #prepare the data files and calculate how many configurations there are to run
    traj_file = args.trajectory[0]
    parallel = args.parallel
    if parallel:
        n_cpus = args.parallel[0]
    num_confs = cal_confs(traj_file)

    #-i will make it only run on a subset of nucleotides.
    #The index file is a space-separated list of particle IDs
    if args.index_file:
        index_file = args.index_file[0]
        with open(index_file, 'r') as f:
            indexes = f.readline().split()
            try:
                indexes = [int(i) for i in indexes]
            except:
                print("ERROR: The index file must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button")
    else: 
        with ErikReader(traj_file) as r:
            indexes = list(range(len(r.read().positions)))

    # load mean structure 
    mean_file = args.mean_structure[0]
    if mean_file.split(".")[-1] == "json":
        with open(mean_file) as file:
            mean_structure = load(file)['g_mean'][indexes]

    elif mean_file.split(".")[-1] == "dat":
        with ErikReader(mean_file) as reader:
            s = reader.read()
            mean_structure = s.positions[indexes]
    print("INFO: mean structure loaded", file=stderr)

    #Calculate centroid, in parallel if available
    if not parallel:
        print("INFO: Computing centroid from the mean of {} configurations using 1 core.".format(num_confs), file=stderr)
        r = ErikReader(traj_file)
        centroid, centroid_a1s, centroid_a3s, centroid_rmsf, centroid_time = compute_centroid(r, mean_structure, indexes, num_confs)

    #If parallel, the trajectory is split into a number of chunks equal to the number of CPUs available.
    #Each of those chunks is then calculated seperatley and the results are compiled .
    if parallel:
        print("INFO: Computing centroid from the mean of {} configurations using {} cores.".format(num_confs, n_cpus), file=stderr)
        candidates = []
        rmsfs = []
        a1s = []
        a3s = []
        ts = []
        out = parallelize_erik_onefile.fire_multiprocess(traj_file, compute_centroid, num_confs, n_cpus, mean_structure, indexes)
        [candidates.append(i[0]) for i in out]
        [rmsfs.append(i[3]) for i in out]
        [a1s.append(i[1]) for i in out]
        [a3s.append(i[2]) for i in out]
        [ts.append(i[4]) for i in out]
        min_id = rmsfs.index(min(rmsfs))
        centroid = candidates[min_id]
        centroid_a1s = a1s[min_id]
        centroid_a3s = a3s[min_id]
        centroid_time = ts[min_id]
        centroid_rmsf = rmsfs[min_id]

    print("INFO: Centroid configuration found at configuration t = {}, RMSF = {}".format(centroid_time, centroid_rmsf), file=stderr)
    
    from oxDNA_analysis_tools.mean2dat import make_dat

    make_dat({'g_mean' : centroid, 'a1_mean': centroid_a1s, 'a3_mean': centroid_a3s}, outfile)

if __name__ == '__main__':
    main()
    
