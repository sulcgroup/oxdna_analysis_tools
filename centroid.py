#!/usr/bin/env python3

from sys import stderr
from Bio.SVDSuperimposer import SVDSuperimposer
from json import loads, dumps
from UTILS.readers import LorenzoReader2, cal_confs
import numpy as np
import argparse
from UTILS import parallelize
from json import load

def compute_centroid(reader, mean_structure, num_confs, start=None, stop=None):
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

    # helper to fetch nucleotide positions 
    fetch_np = lambda conf: np.array([
        n.cm_pos for n in conf._nucleotides 
    ])
    fetch_a1 = lambda conf: np.array([
        n._a1 for n in conf._nucleotides
    ])
    fetch_a3 = lambda conf: np.array([
        n._a3 for n in conf._nucleotides
    ])

    # Use the single-value decomposition method for superimposing configurations 
    sup = SVDSuperimposer()
    lowest_rmsf = 100000 #if you have a larger number than this, we need to talk...
    centroid_candidate = np.zeros_like(mean_structure)
    centroid_a1 = np.zeros_like(mean_structure)
    centroid_a3 = np.zeros_like(mean_structure)

    mysystem = reader._get_system(N_skip = start)
    
    while mysystem != False and confid < stop:
        mysystem.inbox_system()
        # calculate alignment transform
        cur_conf = fetch_np(mysystem)
        cur_conf_a1 = fetch_a1(mysystem)
        cur_conf_a3 = fetch_a3(mysystem)
        sup.set(mean_structure, cur_conf)
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
        
        confid += 1
        mysystem = reader._get_system()

    return centroid_candidate, centroid_a1, centroid_a3, lowest_rmsf

if __name__ == "__main__":
    #handle commandline arguments
    #the positional arguments for this are: 
    # 1. the mean structure from compute_mean.py in json format
    # 2. the trajectory from which to compute the centroid
    # 3. the name of the file to write out the centroid to.  Should be a .dat because oxView uses file extensions
    parser = argparse.ArgumentParser(description="Compute the RMSD of each nucleotide from the mean structure produced by compute_mean.py")
    parser.add_argument('mean_structure', type=str, nargs=1, help="The mean structure .json file from compute_mean.py")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('topology', type=str, nargs=1, help='the topology file associted with the trajectory')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the centroid to')
    args = parser.parse_args()

    #system check
    from config import check_dependencies
    check_dependencies(["python", "Bio", "numpy"])
    from config import set_reference
    INBOXING_REFERENCE_PARTICLE = set_reference()

    #-o names the output file
    if args.output:
        outfile = args.output[0].strip()
    else: 
        outfile = "centroid.dat"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    #prepare the data files and calculate how many configurations there are to run
    top_file  = args.topology[0]
    traj_file = args.trajectory[0]
    parallel = args.parallel
    if parallel:
        n_cpus = args.parallel[0]
    num_confs = cal_confs(traj_file)

    # load mean structure 
    mean_file = args.mean_structure[0]
    if mean_file.split(".")[-1] == "json":
        with open(mean_file) as file:
            mean_structure = load(file)['g_mean']

    elif mean_file.split(".")[-1] == "dat":
        fetch_np = lambda conf: np.array([
            n.cm_pos for n in conf._nucleotides
        ])
        with LorenzoReader2(mean_file, top_file) as reader:
            s = reader._get_system()
            mean_structure = fetch_np(s)
    print("INFO: mean structure loaded", file=stderr)

    #Calculate centroid, in parallel if available
    if not parallel:
        print("INFO: Computing centroid from the mean of {} configurations using 1 core.".format(num_confs), file=stderr)
        r = LorenzoReader2(traj_file,top_file)
        centroid, centroid_a1s, centroid_a3s, rmsf = compute_centroid(r, mean_structure, num_confs)

    #If parallel, the trajectory is split into a number of chunks equal to the number of CPUs available.
    #Each of those chunks is then calculated seperatley and the results are compiled .
    if parallel:
        print("INFO: Computing centroid from the mean of {} configurations using {} cores.".format(num_confs, n_cpus), file=stderr)
        candidates = []
        rmsfs = []
        a1s = []
        a3s = []
        out = parallelize.fire_multiprocess(traj_file, top_file, compute_centroid, num_confs, n_cpus, mean_structure)
        [candidates.append(i[0]) for i in out]
        [rmsfs.append(i[3]) for i in out]
        [a1s.append(i[1]) for i in out]
        [a3s.append(i[2]) for i in out]
        min_id = rmsfs.index(min(rmsfs))
        centroid = candidates[min_id]
        centroid_a1s = a1s[min_id]
        centroid_a3s = a3s[min_id]
    
    from mean2dat import make_dat

    make_dat({'g_mean' : centroid, 'a1_mean': centroid_a1s, 'a3_mean': centroid_a3s}, "centroid.dat")
