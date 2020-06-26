#!/usr/bin/env python3

from sys import stderr
try:
    from Bio.SVDSuperimposer import SVDSuperimposer
except:
    from bio.SVDSuperimposer import SVDSuperimposer
from json import loads, dumps
from UTILS.readers import LorenzoReader2, cal_confs
import numpy as np
import argparse
from UTILS import parallelize

def compute_deviations(reader, mean_structure, indexed_mean_structure, num_confs, start=None, stop=None):
    """
        Computes RMSF of each particle from the mean structure

        Parameters:
            reader (readers.LorenzoReader2): An active reader on the trajectory file to analyze.
            mean_structure (numpy.array): The position of each particle in the mean configuration.  A 3xN array.
            num_confs (int): The number of configurations in the reader.  
            <optional> start (int): The starting configuration ID to begin averaging at.  Used if parallel.
            <optional> stop (int): The configuration ID on which to end the averaging.  Used if parallel.

        Returns:
            deviations (list): Each entry in the list is a numpy.array of the deviations for each particle at a given time.
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

    indexed_fetch_np = lambda conf: np.array([
        n.cm_pos for n in conf._nucleotides if n.index in indexes
    ])

    # Use the single-value decomposition method for superimposing configurations 
    sup = SVDSuperimposer()
    deviations = []

    mysystem = reader._get_system(N_skip = start)
    
    while mysystem != False and confid < stop:
        mysystem.inbox()
        # calculate alignment transform
        cur_conf = fetch_np(mysystem)
        indexed_cur_conf = indexed_fetch_np(mysystem)
        sup.set(indexed_mean_structure, indexed_cur_conf)
        sup.run()
        print("Frame number:",confid, "RMSF:", sup.get_rms())
        # realign frame
        rot, tran = sup.get_rotran()
        # align structures and collect coordinates for each frame 
        # compatible with json 
        deviations.append(
           list(np.linalg.norm(np.einsum('ij, ki -> kj', rot, cur_conf) + tran - mean_structure, axis=1))
        )
        confid += 1
        mysystem = reader._get_system()

    return deviations

if __name__ == "__main__":
    #handle commandline arguments
    #the positional arguments for this are: 
    # 1. the mean structure from compute_mean.py in json format
    # 2. the trajectory from which to compute the deviations
    # 3. the name of the file to write out the deviations to.  Should be a .json because oxView uses file extensions
    parser = argparse.ArgumentParser(description="Compute the RMSD of each nucleotide from the mean structure produced by compute_mean.py")
    parser.add_argument('mean_structure', type=str, nargs=1, help="The mean structure .json file from compute_mean.py")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('topology', type=str, nargs=1, help='the topology file associted with the trajectory')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the deviations json file to')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Compute mean structure of a subset of particles from a space-separated list in the provided file')
    args = parser.parse_args()

    #system check
    from config import check_dependencies
    check_dependencies(["python", "Bio", "numpy"])

    #-o names the output file
    if args.output:
        outfile = args.output[0].strip()
    else: 
        outfile = "devs.json"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    #prepare the data files and calculate how many configurations there are to run
    top_file  = args.topology[0]
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
        with open(top_file, 'r') as f:
            indexes = list(range(int(f.readline().split(' ')[0])))

    # load mean structure 
    mean_structure_file = args.mean_structure[0]
    with open(mean_structure_file) as file:
        mean_data = loads(
            file.read()
        )
    mean_structure = np.array(mean_data["g_mean"])
    indexed_mean_structure = mean_structure[indexes]
    print("INFO: mean structure loaded", file=stderr)

    #Calculate deviations, in parallel if available
    if not parallel:
        print("INFO: Computing deviations from the mean of {} configurations using 1 core.".format(num_confs), file=stderr)
        r = LorenzoReader2(traj_file,top_file)
        deviations = compute_deviations(r, mean_structure, indexed_mean_structure, num_confs)

    #If parallel, the trajectory is split into a number of chunks equal to the number of CPUs available.
    #Each of those chunks is then calculated seperatley and the results are compiled .
    if parallel:
        print("INFO: Computing deviations from the mean of {} configurations using {} cores.".format(num_confs, n_cpus), file=stderr)
        deviations = []
        out = parallelize.fire_multiprocess(traj_file, top_file, compute_deviations, num_confs, n_cpus, mean_structure, indexed_mean_structure)
        [deviations.extend(i) for i in out]

    #compute_deviations() returns the deviation of every particle in every configuration
    #take the mean of the per-configuration deviations to get the RMSF
    means = []
    for i, devs in enumerate(np.transpose(deviations)):
        devs = np.square(devs)
        mean = np.mean(devs)
        means.append(np.sqrt(mean)*0.8518)

    #write the deviations to a json file
    print("INFO: writing deviations to {}".format(outfile), file=stderr)
    with open(outfile,"w") as file:
        file.write(
            dumps({
            "RMSF (nm)" : means
        }))
