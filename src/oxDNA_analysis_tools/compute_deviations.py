#!/usr/bin/env python3

from sys import stderr
try:
    from Bio.SVDSuperimposer import SVDSuperimposer
except:
    from bio.SVDSuperimposer import SVDSuperimposer
import matplotlib.pyplot as plt
from json import loads, dumps
import numpy as np
import argparse
import os
from oxDNA_analysis_tools.UTILS.readers import ErikReader, cal_confs

def compute_deviations(reader, mean_structure, indexed_mean_structure, indexes, num_confs, start=None, stop=None):
    """
        Computes RMSF of each particle from the mean structure

        Parameters:
            reader (readers.ErikReader): An active reader on the trajectory file to analyze.
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

    # Use the single-value decomposition method for superimposing configurations 
    sup = SVDSuperimposer()
    deviations = []
    RMSDs = []

    mysystem = reader.read(n_skip = start)
    
    while mysystem != False and confid < stop:
        mysystem.inbox()
        # calculate alignment transform
        cur_conf = mysystem.positions
        indexed_cur_conf = cur_conf[indexes]
        sup.set(indexed_mean_structure, indexed_cur_conf)
        sup.run()
        print("Frame number:",confid, "Time:", mysystem.time, "RMSD:", sup.get_rms())
        # realign frame
        rot, tran = sup.get_rotran()
        # align structures and collect coordinates for each frame 
        # compatible with json 
        deviations.append(
           list(np.linalg.norm(np.einsum('ij, ki -> kj', rot, cur_conf) + tran - mean_structure, axis=1))
        )
        RMSDs.append(sup.get_rms()*0.8518)
        confid += 1
        mysystem = reader.read()

    return (deviations, RMSDs)

def main():
    #handle commandline arguments
    #the positional arguments for this are: 
    # 1. the mean structure from compute_mean.py in json format
    # 2. the trajectory from which to compute the deviations
    from sys import argv
    print(argv)
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Compute the RMSD of each nucleotide from the mean structure produced by compute_mean.py")
    parser.add_argument('mean_structure', type=str, nargs=1, help="The mean structure .json file from compute_mean.py")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the deviations json file to')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Compute mean structure of a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-r', metavar='rmsd_plot', dest='rmsd_plot', nargs=1, help='The name of the file to save the RMSD plot to.')
    parser.add_argument('-d', metavar='rmsd_data', dest='rmsd_data', nargs=1, help='The name of the file to save the RNSD data in json format.')
    args = parser.parse_args()

    #system check
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "Bio", "numpy", "matplotlib"])

    #-o names the output file
    if args.output:
        outfile = args.output[0].strip()
        if not outfile.split(".")[-1] == 'json':
            outfile += ".json"
    else: 
        outfile = "devs.json"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    #prepare the data files and calculate how many configurations there are to run
    traj_file = args.trajectory[0]
    parallel = args.parallel
    if parallel:
        from oxDNA_analysis_tools.UTILS import parallelize_erik_onefile
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

    #-r names the file to print the RMSD plot to
    if args.rmsd_plot:
        plot_name = args.rmsd_plot[0]
    else:
        plot_name = 'rmsd.png'

    # -d names the file to print the RMSD data to
    if args.rmsd_data:
        data_file = args.rmsd_data[0]

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
        print("INFO: Computing deviations from the mean of {} configurations with an alignment of {} particles using 1 core.".format(num_confs, len(indexed_mean_structure)), file=stderr)
        r = ErikReader(traj_file)
        deviations, RMSDs = compute_deviations(r, mean_structure, indexed_mean_structure, indexes, num_confs)

    #If parallel, the trajectory is split into a number of chunks equal to the number of CPUs available.
    #Each of those chunks is then calculated seperatley and the results are compiled .
    if parallel:
        print("INFO: Computing deviations from the mean of {} configurations with an alignment of {} particles using {} cores.".format(num_confs, len(indexed_mean_structure), n_cpus), file=stderr)
        deviations = []
        RMSDs = []
        out = parallelize_erik_onefile.fire_multiprocess(traj_file, compute_deviations, num_confs, n_cpus, mean_structure, indexed_mean_structure, indexes)
        [deviations.extend(i[0]) for i in out]
        [RMSDs.extend(i[1]) for i in out]

    #compute_deviations() returns the deviation of every particle in every configuration
    #take the mean of the per-configuration deviations to get the RMSF
    rmsfs = np.sqrt(np.mean(np.square(np.array(deviations)), axis=0))*0.8518

    #write the deviations to a json file
    print("INFO: writing deviations to {}".format(outfile), file=stderr)
    with open(outfile,"w") as file:
        file.write(
            dumps({
            "RMSF (nm)" : rmsfs.tolist()
        }))

    #plot RMSDs
    print("INFO: writing RMSD plot to {}".format(plot_name), file=stderr)
    plt.plot(RMSDs)
    plt.axhline(np.mean(RMSDs), color='red')
    plt.xlabel('Configuration')
    plt.ylabel('RMSD (nm)')
    plt.savefig(plot_name)

    #print RMSDs
    print("INFO: writing RMSD data to {}".format(data_file), file=stderr)
    if args.rmsd_data:
        with open(data_file, 'w') as f:
            f.write(dumps({
                "RMSD (nm)" : RMSDs
            }))

if __name__ == '__main__':
    main()
    