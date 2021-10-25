#!/usr/bin/env python3
try:
    from Bio.SVDSuperimposer import SVDSuperimposer
except:
    from bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
from json import loads, dumps
from sys import exit, stderr, argv
from oxDNA_analysis_tools.UTILS.readers import ErikReader, cal_confs
from random import randint
import argparse
import time
import os

from oxDNA_analysis_tools.compute_deviations import compute_deviations
start_t = time.time()

def pick_starting_configuration(traj_file, max_bound, align=None):
    """
        Pick a random conf out of the trajectory file to use as the reference structure.
        
        We assume that that is optimal to align against.  Based on experience, the choice of reference configuration has very little impact on the mean structure

        Parameters: 
            traj_file (string): The name of the trajectory file
            max_bound (int): The reference configuration will be chosen at random from the first max_bound configurations in the trajectory file

        Returns:
            stop_at (int): The configuration ID of the reference configuration
            initial_structure (base.System): The oxDNA system representing the reference configuration.
    """
    with ErikReader(traj_file) as reader:
        if align:
            stop_at = int(align)
        else:
            stop_at = randint(0, max_bound-1)
        print("INFO: We chose {} as reference".format(stop_at), file=stderr)
        initial_structure = reader.read(n_skip=stop_at)
        if not initial_structure:
            print("ERROR: Couldn't read structure at conf num {0}.".format(stop_at), file=stderr)
            exit(1)
        print("INFO: reference structure loaded", file=stderr)
        initial_structure.inbox()
    return stop_at, initial_structure

def normalize(v):
    """
        Return a normalized copy of vector v

        Parameters: 
            v (numpy.array): the vector to be normalized.

        Returns: 
            v / norm(v) (numpy.array).
    """
    norm = np.linalg.norm(v)
    if norm == 0:
        print("WARNING: tried to normalize a vector with 0 length", file=stderr)
        return v
    return v / norm

# helper to prepare a configuration of np.array coordinates
# into smth json is able to serialize
prep_pos_for_json = lambda conf: list(
                        list(p) for p in conf
                    )

def compute_mean (reader, align_conf, indexes, num_confs, start = None, stop = None):
    """
        Computes the mean structure of a trajectory

        Structured to work with the multiprocessing process from UTILS/parallelize.py

        Parameters:
            reader (readers.ErikReader): An active reader on the trajectory file to take the mean of.
            align_conf (numpy.array): The position of each particle in the reference configuration.  A 3xN array.
            num_confs (int): The number of configurations in the reader.  
            <optional> start (int): The starting configuration ID to begin averaging at.  Used if parallel.
            <optional> stop (int): The configuration ID on which to end the averaging.  Used if parallel.
        
        Returns:
            mean_pos_storage (numpy.array): For each particle, the sum of positions in all configurations read.
            mean_a1_storage (numpy.array): For each particle, the sum of a1 orientation vectors in all configuraitons read.
            mean_a3_storage (numpy.array): For each particle, the sum of a3 orientation vectors in all configuraitons read.
            intermediate_mean_structures (list): mean structures computed periodically during the summing to check decoorrelation.
            confid (int): the number of configurations summed for the storage arrays.
    """
    parallel = True
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
        parallel = False
    else: start = int(start)

    mysystem = reader.read(n_skip = start)
    n_nuc = len(mysystem.positions)

    if not parallel:
        #This also computes the mean every num_confs/10 configurations to check decorrelation.
        #Only works when run in serial.
        INTERMEDIATE_EVERY = np.floor(num_confs / 10)
        # storage for the intermediate mean structures

    intermediate_mean_structures = []

    # the class doing the alignment of 2 structures
    sup = SVDSuperimposer()

    mean_pos_storage = np.zeros((n_nuc, 3))
    mean_a1_storage = np.zeros((n_nuc, 3))
    mean_a3_storage = np.zeros((n_nuc, 3))

    # for every conf in the current trajectory we calculate the global mean
    confid = 0

    while mysystem != False and confid < stop:
        mysystem.inbox()
        indexed_cur_conf_pos = mysystem.positions[indexes]

        # calculate alignment
        sup.set(align_conf, indexed_cur_conf_pos)
        sup.run()
        rot, tran = sup.get_rotran()

        #apply alignment
        mysystem.positions = np.einsum('ij, ki -> kj', rot, mysystem.positions) + tran
        mysystem.a1s = np.einsum('ij, ki -> kj', rot, mysystem.a1s)
        mysystem.a3s = np.einsum('ij, ki -> kj', rot, mysystem.a3s)
        mean_pos_storage += mysystem.positions
        mean_a1_storage += mysystem.a1s
        mean_a3_storage += mysystem.a3s

        # print the rmsd of the alignment in case anyone is interested...
        print("Frame:", confid, "Time:", mysystem.time, "RMSF:", sup.get_rms())
        # thats all we do for a frame
        confid += 1
        mysystem = reader.read()

        # We produce 10 intermediate means to check decorrelation.
        # This can't be done neatly in parallel
        if not parallel and confid % INTERMEDIATE_EVERY == 0:
            mp = np.copy(mean_pos_storage)
            mp /= confid
            intermediate_mean_structures.append(
                prep_pos_for_json(mp)
            )
            print("INFO: Calculated intermediate mean for {} ".format(confid))

    return(mean_pos_storage, mean_a1_storage, mean_a3_storage, intermediate_mean_structures, confid)



def main():
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Computes the mean structure of a trajectory file")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the mean structure to')
    parser.add_argument('-f', '--format', metavar='<json/oxDNA/both>', nargs=1, help='Output format for the mean file.  Defaults to json.  Options are \"json\", \"oxdna/oxDNA\", and \"both\"')
    parser.add_argument('-d', '--deviations', metavar='deviation_file', nargs=1, help='Immediatley run compute_deviations.py from the output')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Compute mean structure of a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-a', '--align', metavar='alignment_configuration', nargs=1, help='The id of the configuration to align to, otherwise random')
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "Bio", "numpy"])

    #get file names
    traj_file = args.trajectory[0]
    parallel = args.parallel
    if parallel:
        from oxDNA_analysis_tools.UTILS import parallelize_erik_onefile
        n_cpus = args.parallel[0]

    #-f defines the format of the output file
    outjson = False
    outoxdna = False
    if args.format:
        if "json" in args.format:
            outjson = True
        if "oxDNA" in args.format or "oxdna" in args.format:
            outoxdna = True
        if "both" in args.format:
            outjson = True
            outoxdna = True
        if outjson == outoxdna == False:
            print("ERROR: unrecognized output format\nAccepted formats are \"json\", \"oxDNA/oxdna\", and \"both\"", file=stderr)
            exit(1)
    else:
        print("INFO: No output format specified, defaulting to oxDNA", file=stderr)
        outoxdna = True

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else: 
        if outjson and not outoxdna:
            ext = ".json"
        elif outjson and outoxdna:
            ext = ".json/.dat"
        elif outoxdna and not outjson:
            ext = ".dat"
        outfile = "mean{}".format(ext)
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    #-d will run compute_deviations.py when this script is completed.
    dev_file = None
    if args.deviations:
        dev_file = args.deviations[0]

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


    # The reference configuration which is used to define alignment
    align_conf = []

    #calculate the number of configurations in the trajectory 
    num_confs = cal_confs(traj_file)

    # if we have no align_conf we need to chose one
    # and realign its cms to be @ 0,0,0
    if align_conf == []:
        align = None
        if args.align:
            align = args.align[0]
        align_conf_id, align_poses = pick_starting_configuration(traj_file, num_confs, align)
        # we are just interested in the nucleotide positions
        align_conf = align_poses.positions[indexes]

    #Actually compute mean structure
    if not parallel:
        print("INFO: Computing mean of {} configurations with an alignment of {} particles using 1 core.".format(num_confs, len(align_conf)), file=stderr)
        r = ErikReader(traj_file)
        mean_pos_storage, mean_a1_storage, mean_a3_storage, intermediate_mean_structures, processed_frames = compute_mean(r, align_conf, indexes, num_confs)

    #If parallel, the trajectory is split into a number of chunks equal to the number of CPUs available.
    #Each of those chunks is then calculated seperatley and the result is summed.
    if parallel:
        print("INFO: Computing mean of {} configurations with an alignment of {} particles using {} cores.".format(num_confs, len(align_conf), n_cpus), file=stderr)
        out = parallelize_erik_onefile.fire_multiprocess(traj_file, compute_mean, num_confs, n_cpus, align_conf, indexes)
        mean_pos_storage = np.sum(np.array([i[0] for i in out]), axis=0)
        mean_a1_storage = np.sum(np.array([i[1] for i in out]), axis=0)
        mean_a3_storage = np.sum(np.array([i[2] for i in out]), axis=0)
        intermediate_mean_structures = []
        [intermediate_mean_structures.extend(i[3]) for i in out]
        processed_frames = sum((i[4] for i in out))
    # finished task entry
    print("INFO: processed frames total: {}".format(processed_frames), file=stderr)

    #Convert mean structure to a json file
    mean_file = dumps({
                "i_means" : intermediate_mean_structures,
                "g_mean" : prep_pos_for_json(
                    mean_pos_storage / processed_frames
                ),
                "a1_mean" : prep_pos_for_json(
                   [normalize(v)  for v in (mean_a1_storage / processed_frames)]
                ),
                "a3_mean" : prep_pos_for_json(
                   [normalize(v)  for v in (mean_a3_storage / processed_frames)]
                ),
                "p_frames" : processed_frames,
                "ini_conf":{
                    "conf": prep_pos_for_json(align_conf),
                    "id"  : align_conf_id
                }
            })

    #Save the mean structure to the specified output file.
    if outjson or dev_file:
        #save output as json format
        if outoxdna == True:
            #if making both outputs, automatically set file extensions.
            jsonfile = outfile.split(".")[0]+".json"
        else:
            jsonfile = outfile
        print("INFO: Writing mean configuration to", jsonfile, file=stderr)
        with open(jsonfile, "w") as file:
            file.write(mean_file)
    
    if outoxdna:
        #save output as oxDNA .dat format
        if outjson == True:
            #if making both outputs, automatically set file extensions.
            outname = outfile.split(".")[0]+".dat"
        else:
            outname = outfile
        from oxDNA_analysis_tools.mean2dat import make_dat
        make_dat(loads(mean_file), outname)

    #If requested, run compute_deviations.py using the output from this script.
    if dev_file:
        print("INFO: launching compute_deviations.py", file=stderr)

        #this is probably horrible practice, but to maintain the ability to call things from the command line, I cannot pass arguments between main() calls.
        #so instead we're gonna spoof a global variable to make it look like compute_deviations was called explicitally
        argv.clear()
        argv.extend(['compute_deviations.py', '-o', dev_file, "-r", dev_file.split('.')[0]+"_rmsd.png", "-d", dev_file.split('.')[0]+"_rmsd_data.json"])
        if args.index_file:
            argv.append("-i")
            argv.append(index_file)
        if parallel:
            argv.append("-p") 
            argv.append(str(n_cpus))
        argv.append(jsonfile)
        argv.append(traj_file)

        from oxDNA_analysis_tools import compute_deviations
        from sys import executable
        print(executable)
        print(argv)

        compute_deviations.main()

        #compute_deviations needs the json meanfile, but its not useful for visualization
        #so we dump it
        if not outjson:
            print("INFO: deleting {}".format(jsonfile), file=stderr)
            from os import remove
            remove(jsonfile)

    print(time.time() - start_t)

if __name__ == '__main__':
    main()
    