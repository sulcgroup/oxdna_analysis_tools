import argparse
import os
import time
import numpy as np
from sys import stderr
from collections import namedtuple
from random import randrange
from oxDNA_analysis_tools.align import align
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, inbox, write_conf
from oxDNA_analysis_tools.UTILS.data_structures import Configuration
start_time = time.time()

# A compute context is a single variable (in this case a namedtuple, but could also be a dict or a dataclass)
# which contains the arguments needed for a parallelized function.
# We use this style to avoid having to pass variable numbers of arguments to the parallelizer.
ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "centered_ref_coords",
                                              "indexes"])

# This is the function which computes the sum of particle positions
# This function will be parallelized.                                              
def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*chunk_size, chunk_size, ctx.top_info.nbases)
    confs = (inbox(c, center=True) for c in confs)
    # convert to numpy repr
    np_coords = np.asarray([[c.positions, c.a1s, c.a3s] for c in confs])
    sub_mean = np.zeros(shape=[3,ctx.top_info.nbases,3])
    for c in np_coords:
        sub_mean += align(ctx.centered_ref_coords, c, ctx.indexes)
    
    return sub_mean

# All scripts in oat must have a main with no arguments to work with the command line interface.
def main():

    # A standard way to create and parse command line arguments.
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Computes the mean structure of a trajectory file")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the mean structure to')
    parser.add_argument('-d', '--deviations', metavar='deviation_file', nargs=1, help='Immediatley run compute_deviations.py from the output')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Compute mean structure of a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-a', '--align', metavar='alignment_configuration', nargs=1, help='The id of the configuration to align to, otherwise random')
    args = parser.parse_args()

    # Verify that dependencies are installed and a good version
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    # Get metadata about input files
    traj = args.trajectory[0]
    top_info, traj_info = describe(None, traj)


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

    # -p sets the number of cores to use.  Default is 1.
    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    # alignment requires the ref to be centered at 0
    reference_coords = ref_conf.positions[indexes]
    ref_cms = np.mean(reference_coords, axis=0) # cms prior to centering
    reference_coords = reference_coords - ref_cms

    # The compute context is a single variable containin the arguments needed for the parallelized function. 
    ctx = ComputeContext(
        traj_info, top_info, reference_coords, indexes
    )

    # The callback function is called after each chunk is computed.
    # For computing the mean structure, we just add the sum of particle positions to a running total
    # we'll divide by the number of confs later.
    mean = np.zeros([3, top_info.nbases, 3])
    def callback(i, r):
        nonlocal mean
        mean += r

    # The parallelizer will call the "compute" function with ctx as an argument using "ncpus".
    # The "callback" function is called after each chunk is computed.
    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    # divide the sum by the number of confs to get the mean 
    mean /= traj_info.nconfs
    pos, a1s, a3s = mean

    # renormalize a1 and a3 because weird things happen if they aren't
    a1s = np.array([v/np.linalg.norm(v) for v in a1s])
    a3s = np.array([v/np.linalg.norm(v) for v in a3s])

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else:
        outfile = "mean.dat"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    # Create the mean configuration from the numpy arrays containing the positions and orientations
    # And write it to the outfile
    write_conf(outfile,Configuration(0,ref_conf.box,np.array([0,0,0]), pos, a1s , a3s))
    print("--- %s seconds ---" % (time.time() - start_time))

    # -d runs compute_deviations.py after computing the mean
    if args.deviations:
        from sys import argv
        from oxDNA_analysis_tools import deviations
        dev_file = args.deviations[0]
        print("INFO: Launching compute_deviations")

        #this is probably horrible practice, but to maintain the ability to call things from the command line, I cannot pass arguments between main() calls.
        #so instead we're gonna spoof argv, which is in the global scope, to make it look like deviations was called explicitally
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

        # Call the main method from deviations.
        # It will see the spoofed argv as if it was called from the command line.
        deviations.main()

if __name__ == '__main__':
    main()