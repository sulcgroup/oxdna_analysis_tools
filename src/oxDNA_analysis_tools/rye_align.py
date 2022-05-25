import argparse
import os
import time
import numpy as np
from sys import stderr
from multiprocessing import Pool
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import describe, inbox, conf_to_str
from oxDNA_analysis_tools.UTILS.get_confs import get_confs
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "centered_ref_coords",
                                              "indexes"])

def align(centered_ref_coords, coords, indexes):
    """
    Single-value decomposition-based alignment of configurations

    Parameters
        centered_ref_coords: reference coordinates, centered on [0, 0, 0]
        coords: coordinates to be aligned
        indexes: indexes of the atoms to be aligned

    Returns
        A tuple of the aligned coordinates (coords, a1s, a3s) for the given chunk
    """
    # center on centroid
    av1, reference_coords = np.zeros(3), centered_ref_coords.copy()
    av2 = np.mean(coords[0][indexes], axis=0)
    coords[0] = coords[0] - av2
    # correlation matrix
    a = np.dot(np.transpose(coords[0][indexes]), reference_coords)
    u, _, vt = np.linalg.svd(a)
    rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    # check if we have found a reflection
    if np.linalg.det(rot) < 0:
        vt[2] = -vt[2]
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    tran = av1 - np.dot(av2, rot)
    return  (np.dot(coords[0], rot) + tran,
             np.dot(coords[1], rot),
             np.dot(coords[2], rot))

def compute(ctx:ComputeContext, chunk_size, chunk_id:int):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*chunk_size, chunk_size, ctx.top_info.nbases)
    confs = [inbox(c, center=True) for c in confs]
    # convert to numpy repr
    np_coords = np.asarray([[c.positions, c.a1s, c.a3s] for c in confs])

    # align
    for i, c in enumerate(np_coords):
        c[0], c[1], c[2] = align(ctx.centered_ref_coords, c, ctx.indexes)
        confs[i].positions = c[0]
        confs[i].a1s = c[1]
        confs[i].a3s = c[2]
    #return confs
    out = ''.join([conf_to_str(c) for c in confs])
    return out

def main():

    #handle commandline arguments
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Aligns each frame in a trajectory to the first frame")
    parser.add_argument('top', type=str, nargs=1, help="The trajectory file to align")
    parser.add_argument('traj', type=str, nargs=1, help="The trajectory file to align")
    parser.add_argument('outfile', type=str, nargs=1, help='The name of the new trajectory file to write out')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Align to only a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-r', metavar='reference_structure', dest='reference_structure', nargs=1, help="Align to a provided configuration instead of the first frame.")
    args = parser.parse_args()

    #Parse command line arguments
    top_file = args.top[0]
    traj_file = args.traj[0]
    outfile = args.outfile[0]
    top_info, traj_info = describe(top_file, traj_file)

    #-r will make it align to a provided .dat file instead of the first configuration
    if args.reference_structure:
        #read reference configuration
        _, ref_info = describe(top_file, args.reference_structure[0])
        ref_conf = get_confs(ref_info.idxs, ref_info.path, 0, 1, top_info.nbases)[0]
    else:
        #read the first configuration and use it as the reference configuration for the rest
        ref_conf = get_confs(traj_info.idxs, traj_info.path, 1, 1, top_info.nbases)[0]

    ref_conf = inbox(ref_conf) # Don't need to center now because we're going to after indexing anyway.

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
        indexes = list(range(top_info.nbases))

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    # alignment requires the ref to be centered at 0.  Inboxing did not take the indexing into account.
    reference_coords = ref_conf.positions[indexes]
    ref_cms = np.mean(reference_coords, axis=0) # cms prior to centering
    reference_coords = reference_coords - ref_cms

    # Create a ComputeContext which defines the problem to pass to the worker processes 
    ctx = ComputeContext(
        traj_info, top_info, reference_coords, indexes
    )

    with open(outfile, 'w+') as f:
        def callback(r):
            nonlocal f
            f.write(r)

        oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    print(f"INFO: Wrote aligned trajectory to {outfile}", file=stderr)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()