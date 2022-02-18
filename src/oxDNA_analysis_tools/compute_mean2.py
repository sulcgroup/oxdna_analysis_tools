from Bio.SVDSuperimposer import SVDSuperimposer
from oxDNA_analysis_tools.UTILS.micha_reader import MichaReader
import numpy as np
from os import getenv, path
from oxDNA_analysis_tools.UTILS.micha_reader import partition, flatten, handle_confs
from oxDNA_analysis_tools.p_align import get_n_cpu
from pathos.pools import ProcessPool
from pathos.multiprocessing import cpu_count
from functools import partial
import argparse
import time

start_t = time.time()

#aligner
def align(indexes, ref_conf, mysystem):
    """
    Aligns a single frame to the reference configuration.

    Parameters:
        indexes (list): The indexes of the particles to align.
        ref_conf (base_array): The reference configuration to align to.
        mysystem (base_array): The configuration to align.

    Returns:
        str: The aligned configuration in the format of the original trajectory file.
    """
    sup = SVDSuperimposer()
    #Need to get rid of fix_diffusion artifacts or SVD doesn't work
    mysystem.inbox()
    indexed_cur_conf = mysystem.positions[indexes]

    #Superimpose the configuration to the reference
    sup.set(ref_conf.positions[indexes], indexed_cur_conf)
    sup.run()
    rot, tran = sup.get_rotran()

    #Apply rotation and translation in one step
    mysystem.positions = np.einsum('ij, ki -> kj', rot, mysystem.positions) + tran
    mysystem.a1s = np.einsum('ij, ki -> kj', rot, mysystem.a1s)
    mysystem.a3s = np.einsum('ij, ki -> kj', rot, mysystem.a3s)
    return mysystem #this is the only difference with the one in p_align

def main():

    ncpus = get_n_cpu()
    ntopart = 5
    #prepare to fire multiple processes 
    pool = ProcessPool(nodes=ncpus)

    #handle commandline arguments
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Aligns each frame in a trajectory to the first frame")
    parser.add_argument('top', type=str, nargs=1, help="The trajectory file to align")
    parser.add_argument('traj', type=str, nargs=1, help="The trajectory file to align")
    parser.add_argument('outfile', type=str, nargs=1, help='The name of the new trajectory file to write out')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Align to only a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-r', metavar='reference_structure', dest='reference_structure', nargs=1, help="Align to a provided configuration instead of the first frame.")
    args = parser.parse_args()
    print(args)

        #Parse command line arguments
    top_file = args.top[0]
    traj_file = args.traj[0]
    outfile = args.outfile[0]

    # setup the traj reader for the trajectory file we are aligning
    mr = MichaReader(top_file, traj_file,None,ncpus*ntopart)
    #-r will make it align to a provided .dat file instead of the first configuration
    if args.reference_structure:
        #read reference configuration
        frame = MichaReader(top_file, args.reference_structure[0], None, 1).read()
    else:
        #read the first configuration and use it as the reference configuration for the rest
        frame = mr.read(0)
    #bring reference configuration inbox
    frame.inbox()

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
        indexes = list(range(mr.top_info.bases))


    # compose the final function
    align_frame  = partial(align, indexes, frame) 
    handle_align = partial(handle_confs, mr.top_info,align_frame) # define what to do 

    #figure out how many chunks we are working on 
    conf_chunks = int(len(mr.idxs) / mr.buff_size) + 1

    chunk_mean = []
    chunk_a1 = []
    chunk_a3 = []

    for i, ptr_pos in  enumerate(range(0,len(mr.idxs),mr.buff_size)):
        print("Processing chunk {} of {}".format(i+1, conf_chunks))
        confs_to_parse = mr._get_confs(ptr_pos, mr.buff_size)
        confs_to_parse = partition(confs_to_parse, ntopart)
        parsed_confs =  pool.map(handle_align, confs_to_parse)
        chunk_mean.append(np.sum(np.sum([[c.positions for c in l] for l in parsed_confs], axis=0), axis=0))
        chunk_a1.append(np.sum(np.sum([[c.a1s for c in l] for l in parsed_confs], axis=0), axis=0))
        chunk_a3.append(np.sum(np.sum([[c.a3s for c in l] for l in parsed_confs], axis=0), axis=0))

    chunk_mean = np.sum(chunk_mean, axis=0)  / len(mr.idxs)
    chunk_a1 = np.sum(chunk_a1, axis=0)  / len(mr.idxs)
    chunk_a3 = np.sum(chunk_a3, axis=0)  / len(mr.idxs)

    frame.positions = chunk_mean
    frame.a1s = chunk_a1
    frame.a3s = chunk_a3

    with open('mean.dat', 'w+') as f:
        f.write(frame.conf_to_str())
        
    print(time.time() - start_t)

if __name__ == '__main__':
    main()