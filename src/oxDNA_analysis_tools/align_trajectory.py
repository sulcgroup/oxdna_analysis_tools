#!/usr/bin/env python3

# align_trajectory.py
# Created by: Erik Poppleton
# Date: 8/26/19
# Takes a trajectory and aligns every frame to the first one and writes a new trajectory

try:
    from Bio.SVDSuperimposer import SVDSuperimposer
except:
    from bio.SVDSuperimposer import SVDSuperimposer
from sys import exit
from oxDNA_analysis_tools.UTILS.readers import ErikReader
import numpy as np
import argparse
import os

#aligner
def align_frame(ref_conf, sup, mysystem, indexes):
    """
    Aligns a frame to the reference configuration

    Parameters
        ref_conf (base_array): The reference configuration
        sup (SVDSuperimposer): The superimposer object
        mysystem (base_array): The configuration to align
        indexes (list): The indexes of the particles to align

    Returns
        mysystem (base_array): The aligned configuration
    """
    #Need to get rid of fix_diffusion artifacts or SVD doesn't work
    mysystem.inbox()
    indexed_cur_conf = mysystem.positions[indexes]

    #Superimpose the configuration to the reference
    sup.set(ref_conf, indexed_cur_conf)
    sup.run()
    rot, tran = sup.get_rotran()

    #Apply rotation and translation in one step
    mysystem.positions = np.einsum('ij, ki -> kj', rot, mysystem.positions) + tran
    mysystem.a1s = np.einsum('ij, ki -> kj', rot, mysystem.a1s)
    mysystem.a3s = np.einsum('ij, ki -> kj', rot, mysystem.a3s)

    return mysystem

def main():
    #handle commandline arguments
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Aligns each frame in a trajectory to the first frame")
    parser.add_argument('traj', type=str, nargs=1, help="The trajectory file to align")
    parser.add_argument('outfile', type=str, nargs=1, help='The name of the new trajectory file to write out')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Align to only a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-r', metavar='reference_structure', dest='reference_structure', nargs=1, help="Align to a provided configuration instead of the first frame.")
    args = parser.parse_args()
    
    #run system checks
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy", "Bio"])
    
    #Parse command line arguments
    traj_file = args.traj[0]
    outfile = args.outfile[0]
    sup = SVDSuperimposer()
    
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
    
    #-r will make it align to a provided .dat file instead of the first configuration
    if args.reference_structure:
        #read reference configuration
        r = ErikReader(args.reference_structure[0])
        ref = r.read()
        ref.inbox()
        r = ErikReader(traj_file)
        ref_conf = ref.positions[indexes]
        
        mysystem = align_frame(ref_conf, sup, r.read())
    
    else:
        #read the first configuration and use it as the reference configuration for the rest
        r = ErikReader(traj_file)
        mysystem = r.read()
        mysystem.inbox()
        ref_conf = mysystem.positions[indexes]
        
    #write first configuration to output file
    mysystem.write_new(outfile)
    mysystem = r.read()
    
    
    
    #Read the trajectory one configuration at a time and perform the alignment
    while mysystem != False:
        print("working on t = ", mysystem.time)
        
        mysystem = align_frame(ref_conf, sup, mysystem, indexes)
    
        mysystem.write_append(outfile)
    
        mysystem = r.read()

if __name__ == '__main__':
    main()
