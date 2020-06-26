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
from UTILS.readers import LorenzoReader2
import numpy as np
import argparse
from compute_mean import normalize

#helper function to retrieve positions
fetch_np = lambda conf: np.array([
    n.cm_pos for n in conf._nucleotides 
])

indexed_fetch_np = lambda conf: np.array([
    n.cm_pos for n in conf._nucleotides if n.index in indexes
])

#handle commandline arguments
parser = argparse.ArgumentParser(description="Aligns each frame in a trajectory to the first frame")
parser.add_argument('traj', type=str, nargs=1, help="The trajectory file to align")
parser.add_argument('topology', type=str, nargs=1, help="The topology file corresponding to the trajectory")
parser.add_argument('outfile', type=str, nargs=1, help='The name of the new trajectory file to write out')
parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Align to only a subset of particles from a space-separated list in the provided file')
args = parser.parse_args()

#run system checks
from config import check_dependencies
check_dependencies(["python", "numpy", "Bio"])

#Parse command line arguments
top_file = args.topology[0]
traj_file = args.traj[0]
outfile = args.outfile[0]

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

#read the first configuration and use it as the reference configuration for the rest
r = LorenzoReader2(traj_file, top_file)
ref = r._get_system()
ref.inbox()
ref_conf = indexed_fetch_np(ref)
sup = SVDSuperimposer()

#The topology remains the same so we only write the configuration
ref.print_lorenzo_output(outfile, '/dev/null')
mysystem = r._get_system()

#Read the trajectory one configuration at a time and perform the alignment
while mysystem != False:
    print("working on t = ", mysystem._time)
    #Need to get rid of fix_diffusion artifacts or SVD doesn't work
    mysystem.inbox()
    cur_conf = fetch_np(mysystem)
    indexed_cur_conf = indexed_fetch_np(mysystem)

    #Superimpose the configuration to the reference
    sup.set(ref_conf, indexed_cur_conf)
    sup.run()
    rot, tran = sup.get_rotran()

    #Apply rotation and translation in one step
    cur_conf = np.einsum('ij, ki -> kj', rot, cur_conf) + tran
    
    #Overwrite positions and orientation
    for j,n in enumerate(mysystem._nucleotides):
        n.cm_pos = cur_conf[j]
        n._a1 = normalize(np.dot(n._a1, rot))
        n._a3 = normalize(np.dot(n._a3, rot))

    #print_lorenzo_output will create a new file, print_traj_output extends an existing file
    mysystem.print_traj_output(outfile, '/dev/null')

    mysystem = r._get_system()

