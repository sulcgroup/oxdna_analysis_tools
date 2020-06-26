#!/usr/bin/env python3

#superimpose.py
#Created by: Erik Poppleton
#Date: 2/27/19
#Takes two (or more) configurations and aligns all proceeding ones to the first using the svd superimposer, then spits them out as new dat files.

try:
    from Bio.SVDSuperimposer import SVDSuperimposer
except:
    from bio.SVDSuperimposer import SVDSuperimposer
from sys import exit, stderr
from UTILS.readers import LorenzoReader2
import numpy as np
import argparse
from compute_mean import normalize

fetch_np = lambda conf: np.array([
    n.cm_pos for n in conf._nucleotides 
])

parser = argparse.ArgumentParser(description="superimposes one or more structures sharing a topology to a reference structure")
parser.add_argument('topology', type=str, nargs=1, help="The topology file shared by all the configurations to superimpose")
parser.add_argument('reference', type=str, nargs=1, help="The reference configuration to superimpose to")
parser.add_argument('victims', type=str, nargs='+', help="The configuraitons to superimpose on the reference")
args = parser.parse_args()

#Get the reference files
top_file = args.topology[0]
ref_dat = args.reference[0]

#Create list of configurations to superimpose
to_sup = []
r = LorenzoReader2(ref_dat, top_file)
ref = r._get_system()
ref.inbox()
ref_conf = fetch_np(ref)
for i in args.victims:
    r = LorenzoReader2(i, top_file)
    sys = r._get_system()
    to_sup.append(sys)

sup = SVDSuperimposer()

#Run the biopython superimposer on each configuration and rewrite its configuration file
for i, sys in enumerate(to_sup):
    cur_conf = fetch_np(sys)
    sup.set(ref_conf, cur_conf)
    sup.run()
    rot, tran = sup.get_rotran()
    cur_conf = np.einsum('ij, ki -> kj', rot, cur_conf) + tran
    for j,n in enumerate(sys._nucleotides):
        n.cm_pos = cur_conf[j]
        n._a1 = normalize(np.dot(n._a1, rot))
        n._a3 = normalize(np.dot(n._a3, rot))
    sys.print_lorenzo_output("aligned{}.dat".format(i), "/dev/null")
    print("INFO: Wrote file aligned{}.dat".format(i), file=stderr)