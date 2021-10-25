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
from oxDNA_analysis_tools.UTILS.readers import ErikReader
import numpy as np
import argparse
import os

def main():
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="superimposes one or more structures sharing a topology to a reference structure")
    parser.add_argument('reference', type=str, nargs=1, help="The reference configuration to superimpose to")
    parser.add_argument('victims', type=str, nargs='+', help="The configuraitons to superimpose on the reference")
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Align to only a subset of particles from a space-separated list in the provided file')
    args = parser.parse_args()


    #run system checks
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy", "Bio"])

    #Get the reference files
    ref_dat = args.reference[0]

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
        with ErikReader(ref_dat) as r:
            indexes = list(range(len(r.read().positions)))

    #Create list of configurations to superimpose
    to_sup = []
    r = ErikReader(ref_dat)
    ref = r.read()
    ref.inbox()
    ref_conf = ref.positions[indexes]
    for i in args.victims:
        r = ErikReader(i)
        sys = r.read()
        sys.inbox()
        to_sup.append(sys)

    sup = SVDSuperimposer()

    #Run the biopython superimposer on each configuration and rewrite its configuration file
    for i, sys in enumerate(to_sup):
        indexed_cur_conf = sys.positions[indexes]
        sup.set(ref_conf, indexed_cur_conf)
        sup.run()
        rot, tran = sup.get_rotran()
        sys.positions = np.einsum('ij, ki -> kj', rot, sys.positions) + tran
        sys.a1s = np.einsum('ij, ki -> kj', rot, sys.a1s)
        sys.a3s = np.einsum('ij, ki -> kj', rot, sys.a3s)
        sys.write_new("aligned{}.dat".format(i))
        print("INFO: Wrote file aligned{}.dat".format(i), file=stderr)

if __name__ == '__main__':
    main()