#!/usr/bin/env python

#A utility that prints out H-bond occupancy compared with the design
#Created by: Erik Poppleton
#Date: 11/1/18
#Python3
#Takes a trajectory and a pairs file (indicating the desired design) and prints out in how many configurations the 
#designed pairs are present in.  Also prints out the list of missbonds found in the trajectory and their frequency.

from UTILS.readers import LorenzoReader2, Cal_confs, get_input_parameter
import numpy as np
import sys
from output_bonds import output_bonds
import argparse
from UTILS import parallelize
from os import environ

def bond_analysis(reader, pairs, num_confs, start=None, stop=None):
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)
    
    mysystem = reader._get_system(N_skip = start)
    num_nuc = len(mysystem._nucleotides)

    confid = 0
    tot_bonds = 0.
    tot_missbonds = 0.
    #bond_fraction = np.zeros(len(pairs))
    out_array = np.zeros(num_nuc)
    #missbonds = {}
    while mysystem != False and confid < stop:
        print("\nconf", mysystem._time)
        mysystem.map_nucleotides_to_strands()
        out = output_bonds(inputfile, mysystem)  
        mysystem.read_H_bonds_output_bonds(out)
        #count_bonds = 0
        count_correct_bonds = 0
        count_incorrect_bonds = 0
        for i in range(len(pairs)):
            line = pairs[i]
            a = int(line.split()[0])
            b = int(line.split()[1])
            if b in mysystem._nucleotides[a].interactions:
                #print(a, "and", b, "are correctly bonded")
                count_correct_bonds += 1
                tot_bonds += 1
                #bond_fraction[i] += 1
                out_array[a] += 1
                out_array[b] += 1
            if b not in mysystem._nucleotides[a].interactions and a < b:
                #print(a, "to", b, "is missing")
                if mysystem._nucleotides[a].interactions != []:
                    #m = (a, mysystem._nucleotides[a].interactions[0])
                    #print(a, "(", mysystem._nucleotides[a].get_base(), ") is bonded to", mysystem._nucleotides[a].interactions[0], "(", mysystem._nucleotides[b].get_base(), ")")
                    count_incorrect_bonds += 1
                    #tot_bonds += 1
                    tot_missbonds += 1
                #if m not in missbonds.keys():
                #    missbonds[m] = 1
                #else:
                #    missbonds[m] += 1                   

        print("formed", count_correct_bonds, "out of", len(pairs))
        print("missbonds:", count_incorrect_bonds)
        confid += 1
        mysystem = reader._get_system()
    
    return (tot_bonds, tot_missbonds, out_array, confid)

if __name__ == "__main__":
    #read data from files
    parser = argparse.ArgumentParser(description="Calculate and display the contact map for a structure")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help="The file containing the configurations of which the contact map is needed")
    parser.add_argument('designed_pairs', type=str, nargs=1, help="The file containing the desired nucleotides pairings in the format \n a b\nc d")
    parser.add_argument('output_file', type=str, nargs=1, help="name of the file to save the output json to")
    parser.add_argument('-p', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    
    args = parser.parse_args()
    inputfile = args.inputfile[0]
    dat_file = args.trajectory[0]
    designfile = args.designed_pairs[0]
    outfile = args.output_file[0]
    parallel = args.parallel
    if parallel:
        n_cpus = args.parallel[0]

    top_file = get_input_parameter(inputfile, "topology")
    if "RNA" in get_input_parameter(inputfile, "interaction_type"):
        environ["OXRNA"] = "1"
    else:
        environ["OXRNA"] = "0"
    num_confs = Cal_confs(dat_file, top_file)

    import UTILS.base #this needs to be imported after the model type is set


    with open(designfile, 'r') as file:
        pairs = file.readlines()

    if not parallel:
        r = LorenzoReader2(dat_file,top_file)
        tot_bonds, tot_missbonds, out_array, confid = bond_analysis(r, pairs, num_confs)

    if parallel:
        out = parallelize.fire_multiprocess(dat_file, top_file, bond_analysis, num_confs, n_cpus, pairs)
        tot_bonds = sum((i[0] for i in out))
        tot_missbonds = sum((i[1] for i in out))
        out_array = sum((i[2] for i in out))
        confid = sum((i[3] for i in out))

    print("\nSummary:\navg bonds: {}\navg_missbonds: {}".format(tot_bonds/(int(confid)),tot_missbonds/int(confid)))

    with open(outfile, "w+") as file:
        file.write("{\n\"occupancy\" : [")
        file.write(str(out_array[0]/int(confid)))
        for n in out_array[1:]:
            file.write(", {}".format(n/int(confid)))
        file.write("] \n}") 
    # with open(outfile, "w+") as file:
    #     file.write("Summary:\navg bonds: {0}\nbond\tfrequency\n".format(tot_bonds/(int(confid))))
    #     for j,line in enumerate(pairs):
    #         #if bond_fraction[j] != int(confid):
    #         file.write("({0}, {1})\t{2}\n".format(line.split()[0], line.split()[1], bond_fraction[j]/(int(confid)))) 

#     file.write("missbond\tfrequency\n")
#     for pair in missbonds:
#         file.write("{0}\t{1}\n".format(pair, missbonds[pair]/float(confid)))
