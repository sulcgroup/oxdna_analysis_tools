#!/usr/bin/env python3

#A utility that prints out H-bond occupancy compared with the design
#Created by: Erik Poppleton
#Date: 11/1/18
#Python3
#Takes a trajectory and a pairs file (indicating the desired design) and prints out in how many configurations the 
#designed pairs are present in.  Also prints out the list of missbonds found in the trajectory and their frequency.


import numpy as np
from sys import exit, stderr
import argparse
from os import environ, path
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2, cal_confs, get_input_parameter
from oxDNA_analysis_tools.UTILS import parallelize_lorenzo_onefile
from oxDNA_analysis_tools.output_bonds import output_bonds

def bond_analysis(reader, pairs, inputfile, num_confs, start=None, stop=None):
    """
    Compares a list of desired base pairs with the base pairs present in a trajectory.

    Parameters:
        reader (LorenzoReader2): A reader attached to the trajectory to analyze.
        pairs (list): A list of the designed pairs.  Format is [[a, a'], [b, b']]
        num_confs (int): The total number of configurations in the trajectory.
        <optional> start (int): The starting configuration ID to begin averaging at.  Used if parallel.
        <optional> stop (int): The configuration ID on which to end the averaging.  Used if parallel.

    Returns:
        tot_bonds (float): The average number of bonds found in the trajectory.
        tot_missbonds (float): The average number of missbonded nucleotides in the trajectory.
        out_array (numpy.array): The average occupancy of each designed bond in the trajectory.
        confid (int): The number of configurations analyzed by this function.
    """

    #standard parallelization header
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

    #for every configuration in the trajectory, compare the hydrogen bonds to the design
    while mysystem != False and confid < stop:
        mysystem.map_nucleotides_to_strands()
        out = output_bonds(inputfile, mysystem)  
        if out == '':
            return None, None, [], None #return something to indicate that the configuration is not valid
        mysystem.read_H_bonds_output_bonds(out) #maps the output of output_bonds to the system object
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
        print("Frame:", confid, "Time:", mysystem._time)
        print("formed", count_correct_bonds, "out of", len(pairs))
        print("missbonds:", count_incorrect_bonds)
        confid += 1
        mysystem = reader._get_system()
    
    return (tot_bonds, tot_missbonds, out_array, confid)

def main():
    #read data from files
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Compare the bonds found at each trajectory with the intended design")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help="The trajecotry file to compare against the designed pairs")
    parser.add_argument('designed_pairs', type=str, nargs=1, help="The file containing the desired nucleotides pairings in the format \n a b\nc d")
    parser.add_argument('output_file', type=str, nargs=1, help="name of the file to save the output json overlay to")
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    
    #run system checks
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    args = parser.parse_args()
    inputfile = args.inputfile[0]
    traj_file = args.trajectory[0]
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
    num_confs = cal_confs(traj_file)

    with open(designfile, 'r') as file:
        pairs = file.readlines()

    if not parallel:
        print("INFO: Computing base pairs in {} configurations using 1 core.".format(num_confs), file=stderr)
        r = LorenzoReader2(traj_file,top_file)
        tot_bonds, tot_missbonds, out_array, confid = bond_analysis(r, pairs, inputfile, num_confs)
        try:
            _ = tot_bonds #this will fail if DNAnalysis failed.
        except:
            print("ERROR: DNAnalysis encountered an error and could not analyze the trajectory")
            exit(1)

    if parallel:
        print("INFO: Computing base pairs in {} configurations using {} cores.".format(num_confs, n_cpus), file=stderr)
        out = parallelize_lorenzo_onefile.fire_multiprocess(traj_file, top_file, bond_analysis, num_confs, n_cpus, pairs, inputfile)

        tot_bonds = 0
        tot_missbonds = 0
        out_array = np.zeros(len(open(top_file, 'r').readlines())-1)
        confid = 0
        for i in out:
            if i[0] is not None:
                tot_bonds += i[0]
                tot_missbonds += i[1]
                #out_array += i[2]
                confid += i[3]
            else:
                print("WARNING: Some configurations were invalid and not included in the analysis.  Please check the logs", file=stderr)

            #tot_bonds = sum((i[0] for i in out if i[0] != None))
            #tot_missbonds = sum((i[1] for i in out if i[1] != None))
        out_array = sum((i[2] for i in out if len(i[2]) > 0))
            #confid = sum((i[3] for i in out if i[3] != None))

    print("\nSummary:\navg bonds: {}\navg_missbonds: {}".format(tot_bonds/(int(confid)),tot_missbonds/int(confid)))

    print("INFO: Writing bond occupancy data to {}".format(outfile))
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

if __name__ == '__main__':
    main()
    
