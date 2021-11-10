#!/usr/bin/env python

#Created by: Erik Poppleton
#Date: before 6/21/18
#Python3
#Calculates a vector corresponding to every duplex in a structure.  
#The output from this file is used by duplex_angle_finder.py to make plots of structure flexibility.

from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2, cal_confs, get_input_parameter
import numpy as np
import argparse
from os import environ, path
from sys import stderr
from oxDNA_analysis_tools.UTILS import geom
from oxDNA_analysis_tools.output_bonds import output_bonds
from oxDNA_analysis_tools.UTILS import parallelize_lorenzo_onefile

class Duplex:
    """
        Defines a nucleic acid duplex structure

        Init Parameters:
            index (int): Unique identifier for the duplex
            start_s1 (int): Particle ID of the first particle in the first strand
            end_s1 (int): Particle ID of the last particle in the first strand
            start_s2 (int): Particle ID of the first particle in the complementary strand
            end_s2 (int): Particle ID of the last particle in the complementary strand

        Member Functions:
            get_len()
    """
    def __init__(self, index, start_s1, end_s1, start_s2, end_s2):
        self.index = index
        self.start1 = start_s1
        self.end1 = end_s1
        self.start2 = start_s2
        self.end2 = end_s2

    def get_len(self):
        """
            Returns the length of the duplex
        """
        return self.end1-self.start1

def find_duplex(system):
    """
        Searches a system that has had its hydrogen bonds mapped for duplexes

        Parameters:
            system (base.System): The system to search for duplexes

        Returns:
            duplex_list (list): A list of Duplex objects found in the system
    """
    #create a list of base-paired bases
    paired_bases = []
    for s in system._strands:
        for n in s._nucleotides:
            #print(n.index, n.interactions)
            if len(n.interactions) == 1:
                paired_bases.append(np.array([int(n.index), int(n.interactions[0])]))

    #create a list of duplexes                                             
    duplex_index = 0
    search_start = 0
    duplex_list = []
    while search_start < len(paired_bases)-1:
        start1 = paired_bases[search_start][0]
        end2 = paired_bases[search_start][1]
        if end2 > start1:
            duplex_length = 1
            start_address = search_start
            while start1 + duplex_length == paired_bases[start_address+duplex_length][0] and end2 - duplex_length == paired_bases[start_address+duplex_length][1] and end2-duplex_length > start1+duplex_length:
                duplex_length+=1
                search_start+=1

            if duplex_length > 3:
                duplex_list.append(Duplex(duplex_index, start1,  paired_bases[start_address+duplex_length-1][0], paired_bases[start_address+duplex_length-1][1], end2))
                duplex_index+=1
            search_start+=1
 
        else:
            search_start+=1       

    return duplex_list

def find_angles(reader, inputfile, num_confs, start=None, stop=None):
    """
        Fits a vector to every duplex in each snapshot of a trajectory

        Parameters:
            reader (readers.LorenzoReader2): An active reader on the trajectory file in which to find duplexes.
            num_confs (int): The number of configurations in the reader.  
            <optional> start (int): The starting configuration ID to begin averaging at.  Used if parallel.
            <optional> stop (int): The configuration ID on which to end the averaging.  Used if parallel.

        Returns:
            duplexes_at_step (list): A list containing information about every duplex found in the trajectory.
    """
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)
    confid = 0
    mysystem = reader._get_system(N_skip = start)
                                                                       
    #launches DNAnalysis and prints out a list of all bonds
    duplexes_at_step = []
    while mysystem != False and confid < stop:
        print("Working on  t = {}".format(mysystem._time))
        mysystem.map_nucleotides_to_strands()
        out = output_bonds(inputfile, mysystem)
        if out == '':
            duplex_list = []
            duplexes_at_step.append(duplex_list)
            confid += 1 
            mysystem = reader._get_system()
            continue
        mysystem.read_H_bonds_output_bonds(out)
        duplex_list = find_duplex(mysystem)

        #call geom.get_axis on each duplex
        for i in range(0, len(duplex_list)):
            #print("start1: ", duplex_list[i].start1, "end2: ", duplex_list[i].end2, "end1: ", duplex_list[i].end1, "start2: ", duplex_list[i].start2)
            if environ.get("OXRNA") == "1":
                duplex_list[i].axis, duplex_list[i].final_hel_pos = geom.get_RNA_axis(mysystem, duplex_list[i].start1, duplex_list[i].end1, duplex_list[i].end2, duplex_list[i].start2, only_plane_vector = True)
            else:
                duplex_list[i].axis, duplex_list[i].final_hel_pos = geom.get_DNA_axis(mysystem, duplex_list[i].start1, duplex_list[i].end1, duplex_list[i].end2, duplex_list[i].start2, only_plane_vector = True)
            #print(duplex_list[i].start1, duplex_list[i].end1, duplex_list[i].end2, duplex_list[i].start2)
            duplex_list[i].time = mysystem._time
        duplexes_at_step.append(duplex_list)
        confid += 1 
        mysystem = reader._get_system()

        #for line in duplex_list:
        #    print("[{}, {}, {}, {}, {}, {}],".format(line.final_hel_pos[0], line.final_hel_pos[1], line.final_hel_pos[2], line.axis[0], line.axis[1], line.axis[2]))
    
    return(duplexes_at_step)

def main():
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Fit vectors to every duplex in the structure")
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help="The trajectory file from the simulation")
    parser.add_argument('-o', '--output', metavar='output_file',  type=str, nargs=1, help='name of the file to write the angle list to')
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    #Process command line arguments:
    inputfile = args.inputfile[0]
    traj_file = args.trajectory[0]
    parallel = args.parallel
    if parallel:
        n_cpus = args.parallel[0]

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else: 
        outfile = "angles.txt"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)
        

    #Get relevant parameters from the input file
    top_file = get_input_parameter(inputfile, "topology")
    if "RNA" in get_input_parameter(inputfile, "interaction_type"):
        environ["OXRNA"] = "1"
    else:
        environ["OXRNA"] = "0"

    #Calculate the number of configurations.
    num_confs = cal_confs(traj_file)

    r0 = LorenzoReader2(traj_file, top_file)
    r0._get_system()
    
    #launch find_angle using the appropriate number of threads to find all duplexes.
    if not parallel:
        print("INFO: Fitting duplexes to {} configurations using 1 core.".format(num_confs), file=stderr)
        r = LorenzoReader2(traj_file,top_file)
        duplexes_at_step = find_angles(r, inputfile, num_confs)

    if parallel:
        print("INFO: Fitting duplexes to {} configurations using {} cores.".format(num_confs, n_cpus), file=stderr)
        duplexes_at_step = []
        out = parallelize_lorenzo_onefile.fire_multiprocess(traj_file, top_file, find_angles, num_confs, n_cpus, inputfile)
        [duplexes_at_step.extend(i) for i in out]
    
    if [] in duplexes_at_step:
        print("WARNING: Some configurations were invalid and not included in the analysis.  Please check the log to view the error", file=stderr)

    #print duplexes to a file
    print("INFO: Writing duplex data to {}.  Use duplex_angle_plotter to graph data".format(outfile), file=stderr)
    output = open(outfile, 'w')
    output.write("time\tduplex\tstart1\tend1\tstart2\tend2\taxisX\taxisY\taxisZ\thel_pos\n")
    for i in range (0, len(duplexes_at_step)):
        for j in range(0, len(duplexes_at_step[i])):
            line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t[{},{},{}]\n'.format(duplexes_at_step[i][j].time,duplexes_at_step[i][j].index,duplexes_at_step[i][j].start1,duplexes_at_step[i][j].end1,duplexes_at_step[i][j].start2,duplexes_at_step[i][j].end2,duplexes_at_step[i][j].axis[0],duplexes_at_step[i][j].axis[1],duplexes_at_step[i][j].axis[2],duplexes_at_step[i][j].final_hel_pos[0],duplexes_at_step[i][j].final_hel_pos[1],duplexes_at_step[i][j].final_hel_pos[2])
            output.write(line)
    output.close()

if __name__ == '__main__':
    main()
    
