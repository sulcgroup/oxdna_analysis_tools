import argparse
from sys import exit, stderr
from os import path
from collections import namedtuple
from typing import Tuple
import numpy as np
import oxpy
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_input_parameter

import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "designed_pairs",
                                              "input_file"])

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["analysis_bytes_to_skip"] = str(ctx.traj_info.idxs[chunk_id*chunk_size].offset)
        inp["confs_to_analyse"] = str(chunk_size)
        inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = hb_list \n } \n }'

        if (not inp["use_average_seq"] or inp.get_bool("use_average_seq")) and "RNA" in inp["interaction_type"]:
            print("WARNING: Sequence dependence not set for RNA model, wobble base pairs will be ignored", file=stderr)

        backend = oxpy.analysis.AnalysisBackend(inp)
    
        count_correct_bonds = 0
        count_incorrect_bonds = 0
        tot_bonds = 0
        out_array = np.zeros(ctx.top_info.nbases, dtype=int)
        while backend.read_next_configuration():
            pairs = backend.config_info().get_observable_by_id("my_obs").get_output_string(0).strip().split('\n')
            for p in pairs[1:]:
                p = p.split()
                a = int(p[0])
                b = int(p[1])
                if a in ctx.designed_pairs.keys():
                    if ctx.designed_pairs[a] == b:
                        count_correct_bonds += 1
                        tot_bonds += 1
                        out_array[a] += 1
                        out_array[b] += 1
                    else:
                        count_incorrect_bonds += 1
                        tot_bonds += 1

        return(tot_bonds, count_correct_bonds, count_incorrect_bonds, out_array)


def bond_analysis(traj_info:TrajInfo, top_info:TopInfo, pairs:dict[int, int], inputfile:str, ncpus:int=1) -> Tuple[int, int, int, np.ndarray]:
    '''
        Compare the bond occupancy of a trajectory with a designed structure

        Parameters: 
            traj_info: TrajInfo object containing the trajectory information
            top_info: TopInfo object containing the topology information
            pairs: dict of the designed pairs
            inputfile: the input file used to run the simulation
            ncpus: (optional) number of cores to use

        Returns:
            total_bonds: average number of bonds per configuraiton
            incorrect_bonds: average number of incorrect bonds per configuration
            correct_bonds: average number of correct bonds per configuration
            nt_array: per-nucleotide correct bond occupancy
    '''
    ctx = ComputeContext(traj_info, top_info, pairs, inputfile)

    total_bonds = 0
    correct_bonds = 0
    incorrect_bonds = 0
    nt_array = np.zeros(ctx.top_info.nbases, dtype=int)

    def callback(i, r):
        nonlocal total_bonds, correct_bonds, incorrect_bonds, nt_array
        total_bonds += r[0]
        correct_bonds += r[1]
        incorrect_bonds += r[2]
        nt_array += r[3]

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    total_bonds /= traj_info.nconfs
    correct_bonds /= traj_info.nconfs
    incorrect_bonds /= traj_info.nconfs
    nt_array  = nt_array / traj_info.nconfs

    return(total_bonds, correct_bonds, incorrect_bonds, nt_array)

def main():
    #read data from files
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Compare the bonds found at each trajectory with the intended design")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help="The trajecotry file to compare against the designed pairs")
    parser.add_argument('designed_pairs', type=str, nargs=1, help="The file containing the desired nucleotides pairings in the format \n a b\nc d")
    parser.add_argument('output_file', type=str, nargs=1, help="name of the file to save the output json overlay to")
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    args = parser.parse_args()

    #run system checks
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    inputfile = args.inputfile[0]
    traj_file = args.trajectory[0]
    designfile = args.designed_pairs[0]
    outfile = args.output_file[0]

    top_file = get_input_parameter(inputfile, "topology")
    top_info, traj_info = describe(top_file, traj_file)

    with open(designfile, 'r') as file:
        pairs_txt = file.readlines()

    pairs = {int(p[0]) : int(p[1]) for p in [p.split() for p in pairs_txt]}

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    total_bonds, incorrect_bonds, correct_bonds, nt_array = bond_analysis(traj_info, top_info, pairs, inputfile, ncpus)

    print("\nSummary:\navg bonds: {}\navg_missbonds: {}".format(total_bonds, incorrect_bonds))

    print("INFO: Writing bond occupancy data to {}".format(outfile))

    with open(outfile, "w+") as file:
        file.write("{\n\"occupancy\" : [")
        file.write(str(nt_array[0]))
        for n in nt_array[1:]:
            file.write(", {}".format(n))
        file.write("] \n}")

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()