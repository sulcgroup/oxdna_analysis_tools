import numpy as np
from os import path, getcwd
from sys import stderr
from multiprocessing import Pool
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.RyeReader import no_top_describe
import oxpy

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "input_file",
                                              "visualize",
                                              "conversion_factor",
                                              "ntopart"])

def compute(ctx, chunk_id):
    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["analysis_bytes_to_skip"] = str(ctx.traj_info.idxs[chunk_id*ctx.ntopart].offset)
        inp["confs_to_analyse"] = str(ctx.ntopart)
        inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = pair_energy \n } \n }'

        if (not inp["use_average_seq"] or inp.get_bool("use_average_seq")) and "RNA" in inp["interaction_type"]:
            print("WARNING: Sequence dependence not set for RNA model, wobble base pairs will be ignored", file=stderr)

        backend = oxpy.analysis.AnalysisBackend(inp)

        # The 8 energies are:
        # 0 fene
        # 1 bexc
        # 2 stack
        # 3 nexc
        # 4 hb
        # 5 cr_stack
        # 6 cx_stack
        # 9 Debye-Huckel
        # 7 total
        if ctx.visualize:
            energies = np.zeros((ctx.top_info.nbases, 9))

        while backend.read_next_configuration():
            e_txt = backend.config_info().get_observable_by_id("my_obs").get_output_string(0).strip().split('\n')
            if ctx.visualize:
                for e in e_txt[1:]:
                    if not e[0] == '#':
                        e = e.split()
                        p = int(e[0])
                        q = int(e[1])
                        l = np.array([float(x) for x in e[2:]])*ctx.conversion_factor
                        energies[p] += l
                        energies[q] += l
            else:
                print(e_txt[0])
                for e in e_txt[1:]:
                    if not e[0] == '#':
                        e = e.split()
                        p = int(e[0])
                        q = int(e[1])
                        l = np.array([float(x) for x in e[2:]])*ctx.conversion_factor
                        print("{} {} {} {} {} {} {} {} {} {} {}".format(p, q, l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8]))
                    else: 
                        print(e)
        if ctx.visualize:
            return energies
        else:
            return
                
def main():
    import argparse
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="List all the interactions between nucleotides")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-v', type=str, nargs=1, dest='outfile', help='if you want instead average per-particle energy as a viewer JSON')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-u', '--units', type=str, nargs=1, dest='units', help="(optional) The units of the energy (pNnm or oxDNA)")
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    traj_file = args.trajectory[0]
    inputfile = args.inputfile[0]

    top_info, traj_info  = no_top_describe(traj_file)

    try:
        outfile = args.outfile[0]
        visualize = True
    except:
        visualize = False

    #if path.dirname(inputfile) != getcwd():
    #    sim_directory = path.dirname(inputfile)
    #else:
    #    sim_directory = ""

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    if args.units:
        if args.units[0] == "pNnm":
            units = "pN nm"
            conversion_factor = 41.42
        elif args.units[0] == "oxDNA":
            units = "oxDNA su"
            conversion_factor = 1
        else:
            print("Unrecognized units:", args.units[0], file=stderr)
            exit(1)
    else:
        units = "oxDNA su"
        conversion_factor = 1
        print("INFO: no units specified, assuming oxDNA su", file=stderr)

    # how many confs we want to distribute between the processes
    ntopart = 20
    pool = Pool(ncpus)

    # deduce how many chunks we have to run in parallel
    n_confs  = traj_info.nconfs 
    n_chunks = int(n_confs / ntopart +
                         (1 if n_confs % ntopart else 0))


    ctx = ComputeContext(traj_info, top_info, inputfile, visualize, conversion_factor, ntopart)

    ## Distribute jobs to the worker processes
    print(f"Starting up {ncpus} processes for {n_chunks} chunks")
    results = [pool.apply_async(compute,(ctx,i)) for i in range(n_chunks)]
    print("All spawned, waiting for results")

    energies = np.zeros((ctx.top_info.nbases, 9))

    for i, r in enumerate(results):
        if visualize:
            energies += r.get()
        else:
            r.get()
        print(f"finished {i+1}/{n_chunks}",end="\r")

    if visualize:
        energies /= n_confs
        for i, potential in enumerate(["FENE","bexc", "stack", "nexc", "hb", "cr_stack", "cx_stack", "Debye-Huckel", "Total"]):
            fname = '.'.join(outfile.split('.')[:-1])+"_"+potential+'.json'
            with open(fname, 'w+') as f:
                f.write("{{\n\"{} ({})\" : [".format(potential, units))
                f.write(', '.join([str(x) for x in energies[:,i]]))
                f.write("]\n}")
            print("INFO: Wrote oxView overlay to:", fname)


if __name__ == "__main__":
    main()