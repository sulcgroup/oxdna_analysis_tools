import os
import argparse
from sys import stderr
from collections import namedtuple
from multiprocessing import Pool
from copy import deepcopy
from oxDNA_analysis_tools.UTILS.RyeReader import describe, strand_describe, conf_to_str, get_top_string
from oxDNA_analysis_tools.UTILS.data_structures import Configuration
from oxDNA_analysis_tools.UTILS.get_confs import get_confs



ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "indexes",
                                              "ntopart"])

def compute(ctx:ComputeContext, chunk_id:int):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)
    outstr = [[] for _ in range(len(ctx.indexes))]
    for conf in confs:
        sub_confs = [Configuration(conf.time, conf.box, conf.energy, conf.positions[i], conf.a1s[i], conf.a3s[i]) for i in ctx.indexes]
        for i, sub_conf in enumerate(sub_confs):
            outstr[i].append(conf_to_str(sub_conf))

    return [''.join(out) for out in outstr]

def write_topologies(system, indexes, outfiles):
    top_names = [o+ ".top" for o in outfiles]
    for idx, top_name in zip(indexes, top_names):
        idx = set(idx)
        new_sys = deepcopy(system)
        for s in new_sys:
            if s[0].n3 != None and s[0].id != s[-1].n5:
                print(f"WARNING: Strand {s.id} is circular. Subsetting the trajectory will cut circular strands", file=stderr)
            s.monomers = [n for n in s if n.id in idx]
            
        new_sys.strands = [s for s in new_sys if len(s.monomers) > 0]

        with open(top_name, 'w+') as f:
            f.write(get_top_string(new_sys))

    return top_names


def main():
    #command line arguments
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Extracts parts of a structure into separate trajectories")
    parser.add_argument('trajectory', type=str, nargs=1, help="The trajectory file to subset")
    parser.add_argument('topology', type=str, nargs=1, help="The topology file corresponding to the trajectory")
    parser.add_argument('-i', '--index', metavar='index', action='append', nargs=2, help='A space separated index file and the associated output file name.  This can be called multiple times')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    args = parser.parse_args()

    top_file  = args.topology[0]
    traj_file = args.trajectory[0]
    index_files = [i[0] for i in args.index]
    output_files = [i[1] for i in args.index]
    top_info, traj_info = describe(top_file, traj_file)
    system, _ = strand_describe(top_file)
    indexes = []
    outfiles = []
    for i, o in zip(index_files, output_files):
        with open(i) as f:
            data = f.readline().split()
            try:
                data = sorted([int(i) for i in data])
            except:
                print("ERROR: The index file {} must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button".format(i))
        indexes.append(data)
        outfiles.append(o)

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    # how many confs we want to distribute between the processes
    ntopart = 20
    pool = Pool(ncpus)

    # deduce how many chunks we have to run in parallel
    n_confs  = traj_info.nconfs 
    n_chunks = int(n_confs / ntopart +
                         (1 if n_confs % ntopart else 0))

    # Create a ComputeContext which defines the problem to pass to the worker processes 
    ctx = ComputeContext(
        traj_info, top_info, indexes, ntopart
    )

    ## Distribute jobs to the worker processes
    print(f"Starting up {ncpus} processes for {n_chunks} chunks")
    results = [pool.apply_async(compute,(ctx,i)) for i in range(n_chunks)]
    print("All spawned, waiting for results")

    # Collect the trajectory file contents and write out
    dat_names = [o+ ".dat" for o in outfiles]
    files = [open(f, 'w+') for f in dat_names]
    for i, r in enumerate(results):
        subsets = r.get()
        for f, subset in zip(files, subsets):
            f.write(subset)
        print(f"finished {i+1}/{n_chunks}",end="\r")

    pool.close()
    pool.join()

    for f in files:
        f.close()

    # Write topology files
    top_names = write_topologies(system, indexes, outfiles)

    print("INFO: Wrote trajectories: {}".format(dat_names), file=stderr)
    print("INFO: Wrote topologies: {}".format(top_names), file=stderr)

if __name__ == '__main__':
    main()