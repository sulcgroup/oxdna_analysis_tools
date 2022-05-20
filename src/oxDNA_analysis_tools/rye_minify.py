import argparse
from os import remove, path
from sys import stderr
from collections import namedtuple
from numpy import round
from multiprocessing import Pool
from oxDNA_analysis_tools.UTILS.RyeReader import no_top_describe, conf_to_str
from oxDNA_analysis_tools.UTILS.get_confs import get_confs

import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "d",
                                              "a",
                                              "ntopart"])

def compute(ctx:ComputeContext, chunk_id:int):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)
    
    for conf in confs:
        if ctx.d is not None: #round positions
            conf.positions = round(conf.positions, ctx.d)
            conf.a1s = round(conf.a1s, ctx.d)
            conf.a3s = round(conf.a3s, ctx.d)
        if ctx.a: #discard a vectors
            conf.a1s -= conf.a1s
            conf.a3s -= conf.a3s

    out = ''.join([conf_to_str(c) for c in confs])
    return out

def main():
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Compress given configuration.")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('outfile',    type=str, nargs=1, help='minified file')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-a', action = 'store_true', help='Discard a vectors.')
    parser.add_argument('-d', type=int, nargs=1,  help='Round positions and orientations to the specified number of digits.')

    args = parser.parse_args()

    traj_file = args.trajectory[0]
    out = args.outfile[0]

    top_info, traj_info = no_top_describe(traj_file)

    # -p sets the number of parallel processes
    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    # -d sets the decimals of the output
    if args.d:
        d = args.d[0]
    else:
        d = None

    # -a sets the a vectors to 0
    if args.a:
        a = True
    else:
        a = False
    
    # how many confs we want to distribute between the processes
    ntopart = 20
    pool = Pool(ncpus)

    # deduce how many chunks we have to run in parallel
    n_confs  = traj_info.nconfs 
    n_chunks = int(n_confs / ntopart +
                         (1 if n_confs % ntopart else 0))

    try:
        remove(out)
    except:
        pass

    ctx = ComputeContext(
        traj_info, top_info, d, a, ntopart)

    ## Distribute jobs to the worker processes
    print(f"Starting up {ncpus} processes for {n_chunks} chunks")
    results = [pool.apply_async(compute,(ctx,i)) for i in range(n_chunks)]
    print("All spawned, waiting for results")

    with open(out, 'w+') as f:
        for i,r in enumerate(results):
            chunk = r.get()
            print(f"finished {i+1}/{n_chunks}",end="\r")
            f.write(chunk)

    pool.close()
    pool.join()

    print(f"INFO: Wrote aligned trajectory to {out}", file=stderr)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()