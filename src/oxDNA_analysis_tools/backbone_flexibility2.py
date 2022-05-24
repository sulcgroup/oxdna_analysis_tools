import os
import argparse
from json import dumps
from multiprocessing import Pool
from collections import namedtuple
from sys import stderr
import numpy as np
import matplotlib.pyplot as plt
from oxDNA_analysis_tools.UTILS.RyeReader import describe, strand_describe
from oxDNA_analysis_tools.UTILS.get_confs import get_confs

import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "system",
                                              "ntopart"])

def rad2degree(angle:float):
    """
    Convert radians to degrees

    Parameters:
        angle (float): angle in radians

    Returns:
        angle (float): angle in degrees
    """
    return (angle * 180 / np.pi)

def get_internal_coords():
    pass

def compute(ctx:ComputeContext, chunk_id:int):
    torsions = np.zeros(ctx.top_info.nbases-(2*ctx.top_info.nstrands))
    dihedrals = np.zeros(ctx.top_info.nbases-(3*ctx.top_info.nstrands))

    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)
    for conf in confs:
        for i, strand in enumerate(ctx.system):
            for j, n in enumerate(strand):
                if j == 0:
                    back3 = conf.positions[n.id]
                    continue
                #store second nucleotide
                if j == 1:
                    back2 = conf.positions[n.id]
                    continue
                #store third nucleotide and calculate first torsion
                if j == 2:
                    back1 = conf.positions[n.id]
                    A = back2 - back3
                    B = back2 - back1
                    torsions[n.id-(2*(i+1))] += rad2degree(
                        np.arccos((np.dot(A, B))/(np.linalg.norm(A)*np.linalg.norm(B))))
                    continue

                #actually begin the loop of calculating torsions and dihedrals
                curr = conf.positions[n.id]
                A = back3 - back2
                B = back2 - back1
                C = back1 - curr

                #get torsion angle
                torsions[n.id-(2*(i+1))] += rad2degree(
                    np.arccos((np.dot(B, -C))/(np.linalg.norm(B)*np.linalg.norm(-C))))

                #get dihedral angle
                n1 = np.cross(A, B)
                n2 = np.cross(B, C)
                dihedrals[n.id-(3*(i+1))] += rad2degree(
                    np.arccos(np.linalg.norm(np.dot(n1, n2)) / (np.linalg.norm(n1)*np.linalg.norm(n2))))
                back3 = back2
                back2 = back1
                back1 = curr

    return(torsions, dihedrals)

def main():
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Computes the deviations in the backbone torsion angles")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('topology', type=str, nargs=1, help="The topology file associated with the trajectory file")
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', metavar='output_file', nargs=1, type=str, dest='output', help="(optional) The name of the file to write the graph to")
    parser.add_argument('-d', metavar='data_file', nargs=1, type=str, dest='data', help="(optional) The name of the file to write the data to")
    args = parser.parse_args()

    #run system checks
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy", "matplotlib"])

    top_file  = args.topology[0]
    traj_file = args.trajectory[0]
    top_info, traj_info = describe(top_file, traj_file)

    system, monomers = strand_describe(top_file)

    if args.parallel:
        n_cpus = args.parallel[0]
    else:
        n_cpus = 1

    # how many confs we want to distribute between the processes
    ntopart = 20
    pool = Pool(n_cpus)

    # deduce how many chunks we have to run in parallel
    n_confs  = traj_info.nconfs 
    n_chunks = int(n_confs / ntopart +
                         (1 if n_confs % ntopart else 0))

    # Create a ComputeContext which defines the problem to pass to the worker processes 
    ctx = ComputeContext(
        traj_info, top_info, system, ntopart
    )

    # Allocate memory to store the results
    torsions = np.zeros(ctx.top_info.nbases-(2*top_info.nstrands))
    dihedrals = np.zeros(ctx.top_info.nbases-(3*top_info.nstrands))

    ## Distribute jobs to the worker processes
    print(f"Starting up {n_cpus} processes for {n_chunks} chunks")
    results = [pool.apply_async(compute,(ctx,i)) for i in range(n_chunks)]
    print("All spawned, waiting for results")

    for i,r in enumerate(results):
        t, d = r.get()
        torsions += t
        dihedrals += d
        print(f"finished {i+1}/{n_chunks}",end="\r")

    pool.close()
    pool.join()

    torsions /= n_confs
    dihedrals /= n_confs

    torsions = torsions.tolist()
    dihedrals = dihedrals.tolist()

    if args.output:
        out = args.output[0]
    else:
        out = "ramachandran.png"
        print("INFO: No output file specified, writing to {}".format(out), file=stderr)

    plt.scatter(torsions[top_info.nstrands:], dihedrals)
    plt.xlabel("torsion_angle")
    plt.ylabel("dihedral_angle")
    plt.savefig(out)
    print("INFO: Wrote plot to {}".format(out), file=stderr)

    if args.data:
        out = args.data[0]
        with open(out, "w") as f:
            f.write(dumps({
                "torsions": torsions,
                "dihedrals": dihedrals
            }))
        print("INFO: Wrote angle data to {}".format(out), file=stderr)

    # Maybe some sort of overlay file?  Hard to do since there is 2*nstrands fewer torsions than particles.

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()