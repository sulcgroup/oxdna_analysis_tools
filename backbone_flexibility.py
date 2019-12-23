import numpy as np
from UTILS.readers import LorenzoReader2, Cal_confs, get_input_parameter
from UTILS import parallelize
from os import environ
import argparse
from json import dumps

def rad2degree(angle):
    return (angle * 180 / np.pi)

def get_torsions(reader, num_confs, start = None, stop = None): #YOU ALSO NEED TO RECORD AVERAGE A1, A2, A3 SO YOU CAN BACK CARTESIAN COORDS OUT AT THE END!!!
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)

    confid = 0

    mysystem = reader._get_system(N_skip = start)
    torsions = np.zeros((len(mysystem._nucleotides)-2, stop))

    while mysystem != False and confid < stop:
        for i, strand in enumerate(mysystem._strands):
            for j, n in enumerate(strand._nucleotides):
                if j == 0:
                    back2 = n.cm_pos
                    continue
                if j == 1:
                    back1 = n.cm_pos
                    continue
                curr = n.cm_pos
                A = back1 - back2
                B = back1 - curr
                torsions[n.index-2][confid] = rad2degree(
                    np.arccos((np.dot(A, B))/(np.linalg.norm(A)*np.linalg.norm(B))))
                back2 = back1
                back1 = curr
        confid += 1
        mysystem = reader._get_system()
    
    return(torsions)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computes the deviations in the backbone torsion angles")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('topology', type=str, nargs=1, help="The topology file associated with the trajectory file")
    parser.add_argument('outfile', type=str, nargs=1, help='The file name for the output .json file.')
    parser.add_argument('-p', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    args = parser.parse_args()

    top_file  = args.topology[0]
    traj_file = args.trajectory[0]
    parallel = args.parallel
    if parallel:
        n_cpus = args.parallel[0]

    num_confs = Cal_confs(traj_file, top_file)

    r = LorenzoReader2(traj_file, top_file)

    if not parallel:
        torsions = get_torsions(r, num_confs)

    if parallel:
        out = parallelize.fire_multiprocess(traj_file, top_file, get_torsions, num_confs, n_cpus)
        torsions = np.concatenate([i for i in out], axis=1)

    torsions = np.std(torsions, axis=1).tolist()
    torsions.append(torsions[len(torsions)-1])
    torsions.insert(0, torsions[0])

    print(torsions)

    with open(args.outfile[0],"w") as file:
        file.write(
            dumps({
            "StDev(torsion) (deg)" : torsions
        }))