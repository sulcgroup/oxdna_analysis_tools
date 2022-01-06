#!/usr/bin/env python3

import numpy as np
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2, cal_confs, get_input_parameter
from oxDNA_analysis_tools.UTILS import parallelize_lorenzo_onefile
import os
import argparse
from json import dumps

def rad2degree(angle):
    """
    Convert radians to degrees

    Parameters:
        angle (float): angle in radians

    Returns:
        angle (float): angle in degrees
    """
    return (angle * 180 / np.pi)

def get_internal_coords(reader, num_confs, start = None, stop = None):
    """
    Get the internal coordinates of the nucleic acid backbone

    Parameters:
        reader (LorenzoReader2): the reader object attached to a trajectory file
        num_confs (int): the number of conformations in the trajectory
        start (int): the starting frame to read in
        stop (int): the last frame to read in

    Returns:
        torsions (list[]): the torsion angles for each configuration
        dihedrals (list[]): the dihedral angles for each configuration
    """
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)

    confid = 0

    mysystem = reader._get_system(N_skip = start)
    torsions = np.zeros((len(mysystem._nucleotides)-2, stop))
    dihedrals = np.zeros((len(mysystem._nucleotides)-3, stop))

    while mysystem != False and confid < stop:
        for i, strand in enumerate(mysystem._strands):
            for j, n in enumerate(strand._nucleotides):
                #store first nucleotide
                if j == 0:
                    back3 = n.cm_pos
                    continue
                #store second nucleotide
                if j == 1:
                    back2 = n.cm_pos
                    continue
                #store third nucleotide and calculate first torsion
                if j == 2:
                    back1 = n.cm_pos
                    A = back2 - back3
                    B = back2 - back1
                    torsions[n.index-2][confid] = rad2degree(
                        np.arccos((np.dot(A, B))/(np.linalg.norm(A)*np.linalg.norm(B))))
                    continue

                #actually begin the loop of calculating torsions and dihedrals
                curr = n.cm_pos
                A = back3 - back2
                B = back2 - back1
                C = back1 - curr

                #get torsion angle
                torsions[n.index-2][confid] = rad2degree(
                    np.arccos((np.dot(B, -C))/(np.linalg.norm(B)*np.linalg.norm(-C))))

                #get dihedral angle
                n1 = np.cross(A, B)
                n2 = np.cross(B, C)
                dihedrals[n.index-3][confid] = rad2degree(
                    np.arccos(np.linalg.norm(np.dot(n1, n2)) / (np.linalg.norm(n1)*np.linalg.norm(n2))))
                back3 = back2
                back2 = back1
                back1 = curr
        confid += 1
        mysystem = reader._get_system()
    
    return(torsions, dihedrals)

def main():
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Computes the deviations in the backbone torsion angles")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('topology', type=str, nargs=1, help="The topology file associated with the trajectory file")
    parser.add_argument('outfile', type=str, nargs=1, help='The file name for the output .json file.')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    args = parser.parse_args()

    #run system checks
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    top_file  = args.topology[0]
    traj_file = args.trajectory[0]
    parallel = args.parallel
    if parallel:
        n_cpus = args.parallel[0]

    num_confs = cal_confs(traj_file)

    r = LorenzoReader2(traj_file, top_file)

    if not parallel:
        torsions, dihedrals = get_internal_coords(r, num_confs)

    if parallel:
        out = parallelize_lorenzo_onefile.fire_multiprocess(traj_file, top_file, get_internal_coords, num_confs, n_cpus)
        # Out Dims: 1 Processor, 2 Torsion or Dihedrals, 3 Specific list of torsions listed by conf
        torsions = np.concatenate([out[i][0] for i in range(n_cpus)], axis=1)
        dihedrals = np.concatenate([out[i][1] for i in range(n_cpus)], axis=1)

    torsion_mean = np.mean(torsions, axis=1).tolist()
    dihedral_mean = np.mean(dihedrals, axis=1).tolist()
    #make something akin to a ramachandran plot for DNA origami??
    import matplotlib.pyplot as plt
    plt.scatter(torsion_mean[1:], dihedral_mean)
    plt.xlabel("torsion_angle")
    plt.ylabel("dihedral_angle")
    plt.show()

    torsion_mean.insert(0, torsion_mean[0])
    torsion_mean.insert(0, torsion_mean[0])
    with open(args.outfile[0], "w") as file:
        file.write(dumps({
            "torsion" : torsion_mean
        }))

if __name__ == '__main__':
    main()
            