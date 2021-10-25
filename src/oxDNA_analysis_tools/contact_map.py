#!/usr/bin/env python3

from os import environ, path
import sys
import subprocess
import tempfile
import numpy as np

from oxDNA_analysis_tools.config import set_analysis_path
PROCESSPROGRAM = set_analysis_path()

from oxDNA_analysis_tools.config import check_dependencies
check_dependencies(["numpy"])

def contact_map (inputfile, mysystem, return_full_matrix):
    """
    Computes the distance between every pair of nucleotides and creates a matrix of these distances.

    Parameters:
        inputfile (string): the input file with which the simulation was run,
        mysystem (base.System): a base.py system object containing the system to analyze.
        return_full_matrix (bool): The matrix is symmetric. Return only the lower half or the whole matrix?
    
    Returns:
        distances (numpy.array): The matrix containing pairwise distances between every pair of nucleotides.
    """
    tempfile_obj = tempfile.NamedTemporaryFile()

    #the algorithm to find the distances is written into an oxDNA observable for improved speed on large systems.
    command_for_data =  'analysis_data_output_1 = { \n name = stdout \n print_every = 1 \n col_1 = { \n type=contact_map \n } \n}'
    launchargs = [PROCESSPROGRAM, inputfile ,'trajectory_file='+tempfile_obj.name,command_for_data]
    n = mysystem.N

    #write the system to a tempfile so it can be fed to oxDNA in a subprocess
    mysystem.print_lorenzo_output(tempfile_obj.name,'/dev/null')
    tempfile_obj.flush()
    process = subprocess.run(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    
    #process output
    out = process.stdout.strip()
    err = process.stderr.strip()
    for line in err.split("\n"):
        if "CRITICAL" in line or "ERROR" in line or "Error" in line:
           print (err)
           sys.exit(1)
    distances = np.array(out.split(' '), dtype=float)

    #the distance matrix is diagonal, so the observable only calculates half in a compressed format
    #if the full matrix is needed, this copies it out.
    if return_full_matrix:
        full_matrix = np.empty((n, n))
        for i in range(n):
            for j in range(n):
                if i==j:
                    full_matrix[i][j] = 0
                elif j > i:
                    full_matrix[i][j] = distances[int(i*n - ((i*i)+i)*0.5 + (j-i-1))]
                elif i > j:
                    full_matrix[i][j] = distances[int(j*n - ((j*j)+j)*0.5 + (i-j-1))]
        return(full_matrix)
    else:
        return(distances)

def main():
    import argparse
    import matplotlib.pyplot as plt
    from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2, get_input_parameter

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy", "matplotlib"])

    #get commandline arguments
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Calculate and display the contact map for a structure")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help="The file containing the configurations of which the contact map is needed")
    parser.add_argument('-v', dest='visualize', action='store_const', const=True, default=False, help="should we display the contact map once its calculated? Only recommend if there are few confs.")

    args = parser.parse_args()
    visualize = args.visualize
    inputfile = args.inputfile[0]
    traj_file = args.trajectory[0]

    #process files
    top_file = get_input_parameter(inputfile, "topology")
    if "RNA" in get_input_parameter(inputfile, "interaction_type"):
        environ["OXRNA"] = "1"
    else:
        environ["OXRNA"] = "0"

    #create system object from first configuration in the trajectory
    r = LorenzoReader2(traj_file, top_file)
    system = r._get_system()

    #for every configuration, create a graphical contact map
    while system:
        m = contact_map(inputfile, system, True)
        if visualize:
            fig, ax = plt.subplots()
            a = ax.imshow(m, cmap='viridis', origin='lower')
            ax.set(title = "interaction network",
                ylabel="nucleotide id",
                xlabel="nucleotide id")
            b = fig.colorbar(a, ax=ax)
            b.set_label("distance", rotation = 270)
            plt.show()
        system = r._get_system()

if __name__ == '__main__':
    main()
    
