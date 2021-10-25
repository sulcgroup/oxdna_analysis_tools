#!/usr/bin/env python3

from os import environ, path
import sys
import subprocess
import tempfile
import numpy as np
from functools import partial

from oxDNA_analysis_tools.config import set_analysis_path
PROCESSPROGRAM = set_analysis_path()

def all_vectors (inputfile, mysystem, return_full_matrix):
    tempfile_obj = tempfile.NamedTemporaryFile()
    command_for_data =  'analysis_data_output_1 = { \n name = stdout \n print_every = 1 \n col_1 = { \n type=all_vectors \n } \n}'
    launchargs = [PROCESSPROGRAM, inputfile ,'trajectory_file='+tempfile_obj.name,command_for_data]
    n = mysystem.N

    mysystem.print_lorenzo_output(tempfile_obj.name,'/dev/null')
    tempfile_obj.flush()
    process = subprocess.run(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    out = process.stdout.strip()
    err = process.stderr.strip()
    for line in err.split("\n"):
        if "CRITICAL" in line or "ERROR" in line or "Error" in line:
           print (err)
           sys.exit(1)
    mapper = partial(np.fromstring, sep=' ', dtype=float)
    vectors = np.array(list(map(mapper, out.split('\n'))))

    if return_full_matrix: #there must be a faster way to do this
        full_matrix = np.empty((n, n, 3))
        for i in range(n):
            for j in range(n):
                if i==j:
                    full_matrix[i][j] = np.array([0,0,0])
                elif j > i:
                    full_matrix[i][j] = vectors[int(i*n - ((i*i)+i)*0.5 + (j-i-1))]
                elif i > j:
                    full_matrix[i][j] = vectors[int(j*n - ((j*j)+j)*0.5 + (i-j-1))]
        return(full_matrix)
    else:
        return(vectors)

def main():
    #doesn't actually do anything...
    import argparse
    from UTILS.readers import LorenzoReader2, get_input_parameter
    parser = argparse.ArgumentParser(description="A python wrapper for getting all vectors between nucleotides from a simulation")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help="The file containing the configurations of which the contact map is needed")
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    inputfile = args.inputfile[0]
    traj_file = args.trajectory[0]

    top_file = get_input_parameter(inputfile, "topology")
    if "RNA" in get_input_parameter(inputfile, "interaction_type"):
        environ["OXRNA"] = "1"
    else:
        environ["OXRNA"] = "0"

    import UTILS.base #this needs to be imported after the model type is set

    r = LorenzoReader2(traj_file, top_file)
    system = r._get_system()

    while system:
        m = all_vectors(inputfile, system, True)
        system = r._get_system()

    print("well, it finished...")

if __name__ == '__main__':
    main()
