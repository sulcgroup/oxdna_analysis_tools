#!/usr/bin/env python3
from json import loads
from sys import stderr
import argparse
import numpy as np
import os

def make_dat(mean_info, outfile):
    """
        Converts the json file produced by compute_mean.py into an oxDNA .dat file

        Parameters:
            mean_info (numpy.array): the json file from compute_mean.py parsed as a numpy array.
            outfile (string): the name to save the .dat file to
    """
    mean = np.array(mean_info['g_mean'])
    mmax = mean.max()
    mmin = mean.min()

    #box size is 1.5x the longest dimension of the structure
    box_size = (mmax - mmin).max() * 1.5
    with open(outfile, "w") as file:
        print("INFO: Writing configuration to", outfile, file=stderr)
        file.writelines([
            't = %d\n' % 0,                                       # set time
            'b = %d %d %d\n' % (box_size, box_size, box_size),    # set box size
            'E = 0 0 0\n'])                                       # set energy
        for p, a1, a3 in zip(mean_info["g_mean"], mean_info["a1_mean"] ,mean_info["a3_mean"]):
            file.write(
                "".join([
                        "%f %f %f " % (p[0],p[1],p[2]),
                        "%f %f %f " % (a1[0],a1[1],a1[2]),
                        "%f %f %f " % (a3[0],a3[1],a3[2]),
                        "0 0 0 ",
                        "0 0 0 "
                    ]) + "\n"
            )
def main(): 
    #get commandline arguments 
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Converts a mean structure .json from compute_mean.py to an oxDNA-readable .dat")
    parser.add_argument('mean', type=str, nargs=1, help="A mean structure from compute_mean.py")
    parser.add_argument('output', type=str, nargs=1, help="The name of the output file")
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    #load the mean file, which is in .json format
    with open(args.mean[0], "r") as file:
        mean_info = loads(
            file.read()
        )

    #write the file out in oxDNA format
    outfile = args.output[0]
    make_dat(mean_info, outfile)

if __name__ == '__main__':
    main()