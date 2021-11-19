#!/usr/bin/env python3

# Created by Hao Liu 
# Date 01/22/2019 
# A short script generating force file from the given .dat and .top

from sys import stderr
import os
from oxDNA_analysis_tools.output_bonds import output_bonds
import argparse
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2, get_input_parameter

def main():
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Create an external forces file enforcing the current base-pairing arrangement")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('configuration', type=str, nargs=1, help="The configuration to generate the forces from")
    parser.add_argument('-o', '--output', type=str, nargs=1, help='name of the file to write the forces to. Defaults to forces.txt')
    parser.add_argument('-f', '--pairs', type=str, nargs=1, help='name of the file to write the designed pairs list to')

    args = parser.parse_args()

    #Process command line arguments:
    inputfile = args.inputfile[0]
    conf_file = args.configuration[0]

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else: 
        outfile = "forces.txt"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    if args.pairs:
        pairsfile = args.pairs[0]
    else:
        pairsfile = False     

    #Get relevant parameters from the input file
    top_file = get_input_parameter(inputfile, "topology")

    #get base pairs 
    r = LorenzoReader2(conf_file, top_file)
    mysystem = r._get_system()
    out = output_bonds(inputfile, mysystem)
    out = out.split('\n')

    #Find out the forming bonds series 
    print("INFO: Analyze the output...", file=stderr)
    Bonded = {}
    for i in out:
        if i[0] == '#':
            continue
        splitline = i.split(' ')
        try:
            HB = float(splitline[6])
        except:
            continue
        if HB < -0.001:
            ntid0 = int(splitline[0])
            ntid1 = int(splitline[1])
            if ntid0 not in Bonded:
                Bonded[ntid0] = ntid1
            if ntid1 not in Bonded:
                Bonded[ntid1] = ntid0

    lines = []
    pairlines = []
    mutual_trap_template = '{ \ntype = mutual_trap\nparticle = %d\nstiff = 0.9\nr0 = 1.2\nref_particle = %d\nPBC=1\n}\n'
    for key in sorted(Bonded):
        from_particle_id = key
        to_particle_id = Bonded[key]
        if from_particle_id < to_particle_id:
            if pairsfile:
                pairlines.append("{} {}\n".format(from_particle_id, to_particle_id))
            lines.append(mutual_trap_template % (from_particle_id,to_particle_id))
            lines.append(mutual_trap_template % (to_particle_id,from_particle_id))

    if pairsfile:
        with open(pairsfile, "w") as file:
            file.writelines(pairlines)
            print("INFO: Wrote pairs to {}".format(pairsfile), file=stderr)

    with open(outfile, "w") as file:
            file.writelines(lines)
            print("INFO: Job finished. Wrote forces to {}".format(outfile), file=stderr)

if __name__ == '__main__':
    main()

