#!/usr/bin/env python

#Created by: Erik Poppleton
#Date: 6/29/18
#Python2
#Converts the forces file printed out by tiamat2oxdna to a pairs file containing all designed H-bonds

def main():
    import sys

    if len(sys.argv) < 1:
        print("Usage is: {} infile (outfile)\n\nIf no outfile is given it will print to stdout".format(sys.argv[0]))

    infile = sys.argv[1]

    #if there's no outfile argument, just print to stdout.
    try:
        out = sys.argv[2]
        outfile = open(out, 'w+')

    except:
        outfile = sys.stdout

    #Process the forces file
    with open(infile) as f:
        for line in f:
            if line.startswith("particle"):
                a = line.split()[2]
            if "ref_particle" in line:
                b = line.split()[2]
            if "}" in line:
                if int(a) < int(b):
                    print(a, b, sep = ' ', file=outfile)
                a = -1
                b = -1

if __name__ == '__main__':
    main()
    