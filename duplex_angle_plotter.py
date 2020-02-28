#!/usr/bin/env python3

import numpy as np
from sys import argv, stderr, exit #argv is much better than argparse for how I'm structuring my args
from math import acos, sqrt
from os import environ
import argparse

def rad2degree(angle):
    """
    Convert radians to degrees

    Parameters:
        angle (float): The angle in radians to convert.

    Returns:
        angle (float): The angle converted to degrees.
    """
    return (angle * 180 / np.pi)

def angle_between (axis1, axis2):
    """
    Find the angle between two vectors.

    Parameters:
        axis1 (numpy.array): The first vector.
        axis2 (numpy.array): The second vector.
    
    Returns:
        angle (float): The angle between the vectors in radians.
    """
    return (acos(np.dot(axis1, axis2)/(np.linalg.norm(axis1)*np.linalg.norm(axis2))))

if __name__ == "__main__":
    #Get command line arguments.
    parser = argparse.ArgumentParser(description="Finds the ensemble of distances between any two particles in the system")
    parser.add_argument('-i', metavar=('angle_file', 'particle1', 'particle2'), dest="input", nargs=3, action='append', help='An angle file, particle1, particle2 set.  Can call -i multiple times to plot multiple datasets.')
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The name to save the graph file to')
    parser.add_argument('-f', '--format', metavar='<histogram/trajectory/both>', nargs=1, help='Output format for the graphs.  Defaults to histogram.  Options are \"histogram\", \"trajectory\", and \"both\"')
    args = parser.parse_args()

    try:
        files = [i[0] for i in args.input]
        p1s = [i[1] for i in args.input]
        p2s = [i[2] for i in args.input]
    except Exception as e:
        print("ERROR: Failed to read files")
        print(e)
        parser.print_help()
        exit(1)
    
    n_angles = len(p1s)

    #Make sure that the input is correctly formatted
    if(len(files) != len(p1s) != len(p2s)):
        print("ERROR: bad input arguments\nPlease supply an equal number of trajectory and particle pairs", file=stderr)
        exit(1)

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else: 
        if environ.get('DISPLAY', None) != "":
            print("INFO: No display detected, outputting to \"distance.png\"")
            outfile=False
        else:
            print("INFO: No outfile name provided, defaulting to \"distance.png\"", file=stderr)
            outfile = "distance.png"

    #-f defines which type of graph to produce
    hist = False
    line = False
    if args.format:
        if "histogram" in args.format:
            hist = True
        if "trajectory" in args.format:
            line = True
        if "both" in args.format:
            hist = line = True
        if hist == line == False:
            print("ERROR: unrecognized graph format\nAccepted formats are \"histogram\", \"trajectory\", and \"both\"", file=stderr)
            exit(1)
    else:
        print("INFO: No graph format specified, defaulting to histogram", file=stderr)
        hist = True

    all_angles = []
    means = []
    medians = []
    stdevs = []
    representations = []

    #For each input triplet
    for anglefile, search1, search2 in zip(files, p1s, p2s):

        steps = 1
        last_step = 0
        axis1 = axis2 = np.array([0,0,0])
        angles = []
        found = False

        #the format of the angle file is as follows: (tbh this should be a JSON)
        # 0: time
        # 1: duplex id
        # 2: strand 1 start nucleotide id
        # 3: strand 1 end nucleotide id
        # 4: strand 2 start nucleotide id
        # 5: strand 2 end nucleotide id
        # 6: X-component of the axis vector
        # 7: Y-component of the axis vector
        # 8: Z-component of the axis vector
        # 9: Helix position

        with open(anglefile) as file:
            for l in file.readlines()[1:]: #the first line is a header, so it can be dropped
                try:
                    l = l.split("\t")
                    t = float(l[0])
                except Exception as e:
                    print("ERROR: The following line is incorrectly formatted:")
                    print(l)
                    print("The error was:\n",e)
                    print("skiping the line")
                    continue

                #reset if we're in a new time
                if (t != last_step):
                    found = False
                    steps += 1
                    axis1 = axis2 = np.array([0,0,0])

                #don't need to do anything if both angles were already found for this timestep
                if found:
                    continue

                #look for the nucleotide IDs
                if l[2] == search1 or l[3] == search1:
                    axis1 = np.array([float(l[6]), float(l[7]), float(l[8])])
                if l[2] == search2 or l[3] == search2:
                    axis2 = -1*np.array([float(l[6]), float(l[7]), float(l[8])])

                #once both are found, add them to angle list
                if np.linalg.norm(axis1) > 0 and np.linalg.norm(axis2) > 0:
                    angles.append(rad2degree(angle_between(axis1, axis2))) #add a -90 here if your duplexes in question are antiparallel
                    found = True

                last_step = t

        #compute some statistics
        angles = np.array(angles)
        mean = np.mean(angles)
        median = np.median(angles)
        stdev = np.std(angles)
        representation = len(angles)/steps

        #add to the output data
        all_angles.append(angles)
        means.append(mean)
        medians.append(median)
        stdevs.append(stdev)
        representations.append(representation)

    for i, m in enumerate(means):
        if m > 90:
            all_angles[i] = [180 - a for a in all_angles[i]]
            means[i] = 180 - m
            medians[i] = 180 - medians[i]

    #PUT THE NAMES OF YOUR DATA SERIES HERE
    names = ["1", "2", "3", "4", "5", "6", "7", "8"]
    print("INFO: Name your data series by modifying the \"names\" variable in the script", file=stderr)
    if len(names) < n_angles:
        print("ERROR: Not enough names provided.  There are {} items in the names list and {} data series".format(len(names), n_angles), file=stderr)
        print("INFO: Defaulting to particle IDs as data series names")
        names = ["{}-{}".format(p1, p2) for p1, p2 in zip(p1s, p2s)]

    #print statistical information
    print("name:\t", end='')
    [print("{}\t".format(t), end='') for t in names[:n_angles]]
    print("")

    print("mean:\t", end='')
    [print("{:.2f}\t".format(m), end='') for m in means]
    print("")

    print("stdevs:\t", end='')
    [print("{:.2f}\t".format(s), end='') for s in stdevs]
    print("")

    print("median:\t", end='')
    [print("{:.2f}\t".format(m), end='') for m in medians]
    print("")

    print("freqs:\t", end='')
    [print("{:.2f}\t".format(r), end='') for r in representations]
    print("")

    #make a histogram
    import matplotlib.pyplot as plt
    if outfile and hist == True:
        if line == True:
            out = outfile[:outfile.find(".")]+"_hist"+outfile[outfile.find("."):]
        else:
            out = outfile
    
        bins = np.linspace(0, 180, 60)
    
        artists = []
        for i,alist in enumerate(all_angles):
            a = plt.hist(alist, bins, weights=np.ones(len(alist)) / len(alist),  alpha=0.3, label=names[i], histtype=u'stepfilled', edgecolor='k')
            artists.append(a)
        plt.legend(labels=names)
        plt.xlim((0, 180))
        plt.xlabel("Angle (degrees)")
        plt.ylabel("Normalized frequency")
        if outfile:
            print("INFO: Saving histogram to {}".format(out), file=stderr)
            plt.savefig(out)
        else:
            plt.show()

    #make a trajectory plot
    if outfile and line == True:
        if hist == True:
            plt.clf()
            out = outfile[:outfile.find(".")]+"_traj"+outfile[outfile.find("."):]
        else:
            out = outfile
        
        artists = []
        for i,alist in enumerate(all_angles):
            a = plt.plot(alist)
            artists.append(a)
        plt.legend(labels=names)
        plt.xlabel("Configuration Number")
        plt.ylabel("Angle (degrees)")
        if outfile:
            print("INFO: Saving line plot to {}".format(out), file=stderr)
            plt.savefig(out)
        else:
            plt.show()
