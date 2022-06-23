#!/usr/bin/env python3

from typing import List, Tuple
import numpy as np
from sys import stderr, exit
from os import environ, path
import argparse

def rad2degree(angle:float) -> float:
    """
    Convert radians to degrees

    Parameters:
        angle (float): The angle in radians to convert.

    Returns:
        angle (float): The angle converted to degrees.
    """
    return (angle * 180 / np.pi)

def angle_between (axis1:np.ndarray, axis2:np.ndarray) -> float:
    """
    Find the angle between two vectors.

    Parameters:
        axis1 (numpy.array): The first vector.
        axis2 (numpy.array): The second vector.
    
    Returns:
        angle (float): The angle between the vectors in radians.
    """
    return (np.arccos(np.dot(axis1, axis2)/(np.linalg.norm(axis1)*np.linalg.norm(axis2))))

def get_angle_between(files:List[str], p1s:List[List[int]], p2s:List[List[int]]) -> Tuple[List[List[np.array]], List[List[float]], List[List[float]], List[List[float]], List[List[float]]]:
    """
        Read in a duplex list file and return the angles between specified duplexes.

        Parameters:
            files (List[str]): The list of duplex files to read.
            p1s (List[List[int]]): The list of start1 nucleotide indices for each file.
            p2s (List[List[int]]): The list of start2 nucleotide indices for each file.

        Returns:
            angles (List[List[np.array]]): The list of angles between the specified duplexes.
            means (List[float]): The mean angle between each pair of duplexes.
            medians (List[float]): The median angle between each pair of duplexes.
            stdevs (List[float]): The standard deviation of the angle between each pair of duplexes.
            representations (List[float]): The percentage of confs that have the duplexes.
    """
    all_angles = [[] for _ in files]
    means = [[] for _ in files]
    medians = [[] for _ in files]
    stdevs = [[] for _ in files]
    representations = [[] for _ in files]

    #For each input triplet
    for i, (anglefile, search1, search2) in enumerate(zip(files, p1s, p2s)):

        steps = 0 #counts the number of configurations in the file
        last_step = 0
        count = 0
        all_angles[i] = [[] for _ in p1s[i]]
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
            all_search = search1.copy()
            all_search.extend(search2)
            d = {i : np.array([0, 0, 0]) for i in all_search}
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

                #dump values and reset if we're in a new time (but also skip the first pass)
                if (t != last_step):
                    if (steps != 0):
                        for j, (p1, p2) in enumerate(zip(search1, search2)):
                            if np.linalg.norm(d[p1]) != 0 and np.linalg.norm(d[p2]) != 0:
                                angle = rad2degree(angle_between(d[p1], -1*d[p2])) #add a -90 here if your duplexes in question are antiparallel
                                all_angles[i][j].append(angle)
                            else:
                                all_angles[i][j].append(np.nan)

                    found = False
                    steps += 1
                    d = dict.fromkeys(d, np.array([0, 0, 0]))
                    count = 0 #counts the number of search targets found

                #don't need to do anything if both angles were already found for this timestep
                if found:
                    continue

                #look for the nucleotide IDs.  The -1 on axis 2 assumes you're looking at contiguous duplexes
                for s in all_search:
                    idx = l.index(s, 2, 6) if s in l[2:6] else None
                    if idx:
                        d[s] = np.array([float(l[6]), float(l[7]), float(l[8])])
                        count += 1

                #once all are found, add them to angle list
                if count == len(d):
                    found = True

                last_step = t

        #catch last configuration
        for j, (p1, p2) in enumerate(zip(search1, search2)):
            if np.linalg.norm(d[p1]) != 0 and np.linalg.norm(d[p2]) != 0:
                angle = rad2degree(angle_between(d[p1], -1*d[p2])) #add a -90 here if your duplexes in question are antiparallel
                all_angles[i][j].append(angle)
            else:
                all_angles[i][j].append(np.nan)

        #compute some statistics
        all_angles[i] = [np.array(a) for a in all_angles[i]]
        mean = [np.nanmean(a) for a in all_angles[i]]
        median = [np.nanmedian(a) for a in all_angles[i]]
        stdev = [np.nanstd(a) for a in all_angles[i]]
        representation = [np.count_nonzero(~np.isnan(a))/steps for a in all_angles[i]]

        #add to the output data
        means[i] = mean
        medians[i] = median
        stdevs[i] = stdev
        representations[i] = representation

    return (all_angles, means, medians, stdevs, representations)

def main():
    #Get command line arguments.
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Finds the ensemble of angles between any two duplexes defined by a starting or ending nucleotide in the system")
    parser.add_argument('-i', '--input', metavar='angle_file', dest='input', nargs='+', action='append', help='An angle file from duplex_angle_finder.py and a list of duplex-end particle pairs to compare.  Can call -i multiple times to plot multiple datasets.')
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The name to save the graph file to')
    parser.add_argument('-f', '--format', metavar='<histogram/trajectory/both>', nargs=1, help='Output format for the graphs.  Defaults to histogram.  Options are \"histogram\", \"trajectory\", and \"both\"')
    parser.add_argument('-d', '--data', metavar='data_file', nargs=1, help='If set, the output for the graphs will be dropped as a json to this filename for loading in oxView or your own scripts')
    parser.add_argument('-n', '--names', metavar='names', nargs='+', action='append', help='Names of the data series.  Will default to particle ids if not provided')

    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy", "matplotlib"])

    try:
        files = [i[0] for i in args.input]
        p1s = [i[1::2] for i in args.input]
        p2s = [i[2::2] for i in args.input]
    except Exception as e:
        print("ERROR: Failed to read files")
        print(e)
        parser.print_help()
        exit(1)

    n_angles = sum(len(p) for p in p1s)

    #Make sure that the input is correctly formatted
    if(len(files) != len(p1s) != len(p2s)):
        print("ERROR: bad input arguments\nPlease supply an equal number of trajectory and particle pairs", file=stderr)
        exit(1)

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else: 
        if environ.get('DISPLAY', None) != "":
            print("INFO: No display detected, outputting to \"angle.png\"", file=stderr)
            outfile=False
        else:
            print("INFO: No outfile name provided, defaulting to \"angle.png\"", file=stderr)
            outfile = "angle.png"

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

    all_angles, means, medians, stdevs, representations = get_angle_between(files, p1s, p2s)

    #for i, m in enumerate(means):
    #    if m > 90:
    #        all_angles[i] = [180 - a for a in all_angles[i]]
    #        means[i] = 180 - m
    #        medians[i] = 180 - medians[i]

        # -n sets the names of the data series
    if args.names:
        names = args.names[0]
        if len(names) < n_angles:
            print("WARNING: Names list too short.  There are {} items in names and {} angles were calculated.  Will pad with particle IDs".format(len(names), n_angles), file=stderr)
            for i in range(len(names), n_angles):
                names.append("{}-{}".format([j for sl in p1s for j in sl][i], [j for sl in p2s for j in sl][i]))
        if len(names) > n_angles:
            print("WARNING: Names list too long. There are {} items in names and {} angles were calculated.  Truncating to be the same as distances".format(len(names), n_angles), file=stderr)
            names = names[:n_angles]

    else:
        print("INFO: Defaulting to particle IDs as data series names")
        names = ["{}-{}".format(p1, p2) for p1, p2 in zip([i for sl in p1s for i in sl], [i for sl in p2s for i in sl])]

    # -d will dump the distances as json files for loading with the trajectories in oxView
    if args.data:
        from json import dump
        if len(files) > 1:
            f_names = [path.basename(f) for f in files]
            print("INFO: angle lists from separate trajectories are printed to separate files for oxView compatibility.  Trajectory names will be appended to your provided data file name.", file=stderr)
            file_names = ["{}_{}.json".format(args.data[0].strip('.json'), i) for i,_ in enumerate(f_names)]
        else:
            file_names = [args.data[0].strip('.json')+'.json']
        names_by_traj = [['{}-{}'.format(p1, p2) for p1, p2 in zip(p1l, p2l)] for p1l, p2l in zip(p1s, p2s)]

        for file_name, ns, ang_list in zip(file_names, names_by_traj, all_angles):
            obj = {}
            for n, a in zip(ns, ang_list):
                obj[n] = list(a)
            with open(file_name, 'w+') as f:
                print("INFO: writing data to {}.  This can be opened in oxView using the Order parameter selector".format(file_name))
                dump(obj, f)

    #print statistical information
    print("name:\t", end='')
    [print("{}\t".format(t), end='') for t in names[:n_angles]]
    print("")

    print("mean:\t", end='')
    [print("{:.2f}\t".format(m), end='') for m in [i for sl in means for i in sl]]
    print("")

    print("stdevs:\t", end='')
    [print("{:.2f}\t".format(s), end='') for s in [i for sl in stdevs for i in sl]]
    print("")

    print("median:\t", end='')
    [print("{:.2f}\t".format(m), end='') for m in [i for sl in medians for i in sl]]
    print("")

    print("freqs:\t", end='')
    [print("{:.2f}\t".format(r), end='') for r in [i for sl in representations for i in sl]]
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
        for i,traj_set in enumerate(all_angles):
            for alist in traj_set:
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
        for i,traj_set in enumerate(all_angles):
            for alist in traj_set:
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

if __name__ == '__main__':
    main()
    