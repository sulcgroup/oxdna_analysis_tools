#!/usr/bin/env python

import numpy as np
from sys import exit, stderr
import argparse
import os
import matplotlib.pyplot as plt

#Calculates distance taking PBC into account
def min_image(p1, p2, box):
    """
    Calculates distance between two particles taking PBC into account

    Parameters:
        p1 (np.array): The first particle's position
        p2 (np.array): The second particle's position

    Returns:
        distance (float): The distance between the two particles
    """
    p1 = p1 - (np.floor(p1/box) * box)
    p2 = p2 - (np.floor(p2/box) * box)
    diff = p1 - p2
    diff = diff - (np.round(diff/box)*box)
    return np.linalg.norm(diff)

def get_distances(trajectories, p1s, p2s):
    """
    Calculates specified distances in each of the provided trajectories

    Parameters:
        trajectories (string[]): A list of trajectory file names
        p1s (int[]): A list of particle indexes for the first particle in each distance
        p2s (int[]): A list of particle indexes for the second particle in each distance

    Returns:
        distances (float[][]): A list of distances for each particle pair in each trajectory
    """
    distances = [[] for _ in trajectories]
    for i,trajectory in enumerate(trajectories):
        with open(trajectory, 'r') as traj:
            distances[i] = [[] for _ in p1s[i]]
            l = traj.readline()
            l = traj.readline() #skip the first time line

            #get the box size
            box = np.array([l.split(' ')[2], l.split(' ')[3], l.split(' ')[4]], dtype=float)

            #Read the file line by line and make a dictionary of particle positions
            line_num = 1
            d = {}
            while l:
                #If you're at the right particle, save the COM of the particle to a dictionary
                if line_num-3 in p1s[i] or line_num-3 in p2s[i]:
                    d[line_num-3] = np.array(l.split(' ')[0:3], dtype=float)

                #if its the start of a new configuration, dump the distance data from the last one
                if 't' in l:
                    for j, (p1, p2) in enumerate(zip(p1s[i], p2s[i])):
                        p1 = d[p1]
                        p2 = d[p2]
                        distances[i][j].append(min_image(p1, p2, box)*0.85) #1 oxDNA su = 0.85 nm
                        line_num = 0
                l = traj.readline() #returns false if there's no more conf to load
                line_num += 1
            
            #catch the last configuration
            for j, (p1, p2) in enumerate(zip(p1s[i], p2s[i])):
                p1 = d[p1]
                p2 = d[p2]
                distances[i][j].append(min_image(p1, p2, box)*0.85) #1 oxDNA su = 0.85 nm

    return distances

def main():
    #handle commandline arguments
    #this program has no positional arguments, only flags
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Finds the ensemble of distances between any two particles in the system")
    parser.add_argument('-i', '--input', metavar='input', nargs='+', action='append', help='An input, trajectory, and a list of particle pairs to compare.  Can call -i multiple times to plot multiple datasets.')
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The name to save the graph file to')
    parser.add_argument('-f', '--format', metavar='<histogram/trajectory/both>', nargs=1, help='Output format for the graphs.  Defaults to histogram.  Options are \"histogram\", \"trajectory\", and \"both\"')
    parser.add_argument('-d', '--data', metavar='data_file', nargs=1, help='If set, the output for the graphs will be dropped as a json to this filename for loading in oxView or your own scripts')
    parser.add_argument('-n', '--names', metavar='names', nargs='+', action='append', help='Names of the data series.  Will default to particle ids if not provided')
    parser.add_argument('-c', metavar='cluster', dest='cluster', action='store_const', const=True, default=False, help="Run the clusterer on each configuration's distance?")
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "matplotlib", "numpy"])

    #-i requires 4 or more arguments, the topology file of the structure, the trajectory to analyze, and any number of particle pairs to compute the distance between.
    try:
        input_files = [i[0] for i in args.input]
        trajectories = [i[1] for i in args.input]
        p1s = [i[2::2] for i in args.input]
        p2s = [i[3::2] for i in args.input]
        p1s = [[int(j) for j in i] for i in p1s]
        p2s = [[int(j) for j in i] for i in p2s]

    except Exception as e:
        print("ERROR:", e)
        parser.print_help()
        exit(1)
    
    #get number of distances to calculate
    n_dists = sum([len(l) for l in p1s])

    #Make sure that the input is correctly formatted
    if(len(input_files) != len(trajectories)):
        print("ERROR: bad input arguments\nPlease supply an equal number of input and, trajectory files", file=stderr)
        exit(1)
    if len(p1s) != len(p2s):
        print("ERROR: bad input arguments\nPlease supply an even number of particles", file=stderr)
        exit(1)


    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else: 
        print("INFO: No outfile name provided, defaulting to \"distance.png\"", file=stderr)
        outfile = "distance.png"

    #-f defines which type of graph to produce
    hist = False
    lineplt = False
    if args.format:
        if "histogram" in args.format:
            hist = True
        if "trajectory" in args.format:
            lineplt = True
        if "both" in args.format:
            hist = True
            lineplt = True
        if hist == lineplt == False:
            print("ERROR: unrecognized graph format\nAccepted formats are \"histogram\", \"trajectory\", and \"both\"", file=stderr)
            exit(1)
    else:
        print("INFO: No graph format specified, defaulting to histogram", file=stderr)
        hist = True

    #-c makes it run the clusterer on the output
    cluster = args.cluster
    
    #get the specified distances
    distances = get_distances(trajectories, p1s, p2s)

    # -n sets the names of the data series
    if args.names:
        names = args.names[0]
        if len(names) < n_dists:
            print("WARNING: Names list too short.  There are {} items in names and {} distances were calculated.  Will pad with particle IDs".format(len(names), n_dists), file=stderr)
            for i in range(len(names), len(distances)):
                names.append("{}-{}".format([j for sl in p1s for j in sl][i], [j for sl in p2s for j in sl][i]))
        if len(names) > n_dists:
            print("WARNING: Names list too long. There are {} items in names and {} distances were calculated.  Truncating to be the same as distances".format(len(names), n_dists), file=stderr)
            names = names[:n_dists]

    else:
        print("INFO: Defaulting to particle IDs as data series names")
        names = ["{}-{}".format(p1, p2) for p1, p2 in zip([i for sl in p1s for i in sl], [i for sl in p2s for i in sl])]
    
    # -d will dump the distances as json files for loading with the trajectories in oxView
    if args.data:
        from json import dump
        if len(trajectories) > 1:
            print("INFO: distance lists from separate trajectories are printed to separate files for oxView compatibility.  Trajectory numbers will be appended to your provided data file name.", file=stderr)
            file_names = ["{}_{}.json".format(args.data[0].strip('.json'), i) for i,_ in enumerate(trajectories)]
        else:
            file_names = [args.data[0].strip('.json')+'.json']
        names_by_traj = [['{}-{}'.format(p1, p2) for p1, p2 in zip(p1l, p2l)] for p1l, p2l in zip(p1s, p2s)]
        
        for file_name, ns, dist_list in zip(file_names, names_by_traj, distances):
            obj = {}
            for n, d in zip(ns, dist_list):
                obj[n] = d        
            with open(file_name, 'w+') as f:
                print("INFO: writing data to {}.  This can be opened in oxView using the Order parameter selector".format(file_name))
                dump(obj, f)

    #convert the distance list into numpy arrays because they're easier to work with
    for i, l in enumerate(distances):
        distances[i] = np.array(l)
    
    means = [np.mean(i, axis=1) for i in distances]
    medians = [np.median(i, axis=1) for i in distances]
    stdevs = [np.std(i, axis=1) for i in distances]

    #get some min/max values to make the plots pretty
    lower = min((l.min() for l in distances))
    upper = max((l.max() for l in distances))

    #those horrific list comprehensions unpack lists of lists into a single list
    print("input:\t", end='')
    [print("{}-{}\t".format(p1, p2), end='') for p1, p2 in zip([i for sl in p1s for i in sl], [i for sl in p2s for i in sl])]
    print("")

    print("name:\t", end='')
    [print("{}\t".format(t), end='') for t in names[:n_dists]]
    print("")

    print("mean:\t", end='')
    [print("{:.2f}\t".format(m), end='') for m in [i for sl in means for i in sl]]
    print("")

    print("stdev:\t", end='')
    [print("{:.2f}\t".format(s), end='') for s in [i for sl in stdevs for i in sl]]
    print("")

    print("median:\t", end='')
    [print("{:.2f}\t".format(m), end='') for m in [i for sl in medians for i in sl]]
    print("")

    #make a histogram
    if hist == True:
        if lineplt == True:
            #if making two plots, automatically append the plot type to the output file name
            out = outfile[:outfile.find(".")]+"_hist"+outfile[outfile.find("."):]
        else:
            out = outfile
        bins = np.linspace(np.floor(lower-(lower*0.1)), np.ceil(upper+(upper*0.1)), 60)
        graph_count = 0
        for traj_set in distances:
            for dist_list in traj_set:
                a = plt.hist(dist_list, bins, weights=np.ones(len(dist_list)) / len(dist_list),  alpha=0.5, histtype=u'stepfilled', edgecolor='k', label=names[graph_count])
                graph_count += 1
        plt.xlabel("Distance (nm)")
        plt.ylabel("Normalized frequency")
        plt.legend()
        #plt.show()
        print("INFO: Writing histogram to file {}".format(out), file=stderr)
        plt.savefig("{}".format(out))

    #make a trajectory plot
    if lineplt == True:
        if hist == True:
            #clear the histogram plot
            plt.clf()
            #if making two plots, automatically append the plot type to the output file name
            out = outfile[:outfile.find(".")]+"_traj"+outfile[outfile.find("."):]
        else:
            out = outfile
        graph_count = 0
        for traj_set in distances:
            for dist_list in traj_set:
                a = plt.plot(dist_list, alpha=0.5, label=names[graph_count])
                graph_count += 1
        plt.xlabel("Simulation Steps")
        plt.ylabel("Distance (nm)")
        plt.legend()
        #plt.show()
        print("INFO: Writing trajectory plot to file {}".format(out), file=stderr)
        plt.savefig("{}".format(out))

    if cluster == True:
        if not all([x == trajectories[0] for x in trajectories]):
            print("ERROR: Clustering can only be run on a single trajectory", file=stderr)
            exit(1)
        from oxDNA_analysis_tools.clustering import perform_DBSCAN

        labs = perform_DBSCAN(distances[0].T, len(distances[0][0]), trajectories[0], input_files[0], "euclidean", 12, 8)

if __name__ == '__main__':
    main()
    