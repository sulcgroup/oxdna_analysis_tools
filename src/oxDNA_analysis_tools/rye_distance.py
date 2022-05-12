#!/usr/bin/env python

import numpy as np
from sys import exit, stderr
from multiprocessing import Pool
from collections import namedtuple
import argparse
import os
import matplotlib.pyplot as plt
from oxDNA_analysis_tools.UTILS.RyeReader import describe
from oxDNA_analysis_tools.UTILS.get_confs import get_confs

import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info", 
                                              "p1s",
                                              "p2s",
                                              "ntopart"])

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

def vectorized_min_image(p1s, p2s, box):
    """
    Calculates all mutual distances between two sets of points taking PBC into account
    
    Paramters:
        p1s (np.array): the first set of points (Nx3 array)
        p2s (np.array): the second set of points (Mx3 array)

    returns:
        distances (np.array): the distances between the points (NxM array)
    """

    p1s = p1s - (np.floor(p1s/box) * box)
    p2s = p2s - (np.floor(p2s/box) * box)
    diff = p1s[np.newaxis,:,:] - p2s[:,np.newaxis,:]
    diff = diff - (np.round(diff/box)*box)
    return np.linalg.norm(diff, axis=2)

def compute(ctx:ComputeContext, chunk_id:int):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)
    box = confs[0].box
    distances = np.empty((len(ctx.p1s), len(confs)))

    for i, conf in enumerate(confs):
        distances[:,i] = [min_image(conf.positions[p1], conf.positions[p2], box)* 0.85 for p1, p2 in zip(ctx.p1s, ctx.p2s)]
    
    return distances



def main():
    #handle commandline arguments
    #this program has no positional arguments, only flags
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Finds the ensemble of distances between any two particles in the system")
    parser.add_argument('-i', '--input', metavar='input', nargs='+', action='append', help='A topology, trajectory, and a list of particle pairs to compare.  Can call -i multiple times to plot multiple datasets.')
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The name to save the graph file to')
    parser.add_argument('-f', '--format', metavar='<histogram/trajectory/both>', nargs=1, help='Output format for the graphs.  Defaults to histogram.  Options are \"histogram\", \"trajectory\", and \"both\"')
    parser.add_argument('-d', '--data', metavar='data_file', nargs=1, help='If set, the output for the graphs will be dropped as a json to this filename for loading in oxView or your own scripts')
    parser.add_argument('-n', '--names', metavar='names', nargs='+', action='append', help='Names of the data series.  Will default to particle ids if not provided')
    parser.add_argument('-p', '--parallel', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-c', metavar='cluster', dest='cluster', action='store_const', const=True, default=False, help="Run the clusterer on each configuration's distance?")
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "matplotlib", "numpy"])

    #-i requires 4 or more arguments, the topology file of the structure, the trajectory to analyze, and any number of particle pairs to compute the distance between.
    try:
        topologies = [i[0] for i in args.input]
        trajectories = [i[1] for i in args.input]
        p1ss = [i[2::2] for i in args.input]
        p2ss = [i[3::2] for i in args.input]
        p1ss = [[int(j) for j in i] for i in p1ss]
        p2ss = [[int(j) for j in i] for i in p2ss]

    except Exception as e:
        print("ERROR:", e)
        parser.print_help()
        exit(1)
    
    #get number of distances to calculate
    n_dists = sum([len(l) for l in p1ss])

    #Make sure that the input is correctly formatted
    if(len(topologies) != len(trajectories)):
        print("ERROR: bad input arguments\nPlease supply an equal number of input and, trajectory files", file=stderr)
        exit(1)
    if len(p1ss) != len(p2ss):
        print("ERROR: bad input arguments\nPlease supply an even number of particles", file=stderr)
        exit(1)

    # Get metadata on the inputs
    top_infos = []
    traj_infos = []
    for top, traj in zip(topologies, trajectories):
        top_info, traj_info = describe(top, traj)
        top_infos.append(top_info)
        traj_infos.append(traj_info)

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

    # -p sets the number of cpus to use
    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    # how many confs we want to distribute between the processes
    ntopart = 20
    pool = Pool(ncpus)

    distances = [[] for _ in trajectories]
    for i, (traj_info, top_info, p1s, p2s) in enumerate(zip(traj_infos, top_infos, p1ss, p2ss)):
        # deduce how many chunks we have to run in parallel
        n_confs  = traj_info.nconfs 
        n_chunks = int(n_confs / ntopart +
                    (1 if n_confs % ntopart else 0))
        
        ctx = ComputeContext(traj_info, top_info, p1s, p2s, ntopart)
        distances[i] = [[None]*n_confs for _ in p1s]

        print("INFO: Working on trajectory: {}".format(traj_info.path), file=stderr)

        ## Distribute jobs to the worker processes
        print(f"Starting up {ncpus} processes for {n_chunks} chunks")
        results = [pool.apply_async(compute,(ctx,j)) for j in range(n_chunks)]
        print("All spawned, waiting for results")
 
        for j, r in enumerate(results):
            for k, d in enumerate(r.get()):
                distances[i][k][ntopart*j:ntopart*j+len(d)] = d
            print(f"finished {j+1}/{n_chunks}", end="\r")


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
        print("INFO: Defaulting to particle IDs as data series names", file=stderr)
        names = ["{}-{}".format(p1, p2) for p1, p2 in zip([i for sl in p1ss for i in sl], [i for sl in p2ss for i in sl])]
    
    # -d will dump the distances as json files for loading with the trajectories in oxView
    if args.data:
        from json import dump
        if len(trajectories) > 1:
            print("INFO: distance lists from separate trajectories are printed to separate files for oxView compatibility.  Trajectory numbers will be appended to your provided data file name.", file=stderr)
            file_names = ["{}_{}.json".format(args.data[0].strip('.json'), i) for i,_ in enumerate(trajectories)]
        else:
            file_names = [args.data[0].strip('.json')+'.json']
        names_by_traj = [['{}-{}'.format(p1, p2) for p1, p2 in zip(p1l, p2l)] for p1l, p2l in zip(p1ss, p2ss)]
        
        for file_name, ns, dist_list in zip(file_names, names_by_traj, distances):
            obj = {}
            for n, d in zip(ns, dist_list):
                obj[n] = d        
            with open(file_name, 'w+') as f:
                print("INFO: writing data to {}.  This can be opened in oxView using the Order parameter selector".format(file_name), file=stderr)
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
    [print("{}-{}\t".format(p1, p2), end='') for p1, p2 in zip([i for sl in p1ss for i in sl], [i for sl in p2ss for i in sl])]
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

    #if cluster == True:
    #    if not all([x == trajectories[0] for x in trajectories]):
    #        print("ERROR: Clustering can only be run on a single trajectory", file=stderr)
    #        exit(1)
    #    from oxDNA_analysis_tools.clustering import perform_DBSCAN
#
    #    labs = perform_DBSCAN(distances[0].T, len(distances[0][0]), trajectories[0], input_files[0], "euclidean", 12, 8)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()