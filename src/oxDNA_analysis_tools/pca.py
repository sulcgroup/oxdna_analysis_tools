#!/usr/bin/env python3

#PCA.py
#Written by: Erik Poppleton
#Date: 3/6/2019
#Performs a principal component analysis on a trajectory
#the output JSON can be loaded into the viewer where it will overlay as arrows

import numpy as np 
from oxDNA_analysis_tools.UTILS.readers import ErikReader, cal_confs, get_input_parameter
from sys import exit, stderr
import argparse
from json import load, dumps
try:
    from Bio.SVDSuperimposer import SVDSuperimposer
except:
    from bio.SVDSuperimposer import SVDSuperimposer
from random import randint
from oxDNA_analysis_tools.UTILS import parallelize_erik_onefile
from warnings import catch_warnings, simplefilter
from os import environ, path
from oxDNA_analysis_tools.config import check_dependencies

def make_heatmap(covariance):
    """
        Produces a heatmat plot of the covariance between every particle

        Parameters:
            covariance (numpy.array): The covariance matrix of particle positions

        Displays:
            A matplotlib imshow plot of the covariance
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    a = ax.imshow(covariance, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
       ylabel="nucleotide id*3",
       xlabel="nucleotide id*3")
    b = fig.colorbar(a, ax=ax)
    b.set_label("covariance", rotation = 270)
    plt.savefig("heatmap.png")

def get_cov(reader, align_conf, num_confs, start=None, stop=None):
    """
        Performs principal component analysis on deviations from the mean structure

        Parameters:
            reader (readers.ErikReader): An active reader on the trajectory file to analyze.
            align_conf (numpy.array): The position of each particle in the mean configuration.  A 3xN array.
            num_confs (int): The number of configurations in the reader.  
            <optional> start (int): The starting configuration ID to begin averaging at.  Used if parallel.
            <optional> stop (int): The configuration ID on which to end the averaging.  Used if parallel.

        Returns:
            deviations_marix (numpy.array): The difference in position from the mean for each configuration.
    """
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)
    
    mysystem = reader.read(n_skip = start)
    
    covariation_matrix = np.zeros((len(mysystem.positions)*3, len(mysystem.positions)*3))
    sup = SVDSuperimposer()
    confid = 0

    #for every configuration in the trajectory chunk, align it to the mean and compute positional difference for every particle
    while mysystem != False and confid < stop:
        print("-->", "frame", confid, "time={}".format(mysystem.time))
        mysystem.inbox()
        cur_conf = mysystem.positions
        sup.set(align_conf, cur_conf)
        sup.run()
        rot, tran = sup.get_rotran()
        #equivalent to taking the dot product of the rotation array and every vector in the deviations array
        cur_conf = np.einsum('ij, ki -> kj', rot, cur_conf) + tran
        difference_matrix = (cur_conf - align_conf).flatten()
        covariation_matrix += np.einsum('i,j -> ij', difference_matrix, difference_matrix)

        confid += 1
        mysystem = reader.read()

    return covariation_matrix

def change_basis(reader, align_conf, components, num_confs, start=None, stop=None):
    """
    Transforms each configuration in a trajectory into a point in principal component space

    Parameters:
        reader (readers.ErikReader): An active reader on the trajectory file to analyze.
        align_conf (numpy.array): The position of each particle in the mean configuration.  A 3xN array.
        components (numpy.array): The principal components of the trajectory.  A 3*Nx3*N array.
        num_confs (int): The number of configurations in the reader.  
        <optional> start (int): The starting configuration ID to begin averaging at.  Used if parallel.
        <optional> stop (int): The configuration ID on which to end the averaging.  Used if parallel.

    Returns:
        coordinates (numpy.array): The positions of each frame of the trajectory in principal component space.
    """

    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)
    
    mysystem = reader.read(n_skip = start)

    coordinates = np.empty((stop, len(mysystem.positions)*3))
    coordinates2 = np.empty((stop, len(mysystem.positions)*3))
    sup = SVDSuperimposer()
    confid = 0

    while mysystem != False and confid < stop:
        print("-->", "frame", confid, "time={}".format(mysystem.time))
        mysystem.inbox()
        cur_conf = mysystem.positions
        sup.set(align_conf, cur_conf)
        sup.run()
        rot, tran = sup.get_rotran()
        #equivalent to taking the dot product of the rotation array and every vector in the deviations array
        cur_conf = np.einsum('ij, ki -> kj', rot, cur_conf) + tran
        coordinates[confid] = np.dot(components, cur_conf.flatten())

        confid += 1
        mysystem = reader.read()

    return(coordinates)

def main():
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Calculates a principal component analysis of nucleotide deviations over a trajectory")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('meanfile', type=str, nargs=1, help='The mean structure .json file from compute_mean.py')
    parser.add_argument('outfile', type=str, nargs=1, help='the name of the .json file where the PCA will be written')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")    
    parser.add_argument('-c', metavar='cluster', dest='cluster', action='store_const', const=True, default=False, help="Run the clusterer on each configuration's position in PCA space?")
    args = parser.parse_args()

    check_dependencies(["python", "numpy", "Bio"])

    traj_file = args.trajectory[0]
    inputfile = args.inputfile[0] 
    mean_file = args.meanfile[0]
    outfile = args.outfile[0]
    parallel = args.parallel
    if parallel:
        n_cpus = args.parallel[0]
    #-c makes it run the clusterer on the output
    cluster = args.cluster

    num_confs = cal_confs(traj_file)
    
    if mean_file.split(".")[-1] == "json":
        with open(mean_file) as file:
            align_conf = load(file)['g_mean']

    elif mean_file.split(".")[-1] == "dat" or mean_file.split(".")[-1] == "conf" or mean_file.split(".")[-1] == "oxdna":
        with ErikReader(mean_file) as reader:
            align_conf = reader.read().positions
    else:
        print("ERROR: {} is an unrecognized file type. \nThe mean structure must either be provided as an oxDNA configuration file with the extension .dat, .conf or .oxdna or as the .json file produced by compute_mean.py.", file=stderr)
        exit(1)

    cms = np.mean(align_conf, axis=0) #all structures must have the same center of mass
    align_conf -= cms 
        
    #Compute the deviations
    if not parallel:
        r = ErikReader(traj_file)
        covariation_matrix = get_cov(r, align_conf, num_confs)
    
    if parallel:
        out = parallelize_erik_onefile.fire_multiprocess(traj_file, get_cov, num_confs, n_cpus, align_conf)
        covariation_matrix = np.sum([i for i in out], axis=0)

    covariation_matrix /= (num_confs-1)

    #now that we have the covatiation matrix we're going to use eigendecomposition to get the principal components.
    #make_heatmap(covariance)
    print("INFO: calculating eigenvectors", file=stderr)
    evalues, evectors = np.linalg.eig(covariation_matrix) #these eigenvalues are already sorted
    evectors = evectors.T #vectors come out as the columns of the array
    print("INFO: eigenvectors calculated", file=stderr)
    
    import matplotlib.pyplot as plt
    print("INFO: Saving scree plot to scree.png", file=stderr)
    plt.scatter(range(0, len(evalues)), evalues, s=25)
    plt.xlabel("component")
    plt.ylabel("eigenvalue")
    plt.savefig("scree.png")

    total = sum(evalues)
    running = 0
    i = 0
    while running < 0.9:
        running += (evalues[i] / total)
        i += 1
    
    print("90% of the variance is found in the first {} components".format(i))


    
    #if you want to weight the components by their eigenvectors
    #mul = np.einsum('ij,i->ij',evectors, evalues)
    mul = evectors

    #reconstruct configurations in component space
    #because we donlist't save the difference matrix, this involves running through the whole trajectory again
    if not parallel:
        r = ErikReader(traj_file)
        coordinates = change_basis(r, align_conf, mul, num_confs)
    if parallel:
        out = parallelize_erik_onefile.fire_multiprocess(traj_file, change_basis, num_confs, n_cpus, align_conf, mul)
        coordinates = np.concatenate([i for i in out])

    #make a quick plot from the first three components
    print("INFO: Creating coordinate plot from first three eigenvectors.  Saving to coordinates.png", file=stderr)
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(coordinates[:,0], coordinates[:,1], coordinates[:,2], c='g', s=25)
    plt.savefig("coordinates.png")
    
    #Create an oxView overlays for the first N components
    N = 3
    prep_pos_for_json = lambda conf: list(
                        list(p) for p in conf
                    )
    print("INFO: Change the number of eigenvalues to sum and display by modifying the N variable in the script.  Current value: {}".format(N), file=stderr)
    for i in range(0, N): #how many eigenvalues do you want?
        try:
            if outfile.split(".")[1] != "json":
                raise Exception
            f = outfile.split(".")[0] + str(i) + "." + outfile.split(".")[1]
        except:
            print("ERROR: oxView overlays must have a '.json' extension.  No overlays will be produced", file=stderr)
            break
        out = np.sqrt(evalues[i])*evectors[i]

        with catch_warnings(): #this produces an annoying warning about casting complex values to real values that is not relevant
            simplefilter("ignore")
            output_vectors = out.reshape(int(out.shape[0]/3), 3).astype(float)
        
        with open(f, "w+") as file:
            file.write(dumps({
                "pca" : prep_pos_for_json(output_vectors)
            }))

    #If we're running clustering, feed the linear terms into the clusterer
    if cluster:
        print("INFO: Mapping configurations to component space...", file=stderr)

        #If you want to cluster on only some of the components, uncomment this
        #out = out[:,0:3]

        from oxDNA_analysis_tools.clustering import perform_DBSCAN
        labs = perform_DBSCAN(coordinates, num_confs, traj_file, inputfile, "euclidean", 12, 8)

if __name__ == '__main__':
    main()    
