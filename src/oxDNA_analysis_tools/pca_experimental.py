#!/usr/bin/env python3

#PCA.py
#Written by: Erik Poppleton
#Date: 3/6/2019
#Performs a principal component analysis on a trajectory
#the output JSON can be loaded into the viewer where it will overlay as arrows

import numpy as np 
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2, cal_confs, get_input_parameter
from sys import exit, stderr
import argparse
from json import load, dumps
try:
    from Bio.SVDSuperimposer import SVDSuperimposer
except:
    from bio.SVDSuperimposer import SVDSuperimposer
from random import randint
from oxDNA_analysis_tools.UTILS import parallelize_lorenzo_onefile
from warnings import catch_warnings, simplefilter
from os import environ, path
from oxDNA_analysis_tools.config import check_dependencies
from sklearn.decomposition import PCA

fetch_np = lambda conf: np.array([
    n.cm_pos for n in conf._nucleotides 
])

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

def get_pca(reader, align_conf, num_confs, start=None, stop=None):
    """
        Performs principal component analysis on deviations from the mean structure

        Parameters:
            reader (readers.LorenzoReader2): An active reader on the trajectory file to analyze.
            mean_structure (numpy.array): The position of each particle in the mean configuration.  A 3xN array.
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
    
    mysystem = reader._get_system(N_skip = start)
    
    deviations_matrix = np.empty((stop, (len(align_conf))*3))
    sup = SVDSuperimposer()
    confid = 0

    #for every configuration in the trajectory chunk, align it to the mean and compute positional difference for every particle
    while mysystem != False and confid < stop:
        print("-->", mysystem._time)
        mysystem.inbox()
        cur_conf = fetch_np(mysystem)
        sup.set(align_conf, cur_conf)
        sup.run()
        rot, tran = sup.get_rotran()
        #equivalent to taking the dot product of the rotation array and every vector in the deviations array
        cur_conf = np.einsum('ij, ki -> kj', rot, cur_conf) + tran
        deviations_matrix[confid] = (cur_conf-align_conf).flatten()

        confid += 1
        mysystem = reader._get_system()

    return deviations_matrix

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
    top_file = get_input_parameter(inputfile, "topology")
    if "RNA" in get_input_parameter(inputfile, "interaction_type"):
        environ["OXRNA"] = "1"
    else:
        environ["OXRNA"] = "0"
    import UTILS.base #this needs to be imported after the model type is set
    
    num_confs = cal_confs(traj_file)
    
    if mean_file.split(".")[-1] == "json":
        with open(mean_file) as file:
            align_conf = load(file)['g_mean']

    elif mean_file.split(".")[-1] == "dat":
        fetch_np = lambda conf: np.array([
            n.cm_pos for n in conf._nucleotides
        ])
        with LorenzoReader2(mean_file, top_file) as reader:
            s = reader._get_system()
            align_conf = fetch_np(s)

    cms = np.mean(align_conf, axis=0) #all structures must have the same center of mass
    align_conf -= cms 
        
    #Compute the deviations
    if not parallel:
        r = LorenzoReader2(traj_file,top_file)
        deviations_matrix = get_pca(r, align_conf, num_confs)
    
    if parallel:
        out = parallelize_lorenzo_onefile.fire_multiprocess(traj_file, top_file, get_pca, num_confs, n_cpus, align_conf)
        deviations_matrix = np.concatenate([i for i in out])
    
    #now that we have the deviations matrix we're gonna get the covariance and PCA it
    #note that in the future we might want a switch for covariance vs correlation matrix because correlation (cov/stdev so all diagonals are 1) is better for really floppy structures
    pca = PCA(n_components=3)
    pca.fit(deviations_matrix)
    transformed = pca.transform(deviations_matrix)

    #THIS IS AS FAR AS I GOT

    import matplotlib.pyplot as plt
    print("INFO: Saving scree plot to scree.png", file=stderr)
    plt.scatter(range(0, len(evalues)), evalues, s=25)
    plt.xlabel("component")
    plt.ylabel("eigenvalue")
    plt.savefig("scree.png")

    print("INFO: Creating coordinate plot from first three eigenvectors.  Saving to coordinates.png", file=stderr)
    #if you want to weight the components by their eigenvectors
    #mul = np.einsum('ij,i->ij',evectors[0:3], evalues[0:3])
    mul = evectors

    #reconstruct configurations in component space
    out = np.dot(deviations_matrix, mul).astype(float)

    #make a quick plot from the first three components
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(out[:,0], out[:,1], out[:,2], c='g', s=25)
    plt.savefig("coordinates.png")
    
    #Create an oxView overlay showing the first SUM components
    SUM = 1
    print("INFO: Change the number of eigenvalues to sum and display by modifying the SUM variable in the script.  Current value: {}".format(SUM), file=stderr)
    weighted_sum = np.zeros_like(evectors[0])
    for i in range(0, SUM): #how many eigenvalues do you want?
        weighted_sum += evalues[i]*evectors[i]

    prep_pos_for_json = lambda conf: list(
                            list(p) for p in conf
                        )
    with catch_warnings(): #this produces an annoying warning about casting complex values to real values that is not relevant
        simplefilter("ignore")
        output_vectors = weighted_sum.reshape(int(weighted_sum.shape[0]/3), 3).astype(float)
    with open(outfile, "w+") as file:
        file.write(dumps({
            "pca" : prep_pos_for_json(output_vectors)
        }))

    #If we're running clustering, feed the linear terms into the clusterer
    if cluster:
        print("INFO: Mapping configurations to component space...", file=stderr)

        #If you want to cluster on only some of the components, uncomment this
        #out = out[:,0:3]

        from clustering import perform_DBSCAN
        labs = perform_DBSCAN(out, num_confs, traj_file, inputfile, "euclidean", 12, 8)

if __name__ == '__main__':
    main()    
