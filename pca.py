#!/usr/bin/env python3

#PCA.py
#Written by: Erik Poppleton
#Date: 3/6/2019
#Performs a principal component analysis on a trajectory
#the output JSON can be loaded into the viewere where it will overlay as arrows

import numpy as np 
from UTILS.readers import LorenzoReader2, Cal_confs, get_input_parameter
from sys import exit
import argparse
from json import load, dumps
from Bio.SVDSuperimposer import SVDSuperimposer
from random import randint
from UTILS import parallelize
from warnings import catch_warnings, simplefilter
from os import environ
from compute_mean import compute_cms

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
    plt.show()

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
        #mysystem.inbox_system()
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

def change_basis(reader, align_conf, num_confs, start=None, stop=None):
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)

    linear_terms = np.empty((stop, len(evectors[1])))
    mysystem = reader._get_system(N_skip=start)
    i = 0
    sup = SVDSuperimposer()
    while mysystem != False and i < stop:
        print("-->", mysystem._time)
        mysystem.inbox_system()
        cur_conf = fetch_np(mysystem)
        sup.set(align_conf, cur_conf)
        sup.run()
        rot, tran = sup.get_rotran()
        #equivalent to taking the dot product of the rotation array and every vector in the deviations array
        cur_conf = np.einsum('ij, ki -> kj', rot, cur_conf) + tran
        cur_conf = cur_conf.flatten()
        with catch_warnings(): #this produces an annoying warning about casting complex values to real values that is not relevant
            simplefilter("ignore")
            linear_terms[i] = np.linalg.solve(evectors, cur_conf)
        i += 1
        mysystem = reader._get_system()

    return linear_terms

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates a principal component analysis of nucleotide deviations over a trajectory")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('meanfile', type=str, nargs=1, help='The mean structure .json file from compute_mean.py')
    parser.add_argument('outfile', type=str, nargs=1, help='the name of the .json file where the PCA will be written')
    parser.add_argument('-p', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")    
    
    args = parser.parse_args()
    conf_file = args.trajectory[0]
    inputfile = args.inputfile[0] 
    mean_file = args.meanfile[0]
    outfile = args.outfile[0]
    parallel = args.parallel
    if parallel:
        n_cpus = args.parallel[0]
    top_file = get_input_parameter(inputfile, "topology")
    if "RNA" in get_input_parameter(inputfile, "interaction_type"):
        environ["OXRNA"] = "1"
    else:
        environ["OXRNA"] = "0"
    import UTILS.base #this needs to be imported after the model type is set
    
    num_confs = Cal_confs(conf_file, top_file)
    
    with open(mean_file) as file:
        align_conf = load(file)['g_mean']
    cms = compute_cms(align_conf) #all structures must have the same center of mass
    align_conf -= cms 
        
    if not parallel:
        r = LorenzoReader2(conf_file,top_file)
        deviations_matrix = get_pca(r, align_conf, num_confs)
    
    if parallel:
        out = parallelize.fire_multiprocess(conf_file, top_file, get_pca, num_confs, align_conf)
        deviations_matrix = np.concatenate([i for i in out])
    
    #now that we have the deviations matrix we're gonna get the covariance and PCA it
    #note that in the future we might want a switch for covariance vs correlation matrix because correlation (cov/stdev so all diagonals are 1) is better for really floppy structures
    print("calculating covariance")
    covariance = np.cov(deviations_matrix.T)
    #make_heatmap(covariance)
    print("calculating eigenvectors")
    evalues, evectors = np.linalg.eig(covariance)
    sort = evalues.argsort()[::-1]
    evalues = evalues[sort]
    evectors = evectors[sort].T
    
    import matplotlib.pyplot as plt
    plt.scatter(range(0, len(evalues)), evalues, s=6)
    plt.xlabel("component")
    plt.ylabel("eigen value")
    plt.show()

    mul = np.einsum('ij,i->ij',evectors[0:3], evalues[0:3])
    out = np.matmul(deviations_matrix, mul.T).astype(float)
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(out[:,0], out[:,1], out[:,2], c='g', s=25)
    plt.show()
    
    weighted_sum = np.zeros_like(evectors[0])
    for i in range(0, 1): #how many eigenvalues do you want?
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

    #Now we're going to reconstruct each conf from the eigenvectors and use those weights to cluster the structures
    if not parallel:
        r = LorenzoReader2(conf_file,top_file)
        linear_terms = change_basis(r, align_conf, num_confs)

    if parallel:
        out = parallelize.fire_multiprocess(conf_file, top_file, change_basis, num_confs, n_cpus, align_conf)
        linear_terms = np.concatenate([i for i in out])

    #truncated_terms = linear_terms[:,0:3]

    from clustering import perform_DBSCAN
    labs = perform_DBSCAN(linear_terms, num_confs, conf_file, inputfile)

    