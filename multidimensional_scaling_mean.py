#!/usr/bin/env python3
#Multidimensional_scaling_mean
#Written by: Erik Poppleton
#Date: 2/28/19
#Python 3
#Computes the RMSD of the contact map for a structure.  The average structure is determined
#by Scikit.learn's MDS algorithm, then subtracts the contact map of each individual structure from the man
#This is used to compute a per-nucleotide deviation in the contact map, which can be visualized with oxView

import numpy as np
from UTILS.readers import LorenzoReader2, Cal_confs, get_input_parameter
from sys import exit
from json import dumps, loads
from contact_map import contact_map
import argparse
from UTILS import parallelize
from os import environ

def make_heatmap(contact_map):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    a = ax.imshow(contact_map, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
       ylabel="nucleotide id",
       xlabel="nucleotide id")
    b = fig.colorbar(a, ax=ax)
    b.set_label("distance", rotation = 270)
    plt.show()

def get_mean(reader, num_confs, start=None, stop=None):
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)

    mysystem = reader._get_system(N_skip = start)
    #mysystem.inbox_system()
    cartesian_distances = np.zeros((mysystem.N, mysystem.N))
    confid = 0

    while mysystem != False and confid < stop:
        print("-->", mysystem._time)
        #mysystem.inbox_system()
        cartesian_distances += contact_map(inputfile, mysystem, True)

        confid += 1
        mysystem = reader._get_system()

    return (cartesian_distances)

def get_devs(reader, masked_mean, num_confs, start=None, stop=None):
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)

    mysystem = reader._get_system(N_skip = start)
    #mysystem.inbox_system()
    devs = np.zeros((mysystem.N, mysystem.N))
    confid = 0

    #now that we have a mean structure, we need to compute local deviations...
    while mysystem != False and confid < stop:
        print("-->", confid)
        #mysystem.inbox_system()

        cartesian_distances = contact_map(inputfile, mysystem, True)
        masked_conf = np.ma.masked_array(cartesian_distances, ~(cartesian_distances < cutoff_distance))

        #gonna fill in the masked values with the cutoff distance for now.
        #not sure what the weight on a difference in what is under the distance should be...
        masked_conf = np.ma.filled(masked_conf, cutoff_distance)
        masked_mean = np.ma.filled(masked_mean, cutoff_distance)

        diff = masked_conf - masked_mean
        diff = np.square(diff)
        devs += diff

        confid += 1
        mysystem = reader._get_system()

    return devs

if __name__ == "__main__":
    #2 seems like a good number, at 2.5 you start to see the hard edges caused by end-loops and see some loop interactions
    cutoff_distance = 2.5

    parser = argparse.ArgumentParser(description="Calculate molecular contacts, and assembles an average set of contacts based on MDS")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('meanfile', type=str, nargs=1, help='the name of the .dat file where the mean will be written')
    parser.add_argument('devfile', type=str, nargs=1, help='the name of the .json file where the devs will be written')
    parser.add_argument('-p', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")

    args = parser.parse_args()
    conf_file = args.trajectory[0]
    inputfile = args.inputfile[0]
    meanfile = args.meanfile[0]
    devfile = args.devfile[0]
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

    print("computing mean structure...")

    if not parallel:
        r = LorenzoReader2(conf_file,top_file)
        cartesian_distances = get_mean(r, num_confs)
        mean_distance_map = cartesian_distances * (1/(num_confs))

    if parallel:
        out = parallelize.fire_multiprocess(conf_file, top_file, get_mean, num_confs, n_cpus)
        cartesian_distances = np.sum(np.array([i for i in out]), axis=0)


    mean_distance_map = cartesian_distances * (1/(num_confs))

    r = LorenzoReader2(conf_file,top_file)
    out_conf = r._get_system()

    #make heatmap of the summed distances
    #make_heatmap(mean_distance_map)

    #compute the mean structure using multidimensional scaling on the distances
    output_system = out_conf
    init = np.array([p.cm_pos for p in out_conf._nucleotides])
    from sklearn.manifold import MDS
    #from sklearn.manifold import LocallyLinearEmbedding
    #from megaman.geometry import Geometry
    #from scipy.sparse import csr_matrix
    masked_mean = np.ma.masked_array(mean_distance_map, ~(mean_distance_map < cutoff_distance))
    f = open('test_dist.nmr', 'w+')
    for i, line in enumerate(masked_mean):
        for j, dist in enumerate(line):
            if dist != "--" and dist != 0 and i < j:
                if j%2 == 0:
                    f.write("{}\t{}\t1\t1\t{}\t{}\tn\tn\tn\tn\n".format(i+1, j+1, dist, dist))
                else:
                    f.write("{}\t{}\t1\t1\t{}\t{}\tn\tn\tn\tn\n".format(j+1, i+1, dist, dist))

        
    #super_cutoff_ids = mean_distance_map > cutoff_distance
    #mean_distance_map[super_cutoff_ids] = 0
    #sparse_map = csr_matrix(mean_distance_map)
    mds = MDS(n_components=3, metric=True, max_iter=3000, eps=1e-12, dissimilarity="precomputed", n_jobs=1, n_init=1)
    #lle = LocallyLinearEmbedding(n_neighbors=5, n_components=3, eigen_solver='arpack', max_iter=3000)
    print("fitting local distance data")

    #geom = Geometry()
    #geom = Geometry(adjacency_kwds={'radius':cutoff_distance})#, laplacian_kwds={'scaling_epps':cutoff_distance})
    #geom.set_data_matrix(masked_mean)
    #geom.set_adjacency_matrix(masked_mean)
    #from megaman.embedding import LocallyLinearEmbedding
    #lle = LocallyLinearEmbedding(n_components=3, eigen_solver='arpack', geom=geom)
    #out_coords = lle.fit_transform(masked_mean, input_type='adjacency')

    out_coords = mds.fit_transform(masked_mean)#, init=init)
    #out_coords = lle.fit_transform(masked_mean)
    for i, n in enumerate(output_system._nucleotides):
        n.cm_pos = out_coords[i]
        n._a1 = np.array([0,0,0])
        n._a3 = np.array([0,0,0]) #since the orientation vectors are all 0, this cannot be used in a simulation, but the viewer will handle it

    #Write the mean structure out as a new .dat and .top pair
    output_system.print_lorenzo_output("{}.dat".format(meanfile), "{}.top".format(meanfile))
    print("wrote files: {}.dat, {}.top".format(meanfile, meanfile))

    print("computing per-nucleotide deviations")
    if not parallel:
        r = LorenzoReader2(conf_file,top_file)
        devs = get_devs(r, masked_mean, num_confs)

    if parallel:
        out = parallelize.fire_multiprocess(conf_file, top_file, get_devs, num_confs, n_cpus, masked_mean)
        devs = np.sum(np.array([i for i in out]), axis=0)

    devs = np.ma.masked_array(devs, ~(devs != 0.0)) #mask all the 0s so they don't contribute to the mean
    devs *= (1/num_confs)
    devs = np.mean(devs, axis=0)
    devs = np.sqrt(devs)
    with open(devfile+".json", "w") as file:
        file.write(
            dumps({"contact deviation" : list(devs)})
        )
    print("wrote file {}.json".format(devfile))
