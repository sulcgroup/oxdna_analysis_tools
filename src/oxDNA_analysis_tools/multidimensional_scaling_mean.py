#!/usr/bin/env python3
#Multidimensional_scaling_mean
#Written by: Erik Poppleton
#Date: 2/28/19
#Python 3
#Computes the RMSD of the contact map for a structure.  The average structure is determined
#by Scikit.learn's MDS algorithm, then subtracts the contact map of each individual structure from the man
#This is used to compute a per-nucleotide deviation in the contact map, which can be visualized with oxView

import numpy as np
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2, cal_confs, get_input_parameter
from sys import exit, stderr
from json import dumps, loads
from oxDNA_analysis_tools.contact_map import contact_map
import argparse
from oxDNA_analysis_tools.UTILS import parallelize_lorenzo_onefile
from os import environ, path

def make_heatmap(contact_map):
    """
    Convert a matrix of contact distances to a visual contact map.

    Parameters:
        contact_map (numpy.array): An array of all pairwise distances between nucleotides.
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    a = ax.imshow(contact_map, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
       ylabel="nucleotide id",
       xlabel="nucleotide id")
    b = fig.colorbar(a, ax=ax)
    b.set_label("distance", rotation = 270)
    plt.show()

def get_mean(reader, inputfile, num_confs, start=None, stop=None):
    """
    Computes the mean distance between every pair of nucleotides.

    Parameters:
        reader (readers.LorenzoReader2): An active reader on the trajectory file to process.
        num_confs (int): The number of configurations in the reader.  
        <optional> start (int): The starting configuration ID to begin averaging at.  Used if parallel.
        <optional> stop (int): The configuration ID on which to end the averaging.  Used if parallel.
    
    Returns:
        cartesian_distances (numpy.array): A matrix containing all pairwise distances between nucleotides.
    """
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)

    mysystem = reader._get_system(N_skip = start)
    cartesian_distances = np.zeros((mysystem.N, mysystem.N))
    confid = 0

    while mysystem != False and confid < stop:
        print("Frame:", confid, "Time:", mysystem._time)
        cartesian_distances += contact_map(inputfile, mysystem, True)

        confid += 1
        mysystem = reader._get_system()

    return (cartesian_distances)

def get_devs(reader, masked_mean, inputfile, cutoff_distance, num_confs, start=None, stop=None):
    """
    Computes the RMSD in each particle's distance to other particles in its neighborhood.

    Parameters:
        reader (readers.LorenzoReader2): An active reader on the trajectory file to process.
        num_confs (int): The number of configurations in the reader.
        masked_mean (numpy.ma.masked_array): A masked array containing average distances to other particles in each nucleotides's local neighborhood.
        <optional> start (int): The starting configuration ID to begin averaging at.  Used if parallel.
        <optional> stop (int): The configuration ID on which to end the averaging.  Used if parallel.
    
    Returns:
        devs (np.array): The per-nucleotide RMSD in the distance to neighboring nucleotides.
    """
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)

    mysystem = reader._get_system(N_skip = start)
    devs = np.zeros((mysystem.N, mysystem.N))
    confid = 0

    #now that we have a mean structure, we need to compute local deviations...
    while mysystem != False and confid < stop:
        print("Frame:", confid, "Time:", mysystem._time)

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

def main():
    #at 2.5 you start to see the hard edges caused by end-loops and see some loop interactions
    cutoff_distance = 2.5

    #get commandline arguments
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Calculate molecular contacts, and assembles an average set of contacts based on MDS")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('meanfile', type=str, nargs=1, help='the name of the .dat file where the mean will be written')
    parser.add_argument('devfile', type=str, nargs=1, help='the name of the .json file where the devs will be written')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")

    #process commandline arguments
    args = parser.parse_args()
    traj_file = args.trajectory[0]
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

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    #get the number of configurations in the trajectory
    num_confs = cal_confs(traj_file)

    #Get the mean distance to all other particles
    if not parallel:
        print("INFO: Computing interparticle distances of {} configurations using 1 core.".format(num_confs), file=stderr)
        r = LorenzoReader2(traj_file,top_file)
        cartesian_distances = get_mean(r, inputfile, num_confs)
        mean_distance_map = cartesian_distances * (1/(num_confs))

    if parallel:
        print("INFO: Computing interparticle distances of {} configurations using {} cores.".format(num_confs, n_cpus), file=stderr)
        out = parallelize_lorenzo_onefile.fire_multiprocess(traj_file, top_file, get_mean, num_confs, n_cpus, inputfile)
        cartesian_distances = np.sum(np.array([i for i in out]), axis=0)

    mean_distance_map = cartesian_distances * (1/(num_confs))

    #Making a new configuration file from scratch is hard, so we're just going to read in one and then overwrite the positional information
    r = LorenzoReader2(traj_file,top_file)
    output_system = r._get_system()

    #make heatmap of the summed distances
    #make_heatmap(mean_distance_map)
    
    masked_mean = np.ma.masked_array(mean_distance_map, ~(mean_distance_map < cutoff_distance))
    
    #I tried to use DGSOL to analytically solve this, but origamis were too big
    #f = open('test_dist.nmr', 'w+')
    #for i, line in enumerate(masked_mean):
    #    for j, dist in enumerate(line):
    #        if dist != "--" and dist != 0 and i < j:
    #            if j%2 == 0:
    #                f.write("{}\t{}\t1\t1\t{}\t{}\tn\tn\tn\tn\n".format(i+1, j+1, dist, dist))
    #            else:
    #                f.write("{}\t{}\t1\t1\t{}\t{}\tn\tn\tn\tn\n".format(j+1, i+1, dist, dist))

        
    #super_cutoff_ids = mean_distance_map > cutoff_distance
    #mean_distance_map[super_cutoff_ids] = 0
    #sparse_map = csr_matrix(mean_distance_map)
    
    print("INFO: fitting local distance data", file=stderr)

    #Many embedding algorithms were tried...

    #from sklearn.manifold import LocallyLinearEmbedding
    #from megaman.geometry import Geometry
    #from scipy.sparse import csr_matrix

    #geom = Geometry()
    #geom = Geometry(adjacency_kwds={'radius':cutoff_distance})#, laplacian_kwds={'scaling_epps':cutoff_distance})
    #geom.set_data_matrix(masked_mean)
    #geom.set_adjacency_matrix(masked_mean)
    #from megaman.embedding import LocallyLinearEmbedding
    #lle = LocallyLinearEmbedding(n_neighbors=5, n_components=3, eigen_solver='arpack', max_iter=3000)
    #lle = LocallyLinearEmbedding(n_components=3, eigen_solver='arpack', geom=geom)
    #out_coords = lle.fit_transform(masked_mean, input_type='adjacency')
    #out_coords = lle.fit_transform(masked_mean)
    #init = np.array([p.cm_pos for p in out_conf._nucleotides])

    #Run multidimensional scaling on the average distances to find average positions
    from sklearn.manifold import MDS
    mds = MDS(n_components=3, metric=True, max_iter=3000, eps=1e-12, dissimilarity="precomputed", n_jobs=1, n_init=1)
    out_coords = mds.fit_transform(masked_mean)#, init=init) #this one worked best
    
    #Overwrite the system we made earlier with the coordinates calculated via MDS
    for i, n in enumerate(output_system._nucleotides):
        n.cm_pos = out_coords[i]
        n._a1 = np.array([0,0,0])
        n._a3 = np.array([0,0,0]) #since the orientation vectors are all 0, this cannot be used in a simulation, but the viewer will handle it

    #Write the mean structure out as a new .dat and .top pair
    output_system.print_lorenzo_output("{}.dat".format(meanfile), "{}.top".format(meanfile))
    print("INFO: wrote output files: {}.dat, {}.top".format(meanfile, meanfile), file=stderr)

    #Loop through the trajectory again and calculate deviations from the average distances
    print("INFO: Computing distance deviations of {} configurations using 1 core.".format(num_confs), file=stderr)
    if not parallel:
        r = LorenzoReader2(traj_file,top_file)
        devs = get_devs(r, masked_mean, inputfile, cutoff_distance, num_confs)

    if parallel:
        print("INFO: Computing distance deviations of {} configurations using {} cores.".format(num_confs, n_cpus), file=stderr)
        out = parallelize_lorenzo_onefile.fire_multiprocess(traj_file, top_file, get_devs, num_confs, n_cpus, masked_mean, inputfile, cutoff_distance)
        devs = np.sum(np.array([i for i in out]), axis=0)

    #Dump the deviations to an oxView overlay file
    devs = np.ma.masked_array(devs, ~(devs != 0.0)) #mask all the 0s so they don't contribute to the mean
    devs *= (1/num_confs)
    devs = np.mean(devs, axis=0)
    devs = np.sqrt(devs)
    with open(devfile+".json", "w") as file:
        file.write(
            dumps({"contact deviation" : list(devs)})
        )
    print("INFO: wrote file {}.json".format(devfile), file=stderr)


if __name__ == '__main__':
    main()