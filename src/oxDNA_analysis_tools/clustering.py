#!/usr/bin/env python3

import numpy as np
from os import environ, remove, path
from sys import stderr
import subprocess
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn import metrics
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from json import loads, dump
import codecs
from oxDNA_analysis_tools.output_bonds import output_bonds
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2, get_input_parameter

#Runs principal component analysis of the points to produce orthogonal variation allowing for flattening of data to 3D
def perform_pca(points, out_dims):
    normed_points = points - np.mean(points, axis=0)
    R = np.cov(normed_points, rowvar=False)
    evals, evecs = np.linalg.eigh(R)
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:,idx]
    evals = evals[idx] 
    evecs = evecs[:, :out_dims]
    return np.dot(evecs.T, normed_points.T).T, evals, evecs

def make_heatmap(inputfile, system, filename):
    """
    Creates a heatmap for each cluster centroid

    Parameters: 
        inputfile (str): The input file used to run the simulation.
        system (base.system): The system to make a heatmap from.
        filename (str): The file to write the heatmap to.
    """
    from oxDNA_analysis_tools.contact_map import contact_map
    m = contact_map(inputfile, system, True)
    fig, ax = plt.subplots()
    a = ax.imshow(m, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
        ylabel="nucleotide id",
        xlabel="nucleotide id")
    b = fig.colorbar(a, ax=ax)
    b.set_label("distance", rotation = 270)
    plt.savefig(filename+".png")


def find_element(n, x, array):
    """
    Finds the id of the nth time element x appears in an array.
    """
    c = 0
    for i, j in enumerate(array):
        if (j == x):
            if (c == n):
                return i
            c += 1
    return -1

def split_trajectory(traj_file, inputfile, labs, n_clusters):
    """
    Splits the trajectory into the clustered trajectories

    Parameters:
        traj_file (str): The analyzed trajectory file.
        inputfile (str): The input file used to run the analyzed simulation.
        labs (numpy.array): The cluster each point belongs to.
    """
    top_file = get_input_parameter(inputfile, "topology")

    print ("cluster\tmembers")

    #energies = []
    #H_counts = []

    for cluster in (set(labs)):
        in_cluster = list(labs).count(cluster)

        print ("{}\t{}".format(cluster, in_cluster))

        #energies.append([])
        #H_counts.append([])

        #for making trajectories of each cluster
        try:
            remove("cluster_"+str(cluster)+".dat")
        except: pass

    confid = 0
    r1 = LorenzoReader2(traj_file, top_file)
    system = r1._get_system() 
    
    print ("INFO: splitting trajectory...", file=stderr)
    print ("INFO: Will write cluster trajectories to cluster_<cluster number>.dat", file=stderr)

    while system != False:
        system.print_traj_output("cluster_"+str(labs[confid])+".dat", "/dev/null")

        ###########
        #If you want to get additional information about a cluster, add that code here
        #for example, if you want average energy and hydrogen bonds:
        '''
        energies[labs[confid]].append(0)
        H_counts[labs[confid]].append(0)
        system.map_nucleotides_to_strands()
        out = output_bonds(inputfile, system)

        for line in out.split('\n'):
            if line[0] != '#' and line[0] != '\n':
                line = line.split(" ")
                for m in line[2:9]:
                    energies[labs[confid]][-1] += float(m)
                if float(line[6]) != 0:
                    H_counts[labs[confid]][-1] += 1
        energies[labs[confid]][-1] /= len(system._nucleotides)
        '''
        ############
            
        confid += 1
        system = r1._get_system()

    #This is where you print the information about each cluster
    '''    
    avg_energy = []
    energy_dev = []
    avg_H = []
    H_dev = []
    for e, h in zip(energies, H_counts):
        avg_energy.append(np.mean(e))
        energy_dev.append(np.std(e))
        avg_H.append(np.mean(h))
        H_dev.append(np.std(h))
    print ("cluster\tmembers\tavg_E\tE_dev\tavg_H\tH_dev")
    for cluster in (set(labs)):
        in_cluster = list(labs).count(cluster)
        print ("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(cluster, in_cluster, avg_energy[cluster], energy_dev[cluster], avg_H[cluster], H_dev[cluster]))
    '''

#Calculates the centroid for each cluster
#Creates a .dat and .top file for each centroid
#also makes a heatmap for each cluster and provides some information
def get_centroid(points, metric_name, num_confs, labs, traj_file, inputfile):
    """
    Takes the output from DBSCAN and produces the trajectory and centroid from each cluster.

    Parameters:
        points (numpy.array): The points fed to the clstering algorithm.
        metric_name (str): The type of data the points represent.
        labs (numpy.array): The cluster each point belongs to.
        traj_file (str): The analyzed trajectory file.
        inputfile (str): The input file used to run the analyzed simulation.
    """
    
    print("INFO: splitting clusters...", file=stderr)
    print("INFO: Will write cluster trajectories to traj_<cluster_number>.dat", file=stderr)
    print ("cluster\tn\tavg_E\tE_dev\tavg_H\tH_dev\tcentroid_t")
    for cluster in (set(labs)):
        if metric_name == "precomputed":
            masked = points[labs == cluster]
            in_cluster_id = np.sum(masked, axis = 1).argmin()

        in_cluster = list(labs).count(cluster)
        centroid_id = find_element(in_cluster_id, cluster, labs)
        top_file = get_input_parameter(inputfile, "topology")

        r = LorenzoReader2(traj_file, top_file)
        output = r._get_system(N_skip=centroid_id)
        filename = "centroid"+str(cluster)

        output.print_lorenzo_output(filename+".dat", filename+".top")
        
        make_heatmap(inputfile, output, filename)

#Runs a DBSCAN on the points matrix and creates a 3D plot showing the clusters
#There is code for both animating the plot and just having an interactive 3D plot.  Comment out the one you don't want
def perform_DBSCAN(points, num_confs, traj_file, inputfile, metric_name, eps, min_samples):
    """
    Runs the DBSCAN algorithm using the provided analysis as positions and splits the trajectory into clusters.

    Parameters:
        points (numpy.array): The points fed to the clstering algorithm.
        num_confs (int): The number of configurations in the trajectory.
        traj_file (str): The analyzed trajectory file.
        inputfile (str): The input file used to run the analyzed simulation.
        metric_name (str): The type of data the points represent (usually either "euclidean" or "precomputed").
    
    Returns:
        labels (numpy.array): The clusterID of each configuration in the trajectory.
    """

    #run system checks
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "sklearn", "matplotlib"])
    
    print("INFO: Running DBSCAN...", file=stderr)

    #dump the input as a json file so you can iterate on eps and min_samples
    dump_file = "cluster_data.json"
    print("INFO: Serializing input data to {}".format(dump_file), file=stderr)
    print("INFO: Run just clustering.py with the serialized data to adjust clustering parameters", file=stderr)
    out = [points.tolist(), num_confs, traj_file, inputfile, metric_name]
    dump(out, codecs.open(dump_file, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

    #prepping to show the plot later
    #this only shows the first three dimensions because we assume that this is either PCA data or only a few dimensions anyway

    #components = perform_pca(points, 3)
    dimensions = []
    x = []
    dimensions.append(x)

    if points.shape[1] > 1:
        y = []
        dimensions.append(y)

    if points.shape[1] > 2:
        z = []
        dimensions.append(z)
    
    for i in points:
        for j, dim in enumerate(dimensions):
            dim.append(i[j])

    #DBSCAN parameters:
    #eps: the pairwise distance that configurations below are considered neighbors
    #min_samples: The smallest number of neighboring configurations required to start a cluster
    #metric: If the matrix fed in are points in n-dimensional space, then the metric needs to be "euclidean".
    #        If the matrix is already a square distance matrix, the metrix needs to be "precomputed".
    #the eps and min_samples need to be determined for each input based on the values of the input data
    #If you're making your own multidimensional data, you probably want to normalize your data first.
    print("INFO: Adjust clustering parameters by adding the -e and -m flags to the invocation of this script.", file=stderr)
    print("INFO: Current values: eps={}, min_samples={}".format(eps, min_samples))
    db = DBSCAN(eps=eps, min_samples=min_samples, metric=metric_name).fit(points) 
    labels = db.labels_
    
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print ("Number of clusters:", n_clusters_)

    
    print("INFO: Making cluster plot...")
    if len(dimensions) == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

    plt.xlabel("OP0")
    plt.ylabel("OP1")

    if len(dimensions) == 3:
        ax.set_zlabel("OP2")
        #to show the plot immediatley and interactivley
        '''a = ax.scatter(x, y, z, s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', 7))
        b = fig.colorbar(a, ax=ax)
        plt.show()'''
        
        #to make a video showing a rotating plot
        plot_file = "animated.mp4"
        def init():
            a = ax.scatter(x, y, z, s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', n_clusters_+1))
            fig.colorbar(a, ax=ax)
            return [fig]

        def animate(i):
            ax.view_init(elev=10., azim=i)
            return [fig]

        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=range(360), interval=20, blit=True)
        
        anim.save(plot_file, fps=30, extra_args=['-vcodec', 'libx264'])

    else:
        plot_file = "plot.png"
        if len(dimensions) == 1:
            dimensions.append(np.arange(len(dimensions[0])))
            a = ax.scatter(dimensions[1], dimensions[0], s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', n_clusters_+1))
        else:
            a = ax.scatter(dimensions[0], dimensions[1], s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', n_clusters_+1))
        b = fig.colorbar(a, ax=ax)
        plt.savefig(plot_file)
    print("INFO: Saved cluster plot to {}".format(plot_file), file=stderr)

    if metric_name == "precomputed":
        get_centroid(points, metric_name, num_confs, labels, traj_file, inputfile)

    split_trajectory(traj_file, inputfile, labels, n_clusters_)

    return labels

def main():
    import argparse
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Calculates clusters based on provided order parameters.")
    parser.add_argument('serialized_data', type=str, nargs=1, help="The json-formatted input file")
    parser.add_argument('eps', '-e', type=float, nargs=1, help="The epsilon parameter for DBSCAN")
    parser.add_argument('min_samples', '-m', type=int, nargs=1, help="The min_samples parameter for DBSCAN")
    args = parser.parse_args()
    data_file = args.serialized_data[0]
    if args.eps:
        eps = args.eps[0]
    else:
        eps = 12
    if args.min_samples:
        min_samples = args.min_samples[0]
    else:
        min_samples = 8

    #load a previously serialized dataset
    data = codecs.open(data_file, 'r', encoding='utf-8').read()
    unpacked = loads(data)
    points = np.array(unpacked[0])
    labels = perform_DBSCAN(points, unpacked[1], unpacked[2], unpacked[3], unpacked[4], eps, min_samples)

if __name__ == '__main__':
    main()
    