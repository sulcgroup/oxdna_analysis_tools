#!/usr/bin/env python3

import numpy as np
from os import environ, remove
from sys import stderr
import subprocess
import matplotlib.pyplot as plt
from output_bonds import output_bonds
from UTILS.readers import LorenzoReader2, get_input_parameter
from sklearn.cluster import DBSCAN
from sklearn import metrics
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

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

#creates a heatmap for each cluster centroid
def make_heatmap(inputfile, system, filename):
    from contact_map import contact_map
    m = contact_map(inputfile, system, True)
    fig, ax = plt.subplots()
    a = ax.imshow(m, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
        ylabel="nucleotide id",
        xlabel="nucleotide id")
    b = fig.colorbar(a, ax=ax)
    b.set_label("distance", rotation = 270)
    plt.savefig(filename+".png")


#finds the id of the nth time element x appears in an array
def find_element(n, x, array):
    c = 0
    for i, j in enumerate(array):
        if (j == x):
            if (c == n):
                return i
            c += 1
    return -1


#Calculates the centroid for each cluster
#Creates a .dat and .top file for each centroid
#also makes a heatmap for each cluster and provides some information
def get_centroid(points, num_confs, labs, traj_file, inputfile):
    print("INFO: splitting clusters...", file=stderr)
    print ("cluster\tn\tavg_E\tE_dev\tavg_H\tH_dev\tcentroid_t")
    for cluster in (set(labs)):
        in_cluster = list(labs).count(cluster)
        masked = points[labs == cluster].T
        in_cluster_id = np.sum(masked, axis = 1).argmin()
        centroid_id = find_element(in_cluster_id, cluster, labs)
        top_file = get_input_parameter(inputfile, "topology")
        if "RNA" in get_input_parameter(inputfile, "interaction_type"):
            environ["OXRNA"] = "1"
        else:
            environ["OXRNA"] = "0"
    
        import UTILS.base #this needs to be imported after the model type is set

        r = LorenzoReader2(traj_file, top_file)
        output = r._get_system(N_skip=centroid_id)
        filename = "centroid"+str(cluster)

        output.print_lorenzo_output(filename+".dat", filename+".top")
        
        make_heatmap(inputfile, output, filename)

        #now to get some information about each cluster
        energies = []
        H_counts = []
        confid = 0
        r1 = LorenzoReader2(traj_file, top_file)
        system = r1._get_system()

        #for making trajectories of each cluster
        trajfile_name = "traj_"+str(cluster)+".dat"
        try:
            remove(trajfile_name)
        except: pass
        #traj =  open(trajfile_name, 'w+')
        while system != False:
            if labs[confid] == cluster:
                energies.append(0)
                H_counts.append(0)
                system.print_traj_output(trajfile_name, "/dev/null")
                #system.map_nucleotides_to_strands()
                #out = output_bonds(inputfile, system)
#
                #for line in out.split('\n'):
                #    if line[0] != '#' and line[0] != '\n':
                #        line = line.split(" ")
                #        for m in line[2:9]:
                #            energies[-1] += float(m)
                #        if float(line[6]) != 0:
                #            H_counts[-1] += 1
                #energies[-1] /= len(output._nucleotides)

            confid += 1
            system = r1._get_system()

        avg_energy = np.mean(energies)
        energy_dev = np.std(energies)

        avg_H = np.mean(H_counts)
        H_dev = np.std(H_counts)            
        print ("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{}".format(cluster, in_cluster, avg_energy, energy_dev, avg_H, H_dev, output._time))

#Runs a DBSCAN on the points matrix and creates a 3D plot showing the clusters
#There is code for both animating the plot and just having an interactive 3D plot.  Comment out the one you don't want
def perform_DBSCAN(points, num_confs, traj_file, inputfile, metric_name):
    print("INFO: Running DBSCAN...", file=stderr)
    EPS=3
    MIN_SAMPLES=8

    #prepping to show the plot later
    #components = perform_pca(points, 3)
    dimensions = []
    x = []
    dimensions.append(x)

    if len(points.shape) > 1:
        y = []
        dimensions.append(y)

    if len(points.shape) > 2:
        z = []
        dimensions.append(z)
    
    for i in points:
        i = [i]
        for j, dim in enumerate(dimensions):
            dim.append(i[j])

    if len(points.shape) == 1:
        points = points.reshape(-1, 1)

    #DBSCAN parameters:
    #eps: the pairwise distance that configurations below are considered neighbors
    #min_samples: The smallest number of neighboring configurations required to start a cluster
    #metric: If the matrix fed in are points in n-dimensional space, then the metric needs to be "euclidean".
    #If the matrix is already a square distance matrix, the metrix needs to be "precomputed".
    #the eps and min_samples need to be determined for each structure
    print("INFO: Adjust clustering parameters by modifying the 'EPS' and 'MIN_SAMPLES' values in the script.", file=stderr)
    print("INFO: Current values: eps={}, min_samples={}".format(EPS, MIN_SAMPLES))
    db = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES, metric=metric_name).fit(points) 
    labels = db.labels_
    
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print ("Number of clusters:", n_clusters_)

    

    fig = plt.figure()
    if len(dimensions) == 3:
        ax = fig.gca(projection='3d')
        ax.set_zlabel("component2")
    else:
        ax = fig.add_subplot(1, 1, 1)

    plt.xlabel("component0")
    plt.ylabel("component1")
    
    #to show the plot immediatley and interactivley
    '''a = ax.scatter(x, y, z, s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', 7))
    b = fig.colorbar(a, ax=ax)
    plt.show()
    '''

    if len(dimensions) == 3:
        #to make a video showing a rotating plot
        def init():
            a = ax.scatter(x, y, z, s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', n_clusters_+1))
            fig.colorbar(a, ax=ax)
            return fig,

        def animate(i):
            ax.view_init(elev=10., azim=i)
            return fig,

        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20, blit=True)

        anim.save('animated.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    else:
        if len(dimensions) == 1:
            dimensions.append(np.arange(len(dimensions[0])))
            a = ax.scatter(dimensions[1], dimensions[0], s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', n_clusters_+1))
        else:
            a = ax.scatter(dimensions[0], dimensions[1], s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', n_clusters_+1))
        b = fig.colorbar(a, ax=ax)
        plt.show()

    get_centroid(points, num_confs, labels, traj_file, inputfile)

    return labels