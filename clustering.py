#!/usr/bin/env python3

import numpy as np
from os import environ, remove
import subprocess
import matplotlib.pyplot as plt
from output_bonds import output_bonds
from UTILS.readers import LorenzoReader2, get_input_parameter
from sklearn.cluster import DBSCAN
from sklearn import metrics
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

#Runs principal component analysis of the differences to produce orthogonal variation allowing for flattening of data to 3D
def perform_pca(differences, out_dims):
    normed_differences = differences - np.mean(differences, axis=0)
    R = np.cov(normed_differences, rowvar=False)
    evals, evecs = np.linalg.eigh(R)
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:,idx]
    evals = evals[idx] 
    evecs = evecs[:, :out_dims]
    return np.dot(evecs.T, normed_differences.T).T, evals, evecs

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
def get_centroid(differences, num_confs, labs, conf_file, inputfile):
    print("splitting clusters...")
    print ("cluster\tn\tavg_E\tE_dev\tavg_H\tH_dev\tcentroid_t")
    for cluster in (set(labs)):
        in_cluster = list(labs).count(cluster)
        masked = differences[labs == cluster].T
        in_cluster_id = np.sum(masked, axis = 1).argmin()
        centroid_id = find_element(in_cluster_id, cluster, labs)
        top_file = get_input_parameter(inputfile, "topology")
        if "RNA" in get_input_parameter(inputfile, "interaction_type"):
            environ["OXRNA"] = "1"
        else:
            environ["OXRNA"] = "0"
    
        import UTILS.base #this needs to be imported after the model type is set

        r = LorenzoReader2(conf_file, top_file)
        output = r._get_system(N_skip=centroid_id)
        filename = "centroid"+str(cluster)

        output.print_lorenzo_output(filename+".dat", filename+".top")
        
        make_heatmap(inputfile, output, filename)

        #now to get some information about each cluster
        energies = []
        H_counts = []
        confid = 0
        r1 = LorenzoReader2(conf_file, top_file)
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

#Runs a DBSCAN on the differences matrix and creates a 3D plot showing the clusters
#There is code for both animating the plot and just having an interactive 3D plot.  Comment out the one you don't want
def perform_DBSCAN(differences, num_confs, conf_file, inputfile):
    print("Running DBSCAN...")
    if (len(differences.shape) == 1):
        differences = np.abs(differences - differences.reshape(-1, 1))
    db = DBSCAN(eps=12, min_samples=10).fit(differences) #the eps and min_samples need to be determined for each structure
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print ("number of clusters:", n_clusters_)

    get_centroid(differences, num_confs, labels, conf_file, inputfile)
    #components = perform_pca(differences, 3)

    x = []
    y = []
    z = []
    for i in differences:
        x.append(i[0])
        y.append(i[1])
        z.append(i[2])

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.xlabel("component0")
    plt.ylabel("component1")
    ax.set_zlabel("component2")
    #to show the plot immediatley and interactivley
    '''a = ax.scatter(x, y, z, s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', 7))
    b = fig.colorbar(a, ax=ax)
    plt.show()
    '''
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

    return labels