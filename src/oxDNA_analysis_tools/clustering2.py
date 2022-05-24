import numpy as np
import argparse
from os import environ, remove, path
from sys import stderr
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn import metrics
from matplotlib import animation
from json import dump, load
from collections import namedtuple
from oxDNA_analysis_tools.config import check_dependencies
from oxDNA_analysis_tools.UTILS.data_structures import TrajInfo, TopInfo
from oxDNA_analysis_tools.UTILS.RyeReader import linear_read, conf_to_str

#not tested, probably works.

def split_trajectory(traj_info, top_info, labs):
    """
    Splits the trajectory into the clustered trajectories

    Parameters:
        traj_info (TrajInfo): Metadata on the trajectory file
        top_info (TopInfo): Metadata on the topology file
        labs (numpy.array): The cluster each configuration belongs to.
    """

    #How many are in each cluster?
    print ("cluster\tmembers")

    slabs = set(labs)

    for cluster in slabs:
        in_cluster = list(labs).count(cluster)
        print ("{}\t{}".format(cluster, in_cluster))

        #Clear old trajectory files
        try:
            remove("cluster_"+str(cluster)+".dat")
        except: pass

    print ("INFO: splitting trajectory...", file=stderr)
    print ("INFO: Trajectories for each cluster will be written to cluster_<cluster number>.dat", file=stderr)

    files = [open("cluster_"+str(cluster)+".dat", 'w+') for cluster in slabs]
    i = 0

    for chunk in linear_read(traj_info, top_info, ntopart=20):
        for conf in chunk:
            files[labs[i]].write(conf_to_str(conf))
            i += 1

def perform_DBSCAN(traj_info:TrajInfo, top_info:TopInfo, op:np.array, metric:str, eps:float, min_samples:int):
    check_dependencies(["python", "sklearn", "matplotlib"])
    
    #dump the input as a json file so you can iterate on eps and min_samples
    dump_file = "cluster_data.json"
    print("INFO: Serializing input data to {}".format(dump_file), file=stderr)
    print("INFO: Run  `oat clustering {} -e<eps> -m<min_samples>`  to adjust clustering parameters".format(dump_file), file=stderr)
    out = [op.tolist(), traj_info, top_info, metric]
    dump(out, open(dump_file, 'w+'), separators=(',', ':'), sort_keys=True, indent=4)

    # Prepping a plot of the first 3 dimensions of the provided op
    dimensions = []
    x = []
    dimensions.append(x)

    if op.shape[1] > 1:
        y = []
        dimensions.append(y)

    if op.shape[1] > 2:
        z = []
        dimensions.append(z)
    
    for i in op:
        for j, dim in enumerate(dimensions):
            dim.append(i[j])
    
    print("INFO: Running DBSCAN...", file=stderr)

    #DBSCAN parameters:
    #eps: the pairwise distance that configurations below are considered neighbors
    #min_samples: The smallest number of neighboring configurations required to start a cluster
    #metric: If the matrix fed in are points in n-dimensional space, then the metric needs to be "euclidean".
    #        If the matrix is already a square distance matrix, the metrix needs to be "precomputed".
    #the eps and min_samples need to be determined for each input based on the values of the input data
    #If you're making your own multidimensional data, you probably want to normalize your data first.
    print("INFO: Adjust clustering parameters by adding the -e and -m flags to the invocation of this script.", file=stderr)
    print("INFO: Current values: eps={}, min_samples={}".format(eps, min_samples))
    db = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(op) 
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

    split_trajectory(traj_info, top_info, labels, n_clusters_)

    return labels

def main():
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
    with open(data_file, 'r') as f:
        data = load(f)
    points = np.array(data[0])
    labels = perform_DBSCAN(points, data[1], data[2], data[3], data[4], eps, min_samples)

if __name__ == '__main__':
    main()