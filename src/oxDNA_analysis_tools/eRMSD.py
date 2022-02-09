#!/usr/bin/env python3

#A contact-map based clustering and flexibility quantification algorithm
#Created by: Erik Poppleton
#Date: 11/13/18
#Python2
#Takes topology and trajectory and calculates contact mapping between all nucleotides and calculates deviation score
#between snapshots

import numpy as np
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2, cal_confs, get_input_parameter
from os import environ, path
import subprocess
import pickle
import matplotlib.pyplot as plt
from oxDNA_analysis_tools.UTILS.all_vectors import all_vectors
import argparse
from oxDNA_analysis_tools.clustering import perform_DBSCAN
from oxDNA_analysis_tools.UTILS import parallelize_lorenzo_onefile

#a matrix of vectors in local cylindrical coordinates
def calc_matrix(system, inputfile):
    interaction_sites = np.empty([len(system._nucleotides),4, 3]) #length of nucleotide list, 3-element position, 3-element a1 vector, 3-element a2 vector 3-element a3 vector
    for i,n in enumerate(system._nucleotides):
        interaction_sites[i] = np.array([n.cm_pos, n._a1, n._a2, n._a3]) #note that this uses the CMpos, not the interaction site
        #Because I don't know whether you're feeding me DNA or RNA without asking

    #calculate the interaction network    
    #from warnings import catch_warnings, simplefilter
    Vvec = all_vectors(inputfile, system, True)
    Vtransform_matrix = np.linalg.inv(np.array(interaction_sites[:, 1:].transpose(0, 2, 1))) 
    Vvec = np.einsum('ijk,ilk -> ilj', Vtransform_matrix, Vvec) #magic to apply the transform matrix to each vector in Vvec
    #in the bussi paper, they scale the structures to make the stacking and H-bonding the same length
    #I'm not going to do that to save us a computation, also because oxDNA has very similar H-bonding
    #and stacking lengths.  Using energies would probably be a more accurate method to get differences.
    #Vvec = Vvec*np.array([])
    #in the paper they also convert to cyllindrical coordinates to make figure 1.  To reproduce that I also convert to cyllindrical coords, but its not necessary for the eRMSD calculation
    #Vnetwork_matrix[:, :, 0] = np.linalg.norm(np.array([Vvec[:, :, 0], Vvec[:, :, 1], 0])) #magnitude of a1-a2 projection
    #with catch_warnings(): #dividing the whole matrix produces a divide by 0 warning on the diagonal, which is annoying
    #    simplefilter("ignore")
    #    Vnetwork_matrix[:, :, 1] = np.arctan(Vvec[:, :, 1]/Vvec[:, :, 0]) #angle of a1-a2 projection
    #Vnetwork_matrix[:, :, 2] = Vvec[:, :, 2] #a3 projection
    
    #Leaving this in here for posterity, it does the exact same thing as the above block, but ~80x slower
    '''start = time.clock()
    for i,val1 in enumerate(interaction_sites): #reference base
        transform_matrix = np.array([val1[1], val1[2], val1[3]]).T
        for j,val2 in enumerate(interaction_sites): #other base
            vec = val2[0] - val1[0]
            #if i == 0 and j == 1: print transform_matrix
            if i != j:
                transformed_vec = np.linalg.solve(transform_matrix, vec)
                proj = [transformed_vec[0], transformed_vec[1], 0]
                p = np.linalg.norm(proj)
                theta = np.arctan((transformed_vec[1] / transformed_vec[0]))
                z = transformed_vec[2]
            else: 
                theta = 0.
                p = 0.
                z = 0.
            network_matrix[i][j] = np.array([p, theta, z]) #cylindrical coords
    print network_matrix[0, 1]
    print time.clock() - start'''
    return(Vvec)

#part of the eRMSD calculation
def calc_gvec(r):
    #if we actually used cyllindrical coords, this line would be necessary
    #r = np.array([connection[0] * np.cos(connection[1]), connection[0] * np.sin(connection[1]), connection[2]]) #this is cartesian wrt the a-vectors
    r_norm = np.linalg.norm(r)
    y = 2.9 #cutoff distance
    if r_norm < y and r_norm != 0:
        G = np.array([np.sin(y * r_norm) * (r[0]/r_norm), np.sin(y * r_norm) * (r[1]/r_norm), np.sin(y * r_norm) * (r[2]/r_norm), 1+np.cos(y * r_norm)])
    else:
        G = np.array([0., 0., 0., 0.])
    return G

#calculate eRMSD, for detailed information see: 
# Bottaro, S., Di Palma, F., & Bussi, G. (2014). 
# The role of nucleobase interactions in RNA structure and dynamics. 
def calc_eRMSD(mat1, mat2):
    SD = 0
    for ie, i in enumerate(mat1):
        for je,  j in enumerate(i):
            G1 = calc_gvec(j)
            G2 = calc_gvec(mat2[ie][je])
            SD += np.square(np.linalg.norm(G1-G2))
    SD /= len(i) #now its MSD
    SD = np.sqrt(SD) #and now its RMSD
    return SD

#reads in every combination of configurations and calls calc_eRMSD()
#only fills in half the matrix since eRMSDs are symmetric, though the connection matrix is not
def get_eRMSDs(r1, r2, inputfile, traj_file, top_file, num_confs, start=None, stop=None):
    if stop is None:
        stop = num_confs
    else: stop = int(stop)
    if start is None:
        start = 0
    else: start = int(start)
    confid = 0

    system1 = r1._get_system(N_skip=start)
    system2 = r2._get_system(N_skip=start+1)
    eRMSDs = np.zeros((num_confs, num_confs))
    i = start
    j = start+1
    while system1 != False and confid < stop: 
        print ("working on configuration", i, "t =", system1._time)
        system1.inbox()
        mat1 = calc_matrix(system1, inputfile)
        while system2:
            print ("working on configuration", i, "compared to", j)
            system2.inbox()
            mat2 = calc_matrix(system2, inputfile)
            eRMSDs[i][j] = calc_eRMSD(mat1, mat2)
            system2 = r2._get_system()
            j+=1

        i += 1
        j = i+1
        confid += 1
        system1 = r1._get_system()
        r2 = LorenzoReader2(traj_file, top_file)
        system2 = r2._get_system(N_skip=j)

    return(eRMSDs)

##############################################################################################################
#Welcome to the meat
#First, we read arguments and calculate the eRMSD between all structures
#This will not produce any outputs
def main():
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Calculate differences between structures and automatically apply DBSCAN to retrieve clusters")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy", "matplotlib"])

    traj_file = args.trajectory[0]
    inputfile = args.inputfile[0] 
    parallel = args.parallel
    if parallel:
        n_cpus = args.parallel[0]

    top_file = get_input_parameter(inputfile, "topology")
    if "RNA" in get_input_parameter(inputfile, "interaction_type"):
        environ["OXRNA"] = "1"
    else:
        environ["OXRNA"] = "0"
    num_confs = cal_confs(traj_file)
    import UTILS.base #this needs to be imported after the model type is set


    r2 = LorenzoReader2(traj_file, top_file)

    #how do you want to get your eRMSDs?  Do you need to do the time-consuming calculation or is it done and you have a pickle?
    if not parallel:
        r1 = LorenzoReader2(traj_file, top_file)

        eRMSDs = get_eRMSDs(r1, r2, inputfile, traj_file, top_file, num_confs)
    if parallel:
        out = parallelize_lorenzo_onefile.fire_multiprocess(traj_file, top_file, get_eRMSDs, num_confs, n_cpus, r2, inputfile, traj_file, top_file, matrix=True)
        eRMSDs = np.sum((i for i in out), axis=0)
    #eRMSDs = pickle.load(open('tmp_eRMSDs', 'rb'))

    #the eRMSD matrix is actually only half a matrix
    for ni,i in enumerate(eRMSDs):
        for nj,j in enumerate(i):
            eRMSDs[nj][ni] = j
            if ni == nj:
                eRMSDs[ni][nj] = 0

    #since calculating the eRMSDs are so time-consuming to calculate we're gonna pickle it to iterate the DBSCAN later.
    with open("tmp_eRMSDs", "wb") as file:
        pickle.dump(eRMSDs, file)

    ###############################################################################################################
    #Next, we're going to perform a DBSCAN on that matrix of eRMSDs to find clusters of similar structures
    perform_DBSCAN(eRMSDs, num_confs, traj_file, inputfile, "precomputed", 12, 8)

    ##############################################################################################################

if __name__ == '__main__':
    main()
    