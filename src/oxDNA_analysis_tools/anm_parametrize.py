from sys import path
import numpy as np
import os
import json

from oxDNA_analysis_tools.UTILS.readers import ErikReader
import argparse

# calculate the position of the super particle
def mean_pos(conf,idxs):
    return np.sum(conf.positions[idxs], axis=0) / idxs.shape[0]
# calculate the particle positions
def get_mean_positions(conf, particles_array):
    return np.array(
        list(map(lambda idx: mean_pos(conf,idx),particles_array))
    )
# calculate deviation from mean position
def calculate_deviations(positions, reference_configuration):
    d = np.subtract(positions, reference_configuration) # Positions subtracted (dx, dy, dz)
    devs = np.sqrt(np.sum(np.square(d), axis=1)) # sqrt(dx**2 + dy**2 + dz**2)
    #print(devs.shape)
    return devs

def main():
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="compute par file for DNA-ANM model")
    parser.add_argument('index_file', type=str, help="Index file describing what Bases belong to which Super particle.")
    parser.add_argument('mean_file', type=str,  help="Reference configuratio.")
    parser.add_argument('trajectory', type=str,  help="Trajectory to evaluate.")
    parser.add_argument('out_file', type=str, help="Output par file name.")
    args = parser.parse_args()

    #if not "out_file" in args.format:
    #    exit(1)
    particles_array = []
    # import index file
    with open(args.index_file) as file:
        index_lines = file.readlines()
    # parse the file
    for line in index_lines:
        particles_array.append(
                np.array(line.split(" "),dtype=np.int)
        )
    # parsed index file
    particles_array = np.array(particles_array,dtype=object)

    # import mean configuration
    r = ErikReader(args.mean_file)
    ref = r.read()
    ref.inbox()

    # calculate reference conf particle positions
    reference_configuration = get_mean_positions(ref, particles_array)

    # to collect the distance data of the superparticles 
    trajectory_devs = []
    # go over the trajectory to compute the distances array
    reader = ErikReader(args.trajectory)
    system = reader.read()
    i=0
    while system:
        if(i%10 == 0):
            print(i)
        system.inbox()
        # the position of the superparticles
        positions_of_current_conf = get_mean_positions(system)
        # Calculate Deviations of current conf from mean
        current_devs = calculate_deviations(positions_of_current_conf, reference_configuration)
        # Add to Storage Array
        trajectory_devs.append(current_devs)
        system = reader.read()
        i+=1

    # make it into numpy as that is easier to compute with
    trajectory_devs = np.array(trajectory_devs)

    devs_sqrd = np.square(trajectory_devs) # Square all devs

    mean_devs_sqrd = np.mean(devs_sqrd, axis=0) # Mean of Squared Devs
    print(mean_devs_sqrd.shape)

    #deviations are in sim units, convert to nm in next step

    # Convert from sim units to nm
    rmsf = np.multiply(np.sqrt(mean_devs_sqrd),0.8518)

    #Format for json output
    dict = {"RMSF (nm)": rmsf.tolist()}

    with open(args.out_file, 'w') as f:
        json.dump(dict, f)

if __name__ == '__main__':
    main()
    