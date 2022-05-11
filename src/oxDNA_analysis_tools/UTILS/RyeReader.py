import numpy as np
from pickle import loads, dumps
from io import StringIO
from os.path import exists
import os
from oxDNA_analysis_tools.UTILS.data_structures import *

def Chunker(file, fsize, size=1000000):
    current_chunk = 0  
    while True:
        b = file.read(size)
        if not b: break
        yield Chunk(b,current_chunk*size, current_chunk * size + size > fsize, fsize)
        current_chunk+=1

#calculates the length of a trajectory file
def _index(traj_file): 
    def find_all(a_str, sub):
        #https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
        start = 0
        idxs = []
        while True:
            start = a_str.find(sub, start)
            if start == -1: return idxs
            idxs.append(start)
            start += len(sub) 
    # our trajectory files can be splitted by the occurence of t
    val = b"t" 
    conf_starts = []
    fsize = os.stat(traj_file).st_size 
    for chunk in Chunker(open(traj_file, 'rb'), fsize):
        idxs = np.array(find_all(chunk.block,val)) + chunk.offset # find all offsets of t in file
        conf_starts.extend(idxs)
    conf_starts = np.array(conf_starts)
    #generate a list of confs for all but the last one
    idxs = [ConfInfo(conf_starts[i], conf_starts[i+1] - conf_starts[i],i) 
                                            for i in range(len(conf_starts)-1)]
    #handle last offset
    idxs.append(ConfInfo(conf_starts[-1], fsize - conf_starts[-1], len(conf_starts)-1))
    return idxs


def get_traj_info(traj):
    #if idxs is None: # handle case when we have no indexes provided
    if not(exists(traj+".pyidx")):
        idxs = _index(traj) # no index created yet
        with open(traj+".pyidx","wb") as file:
            file.write(dumps(idxs)) # save it
    else:
        #we can load the index file
        with open(traj+".pyidx","rb") as file:
            idxs = loads(file.read())
    return TrajInfo(traj,len(idxs),idxs)

def inbox(conf, center=False):
    """
        Modify the positions attribute such that all positions are inside the box.
    """
    def realMod (n, m):
        return(((n % m) + m) % m)
    def coord_in_box(p):
        p = realMod(p, conf.box)
        return(p)
    def calc_PBC_COM(conf):
        angle = (conf.positions * 2 * np.pi) / conf.box
        cm = np.array([[np.sum(np.cos(angle[:,0])), np.sum(np.sin(angle[:,0]))], 
        [np.sum(np.cos(angle[:,1])), np.sum(np.sin(angle[:,1]))], 
        [np.sum(np.cos(angle[:,2])), np.sum(np.sin(angle[:,2]))]]) / len(angle)
        return conf.box / (2 * np.pi) * (np.arctan2(-cm[:,1], -cm[:,0]) + np.pi)
    target = np.array([conf.box[0] / 2, conf.box[1] / 2, conf.box[2] / 2])
    cms = calc_PBC_COM(conf)
    positions = conf.positions + target - cms   
    new_poses = coord_in_box(positions)
    positions += (new_poses - conf.positions)
    if center:
        cms = np.mean(positions, axis=0)
        positions -= cms
    return Configuration(
        conf.time, conf.box, conf.energy,
        positions, conf.a1s, conf.a3s
    )

def write_conf(path,conf, append=False):
    """
        write the conf to a file
    """
    out = []
    out.append('t = {}'.format(int(conf.time)))
    out.append('b = {}'.format(' '.join(conf.box.astype(str))))
    out.append('E = {}'.format(' '.join(conf.energy.astype(str))))
    for p, a1, a3 in zip(conf.positions, conf.a1s, conf.a3s):
        out.append('{} {} {} 0 0 0 0 0 0'.format(' '.join(p.astype(str)), ' '.join(a1.astype(str)), ' '.join(a3.astype(str))))
    
    mode = 'a' if append else 'w'
    with open(path,mode) as f:
        f.write("\n".join(out))

def conf_to_str(conf):
    """
    Write configuration as a string

    Parameters
    ----------
    conf (Configuration) : The configuration to write

    Returns
    -------
    (str) : The configuration as a string
    """
    # When writing a configuration to a file, the conversion from ndarray to string is the slowest part
    # This horrific list comp is the best solution we found
    header = f't = {int(conf.time)}\nb = {" ".join(conf.box.astype(str))}\nE = {" ".join(conf.energy.astype(str))}\n'
    return(''.join([header, ''.join([('{} {} {} 0 0 0 0 0 0\n'.format(' '.join(p.astype(str)), ' '.join(a1.astype(str)), ' '.join(a3.astype(str)))) for p, a1, a3 in zip(conf.positions, conf.a1s, conf.a3s)])]))

def get_top_info(top):
    """
        bare bones of topology info
    """
    with open(top) as f:
        my_top_info = f.readline().split(' ')
        if len(my_top_info)  == 2:
            nbases, nstrands = my_top_info
        elif len(my_top_info) == 5:
            nbases, nstrands, ndna, nres, ndnastrands = my_top_info
    return TopInfo(int(nbases), int(nstrands))

def describe(top, traj):
    """
        retrieve top and traj info for a provided pair
    """
    return get_top_info(top), get_traj_info(traj)

def no_top_describe(traj):
    """
        Retrieve top and traj info without providing a topology. 

        Note that the resulting top_info will have 0 strands because that information cannot be found in the trajectory. 
    """

    with open(traj) as f:
        l = ''
        # dump the header
        for _ in range(4):
            l = f.readline()
        n_bases = 0
        while (l[0] != 't'):
            n_bases += 1
            l = f.readline()
            if l == '':
                break

    return TopInfo(int(n_bases), 0), get_traj_info(traj)

def strand_describe(top):
    """
        Retrieve all information from topology file mapping nucleotides to strands.
        
        Parameters:
            top (str) : path to topology file

        Returns:
            system (System) : system object 
            monomers (list of Monomer) : list of monomers
    """
    get_neighbor = lambda x: monomers[x] if x != -1 else None

    with open (top) as f:
        l = f.readline().split()
        nmonomers = int(l[0])
        nstrands = int(l[1])

        system = System([None] * nstrands)
        monomers = [Monomer(i, None, None, None, None) for i in range(nmonomers)]

        ls = f.readlines()

        l = ls[0].split()
        curr = int(l[0])
        mid = 0
        s_start = 0
        s = Strand(curr)
        monomers[mid].type = l[1]
        monomers[mid].strand = s
        monomers[mid].n3 = get_neighbor(int(l[2]))
        monomers[mid].n5 = get_neighbor(int(l[3]))
        l = ls[1].split()
        mid += 1
        while l:
            if int(l[0]) != curr:
                s.monomers = monomers[s_start:mid]
                system[curr-1] = s #this is going to do something weird with proteins
                curr = int(l[0])
                s = Strand(curr)
                s_start = mid
            
            monomers[mid].type = l[1]
            monomers[mid].strand = s
            monomers[mid].n3 = get_neighbor(int(l[2]))
            monomers[mid].n5 = get_neighbor(int(l[3]))

            mid += 1
            try:
                l = ls[mid].split()
            except IndexError:
                break  

    s.monomers = monomers[s_start:mid]
    system[curr-1] = s #this is going to do something weird with proteins

    return system, monomers

