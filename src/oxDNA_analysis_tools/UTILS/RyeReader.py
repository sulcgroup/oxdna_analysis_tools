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

def get_confs(traj_info, start, nconfs):
    """
        low level function getting the unparsed configurations
    """
    #conf_count = len(traj_info.idxs)
    if(start+nconfs >= traj_info.nconfs): # make sure we stay in bounds 
        nconfs = traj_info.nconfs - start
    # go to start position 
    with open(traj_info.path) as traj_file:
        traj_file.seek(traj_info.idxs[start].offset)
        # figure out how big of a chunk we want to read 
        size = sum([traj_info.idxs[i].size for i in range(start,start+nconfs)])
        # read chunk and prepare to split
        chunk = StringIO(traj_file.read(size)) # work with the string like a file 
        return [chunk.read(traj_info.idxs[i].size)
                            for i in range(start,start+nconfs)] 

def inbox(conf):
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
    center = calc_PBC_COM(conf)
    positions = conf.positions + target - center   
    new_poses = coord_in_box(positions)
    positions += (new_poses - conf.positions)
    return Configuration(
        conf.time, conf.box, conf.energy,
        positions, conf.a1s, conf.a3s
    )

def write_conf(path,conf):
    """
        write the conf to a file
    """
    out = []
    out.append('t = {}'.format(int(conf.time)))
    out.append('b = {}'.format(' '.join(conf.box.astype(str))))
    out.append('E = {}'.format(' '.join(conf.energy.astype(str))))
    for p, a1, a3 in zip(conf.positions, conf.a1s, conf.a3s):
        out.append('{} {} {} 0.0 0.0 0.0 0.0 0.0 0.0'.format(' '.join(p.astype(str)), ' '.join(a1.astype(str)), ' '.join(a3.astype(str))))
    with open(path,"w") as f:
        f.write("\n".join(out))
 
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

