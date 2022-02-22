import numpy as np
from os import getenv, path
from multiprocessing import cpu_count, Pool
from collections import namedtuple
from io import StringIO
from pickle import loads, dumps
from os.path import exists
import os

import time
start_time = time.time()



#top, traj = "./hinge_correct_seq.top","./100confs.dat"
top, traj = r"/mnt/g/hinge1/hinge_correct_seq.top",r"/mnt/g/hinge1/aligned.dat"


#figure out how many chunks we are working on 

def get_n_cpu():
    try:
        available_cpus = int(getenv('SLURM_NTASKS'))
    except:
        available_cpus = cpu_count()
    return available_cpus


def align(centered_ref_coords, cms_ref_cords ,coords):
    # center on centroid
    av1, reference_coords = cms_ref_cords, centered_ref_coords
    av2 = np.mean(coords[0], axis=0)
    coords[0]= coords[0] - av2
    # correlation matrix
    a = np.dot(np.transpose(coords[0]), reference_coords)
    u, _, vt = np.linalg.svd(a)
    rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    # check if we have found a reflection
    if np.linalg.det(rot) < 0:
        vt[2] = -vt[2]
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    tran = av1 - np.dot(av2, rot)
    return  (np.dot(coords[0], rot) + tran,
             np.dot(coords[1], rot),
             np.dot(coords[2], rot))


BaseArray = namedtuple("BaseArray", ["time","box", "energy", "positions", "a1s", "a3s"])
#def parse_conf(lines,nbases):
#    lines = lines.split('\n')
#    parameters =  np.array([line.split(maxsplit=9)[0:9] for line in lines[3:3+nbases]],dtype=float)
#    conf = BaseArray(
#        float(lines[0][lines[0].index("=")+1:]),
#        np.array(lines[1].split("=")[1].split(), dtype=float),
#        np.array(lines[2].split("=")[1].split(), dtype=float),
#        parameters[:,0:3],
#        parameters[:,3:6],
#        parameters[:,6:9],
#    )
#    return conf
from parse_conf import parse_conf

def get_confs(idxs,traj_path, start, nconfs):
    conf_count = len(idxs)
    if(start+nconfs >= conf_count): # make sure we stay in bounds 
        nconfs = conf_count - start
    # go to start position 
    with open(traj_path, 'rb') as traj_file:
        traj_file.seek(idxs[start].offset)
        # figure out how big of a chunk we want to read 
        size = sum([idxs[i].size for i in range(start,start+nconfs)])
        # read chunk and prepare to split
        chunk = traj_file.read(size)
        return chunk.split(b't')


def compute(traj,nbases, idxs, centered_ref_coords, cms_ref_cords, ntopart, ptr):
    #reader = MichaReader(top,traj)
    confs = get_confs(idxs,traj, ptr*ntopart,ntopart)
    #TODO: NEEDS INBOX
    confs = [inbox(parse_conf(c,nbases)) for c in confs[1:]]
    # convert to numpy repr
    aligned_coords = np.asarray([[align_conf.positions, align_conf.a1s, align_conf.a3s] 
                                                                  for align_conf in confs])
    sub_mean = np.zeros(shape=[3,len(confs[0].positions),3])
    for c in aligned_coords:
        sub_mean += align(centered_ref_coords,cms_ref_cords, c)
    #aligned_coords = [align(centered_ref_coords,cms_ref_cords, c) for c in aligned_coords]
    #sub_mean = np.sum(aligned_coords, axis=0)

    return sub_mean


def write_conf(path,conf):
    out = []
    out.append('t = {}'.format(int(conf.time)))
    out.append('b = {}'.format(' '.join(conf.box.astype(str))))
    out.append('E = {}'.format(' '.join(conf.energy.astype(str))))
    for p, a1, a3 in zip(conf.positions, conf.a1s, conf.a3s):
        out.append('{} {} {} 0.0 0.0 0.0 0.0 0.0 0.0'.format(' '.join(p.astype(str)), ' '.join(a1.astype(str)), ' '.join(a3.astype(str))))
    with open(path,"w") as f:
        f.write("\n".join(out))





# chunk file into blocks of given size
Chunk = namedtuple('Chunk', ['block','offset', 'is_last','file_size'])
def blocks(file, fsize, size=1000000):
    current_chunk = 0  
    while True:
        b = file.read(size)
        if not b: break
        yield Chunk(b,current_chunk*size, current_chunk * size + size > fsize, fsize)
        current_chunk+=1

def find_all(a_str, sub):
    #https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
    start = 0
    idxs = []
    while True:
        start = a_str.find(sub, start)
        if start == -1: return idxs
        idxs.append(start)
        start += len(sub) # use start += 1 to find overlapping matches

#calculates the length of a trajectory file
ConfInfo = namedtuple('ConfInfo', ['offset', 'size','id'])
def index(traj_file): 
    val = b"t" 
    conf_starts = []
    fsize = os.stat(traj_file).st_size 
    for chunk in blocks(open(traj_file, 'rb'), fsize, size=1000000):#5500000):
        idxs = np.array(find_all(chunk.block,val)) + chunk.offset # find all offsets of t in file
        conf_starts.extend(idxs)

    conf_starts = np.array(conf_starts)
    #generate a list of confs for all but the last one
    idxs = [ConfInfo(conf_starts[i], conf_starts[i+1] - conf_starts[i],i) 
                                            for i in range(len(conf_starts)-1)]
    #handle last offset
    idxs.append(ConfInfo(conf_starts[-1], fsize - conf_starts[-1], len(conf_starts)-1))
    return idxs


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
    return BaseArray(
        conf.time, conf.box, conf.energy,
        positions, conf.a1s, conf.a3s
    )


#if idxs is None: # handle case when we have no indexes provided
if not(exists(traj+".pyidx")):
    idxs = index(traj) # no index created yet
    with open(traj+".pyidx","wb") as file:
        file.write(dumps(idxs)) # save it
else:
    #we can load the index file
    with open(traj+".pyidx","rb") as file:
        idxs = loads(file.read())

#TODO working only on 10
#idxs = idxs[:10] #first 10

with open(top) as f:
    my_top_info = f.readline().split(' ')
    if len(my_top_info)  == 2:
        nbases, nstrands = my_top_info
    elif len(my_top_info) == 5:
        nbases, nstrands, ndna, nres, ndnastrands = my_top_info
    nbases = int(nbases)

ref_conf = parse_conf(
    get_confs(idxs, traj, 0, 1)[1],
    nbases
)
#TODO: NEEDS INBOX
ref_conf = inbox(ref_conf)

ncpus = get_n_cpu()
ntopart = 50
#prepare to fire multiple processes 
pool = Pool(ncpus)

n_confs  = len(idxs) 
n_chunks = int(n_confs / ntopart +
                     (1 if n_confs % ntopart else 0))

reference_coords = ref_conf.positions
av1 = np.mean(reference_coords, axis=0)
reference_coords = reference_coords - av1



print("starting up processes")
results = [pool.apply_async(compute,(traj,len(ref_conf.positions), idxs,reference_coords, av1, ntopart,i)) 
                        for i in range(n_chunks)]
print("all spawned")


# get the results
acc = np.zeros([3, nbases, 3])
for i,r in enumerate(results):
    print(f"finished {i+1}/{n_chunks}",end="\r")
    acc += r.get()
print()
# compute the mean 
acc /= n_confs
pos, a1s, a3s = acc

# renormalize
a1s = np.array([v/np.linalg.norm(v) for v in a1s])
a3s = np.array([v/np.linalg.norm(v) for v in a3s])


write_conf(r"/mnt/g/hinge1/mean_m.dat",BaseArray(0,ref_conf.box,np.array([0,0,0]), pos, a1s , a3s))
print("--- %s seconds ---" % (time.time() - start_time))