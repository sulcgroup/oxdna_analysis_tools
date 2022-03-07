from collections import namedtuple
from UTILS.RyeReader import Configuration, describe, inbox, write_conf
from multiprocessing import Pool
import numpy as np
import time
from oxDNA_analysis_tools.UTILS.get_confs import get_confs
start_time = time.time()

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



ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "centered_ref_coords",
                                              "cms_ref_coords",
                                              "ntopart"])
def compute(ctx:ComputeContext,chunk_id):
    confs = get_confs(ctx.traj_info.idxs, ctx.traj_info.path, chunk_id*ctx.ntopart, ctx.ntopart, ctx.top_info.nbases)
    confs = (inbox(c) for c in confs)
    # convert to numpy repr
    aligned_coords = np.asarray([[align_conf.positions, align_conf.a1s, align_conf.a3s] 
                                                                  for align_conf in confs])
    
    sub_mean = np.zeros(shape=[3,ctx.top_info.nbases,3])
    for c in aligned_coords:
        sub_mean += align(ctx.centered_ref_coords,ctx.cms_ref_coords, c)
    

    return sub_mean

#top, traj = r"/mnt/c/Users/mmatthi3/test2/hinge_correct_seq.top", r"/mnt/c/Users/mmatthi3/test2/aligned.dat"
top, traj = "./hinge_correct_seq.top","./100confs.dat"
top_info, traj_info = describe(top, traj)

# fetch reference conf
ref_conf = get_confs(traj_info.idxs, traj, 0, 1, top_info.nbases)[0]
ref_conf = inbox(ref_conf)

# figure out how much resorces we have
ncpus = 5#get_n_cpu()
# how many confs we want to distribute between the processes
ntopart = 20
#prepare to fire multiple processes 
pool = Pool(ncpus)

# deduce how many chunks we have to run in parallel
n_confs  = traj_info.nconfs 
n_chunks = int(n_confs / ntopart +
                     (1 if n_confs % ntopart else 0))

# alignment requires the ref to be centered at 0
reference_coords = ref_conf.positions
ref_cms = np.mean(reference_coords, axis=0) # cms prior to centering
reference_coords = reference_coords - ref_cms


ctx = ComputeContext(
    traj_info, top_info, reference_coords, ref_cms, ntopart
)
print(f"Starting up {ncpus} processes for {n_chunks} chunks")
results = [pool.apply_async(compute,(ctx,i)) for i in range(n_chunks)]
print("All spawned")


# get the results
acc = np.zeros([3, top_info.nbases, 3])
for i,r in enumerate(results):
    print(f"finished {i+1}/{n_chunks}",end="\r")
    acc += r.get()
pool.close()
pool.join()

print()
# compute the mean 
acc /= n_confs
pos, a1s, a3s = acc

# renormalize
a1s = np.array([v/np.linalg.norm(v) for v in a1s])
a3s = np.array([v/np.linalg.norm(v) for v in a3s])


write_conf("./mean_r.dat",Configuration(0,ref_conf.box,np.array([0,0,0]), pos, a1s , a3s))
print("--- %s seconds ---" % (time.time() - start_time))