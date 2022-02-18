from collections import namedtuple
from UTILS.RyeReader import Configuration, describe, get_confs, parse_conf, inbox, get_n_cpu, write_conf
from multiprocessing import Pool
import numpy as np
import time
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
def compute(ctx,ptr):
    confs = get_confs(ctx.traj_info, ptr*ctx.ntopart,ctx.ntopart)
    confs = (inbox(parse_conf(ctx.top_info, c)) for c in confs)
    # convert to numpy repr
    aligned_coords = np.asarray([[align_conf.positions, align_conf.a1s, align_conf.a3s] 
                                                                  for align_conf in confs])
    aligned_coords = [align(ctx.centered_ref_coords,ctx.cms_ref_coords, c) for c in aligned_coords]
    sub_mean = np.sum(aligned_coords, axis=0)
    return sub_mean


top, traj = r"/mnt/g/hinge1/hinge_correct_seq.top",r"/mnt/g/hinge1/aligned.dat"
top_info, traj_info = describe(top, traj)


# fetch refference conf
ref_conf = parse_conf(top_info, get_confs(traj_info,0,1)[0])
ref_conf = inbox(ref_conf)

# figure out how much resorces we have
ncpus = get_n_cpu()
# how many confs we want to distribute between the processes
ntopart = 50
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
print("starting up processes")
results = [pool.apply_async(compute,(ctx,i)) 
                        for i in range(n_chunks)]
print("all spawned")


# get the results
acc = np.zeros([3, top_info.nbases, 3])
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


write_conf("/mnt/g/hinge1/mean_m.dat",Configuration(0,ref_conf.box,np.array([0,0,0]), pos, a1s , a3s))
print("--- %s seconds ---" % (time.time() - start_time))