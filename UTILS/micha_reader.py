from collections import namedtuple
import os
import numpy as np
from base_array import base_array
from copy import deepcopy

# chunk file into blocks of given size
Chunk = namedtuple('Chunk', ['block','offset', 'is_last','file_size'])
def blocks(file, size=1000000):
    current_chunk = 0 
    fsize = os.stat(traj_file).st_size
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
    counter = 0
    for chunk in blocks(open(traj_file, 'rb'), size=5000000):
        idxs = np.array(find_all(chunk.block,val)) + chunk.offset # find all offsets of t in file
        conf_starts.extend(idxs)

    conf_starts = np.array(conf_starts)
    #generate a list of confs for all but the last one
    idxs = [ConfInfo(conf_starts[i], conf_starts[i+1] - conf_starts[i],i) 
                                            for i in range(len(conf_starts)-1)]
    #handle last offset
    idxs.append(ConfInfo(conf_starts[-1], size - conf_starts[-1], len(conf_starts)-1))
    return idxs

TopInfo = namedtuple('TopInfo', ['bases', 'strands'])
class MichaReader:
    def __init__(self, top, traj_file, idxs = None):
        # setup configuration index
        if idxs is None: # handle case when we have no indexes provided
            self.idxs = index(traj_file)
        else:
            self.idxs = idxs
        # store the number of configurations
        self.conf_count = len(self.idxs)

        #get the topology file info 
        with open(top) as f:
            nbases, nstrands = f.readline().split(' ')
            self.top_info = TopInfo(int(nbases), int(nstrands))
        
        #open trajectory file
        self.traj_file = open(traj_file, 'r')

        #setup dummy conf 
        self._conf = base_array(
            0, np.zeros(3), 0,
            np.zeros([self.top_info.bases, 3], dtype=float),
            np.zeros([self.top_info.bases, 3], dtype=float),
            np.zeros([self.top_info.bases, 3], dtype=float),
        )
        self.state = -1
    
    def _get_conf(self, idx):
        if idx >= self.conf_count:
            raise Exception("Invalid configuration index")

        self.traj_file.seek(self.idxs[idx].offset)
        conf = self.traj_file.read(self.idxs[idx].size)
        return conf

    def _parse_conf(self,idx):
        lines = self._get_conf(idx).split('\n') # adds an extra empty one at the end
        if len(lines) -4 != self.top_info.bases:
            raise Exception("Invalid configuration.")
        
        # populate our dummy conf by data from the conf
        self._conf.time = float(lines[0][lines[0].index("=")+1:])
        self._conf.box = np.array(lines[1].split("=")[1].split(), dtype=float)
        self._conf.energy = np.array(lines[2].split("=")[1].split(), dtype=float)

        # parse out the pos and a's 
        for i in range(self.top_info.bases):
            line = lines[3+i].split()
            self._conf.positions[i] = np.array(line[0:3], dtype=float)
            self._conf.a1s[i] = np.array(line[3:6], dtype=float)
            self._conf.a3s[i] = np.array(line[6:9], dtype=float)
        return deepcopy(self._conf)
    
    def read(self,idx=None):
        if(idx):
            self.state=idx
        else:
            self.state+=1
        if self.state >= self.conf_count:
            return None
        return self._parse_conf(self.state)