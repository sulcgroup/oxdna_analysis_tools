from collections import namedtuple
import os
import numpy as np
from oxDNA_analysis_tools.UTILS.base_array import base_array
from copy import deepcopy
from io import StringIO
from os.path import exists
from pickle import loads, dumps

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
    counter = 0
    fsize = os.stat(traj_file).st_size
    for chunk in blocks(open(traj_file, 'rb'), fsize, size=1000000):
        idxs = np.array(find_all(chunk.block,val)) + chunk.offset # find all offsets of t in file
        conf_starts.extend(idxs)

    conf_starts = np.array(conf_starts)
    #generate a list of confs for all but the last one
    idxs = [ConfInfo(conf_starts[i], conf_starts[i+1] - conf_starts[i],i) 
                                            for i in range(len(conf_starts)-1)]
    #handle last offset
    idxs.append(ConfInfo(conf_starts[-1], fsize - conf_starts[-1], len(conf_starts)-1))
    return idxs


def parse_conf(lines,nbases):
    lines = lines.split('\n')
    if lines[-1] == '' and len(lines) -4 != nbases:
        raise Exception("Incorrect number of bases in topology file")
    elif lines[-1] != '' and len(lines) -3 != nbases:
        raise Exception("Incorrect number of bases in topology file")   
    #setup dummy conf 
    conf = base_array(
        0, np.zeros(3), 0,
        np.zeros([nbases, 3], dtype=float),
        np.zeros([nbases, 3], dtype=float),
        np.zeros([nbases, 3], dtype=float),
    )
    # populate our dummy conf by data from the conf
    conf.time = float(lines[0][lines[0].index("=")+1:])
    conf.box = np.array(lines[1].split("=")[1].split(), dtype=float)
    conf.energy = np.array(lines[2].split("=")[1].split(), dtype=float)
    # parse out the pos and a's 
    for i in range(nbases):
        line = lines[3+i].split()
        conf.positions[i] = np.array(line[0:3], dtype=float)
        conf.a1s[i] = np.array(line[3:6], dtype=float)
        conf.a3s[i] = np.array(line[6:9], dtype=float)
    return conf

def partition(confs, n):
    return [confs[i:i+n] for i in range(0, len(confs), n)]

def flatten(t): # more versions @ https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-a-list-of-lists
    return [item for sublist in t for item in sublist]

def handle_confs(ti, process, confs):
    return [process(parse_conf(x,ti.bases)) for x in confs]

TopInfo = namedtuple('TopInfo', ['bases', 'strands'])
class MichaReader:
    def __init__(self, top, traj_file, idxs = None, buff_size = 50):
        # setup configuration index
        if idxs is None: # handle case when we have no indexes provided
            if not(exists(traj_file+".pyidx")):
                self.idxs = index(traj_file) # no index created yet
                with open(traj_file+".pyidx","wb") as file:
                    file.write(dumps(self.idxs)) # save it
            else:
                #we can load the index file
                with open(traj_file+".pyidx","rb") as file:
                    self.idxs = loads(file.read())
        else: # index provided
            self.idxs = idxs
        # store the number of configurations
        self.conf_count = len(self.idxs)

        #get the topology file info 
        with open(top) as f:
            my_top_info = f.readline().split(' ')
            if len(my_top_info)  == 2:
                nbases, nstrands = my_top_info
            elif len(my_top_info) == 5:
                nbases, nstrands, ndna, nres, ndnastrands = my_top_info
            self.top_info = TopInfo(int(nbases), int(nstrands))

        
        #open trajectory file
        self.traj_file = open(traj_file, 'r')

        #information to track the file position and the current confs we handle
        self.ptr = 0
        self.buff_size = buff_size
        self.buff = []
    
    def _get_conf(self, idx):
        if idx > self.conf_count:
            raise Exception("Invalid configuration index")
        self.traj_file.seek(self.idxs[idx].offset)
        conf = self.traj_file.read(self.idxs[idx].size)
        return conf
    
    def _get_confs(self, start, nconfs):
        if(start+nconfs >= self.conf_count): # make sure we stay in bounds 
            nconfs = self.conf_count - start
        # go to start position 
        self.traj_file.seek(self.idxs[start].offset)
        # figure out how big of a chunk we want to read 
        size = sum([self.idxs[i].size for i in range(start,start+nconfs)])
        # read chunk and prepare to split
        chunk = StringIO(self.traj_file.read(size)) # work with the string like a file 
        return [chunk.read(self.idxs[i].size)
                            for i in range(start,start+nconfs)] 

    def _parse_conf(self,lines):
        return parse_conf(lines, self.top_info.bases)
    
    def read(self,idx=None):
        if(idx):
            if idx >= self.conf_count: return None
            lines = self._get_conf(idx)
            return self._parse_conf(lines)
        if(not self.buff and self.ptr < self.conf_count): # no confs in the buff - try get some 
            self.buff.extend( 
                self._get_confs(self.ptr, self.buff_size)[::-1] # save reversed list for speed with pop
            )
            self.ptr += self.buff_size
        if(self.buff): #handle the conf 
            return self._parse_conf(
                self.buff.pop()  # pops the last element out (as it's faster than using pop(0) which is O(n))
            )
        return None
      

    def __del__(self):
        self.traj_file.close()
