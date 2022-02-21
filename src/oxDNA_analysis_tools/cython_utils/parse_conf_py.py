#cimport numpy as np
#from ..m_mean import BaseArray
#
#cdef tuple parse_conf(str lines, int nbases):
#    lines = lines.split('\n')

import cython
from collections import namedtuple
import numpy as np

BaseArray = namedtuple("BaseArray", ["time","box", "energy", "positions", "a1s", "a3s"])
def parse_conf(lines,nbases):
    lines = lines.split('\n')
    parameters =  np.array([line.split(maxsplit=9)[0:9] for line in lines[3:3+nbases]],dtype=float)
    conf = BaseArray(
        float(lines[0][lines[0].index("=")+1:]),
        np.array(lines[1].split("=")[1].split(), dtype=float),
        np.array(lines[2].split("=")[1].split(), dtype=float),
        parameters[:,0:3],
        parameters[:,3:6],
        parameters[:,6:9],
    )
    return conf