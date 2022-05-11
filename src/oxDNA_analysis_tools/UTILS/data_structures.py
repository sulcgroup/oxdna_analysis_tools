from __future__ import annotations
from dataclasses import dataclass
import numpy as np

# Data class to hold information about a chunk of a trajectory
@dataclass(slots=True)
class Chunk:
    block : str
    offset : int
    is_last : bool
    file_size : int

# Data class to hold metadata about an individual configuration
@dataclass(slots=True)
class ConfInfo:
    offset : int
    size : int
    id : int

# Data class to hold information about a whole trajectory file
@dataclass(slots=True)
class TrajInfo:
    path : str
    nconfs : int
    idxs : ConfInfo

# Data class which actually contains configuration information
@dataclass(slots=True)
class Configuration:
    time : int
    box : np.array
    energy : np.array
    positions : np.array
    a1s : np.array
    a3s : np.array

# Data class which contains topology information
@dataclass(slots=True)
class TopInfo:
    nbases : int
    nstrands : int

class System:
    __slots__ = ('strands')
    strands : list

    def __init__(self, strands = []):
        self.strands = strands

    def __getitem__(self, key):
        return self.strands[key]

    def __setitem__(self, key, value):
        self.strands[key] = value

    def __iter__(self):
        return (s for s in self.strands)

    def append(self, strand):
        self.strands.append(strand)

class Strand:
    __slots__ = ('id', 'monomers')
    id : int
    monomers : list

    def __init__(self, id):
        self.id = id

    def __getitem__(self, key):
        return self.monomers[key]

    def __setitem__(self, key, value):
        self.monomers[key] = value

    def __iter__(self):
        return (m for m in self.monomers)

@dataclass(slots=True)
class Monomer:
    id : int
    type : str
    strand : Strand
    n3 : Monomer
    n5 : Monomer