from dataclasses import dataclass
from __future__ import annotations
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

@dataclass(slots=True)
class System:
    top_info : TopInfo
    strands : list

    def __getitem__(self, key):
        return self.strands[key]

    def __setitem__(self, key, value):
        self.strands[key] = value

    def __iter__(self):
        return (s for s in self.strands)

@dataclass(slots=True)
class Strand:
    nbases : int
    bases : list

    def __getitem__(self, key):
        return self.bases[key]

    def __setitem__(self, key, value):
        self.bases[key] = value

    def __iter__(self):
        return (b for b in self.bases)

@dataclass(slots=True)
class Base:
    id : int
    type : str
    strand : Strand
    n3 : Base
    n5 : Base