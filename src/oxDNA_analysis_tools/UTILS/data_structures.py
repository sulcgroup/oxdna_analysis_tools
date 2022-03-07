from dataclasses import dataclass
import numpy as np

# Data class to hold information about a chunk of a trajectory
@dataclass
class Chunk:
    block : str
    offset : int
    is_last : bool
    file_size : int

# Data class to hold metadata about an individual configuration
@dataclass
class ConfInfo:
    offset : int
    size : int
    id : int

# Data class to hold information about a whole trajectory file
@dataclass
class TrajInfo:
    path : str
    nconfs : int
    idxs : ConfInfo

# Data class which actually contains configuration information
@dataclass
class Configuration:
    time : int
    box : np.array
    energy : np.array
    positions : np.array
    a1s : np.array
    a3s : np.array

# Data class which contains topology information
@dataclass
class TopInfo:
    nbases : int
    nstrands : int