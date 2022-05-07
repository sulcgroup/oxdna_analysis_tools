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