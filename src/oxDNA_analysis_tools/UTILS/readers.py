import os.path
import sys
import numpy as np
from oxDNA_analysis_tools.UTILS import base
from oxDNA_analysis_tools.UTILS.base_array import base_array

#helper for cal_confs
def blocks(file, size=1000000):
    """
    An iterator that returns chunks of the given file

    Parameters:
        file (file): An open file handler to chunk.
        <optional> size (int): The size of the chunks.  Defaults to 1 MB

    Yields:
        Chunks of the provided file.
    """
    while True:
        b = file.read(size)
        if not b: break
        yield b

#calculates the length of a trajectory file
def cal_confs(traj_file):
    """
    Calculates the number of configurations in a trajectory
    """
    with open(traj_file, "rb") as f:
        return (sum(bl.count(b"t") for bl in blocks(f)))

#gets the value out of an oxDNA input file
def get_input_parameter(input_file, parameter):
    """
    Gets the value of a parameter in an oxDNA input file
    Parameters:
        input_file (str): The path to the input file
        parameter (str): The parameter you want to get the value of

    Returns:
        value (str): The value of the parameter
    """
    fin = open(input_file)
    value = ''
    for line in fin:
        line = line.lstrip()
        if not line.startswith('#'):
            if parameter in line:
                value = line.split('=')[1].replace(' ','').replace('\n','')
    fin.close()
    if value == '':
        print("ERROR: Key {} not found in input file {}".format(parameter, input_file))
    return value

#The Python3 successor to the trusty LorenzoReader
#Reads in oxDNA trajectory files one configuration at a time to avoid memory overflow
#_get_system() returns a System object as defined in base.py
class LorenzoReader2:
    """
    Reads oxDNA trajectory files and creates System objects from individual configurations

    The System object is a descriptive representation of DNA which is made of strands which are made of nucleotides.  It has many built-in methods defined in base.py.
    """
    def __enter__(self):
        return self

    def __exit__(self,  exc_type, exc_val, exc_tb):
        self.__del__()

    def __init__(self, configuration, topology, check_overlap=False):
        # NOTE: on paths with symlinks kills the output stream 
        # so script continues to run ... probably  better to crash 
        if not os.path.isfile(configuration):
            print("Configuration file '{}' is not readable".format(configuration))
        if not os.path.isfile(topology):
            print("Topology file '{}' is not readable".format(topology))
        self._check_overlap = check_overlap
        self._conf = open(configuration, "r")
        with open(topology, "r") as file:
            file.readline()
            self._top_lines = file.readlines()
            self.particle_count = len(self._top_lines)

    def __del__(self):
        try: 
            if self._conf: 
                self._conf.close()
        except:
            print("ERROR: The reader could not load the provided configuration/topology pair")

    def __iter__(self):
        return self

    def __next__(self):
        s = self._get_system()
        if s:
            s.inbox_system()
            return s
        else:
            raise StopIteration
  
    def _read(self, only_strand_ends=False, skip=False):
        """
        Read a single configuration from the trajectory file

        Parameters:
            <optional> only_strand_ends (bool): If True, only the first and last nucleotide of a strand will be read.  Defaults to False.
            <optional> skip (bool): If True, the configuration will be skipped. This is how you track to a specific configuration.  Defaults to False.
        
        Returns:
            system (System): The System object representing the configuration
        """
        timeline = self._conf.readline()
        time = 0
        if  len(timeline) == 0:
            return False
        else:
            time = float(timeline.split()[2])
        
        box = np.array([float(x) for x in self._conf.readline().split()[2:]])
        [E_tot, E_pot, E_kin] = [float(x) for x in self._conf.readline().split()[2:5]]

        if skip: #this could be slightly optimized by moving it up, but then it wouldn't be same length as toplines
            for tl in self._top_lines:
                self._conf.readline()
            return False

        system = base.System(box, time=time, E_pot=E_pot, E_kin=E_kin)
        base.Nucleotide.index = 0
        base.Strand.index = 0

        s = False
        strandid_current = 0
        for tl in self._top_lines:
            tls = tl.split()
            n3 = int(tls[2])
            n5 = int(tls[3])
            strandid = int(tls[0])
            if (len (tls[1]) == 1):
                try:
                    if strandid > 0:
                        b = base.base_to_number[tls[1]]
                    if strandid < 0 :
                        b = base.aa_to_number[tls[1]]
                except:
                    b = 0
                bb = b
            else:
                try:
                    tmp = int (tls[1])
                except:
                    raise ValueError ("problems in topology file with specific base pairing")

                if tmp > 0:
                    b = tmp % 4
                else:
                    b = (3 - ((3 - tmp) % 4))
                bb = tmp

            if strandid != strandid_current:
                # check for circular strand
                if n3 != -1:
                    iscircular = True
                else:
                    iscircular = False

                if s:
                    system.add_strand(s, self._check_overlap)
                if strandid >= 0:
                    s = base.Strand()
                else:
                    s = base.Peptide()
                if iscircular:
                    s.make_circular()
                strandid_current = strandid

            ls = self._conf.readline().split()
            cm = [float(x) for x in ls[0:3]]
            a1 = [float(x) for x in ls[3:6]]
            a3 = [float(x) for x in ls[6:9]]
            v = [float(x) for x in ls[9:12]]
            L = [float(x) for x in ls[12:15]]
            if not only_strand_ends or n3 == -1 or n5 == -1:
                try:
                    s.add_nucleotide(base.Nucleotide(cm, a1, a3, b, bb, v, L, n3))
                except Exception as e:
                    print("ERROR: Reader died while reading configuration with t = {}.\n\nError message:\n{}".format(time, e))
                    return False

        system.add_strand(s, self._check_overlap)

        return system


    # if only_strand_ends == True then a strand will contain only the first and the last nucleotide
    # useful for some analysis like csd for coaxial interactions
    def _get_system(self, only_strand_ends=False, N_skip=0):
        for _ in range(N_skip):
            self._read(skip=True)

        return self._read(only_strand_ends=only_strand_ends, skip=False)


class ErikReader:
    """
    A bare-bones trajectory reader which only reads positions and rotations, ignoring connectivities. Produces a base_array object.
    """
    def __enter__(self):
        return self

    def __exit__(self,  exc_type, exc_val, exc_tb):
        self.__del__()

    def __del__(self):
        try: 
            if self._conf: 
                self._conf.close()
        except:
            print("ERROR: The reader could not load the provided configuration")
            sys.exit(1)

    def __init__(self, configuration):
        if not os.path.isfile(configuration):
            print("Configuration file '{}' is not readable".format(configuration))
            sys.exit(1)
        self._conf = open(configuration, "r")
        self._time = None
        self._box = np.zeros(3)
        self._conf_energy = np.zeros(3)
        self._len = 0


    def _read_first(self):
        """
        Special reader that handles the first line and allocates the memory for storing the system
        """
        self._line = self._conf.readline().split()
        if len(self._line) == 0:
            print("ERROR: the configuration file is empty")
        else:
            self._time = float(self._line[2])

        self._line = self._conf.readline().split()
        self._box = np.array(self._line[2:], dtype=float)
        
        self._line = self._conf.readline().split()
        self._conf_energy = np.array(self._line[2:], dtype=float) # keep the energy info (cause it's interesting for remd)

        #we make lists here because we don't know how big the configuration is (that is in the topology)
        positions = []
        a1s = []
        a3s = []
        self._line = self._conf.readline().split()
        while self._line[0] != "t":
            positions.append(self._line[0:3])
            a1s.append(self._line[3:6])
            a3s.append(self._line[6:9])
            self._line = self._conf.readline().split()
            if self._line == []:
                break

        self._configuration = base_array(self._time, self._box, self._conf_energy, np.array(positions, dtype=float), np.array(a1s, dtype=float), np.array(a3s, dtype=float))
        self._len = len(positions)

        return (self._configuration)
    

    def read(self, n_skip=0):
        """
        Read the next configuration in the trajectory.

        Parameters:
            <optional> n_skip (int): skip the next n configurations before returning

        Returns:
            system (base_array): a set of numpy arrays containing positions and orientations of each nucleotides
        """
        if not self._time:
            if n_skip == 0:
                return self._read_first()
            else:
                self._read_first()
                n_skip -= 1

        if n_skip > 0:
            for i in range((n_skip * (self._len + 3)) -1):
                self._conf.readline()
            self._line = self._conf.readline().split() # get the time 

        if  len(self._line) == 0:
            return False
        else:
            self._time = float(self._line[2])
            self._configuration.time = self._time

        self._conf.readline() # we already know the box

        #self._conf.readline() # still don't care about the energy
        self._line = self._conf.readline().split()
        self._conf_energy = np.array(self._line[2:], dtype=float) # keep the energy info (cause it's interesting for remd)
        self._configuration.energy = self._conf_energy

        self._line = self._conf.readline().split() #first line of the configuration
        for i in range(self._len):
            if (len(self._line) == 0 or "t" in self._line) and i < self._len:
                print("ERROR: Reader encountered a partial configuration with only {} lines ({} expected)".format(i, self._len))
            self._configuration.positions[i] = np.array(self._line[0:3], dtype=float)
            self._configuration.a1s[i] = np.array(self._line[3:6], dtype=float)
            self._configuration.a3s[i] = np.array(self._line[6:9], dtype=float)
            self._line = self._conf.readline().split()
        return (self._configuration)
    



        
