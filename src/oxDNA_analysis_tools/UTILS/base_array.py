import numpy as np

"""
A class that contains a series of numpy arrays corresponding to the positions and orientations of all particles in an oxDNA system.  This does not have topology information.
"""
class base_array(object):
    """
    Create a base array that contains oxDNA configuration information

    Parameters
    ----------
    time : int
        The timestamp of the simulation.
    box : numpy.array(3)
        The box size of the simulation.
    positions : numpy.array((n, 3))
        The x, y, z particle positions of n particles in the system.
    a1s : numpy.array((n, 3))
        The x, y, z orientation of the a1 (pairing) vectors of n particles in the system.
    a3s: numpy.array((n, 3))
        The x, y, z orientation of the a3 (stacking) vectors of n particles in the system.

    Attributes
    ----------
    time : int
        The timestamp of the simulation.
    box : numpy.array(3)
        The box size of the simulation.
    positions : numpy.array((n, 3))
        The x, y, z particle positions of n particles in the system.
    a1s : numpy.array((n, 3))
        The x, y, z orientation of the a1 (pairing) vectors of n particles in the system.
    a3s: numpy.array((n, 3))
        The x, y, z orientation of the a3 (stacking) vectors of n particles in the system.
    
    Returns
    -------
    base_array 
        An object containing the oxDNA configuration information

    """ 
    def __init__(self, time, box, energy, positions, a1s, a3s):
        self.time = time
        self.box = box
        self.energy = energy
        self.positions = positions
        self.a1s = a1s
        self.a3s = a3s
    
    def inbox(self):
        """
        Modify the positions attribute such that all positions are inside the box.
        """
        def realMod (n, m):
            return(((n % m) + m) % m)

        def coord_in_box(p):
            p = realMod(p, self.box)
            return(p)

        def calc_PBC_COM(self):
            angle = (self.positions * 2 * np.pi) / self.box
            cm = np.array([[np.sum(np.cos(angle[:,0])), np.sum(np.sin(angle[:,0]))], 
            [np.sum(np.cos(angle[:,1])), np.sum(np.sin(angle[:,1]))], 
            [np.sum(np.cos(angle[:,2])), np.sum(np.sin(angle[:,2]))]]) / len(angle)
            return self.box / (2 * np.pi) * (np.arctan2(-cm[:,1], -cm[:,0]) + np.pi)

        target = np.array([self.box[0] / 2, self.box[1] / 2, self.box[2] / 2])
        center = calc_PBC_COM(self)

        self.positions += (target - center)
        old_poses = self.positions.copy()
        new_poses = coord_in_box(self.positions.copy())
        self.positions += (new_poses - old_poses)

    def conf_to_str(self):
        """
        Write array to an oxDNA configuration file

        Parameters
        ----------
        f : _io.TextIOWrapper
            An open Python file handler.
        """
        out = []
        out.append('t = {}\n'.format(int(self.time)))
        out.append('b = {}\n'.format(' '.join(self.box.astype(str))))
        out.append('E = {}\n'.format(' '.join(self.energy.astype(str))))
        for p, a1, a3 in zip(self.positions, self.a1s, self.a3s):
            out.append('{} {} {} 0.0 0.0 0.0 0.0 0.0 0.0\n'.format(' '.join(p.astype(str)), ' '.join(a1.astype(str)), ' '.join(a3.astype(str))))
        return ''.join(out)

    def _write_configuration(self, f):
        """
        Write a configuration line to a file

        Parameters
        ----------
        f: _io.TextIOWrapper
            An open Python file handler to write to.
        """
        f.write('t = {}\n'.format(int(self.time)))
        f.write('b = {}\n'.format(' '.join(self.box.astype(str))))
        f.write('E = {}\n'.format(' '.join(self.energy.astype(str))))
        for p, a1, a3 in zip(self.positions, self.a1s, self.a3s):
            f.write('{} {} {} 0.0 0.0 0.0 0.0 0.0 0.0\n'.format(' '.join(p.astype(str)), ' '.join(a1.astype(str)), ' '.join(a3.astype(str))))

    def write_new(self, filename):
        """
        Starts writing a new configuration file.  Will overwrite existing content.
        
        Parameters
        ----------
        filename: str
            The filename to write out to.  Should end in .dat, .conf or .oxdna for oxView compatibility
        """
        with open(filename, 'w') as f:
            self._write_configuration(f)

    def write_append(self, filename):
        """
        Appends to an existing trajectory file.

        Parameters
        ----------
        filename: str
            The filename to write out to.  Should end in .dat, .conf or .oxdna for oxView compatibility
        """
        with open(filename, 'a') as f:
            self._write_configuration(f)