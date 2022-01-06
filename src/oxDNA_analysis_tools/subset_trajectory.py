from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2
from oxDNA_analysis_tools.UTILS import base
import argparse
from sys import stderr
import os

def subset_traj(system, ids):
    """
    Extracts specified parts of a configuration into separate configuration files

    Parameters:
        system (System): the system object representing a single configuration
        ids (list[]): a list of lists of particle ids to extract

    Returns:
        new_systems (list[System]): a list of systems representing the extracted configurations
    """
    new_systems = [base.System(system._box) for _ in ids]
    for s in new_systems:
        s._time = system._time
    system.inbox()
    for s in system._strands:
        tmp_strands = [None for _ in ids]
        for n in s._nucleotides:
            for i, idx in enumerate(ids):
                if n.index in idx:
                    if not tmp_strands[i]:
                        if s.index >= 0:
                            tmp_strands[i] = base.Strand()
                        elif s.index < 0:
                            tmp_strands[i] = base.Peptide()

                    tmp_strands[i].add_nucleotide(n)

        for st, sys in zip(tmp_strands, new_systems):
            if st:
                sys.add_strand(st, check_overlap=False)

    return new_systems

def main():
    #command line arguments
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Extracts parts of a structure into separate trajectories")
    parser.add_argument('trajectory', type=str, nargs=1, help="The trajectory file to subset")
    parser.add_argument('topology', type=str, nargs=1, help="The topology file corresponding to the trajectory")
    parser.add_argument('-i', '--index', metavar='index', action='append', nargs=2, help='A space separated index file and the associated output file name.  This can be called multiple times')
    args = parser.parse_args()

    top_file  = args.topology[0]
    traj_file = args.trajectory[0]
    index_files = [i[0] for i in args.index]
    output_files = [i[1] for i in args.index]
    ids = []
    for i in index_files:
        with open(i) as f:
            data = f.readline().split()
            try:
                data = [int(i) for i in data]
            except:
                print("ERROR: The index file {} must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button".format(i))
        ids.append(data)

    r = LorenzoReader2(traj_file, top_file)
    system = r._get_system()

    new_systems = subset_traj(system, ids)

    for sys, o in zip(new_systems, output_files):
        print("INFO: writing subset to {}".format(o), file=stderr)
        sys.print_lorenzo_output(o+".dat", o+".top")

    system = r._get_system()
    while system:
        print("INFO: working on t = {}".format(system._time), file=stderr)
        new_systems = subset_traj(system, ids)

        for sys, o in zip(new_systems, output_files):
            sys.print_traj_output(o+".dat", '/dev/null')

        system = r._get_system()

if __name__ == '__main__':
    main()