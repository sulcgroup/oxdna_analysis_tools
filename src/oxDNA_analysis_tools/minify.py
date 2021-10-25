from oxDNA_analysis_tools.UTILS.readers import ErikReader, cal_confs
import argparse
from os import remove, path
from numpy import round


def main():
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Compress given configuration.")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('outfile',    type=str, nargs=1, help='minified file')

    parser.add_argument('-a', action = 'store_true', help='Discard a vectors.')
    parser.add_argument('-p', type=int, nargs=1,  help='Round positions and orientations to the specified number of digits.')

    args = parser.parse_args()

    traj_file = args.trajectory[0]
    out = args.outfile[0]
    # get the number of configurations
    n_confs = cal_confs(traj_file)

    try: # make sure there is no out file
        remove(out)
    except:
        pass
    
    with ErikReader(traj_file) as reader:
        for i in range(n_confs):
            print(i+1,":",n_confs)
            # Erik reader ignores velocities 
            system = reader.read()
            if args.p: # round positions
                system.positions = round(system.positions, args.p[0])
                system.a1s = round(system.a1s, args.p[0])
                system.a3s = round(system.a3s, args.p[0])
            if args.a: # discard a vectors
                system.a1s -= system.a1s 
                system.a3s -= system.a3s
            # output conf
            system.write_append(out)

if __name__ == '__main__':
    main()