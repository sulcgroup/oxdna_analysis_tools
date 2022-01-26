from copyreg import pickle
from oxDNA_analysis_tools.UTILS.micha_reader import MichaReader
import argparse
from os import getenv, path
from pickle import dumps


def main():
    #handle commandline arguments
    parser = argparse.ArgumentParser(prog = path.basename(__file__), description="Creates an index file of the provided trajectory.")
    parser.add_argument('top', type=str, nargs=1, help="The trajectory file to align")
    parser.add_argument('traj', type=str, nargs=1, help="The trajectory file to align")
    args = parser.parse_args()

    # Parse command line arguments
    top_file = args.top[0]
    traj_file = args.traj[0]

    print("indexing provided trajectory")
    # initiate the reader
    reader = MichaReader(top_file, traj_file)
    print("done")


if __name__ == '__main__':
    main()