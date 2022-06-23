#!/usr/bin/env python

from os import path
from sys import exit, stderr
from typing import List
import argparse

#checking dependencies to make sure everything is correct
def check_dependencies(to_check:List[str]):
    flag = False
    dependencies = {
        "numpy": 1.14,
        "matplotlib": 3.0,
        "Bio": 1.73,
        "sklearn": 0.21,
        "oxpy": 3.2,
    }
    real_names = {
        "numpy": "Numpy",
        "matplotlib": "MatPlotLib",
        "Bio": "BioPython",
        "sklearn": "SciKit-Learn",
        "oxpy": "oxpy"
    }
    websites = {
        "numpy": "numpy.org", 
        "matplotlib": "matplotlib.org",
        "Bio": "biopython.org",
        "sklearn": "scikit-learn.org",
        "oxpy": ""
    }

    #get version of this package
    oat = __import__('oxDNA_analysis_tools')
    print("INFO: oxDNA_analysis_tools version: {}".format(oat.__version__), file=stderr)
    print("INFO: running config.py installed at: ", path.realpath(__file__), file=stderr)


    #check python version
    if "python" in to_check:
        from sys import version_info
        ver = '.'.join([str(i) for i in version_info[0:2]])
        print("INFO: Python version: {}".format('.'.join([str(i) for i in version_info[0:3]])), file=stderr)
        if version_info < (3, 6):
            flag = True
            print("WARNING: Some scripts will not run with Python versions earler than 3.6.  You have {}, please update your environment", file=stderr)

    #check packages
    for package in to_check:
        if package == "python": continue
        try:
            mod = __import__(package)
            print("INFO: Package {} found. Version: {}".format(real_names[package], mod.__version__), file=stderr)
        except:
            if package == "Bio":
                try:
                    mod = __import__('bio')
                    print("INFO: Package {} found. Version: {}".format(real_names[package], mod.__version__), file=stderr)
                except:
                    flag = True
                    print("ERROR: Unable to find package {}.  Please check your environment or follow the installation instructions at {}".format(real_names[package], websites[package]), file=stderr)
                    continue
            else:
                flag = True
                print("ERROR: Unable to find package {}.  Please check your environment or follow the installation instructions at {}".format(real_names[package], websites[package]), file=stderr)
                continue
        ver = float('.'.join(mod.__version__.split(".")[0:2]))
        if ver < dependencies[package]:
            flag = True
            print("WARNING: Your version for package {} is {}.  This tool was tested using {}.  You may need to update your environment".format(real_names[package], ver, dependencies[package]), file=stderr)

    if flag:
        print("WARNING: Some packages need to be installed/updated.", file=stderr)
    else:
        print("INFO: No dependency issues found.", file=stderr)

    return flag

def set_chunk_size(chunk_size:int):
    with open(path.realpath(__file__).strip('config.py')+"UTILS/chunksize.py", 'w') as f:
        f.write("CHUNKSIZE = "+str(chunk_size))

def get_chunk_size():
    try:
        from oxDNA_analysis_tools.UTILS.chunksize import CHUNKSIZE
        print("INFO: Analyses will be computed in chunks of {} configurations at a time".format(CHUNKSIZE), file=stderr)
        print("INFO: You can modify this number by running oat config -n <number>, which will be persistent between analyses.", file=stderr)
    except:
        raise Exception("Unable to read chunksize from file. UTILS/chunksize.py should contain a line like CHUNKSIZE = 100")

def main():
    parser = argparse.ArgumentParser(description='Configure oxDNA_analysis_tools')
    parser.add_argument('-n', '--chunk_size', type=int, help='Number of configurations to per chunk.  Persistent across analyses.')
    args = parser.parse_args()
    if args.chunk_size:
        set_chunk_size(args.chunk_size)
        print("INFO: future analyses will calculate in blocks of {} confs at a time".format(args.chunk_size), file=stderr)

    check_dependencies(["python", "numpy", "matplotlib", "Bio", "sklearn", "oxpy"])

    print()
    get_chunk_size()

if __name__ == '__main__':
    main()
    
