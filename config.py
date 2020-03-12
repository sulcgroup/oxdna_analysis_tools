#!/usr/bin/env python

from os import path
from sys import exit, stderr

#The inboxing method used here moves everything relative to a reference particle because of oxDNA's fix_diffusion artifacts
#If the largest dimension of the structure is further away from the reference particle than half the box size the structure will be broken by the PBC
#This particle needs to be in the first strand to avoid having to map nucleotides to strands, which would significantly slow down many scripts
def set_reference():
	INBOXING_REFERENCE_PARTICLE = 0
	print("INFO: Particle {} will be used as inboxing reference".format(INBOXING_REFERENCE_PARTICLE), file=stderr)
	print("INFO: The reference should be in the middle of the structure and can be set in config.py", file=stderr)
	return INBOXING_REFERENCE_PARTICLE

#set the path to your compiled copy of DNAnalysis here
def set_analysis_path():
	PROCESSPROGRAM = '/home/erik/Simulations/oxdna-code/oxDNA/bin/DNAnalysis'

	if not path.isfile(PROCESSPROGRAM):
		print ("ERROR: Cannot execute DNAnalysis program. Please edit config.py to point to your compiled DNAnalysis. Current target:", PROCESSPROGRAM, file=stderr)
		exit()
	return PROCESSPROGRAM

#checking dependencies to make sure everything is correct
def check_dependencies(to_check):
	flag = False
	dependencies = {
		"numpy": 1.14,
		"matplotlib": 3.0,
		"Bio": 1.73,
		"sklearn": 0.21,
		"pathos": 0.2,
	}
	real_names = {
		"numpy": "Numpy",
		"matplotlib": "MatPlotLib",
		"Bio": "BioPython",
		"sklearn": "SciKit-Learn",
		"pathos": "Pathos"
	}
	websites = {
		"numpy": "numpy.org", 
		"matplotlib": "matplotlib.org",
		"Bio": "biopython.org",
		"sklearn": "scikit-learn.org",
		"pathos": "pypi.org/project/pathos/"
	}

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
			flag = True
			print("ERROR: Unable to find package {}.  Please check your environment or follow the installation instructions at {}".format(real_names[package], websites[package]), file=stderr)
		ver = float('.'.join(mod.__version__.split(".")[0:2]))
		if ver < dependencies[package]:
			flag = True
			print("WARNING: Your version for package {} is {}.  This tool was tested using {}.  You may need to update your environment".format(real_names[package], ver, dependencies[package]), file=stderr)

	if flag:
		print("WARNING: Some packages need to be installed/updated.", file=stderr)
	else:
		print("INFO: No dependency issues found.", file=stderr)

	return flag

if __name__ == "__main__":
	check_dependencies(["python", "numpy", "matplotlib", "Bio", "sklearn", "pathos"])
	p = set_analysis_path()
	print("INFO: DNAnalysis found at:", p, file=stderr)