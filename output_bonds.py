#!/usr/bin/env python
# from original oxDNA repository
#A utility that prints out the number of hydrogen bonds between different strands in the system

import numpy as np
from os import environ, path
from sys import stderr, exit
from UTILS.readers import LorenzoReader2, get_input_parameter
import subprocess
import tempfile

from config import set_analysis_path
PROCESSPROGRAM = set_analysis_path()

command_for_data =  'analysis_data_output_1={ \n name = stdout \n print_every = 1 \n col_1 = { \n type=pair_energy \n} \n}'

def output_bonds (inputfile, system):
	tempfile_obj = tempfile.NamedTemporaryFile()
	system.print_lorenzo_output(tempfile_obj.name,'/dev/null')
	tempfile_obj.flush()

	launchargs = [PROCESSPROGRAM,inputfile ,'trajectory_file='+tempfile_obj.name,command_for_data]

	myinput = subprocess.run(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	out = myinput.stdout.strip()
	err = myinput.stderr.strip()
	for line in err.split('\n'):
		if "CRITICAL" in line or "ERROR" in line:
			print(err, file=stderr)
			exit(1)
	return out

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="List all the interactions between nucleotides")
	parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
	parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
	parser.add_argument('-v', type=str, nargs=1, dest='outfile', help='if you want instead average per-particle energy as a viewer JSON')
	args = parser.parse_args()

	from config import check_dependencies
	check_dependencies(["python", "numpy"])

	traj_file = args.trajectory[0]
	inputfile = args.inputfile[0]
	try:
		outfile = args.outfile[0]
		visualize = True
	except:
		visualize = False

	top_file = get_input_parameter(inputfile, "topology")
	if "RNA" in get_input_parameter(inputfile, "interaction_type"):
		environ["OXRNA"] = "1"
	else:
		environ["OXRNA"] = "0"
	import UTILS.base #this needs to be imported after the model type is set


	myreader = LorenzoReader2(traj_file,top_file)
	mysystem = myreader._get_system()

	energies = np.zeros(mysystem.N)
	count = 0

	while mysystem != False:
		out = output_bonds(inputfile, mysystem)
		if visualize:
			for line in out.split('\n'):
				if not line.startswith('#'):
					line = [float(l) for l in line.split(' ')]
					energies[int(line[0])] += sum(line[2:])
					energies[int(line[1])] += sum(line[2:])
		else:
			print(out)

		count += 1
		mysystem = myreader._get_system()

	if visualize:
		energies *= (41.42/count)
		with open(outfile, "w+") as file:
			file.write("{\n\"Energy (pN nm)\" : [")
			file.write(str(energies[0]))
			for n in energies[1:]:
				file.write(", {}".format(n))
			file.write("] \n}")




