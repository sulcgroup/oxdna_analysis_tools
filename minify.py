from UTILS.readers import ErikReader, cal_confs
from sys import argv, exit
from os import remove

if len(argv)<3:
	print("Discards velocity information from the given trajectory.")
	print("Usage is >minify.py traj_file out_file")
	exit()

# input files
traj_file,out = argv[1],argv[2]

# make sure there is no outfile
try:
	remove(out)
except:
	pass
# get the number of configurations
n_confs = cal_confs(traj_file)


with ErikReader(traj_file) as reader:
	for i in range(n_confs):
		print(i+1,":",n_confs)
		system = reader.read()
		system.write_append(out)