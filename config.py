from os import path
from sys import exit, stderr

INBOXING_REFERENCE_PARTICLE = 0

PROCESSPROGRAM = '/home/erik/Simulations/oxdna-code/oxDNA/bin/DNAnalysis'

if not path.isfile(PROCESSPROGRAM):
	print ("ERROR: Cannot execute DNAnalysis program. Please edit config.py to point to your compiled DNAnalysis. Current target:", PROCESSPROGRAM, file=stderr)
	exit()
