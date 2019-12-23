#!/usr/bin/env python

#Created by: Erik Poppleton
#Date: 6/29/18
#Python2
#Converts the forces file printed out by tiamat2oxdna to a pairs file containing all designed H-bonds

from __future__ import print_function
import sys

infile = sys.argv[1]
with open(infile) as f:
    for line in f:
        if line.startswith("particle"):
            a = line.split()[2]
        if "ref_particle" in line:
            b = line.split()[2]
        if "}" in line:
            if int(a) < int(b):
                print(a, b, sep = ' ')
            a = -1
            b = -1