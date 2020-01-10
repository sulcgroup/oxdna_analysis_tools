#!/usr/bin/env python3
import sys
import subprocess

# Created by Hao Liu 
# Date 01/22/2019 
# A short script generating force file from the given .dat and .top

'''
***********************************************************************************************************************************
#If you want the bash command to get sorted bonding pairs file:
python ~/software/oxDNA/oxdna-code/oxDNA/UTILS/output_bonds.py input_DNA_left40 right40.dat | gawk '{if ($7 < -0.01){print $1 " " $2}}' > right40outputbonds.txt
cat right10outputbonds.txt | gawk '{if($1 > max){max = $1} arr[$1]=$0} END{for(i = 1; i<=max; i++){print arr[i]}}' > sortedright10.txt
***********************************************************************************************************************************
'''

if len(sys.argv) < 4:
    print("The usage is %s , input conf topology outputfilename"%(sys.argv[0]))
    sys.exit()

#PROCESSDIR = '/home/epopplet/oxDNA/oxdna-code/oxDNA/UTILS/'
PROCESSDIR = '/home/erik/Simulations/oxdna-code/oxDNA/UTILS/'

#read data from files
inputfile = sys.argv[1]
traj_file = sys.argv[2]
topfile = sys.argv[3]
output = sys.argv[4]

#Lanch output_bonds.py 
#Theoretically "run" method is faster than "popen" method. 
print("Launch output_bonds.py ...")
launchargs = ['python2', PROCESSDIR + 'output_bonds.py',inputfile,traj_file]
myinput = subprocess.run(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
out = myinput.stdout.split("\n")
err = myinput.stderr.split("\n")

#Error reporting 
for line in err:
    if "CRITICAL" in line or "ERROR" in line or "python:" in line:
        print (line)
        sys.exit(1) 

#Find out the forming bonds series 
print("Analyze the output...")
Bonded = {}
for i in out:
    splitline = i.strip("\n").split()
    try:
        HB = float(splitline[6])
    except:
        continue
    if HB < -0.001:
        ntid0 = int(splitline[0])
        ntid1 = int(splitline[1])
        if ntid0 not in Bonded:
            Bonded[ntid0] = ntid1
        if ntid1 not in Bonded:
            Bonded[ntid1] = ntid0

lines = []
mutual_trap_template = '{ \ntype = mutual_trap\nparticle = %d\nstiff = 0.9\nr0 = 1.2\nref_particle = %d\nPBC=1\n}\n'
for key in sorted(Bonded):
        print(key, Bonded[key])
        from_particle_id = key
        to_particle_id = Bonded[key]
        if from_particle_id >= 0:
                lines.append(mutual_trap_template % (from_particle_id,to_particle_id))
                lines.append(mutual_trap_template % (to_particle_id,from_particle_id))

with open(output, "w") as file:
        file.writelines(lines)
        print("Job finished. Force file created")
        


