
#This script should not assume that center of the helix is on the line between two bonded bases
#Created by: Petr Sulc
import oxDNA_analysis_tools.UTILS.base
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import math

prune = 1


def my_norm(vec):
	return np.sqrt(np.dot(vec,vec))

def get_turn_per_bp(s,nucid,hel_vector):
    #this function returns the distance between midpoints of hb sites
    box = s._box
    i1A = nucid
    i1B = len(s._strands[1]._nucleotides)-1-nucid
    firstA = s._strands[0]._nucleotides[i1A]
    firstB = s._strands[1]._nucleotides[i1B] 	
    secondA = s._strands[0]._nucleotides[i1A+1]
    secondB = s._strands[1]._nucleotides[i1B-1]
    first_midpos = (firstA.get_pos_back() - firstB.get_pos_back())   
    second_midpos = (secondA.get_pos_back() - secondB.get_pos_back())   
    first_midpos /= math.sqrt(np.dot(first_midpos,first_midpos))
    second_midpos /= math.sqrt(np.dot(second_midpos,second_midpos))
    first_midpos = first_midpos - np.dot(hel_vector,first_midpos) * hel_vector
    second_midpos = second_midpos - np.dot(hel_vector,second_midpos) * hel_vector
    first_midpos /= math.sqrt(np.dot(first_midpos,first_midpos))
    second_midpos /= math.sqrt(np.dot(second_midpos,second_midpos))
 
    #angle = math.acos(np.dot(first_midpos,hel_vector)) - math.acos(np.dot(second_midpos,hel_vector))
    angle = math.acos(np.dot(first_midpos,second_midpos))
    return angle 

def fit_plane(points):
	#fits plane through points, return normal to the plane:
	rc= np.array(np.zeros(3))
	A = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])
	for point in points:
		rc += point
	rc /= len(points)
	#print(rc)

	for point in points:
		point  = point - rc
		for i in range(3):
			for j in range(3):
				A[i][j] += point[i]*point[j]

	#this could be made faster by using a symmetric eigensolver
	vals, vecs = np.linalg.eigh(A) 		
	#print vals, vecs		
	return vecs[:,0]

def get_local_axis(s,base_id,local_length=4):
	#i = base_id
	posBs = []
	posAs = []
	vectorsA = []
	back_poses = []
	for i in range(base_id,base_id+local_length):	
		nucA = s._nucleotides[i]
		nucB= s._nucleotides[s._N - i-1]
		nucAc = s._nucleotides[i+1]
		nucBc = s._nucleotides[s._N - i- 1 -1]
		midpoint_A = 0.5*(nucA.get_pos_stack() + nucB.get_pos_stack())	
		midpoint_B = 0.5*(nucBc.get_pos_stack() + nucAc.get_pos_stack())
		vectorsA.append(-midpoint_A + midpoint_B)
		posAs.append(nucA.get_pos_back())
		posBs.append(nucB.get_pos_back())
		posAs.append(nucAc.get_pos_back())
		posBs.append(nucBc.get_pos_back())

		back_poses.append(nucA.get_pos_back() - nucAc.get_pos_back() )
		back_poses.append(nucB.get_pos_back() - nucBc.get_pos_back() )
		back_poses.append(midpoint_A - midpoint_B)	
		guess =  (-midpoint_A + midpoint_B) #this is the vector pointing from midpoint of the first bp to the last bp, an estimate of the vector
		guess = guess / my_norm(guess)
		plane_vector = fit_plane(back_poses)
		
		if(np.dot(guess,plane_vector) < 0):
			#print 'Debug Warning, plane vector pointing in opposite direction'
			plane_vector = -1. * plane_vector
	hel_pos = []
	
	#project to the plane
	apos = posAs[0]
	bpos = posBs[0]
	apos = apos - np.dot(apos,plane_vector) * plane_vector
	bpos = bpos - np.dot(bpos,plane_vector) * plane_vector
	bp_vec = -apos + bpos
	midpointA = 0.5 * (apos + bpos)
	perpendicular_vecA = np.cross(bp_vec,plane_vector)
	perpendicular_vecA /= my_norm(perpendicular_vecA)
	apos = posAs[1]
	bpos = posBs[1]
	apos = apos - np.dot(apos,plane_vector) * plane_vector
	bpos = bpos - np.dot(bpos,plane_vector) * plane_vector
	bp_vec = -apos + bpos
	midpointB = 0.5 * (apos + bpos)
	perpendicular_vecB = np.cross(bp_vec,plane_vector)
	perpendicular_vecB /= my_norm(perpendicular_vecB)
	mat = np.array([perpendicular_vecA, -perpendicular_vecB ]).transpose()
	y = midpointB - midpointA
	t,c = np.linalg.lstsq(mat,y)[0]
	if my_norm(midpointA + t*perpendicular_vecA  -  (midpointB + c *perpendicular_vecB)) > 1.e-6 :
		print('Debug Error in finding common intersection point',midpointA + t*perpendicular_vecA , mipointB + c *perpendicular_vecB)
	hel_position = midpointA + t*perpendicular_vecA	
	hel_pos.append(hel_position)	

	print('Position of bp',i,' is ',hel_position)
		
	final_hel_pos = np.zeros(3)
	for pos in hel_pos:
		final_hel_pos += pos
	final_hel_pos /= len(hel_pos)	
					

	return plane_vector,final_hel_pos


def get_inclination(s,plane_vector,i):
	posa = s._nucleotides[i].get_pos_back()
	posb = s._nucleotides[s._N - 1 - i].get_pos_back()
	a1a = s._nucleotides[i]._a1
	a1b = -s._nucleotides[s._N - 1 - i]._a1

	inclination_vector = posa - posb
	inclination_vector /= my_norm(inclination_vector)
	angle = np.rad2deg(math.acos( np.dot(inclination_vector,plane_vector ) ))
	angleA = np.rad2deg(math.acos( np.dot(a1a,plane_vector ) ))
	angleB = np.rad2deg(math.acos( np.dot(a1b,plane_vector ) ))
	angle = 0.5*(angleA + angleB)
	return angle

	
def get_data_with_local_axis(s,first_base=0,last_base=-1,only_plane_vector=False):
	if last_base == -1:
		last_base = s._N/2 - 1
	i = first_base
	nucA = s._nucleotides[i]
	nucB= s._nucleotides[s._N - i-1]
	midA = (nucA.get_pos_back() + nucB.get_pos_back() ) / 2.
	i = last_base
	nucAc = s._nucleotides[i]
	nucBc = s._nucleotides[s._N - i -1]
	midAc = (nucAc.get_pos_back() + nucBc.get_pos_back() ) / 2.
	guess =  (-midA + midAc) #this is the vector pointing from midpoint of the first bp to the last bp, an estimate of the vector
	guess = guess / my_norm(guess)

	vectorsA = []
	posAs = []
	posBs = []
	back_poses = []	
	for i in range(first_base,last_base):
		nucA = s._nucleotides[i]
		nucB= s._nucleotides[s._N - i-1]
		nucAc = s._nucleotides[i+1]
		nucBc = s._nucleotides[s._N - i- 1 -1]
		midpoint_A = 0.5*(nucA.get_pos_back() + nucB.get_pos_back())	
		midpoint_B = 0.5*(nucBc.get_pos_back() + nucAc.get_pos_back())
		vectorsA.append(-midpoint_A + midpoint_B)
		posAs.append(nucA.get_pos_back())
		posBs.append(nucB.get_pos_back())

		back_poses.append(nucA.get_pos_back() - nucAc.get_pos_back() )
		back_poses.append(nucB.get_pos_back() - nucBc.get_pos_back() )
		
		if(i == last_base - 1):
			posAs.append(nucAc.get_pos_back())
			posBs.append(nucBc.get_pos_back())
	
	#now we find the point where the helical vecotr originates

	#now, we can do some statistics: inclination, angle_per_bp, width, rise_per_bp
	inclinations = []
	angles = []
	widths = []
	rises = []
	dislocations = []

	for i in range(len(posAs)-1):
		#print 'Debug, len is ',len(posAs), 'i is',i
		plane_vector,final_hel_pos = get_local_axis(s,i)
		#print 'Debug: local plane vector is : ',plane_vector
		vec_bp = -posAs[i] + posBs[i]
		vec_bp_norm = -vec_bp / my_norm(vec_bp)	
		inclinations.append(90-get_inclination(s,plane_vector,i))
		#inclinations.append(90. - np.rad2deg(math.acos(np.dot(vec_bp_norm,plane_vector))))
		
		#now lift the plane so that it is in the plane of the bp
		midpointA = 0.5 * (posAs[i] + posBs[i])
		midpointB = 0.5 * (posAs[i+1] + posBs[i+1])
		dA = np.dot(plane_vector,midpointA)
		dB = np.dot(plane_vector,midpointB)
		rises.append(dB - dA)
			
		
		midpointA_projection = midpointA - plane_vector * np.dot(midpointA,plane_vector)
		midpointB_projection = midpointB - plane_vector * np.dot(midpointB,plane_vector)
		avec = midpointA_projection - final_hel_pos
		bvec = midpointB_projection - final_hel_pos
		dislocations.append(my_norm(avec))
		dislocations.append(my_norm(bvec))
		#avec = avec / my_norm(avec)
		#bvec = bvec /  my_norm(bvec)
		#angle = np.rad2deg(math.acos(np.dot(avec,bvec)))
		#vec_bp_2 = -posAs[i+1] + posBs[i+1]
		#vec_bp_2 /= my_norm(vec_bp_2)
		#angle = 180. - np.rad2deg(math.acos(np.dot(vec_bp_norm,vec_bp_2)))
		vec_bp_2 = -posAs[i+1] + posBs[i+1]
		vec_bp_2 /= my_norm(vec_bp_2)
		#do projection to the helical axis
		vec_bp_norm = vec_bp_norm - plane_vector * np.dot(vec_bp_norm,plane_vector)
		vec_bp_2 = vec_bp_2 - plane_vector * np.dot(vec_bp_2,plane_vector)
		vec_bp_2 /= my_norm(vec_bp_2)
		vec_bp_norm /= my_norm(vec_bp_norm)
		
		angle = 180. - np.rad2deg(math.acos(np.dot(vec_bp_norm,vec_bp_2)))
		 
		#angle = np.rad2deg(get_turn_per_bp(s,i,plane_vector))
		angles.append(angle)
		
		posA_projection = posAs[i] - plane_vector * np.dot(posAs[i],plane_vector)
		posB_projection = posBs[i] - plane_vector * np.dot(posBs[i],plane_vector)
		wA = my_norm(posA_projection - final_hel_pos)
		wB = my_norm(posB_projection - final_hel_pos)
		widths.append(wA)
		widths.append(wB)	

	return plane_vector, final_hel_pos, inclinations, angles, widths, rises, dislocations	

#needs to be modified so it can handle strand id and first and last base in duplex
def get_RNA_axis(s,first_base, last_base, cfirst_base, clast_base, only_plane_vector=False):

	nucA = s._nucleotides[first_base]
	nucB= s._nucleotides[cfirst_base]
	midA = (nucA.get_pos_back() + nucB.get_pos_back() ) / 2.
	nucAc = s._nucleotides[last_base]
	nucBc = s._nucleotides[clast_base]
	midAc = (nucAc.get_pos_back() + nucBc.get_pos_back() ) / 2.
	guess =  (-midA + midAc) #this is the vector pointing from midpoint of the first bp to the last bp, an estimate of the vector
	guess = guess / my_norm(guess)

	vectorsA = []
	posAs = []
	posBs = []
	back_poses = []	
	
	for i in range(0, last_base-first_base):
		nucA = s._nucleotides[first_base+i]
		nucB= s._nucleotides[cfirst_base-i]
		nucAc = s._nucleotides[first_base+i+1]
		#print(nucA.index, nucAc.index)
		nucBc = s._nucleotides[cfirst_base-i-1]
		midpoint_A = 0.5*(nucA.get_pos_back() + nucB.get_pos_back())	
		midpoint_B = 0.5*(nucBc.get_pos_back() + nucAc.get_pos_back())
		vectorsA.append(-midpoint_A + midpoint_B)
		posAs.append(nucA.get_pos_back())
		posBs.append(nucB.get_pos_back())
		if(i == last_base - 1):
			posAs.append(nucAc.get_pos_back())
			posBs.append(nucBc.get_pos_back())
		back_poses.append(-nucA.get_pos_back() + nucAc.get_pos_back() )
		back_poses.append(-nucB.get_pos_back() + nucBc.get_pos_back() )

	#now we try to fit a plane through all these points
	plane_vector = fit_plane(back_poses)
	#print(plane_vector)
	#plane_vector = fit_plane(vectorsA)

	if(np.dot(guess/my_norm(guess),plane_vector) < 0):
		#print 'Warning, plane vector pointing in opposite direction'
		plane_vector = -1. * plane_vector
	if(np.rad2deg(math.acos(np.dot(plane_vector,guess))) > 20):
		print ('Warning, guess vector and plane vector have angles:', np.rad2deg(math.acos(np.dot(guess/my_norm(guess),plane_vector))))

	#now we find the point where the helical vecotr originates
	hel_pos = []
	for i in range(len(posAs)-1):
		#project to the plane
		apos = posAs[i]
		bpos = posBs[i]
		apos = apos - np.dot(apos,plane_vector) * plane_vector
		bpos = bpos - np.dot(bpos,plane_vector) * plane_vector
		bp_vec = -apos + bpos
		if np.linalg.norm(bp_vec) == 0:
			    continue
		midpointA = 0.5 * (apos + bpos)
		perpendicular_vecA = np.cross(bp_vec,plane_vector)
		perpendicular_vecA /= my_norm(perpendicular_vecA)
		apos = posAs[i+1]
		bpos = posBs[i+1]
		apos = apos - np.dot(apos,plane_vector) * plane_vector
		bpos = bpos - np.dot(bpos,plane_vector) * plane_vector
		bp_vec = -apos + bpos
		if np.linalg.norm(bp_vec) == 0:
			continue
		midpointB = 0.5 * (apos + bpos)
		perpendicular_vecB = np.cross(bp_vec,plane_vector)
		perpendicular_vecB /= my_norm(perpendicular_vecB)
		mat = np.array([perpendicular_vecA, -perpendicular_vecB ]).transpose()
		y = midpointB - midpointA
		t,c = np.linalg.lstsq(mat,y)[0]
		if my_norm(midpointA + t*perpendicular_vecA  -  (midpointB + c *perpendicular_vecB)) > 1.e-6 :
			print ('Error in finding common intersection point',midpointA + t*perpendicular_vecA , mipointB + c *perpendicular_vecB)
		hel_position = midpointA + t*perpendicular_vecA	
		#print ('Hel position from BP',i,'is ',hel_position)
		hel_pos.append(hel_position)	
				
	final_hel_pos = np.zeros(3)
	for pos in hel_pos:
		final_hel_pos += pos
	final_hel_pos /= len(hel_pos)				

	if only_plane_vector:
		return plane_vector, final_hel_pos#, back_poses, guess 
		
	#now, we can do some statistics: inclination, angle_per_bp, width, rise_per_bp
	inclinations = []
	angles = []
	widths = []
	rises = []
	dislocations = []

	for i in range(len(posAs)-1):
		vec_bp = -posAs[i] + posBs[i]
		vec_bp_norm = -vec_bp / my_norm(vec_bp)	
		
		inclinations.append(90-get_inclination(s,plane_vector,i))
		#inclinations.append(90. - np.rad2deg(math.acos(np.dot(vec_bp_norm,plane_vector))))
		
		#now lift the plane so that it is in the plane of the bp
		midpointA = 0.5 * (posAs[i] + posBs[i])
		midpointB = 0.5 * (posAs[i+1] + posBs[i+1])
		dA = np.dot(plane_vector,midpointA)
		dB = np.dot(plane_vector,midpointB)
		rises.append(dB - dA)
		#print 'MidpointA', midpointA, 'MidpointB', midpointB, 'Difference', dB - dA
			
		hel_pointA = final_hel_pos + dA * plane_vector
		hel_pointB = final_hel_pos + dB * plane_vector	
		
		midpointA_projection = midpointA - plane_vector * np.dot(midpointA,plane_vector)
		midpointB_projection = midpointB - plane_vector * np.dot(midpointB,plane_vector)
		avec = midpointA_projection - final_hel_pos
		bvec = midpointB_projection - final_hel_pos
		dislocations.append(my_norm(avec))
		dislocations.append(my_norm(bvec))
		#avec = avec / my_norm(avec)
		#bvec = bvec /  my_norm(bvec)
		#angle =  np.rad2deg(math.acos(np.dot(avec,bvec)))
		#angle =  np.rad2deg(math.acos(np.dot(avec,bvec)))

		vec_bp_2 = -posAs[i+1] + posBs[i+1]
		vec_bp_2 /= my_norm(vec_bp_2)
		#do projection to the helical axis
		vec_bp_norm = vec_bp_norm - plane_vector * np.dot(vec_bp_norm,plane_vector)
		vec_bp_2 = vec_bp_2 - plane_vector * np.dot(vec_bp_2,plane_vector)
		vec_bp_2 /= my_norm(vec_bp_2)
		vec_bp_norm /= my_norm(vec_bp_norm)
		
		angle = 180. - np.rad2deg(math.acos(np.dot(vec_bp_norm,vec_bp_2)))
		 
		#angle = np.rad2deg(get_turn_per_bp(s,i,plane_vector))
		angles.append(angle)
		
		posA_projection = posAs[i] - plane_vector * np.dot(posAs[i],plane_vector)
		posB_projection = posBs[i] - plane_vector * np.dot(posBs[i],plane_vector)
		wA = my_norm(posA_projection - final_hel_pos)
		wB = my_norm(posB_projection - final_hel_pos)
		widths.append(wA)
		widths.append(wB)	

	return plane_vector, final_hel_pos, inclinations, angles, widths, rises, dislocations	

def get_DNA_axis (s,first_base, last_base, cfirst_base, clast_base, only_plane_vector=False):
	vec = np.empty((last_base-first_base, 3))
	for i in range(last_base-first_base):
		nucA = s._nucleotides[first_base+i]
		nucAc= s._nucleotides[cfirst_base-i]
		midA = (nucA.get_pos_base() + nucAc.get_pos_base() ) / 2.
		nucB = s._nucleotides[first_base+i+1]
		nucBc = s._nucleotides[cfirst_base-i-1] 
		midB = (nucB.get_pos_base() + nucBc.get_pos_base()) / 2.
		vec[i] = (midB - midA)/(np.linalg.norm(midB - midA))

	#print(s._nucleotides[cfirst_base].index, s._nucleotides[cfirst_base].get_pos_base())

	pos = (s._nucleotides[first_base].get_pos_base() + s._nucleotides[cfirst_base].get_pos_base() ) / 2.
	vector = np.mean(vec, axis=0)
	return(vector, pos)

#-----------------

def get_bb_dist(s,nucid): 
    #this function measures the distance between backbone sites of the opposite nucleotides
    box = s._box
    i1A = nucid
    i1B = len(s._strands[1]._nucleotides)-1-nucid
    firstA = s._strands[0]._nucleotides[i1A]
    firstB = s._strands[1]._nucleotides[i1B] 	
    first_midpos = (firstA.get_pos_back() - firstB.get_pos_back())   
    return math.sqrt(np.dot(first_midpos,first_midpos)) 

def get_back_back_distance(s,nucid):
    box = s._box
    i1A = nucid
    i1B = len(s._strands[1]._nucleotides)-1-nucid
    firstA = s._strands[0]._nucleotides[i1A]
    firstB = s._strands[1]._nucleotides[i1B] 	
    secondA = s._strands[0]._nucleotides[i1A+1]
    secondB = s._strands[1]._nucleotides[i1B-1]
    first_midpos = (firstA.get_pos_back() - secondA.get_pos_back())   
    second_midpos = (firstB.get_pos_back() - secondB.get_pos_back())
    return math.sqrt( np.dot(first_midpos,first_midpos)), math.sqrt(np.dot(second_midpos,second_midpos))

def get_end_j(s,nid1,nid2):
    i1A = nid1
    i2A = nid2
    i1B = len(s._strands[1]._nucleotides) - nid1 - 1
    i2B = len(s._strands[1]._nucleotides) - nid2 - 1 
    
    firstA = s._strands[0]._nucleotides[i1A]
    firstB = s._strands[1]._nucleotides[i1B] 	
    lastA = s._strands[0]._nucleotides[i2A]
    lastB = s._strands[1]._nucleotides[i2B]
   
    first_midpos = (firstA.get_pos_base() + firstB.get_pos_base()) / 2.0  
    last_midpos = (lastA.get_pos_base() + lastB.get_pos_base()) / 2.0  
  
    box = s._box  
    r0N = last_midpos - first_midpos 
    r0N -= box * np.rint (r0N / box)
    return math.sqrt(np.dot(r0N,r0N)) / float((nid2 - nid1))



'''
if len(sys.argv) < 3:
    base.Logger.log("Usage is %s configuration topology [offset=4] [prune=1]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()

offset = 4
prune = 1
if len(sys.argv) >= 4:
    offset = int(sys.argv[3])

if len(sys.argv) >= 5:
    prune = int(sys.argv[4])
    


i1A = offset
i1B = len(s._strands[1]._nucleotides) - offset - 1

print  "#Nucleotides", i1A, i1B


L2 = 0.
l0 = 0.
Ll0 = 0.
niter = 0
rises = []
read_confs = 1
end_rises = []
angles = []
helix_widths = []
incsA = []
incsB = []

bbAs = []
bbBs = []


all_inclinations = []
all_widths = []
all_rises = []
all_dislocations = []
all_angles = []

read_confs = 0

while s:
    if(read_confs % prune != 0):
	read_confs += 1
    	s = l.get_system()
	continue     
    end_rises.append(get_end_j(s,i1A,i1B))
 
    box = s._box  
    plane_vector, final_hel_pos, inclinations, angles, widths, rises, dislocations = get_axis(s,first_base=offset, last_base=len(s._strands[0]._nucleotides)-offset-2,only_plane_vector=False)
    all_inclinations = all_inclinations + inclinations
    all_angles = all_angles + angles
    all_widths = all_widths + widths
    all_rises = all_rises + rises
    all_dislocations = all_dislocations + dislocations
    #print dislocations
    print angles
    print inclinations
    print '#Plane vector:',plane_vector
    print '#Position of helix',final_hel_pos
    for j in range (offset,len(s._strands[0]._nucleotides) - offset - 2):
	jA = j
	
	bbA,bbB = get_back_back_distance(s,jA)
	bbAs.append(bbA)
	bbBs.append(bbB)

    s = l.get_system()
    niter += 1


print '# ',sys.argv[1]
print '# read configurations: ', niter
print 'Rise per bp: %f (+/- %f), Rise (from end-end-dist): %f, Angle (per bp): %f (+/- %f) = pitch %f (+/- %f), width %f (+/- %f), inclination : %f (+- %f), back-back distances:A %f (+- %f), B: %f (+- %f), dislocations: %f (+- %f) \n' % ( np.mean(all_rises)  ,np.std(all_rises) ,np.mean(end_rises),np.mean(all_angles),np.std(all_angles),360. / np.mean(all_angles) , np.std(360. / np.array(all_angles)), np.mean(all_widths), np.std(all_widths) ,(np.mean(all_inclinations)),(np.std(all_inclinations)),  np.mean(bbAs),np.std(bbAs),np.mean(bbBs),np.std(bbBs), np.mean(all_dislocations), np.std(all_dislocations))
'''

