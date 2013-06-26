#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# check which person the seed has the highest similarity

import os, glob, subprocess, random, operator, time
from optparse import OptionParser


file_path = "/home/guoxing/morehouse_new_hap/fils/"
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--seed", type="string", dest="seedFile",help = "Input seed File Name", default="null")
parser.add_option("-r", "--ref", type="string", dest="refFile",help = "Input ref File Name", default="null")

(options, args) = parser.parse_args()

seed_file = options.seedFile
ref_file = options.refFile
seed_dict = {}
ref_dict = {}

# A from Father, B from Mother
class seed:
	def __init__(me, rsID, position, A):
		me.rsID = rsID
		me.position = position
		me.A = A
		
# ref class
class ref:
	def __init__(me, name, index, same2seed, similarity, snp_dict):
		me.name = name
		me.index = index
		me.same2seed = same2seed
		me.similarity = similarity
		me.snp_dict = snp_dict		
		
print "seed file: ", seed_file
print "ref file: ", ref_file

inputFile_seed = open(currentPath + seed_file, "r")
inputFile_ref = open(currentPath + ref_file, "r")
data_record_file_name = seed_file + "_record"
data_record_file = open(currentPath + data_record_file_name, "w")

total_seed_number = 0

for line in inputFile_seed:
	if not line.startswith("rsID"):
		total_seed_number += 1
		elements = line.strip().split()
		rsID = elements[0].strip()
		position = elements[1].strip()
		A = elements[2].strip() # compare A
		#A = elements[3].strip() # compare B
		seed_dict[position] = seed(rsID, position, A)
inputFile_seed.close()

print len(seed_dict)
for line in inputFile_ref:
	if line.startswith("rsID"):
		elements = line.strip().split()
		i = 0
		while i < len(elements):
			name = elements[i].strip()
			if name.startswith("NA") and name not in ref_dict:			
				index = i
				same2seed = 0
				similarity = 0
				snp_dict = {}			
				ref_dict[i] = ref(name, index, same2seed, similarity, snp_dict)
			else:
				print "duplicate NA name: ", name
			i += 1
	if not line.startswith("rsID"):
		elements = line.strip().split()
		try:
			position = elements[1].strip()
			if position in seed_dict:
				seed_allele = seed_dict[position].A
				for index, ref in ref_dict.iteritems():
					ref_allele = elements[index].strip()
				#	print ref_allele
					ref.snp_dict[position] = ref_allele
					if seed_allele == ref_allele:
				#		print seed_allele
						ref.same2seed += 1
		except ValueError:
			print "error in", line		

inputFile_ref.close()

#max_similarity = 0

print >> data_record_file, "seed file: ", seed_file
print >> data_record_file, "ref file: ", ref_file
print >> data_record_file, "total seed number: ", total_seed_number
print >> data_record_file, "total NA number: ", len(ref_dict)

for index, ref in ref_dict.iteritems():
	ref.similarity = float(ref.same2seed)/float(total_seed_number)
	print >> data_record_file, ref.name, ref.same2seed, ref.similarity


print "total seed number: ", total_seed_number
print "total NA number: ", len(ref_dict)
data_record_file.close()


