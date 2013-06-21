#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# To correct and extend the seed

import os, glob, subprocess, random, operator, time, math
from optparse import OptionParser

program_path = "/home/guoxing/tool/morehouse"
file_path = "/home/guoxing/disk2/solid/common_files/"

currentPath = os.getcwd() + '/'


def wccount(filename):
    out = subprocess.Popen(['wc', '-l', filename],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

class seed:
	def __init__(me, rsID, position, allele_ori):
		me.rsID = rsID
		me.position = position
		me.allele_ori = allele_ori
		me.allele_new = ""
		me.allele_new_percentage = 0
		me.allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}

def print_f(var):
	print '\"' + var

def load_seed_data(file_name):
	snp_dict = {}
	global title_haplotype	# save the title line of hap file for latter use
	global total_seed_number
	total_seed_number = 0
	input_file = open(currentPath + file_name, "r")	
	for line in input_file:
		if line.startswith("rsID"):
			title_haplotype = line.strip()
		else:
			total_seed_number += 1
			elements = line.strip().split()
			rsID = elements[0].strip()
			position = int(elements[1].strip())
			allele_ori = elements[2].strip()
			#allele_new = ""
			#allele_new_percentage = 0
			#allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
			snp_dict[position] = seed(rsID, position, allele_ori)
	input_file.close()
	return snp_dict

# start time
start = time.time()		
		
# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--seed", type="string", dest="seedFile",help = "Input seed File Name", default="null")


parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
parser.add_option("-d", "--threshold", type="string", dest="threshold",help = "Input the depth threshold", default="3")
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")

(options, args) = parser.parse_args()

seed_file = options.seedFile
sam_file = options.samFile
depth_threshold = int(options.threshold)
chr_name = options.chrName

print "seed_file is :", seed_file
#total_seed_number = int(wccount(seed_file)) - 1

seed_file_name = seed_file[:seed_file.find('.')].strip()
seed_dict = load_seed_data(seed_file)
print len(seed_dict)

print "title_haplotype", title_haplotype
print "total_seed_number: ", total_seed_number

def seed_spliter(seed_dict, number_of_subfile):
	global seed_sorted_list
	seed_sorted_list = [x for x in seed_dict.iteritems()] 
	seed_sorted_list.sort(key=lambda x: x[0]) # sort by key
	seed_number_ceilling = int(math.ceil(float(total_seed_number)/100)*100)
	print "seed_number_ceilling: ", seed_number_ceilling
	seed_in_each_subfile = seed_number_ceilling/number_of_subfile
	print "seed_in_each_subfile: ", seed_in_each_subfile
	seed_in_last_subfile = int(math.fmod(total_seed_number, seed_in_each_subfile))
	print "seed_in_last_subfile: ", seed_in_last_subfile
	
	file_number = 0
	while file_number < number_of_subfile:
		line_number = 0
		output_subfile_name = seed_file_name + "_" + str(file_number) + ".txt"
		print "output: ", output_subfile_name
		output_subfile = open(currentPath + output_subfile_name, "w")
		print >> output_subfile, title_haplotype
		i = 0
		sub_seed_list = seed_sorted_list[:seed_in_each_subfile*file_number] + seed_sorted_list[seed_in_each_subfile*(file_number+1):]
		print len(sub_seed_list)
		for seed in sub_seed_list:
			line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
			print >> output_subfile, line	
		file_number += 1
		output_subfile.close()

	
	
	

number_of_subfile = 3
seed_spliter(seed_dict, number_of_subfile)
