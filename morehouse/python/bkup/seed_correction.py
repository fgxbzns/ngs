#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# To correct and extend the seed

import os, glob, subprocess, random, operator, time, math
from optparse import OptionParser

program_path = "/home/guoxing/disk2/ngs/morehouse/python/"
file_path = "/home/guoxing/disk2/solid/common_files/"

currentPath = os.getcwd() + '/'


def wccount(filename):
    out = subprocess.Popen(['wc', '-l', filename],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

class seed:
	def __init__(me, rsID, position, allele_ori, allele_new, allele_dict):
		me.rsID = rsID
		me.position = position
		me.allele_ori = allele_ori
		me.allele_new = allele_new
		me.allele_new_percentage = 0
		me.allele_dict = allele_dict

def keywithmaxval(dict):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value """  
     v=list(dict.values())
     k=list(dict.keys())
     return k[v.index(max(v))]

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
			allele_new = ""
			#allele_new_percentage = 0
			allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
			snp_dict[position] = seed(rsID, position, allele_ori, allele_new, allele_dict)
	input_file.close()
	return snp_dict

def load_hifi_result(file_name, hifi_dict):
	#hifi_dict = {}
	#for position, seed in seed_dict.iteritems():
	#	hifi_dict[position] = seed
	total_hifi_snp_number = 0
	input_file = open(currentPath + file_name, "r")	
	for line in input_file:
		if line.startswith("rsID"):
			title_haplotype = line.strip()
		else:
			total_hifi_snp_number += 1
			elements = line.strip().split()
			rsID = elements[0].strip()
			try:
				position = int(elements[1].strip())
				if position not in hifi_dict:
					allele_ori = ""	
					allele_new = ""
					allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
					hifi_dict[position] = seed(rsID, position, allele_ori, allele_new, allele_dict)							
				allele_new = elements[2].strip()	# need to modify this if seed is from second snp column, mother side
				hifi_dict[position].allele_new = allele_new
				if allele_new == 'A':
					hifi_dict[position].allele_dict['A'] += 1
				elif allele_new == 'T':
					hifi_dict[position].allele_dict['T'] += 1
				elif allele_new == 'C':
					hifi_dict[position].allele_dict['C'] += 1
				elif allele_new == 'G':
					hifi_dict[position].allele_dict['G'] += 1
			except:
				#print line.strip()
				pass
	input_file.close()
	print "total hifi snp number in : ", file_name, total_hifi_snp_number 
	return hifi_dict


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
		sub_seed_list = seed_sorted_list[:seed_in_each_subfile*file_number] \
						+ seed_sorted_list[seed_in_each_subfile*(file_number+1):(len(seed_sorted_list)-1)]
		sub_seed_list.append(seed_sorted_list[-1])				
		print len(sub_seed_list)
		for seed in sub_seed_list:
			line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
			print >> output_subfile, line	
		file_number += 1
		output_subfile.close()
		"""
		cp = "cp " + output_subfile_name + " haplotype.txt"
		cp_process = subprocess.Popen(cp, shell=True)
		cp_process.wait()
		
		hifi = program_path + "hifi &"
		hifi_process = subprocess.Popen(hifi, shell=True)
		hifi_process.wait()
		
		cp = "cp imputedhaplotype.txt imputedhaplotype_" + str(file_number) + ".txt"
		cp_process = subprocess.Popen(cp, shell=True)
		cp_process.wait()
		"""

number_of_subfile = 5
seed_spliter(seed_dict, number_of_subfile)
"""
hifi_dict = {}
for position, seeds in seed_dict.iteritems():
	hifi_dict[position] = seeds
print len(hifi_dict)
	
file_number = 0
while file_number < number_of_subfile:	
	input_subfile_name = "imputedhaplotype_0" + str(file_number) + ".txt"
	hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)
	print file_number," :", len(hifi_dict)
	file_number += 1

seed_new_file = open(currentPath + "haplotype_new.txt", "w")
print >> seed_new_file, title_haplotype
num_100 = 0

hifi_sorted_list = [x for x in hifi_dict.iteritems()] 
hifi_sorted_list.sort(key=lambda x: x[0]) # sort by key
for snp in hifi_sorted_list:
	seeds = snp[1]
	max_base = keywithmaxval(seeds.allele_dict)
	max_value = seeds.allele_dict[max_base]
	seeds.allele_new_percentage = float(max_value)/float(number_of_subfile)
	if int(seeds.allele_new_percentage*100) == 100:
		num_100 += 1
		seeds.allele_new = max_base
		line = seeds.rsID + "\t" + str(seeds.position) + "\t" + max_base
		print >> seed_new_file, line

seed_new_file.close()
	
print "num_100", num_100

"""



