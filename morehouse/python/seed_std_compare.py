#!/usr/bin/python
#######################################################################################
# Compare seed and std hap, to check purity of seed
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import file_path, data_record_path, currentPath
from tools import sort_dict_by_key, load_raw_data

def seed_std_compare(seed_input_file, chr_name):
	std_file_name = "ASW_"+chr_name+"_child_hap_refed.txt"
	std_input_file = file_path + std_file_name
	snp_hap_std_dict = load_raw_data(std_input_file)[1]
	
	snp_hap_seed_dict = load_raw_data(seed_input_file)[1]
	seed_file_name = seed_input_file[:seed_input_file.find('.')].strip()
	print "snp_hap_std_total_number", len(snp_hap_std_dict)
	print "snp_hap_seed_total_number", len(snp_hap_seed_dict)
	
	seed_same_to_A = 0
	seed_same_to_B = 0
	seed_same_to_AB = 0
	seed_X = 0
	seed_N = 0
	seed_not_AB = 0
	
	same_to_B_dict = {}

	compare_output_file_name = seed_file_name + "_" + chr_name + "_compare.txt"
	compare_output_file = open(currentPath + compare_output_file_name, "w")
	print >> compare_output_file, "seed file: ", seed_input_file
	print >> compare_output_file, "std hap file: ", std_input_file

	snp_hap_seed_dict_sorted_list = sort_dict_by_key(snp_hap_seed_dict)

	for snp in snp_hap_seed_dict_sorted_list:
		position = snp[0]	
		# check hifi seeds, these position need to be in std too
		if position in snp_hap_std_dict:
			elements_std = snp_hap_std_dict[position]
			std_A = elements_std[2].strip()
			std_B = elements_std[3].strip()
			elements_seed = snp_hap_seed_dict[position]
			seed_A = elements_seed[2].strip()
			#seed_A = elements_seed[3].strip() # for solid data 4 and 6, chr from mother
			if seed_A == std_A:
				if seed_A == std_B:
					seed_same_to_AB += 1
				else:
					seed_same_to_A += 1			
			elif seed_A == std_B:
				seed_same_to_B += 1
				same_to_B_dict[position] = elements_seed
			elif std_A == "X" or std_B == "X":
				seed_X += 1
			elif std_A == "N" or std_B == "N":
				seed_N += 1
			else:
				seed_not_AB += 1
	A_in_hetero = format((float(seed_same_to_A)/float(seed_same_to_A + seed_same_to_B))*100, "0.2f")
	B_in_hetero = format((float(seed_same_to_B)/float(seed_same_to_A + seed_same_to_B))*100, "0.2f")
	
	print "seed_same_to_A", seed_same_to_A
	print "seed_same_to_B", seed_same_to_B
	print "A_in_hetero", A_in_hetero, "%"
	print "B_in_hetero", B_in_hetero, "%"
	print "homo seed", seed_same_to_AB
	print "seed_not_AB", seed_not_AB

	print >> compare_output_file, "seed_same_to_A", seed_same_to_A
	print >> compare_output_file, "seed_same_to_B", seed_same_to_B
	print >> compare_output_file, "A_in_hetero", A_in_hetero, "%"
	print >> compare_output_file, "B_in_hetero", B_in_hetero, "%"
	print >> compare_output_file, "homo seed", seed_same_to_AB
	print >> compare_output_file, "seed_not_AB", seed_not_AB	
	compare_output_file.close()
	return same_to_B_dict

"""
def usage():
	print "%s [seed_file] [chr]" % sys.argv[0]
if __name__=='__main__':
	if len(sys.argv)!=3:
		usage()
		sys.exit(1)
	seed_input_file=sys.argv[1]
	chr_name=sys.argv[2]
	seed_std_compare(seed_input_file, chr_name)
"""

def get_args():
	desc="Compare seed and std hap, to check purity of seed"

	usage = "seed_std_compare -i seed_file -c chr#" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
	parser.add_option("-i", "--seed", type="string", dest="hifiSeed",help = "Input seed file Name", default="null")
	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.hifiSeed == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__=='__main__':
	options = get_args()
	chr_name = options.chrName
	seed_input_file = options.hifiSeed	
	seed_std_compare(seed_input_file, chr_name)
