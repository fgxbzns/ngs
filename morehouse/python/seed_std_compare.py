#!/usr/bin/python
#######################################################################################
# Compare seed and std hap, to check purity of seed
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *
from calculate_maf import calculate_maf

def output_seed(file_name, seed_title_info, snp_hap_seed_dict, output_seed_dict):
	seed_new_file = open(currentPath + file_name, "w")
	print >> seed_new_file, seed_title_info
	output_pos = output_seed_dict.keys()
	output_pos.sort()
	for pos in output_pos:
		print >> seed_new_file, list_to_line(snp_hap_seed_dict[pos])
	seed_new_file.close()

def seed_std_compare(seed_input_file, chr_name):
	std_file_name = "ASW_"+chr_name+"_child_hap_refed.txt"	# for solid data
	#std_file_name = "NA12878_hap_new_refed.txt"	# simulation data hg18 chr6

	std_input_file = file_path + std_file_name
	snp_hap_std_dict = load_raw_data(std_input_file, raw_data_format)[1]
	
	seed_title_info, snp_hap_seed_dict = load_raw_data(seed_input_file, raw_data_format)
	seed_file_name = seed_input_file[:seed_input_file.find('.')].strip()
	print "snp_hap_std_total_number", len(snp_hap_std_dict)
	print "snp_hap_seed_total_number", len(snp_hap_seed_dict)
	
	seed_same_to_A = 0
	seed_same_to_B = 0
	seed_same_to_AB = 0
	seed_X = 0
	seed_N = 0
	seed_not_AB = 0
	
	same_to_A_dict = {}
	same_to_B_dict = {}

	compare_output_file_name = seed_file_name + "_" + chr_name + "_compare.txt"
	compare_output_file = open(compare_output_file_name, "w")
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
			if seed_A == std_A :
				if seed_A == std_B:
					seed_same_to_AB += 1
				elif std_B != "N":
					seed_same_to_A += 1	
					same_to_A_dict[position] = elements_seed		
			elif seed_A == std_B and std_A != "N" and std_A != "X":
				seed_same_to_B += 1
				same_to_B_dict[position] = elements_seed
			elif std_A == "X" or std_B == "X":
				seed_X += 1
			elif std_A == "N" or std_B == "N":
				seed_N += 1
			else:
				seed_not_AB += 1
	
	
	A_in_hetero = 0
	B_in_hetero = 0
	A_in_hetero = format((float(seed_same_to_A)/float(seed_same_to_A + seed_same_to_B))*100, "0.2f")
	B_in_hetero = format((float(seed_same_to_B)/float(seed_same_to_A + seed_same_to_B))*100, "0.2f")
	
	print "total hetero", seed_same_to_A + seed_same_to_B
	print "seed_same_to_A", seed_same_to_A
	print "seed_same_to_B", seed_same_to_B
	print "A_in_hetero", A_in_hetero, "%"
	print "B_in_hetero", B_in_hetero, "%"
	print "homo seed", seed_same_to_AB
	print "seed_not_AB", seed_not_AB
	print "seed_X", seed_X
	print "seed_N", seed_N

	print >> compare_output_file, "seed_same_to_A", seed_same_to_A
	print >> compare_output_file, "seed_same_to_B", seed_same_to_B
	print >> compare_output_file, "A_in_hetero", A_in_hetero, "%"
	print >> compare_output_file, "B_in_hetero", B_in_hetero, "%"
	print >> compare_output_file, "homo seed", seed_same_to_AB
	print >> compare_output_file, "seed_not_AB", seed_not_AB	
	compare_output_file.close()
	
	#print same_to_B_dict
	hetero_seed_dict = dict_add(same_to_A_dict, same_to_B_dict)
	hetero_seed_sorted_list = sort_dict_by_key(hetero_seed_dict)
	
	distribution = 0
	hetero_seed_sorted_list.reverse()
	for i in range(len(hetero_seed_sorted_list)-1):
		distribution = hetero_seed_sorted_list[i][0] - hetero_seed_sorted_list[i+1][0]
	print "hetero distribution", distribution
	print "hetero distribution mean", round(float(distribution)/len(hetero_seed_sorted_list), 3)
	print "hetero distribution mean over chr", round(float(distribution)/len(hetero_seed_sorted_list)/chr_length_dict[chr_name], 10)
		
	distribution = 0
	snp_hap_seed_dict_sorted_list.reverse()
	for i in range(len(snp_hap_seed_dict_sorted_list)-1):
		distribution = snp_hap_seed_dict_sorted_list[i][0] - snp_hap_seed_dict_sorted_list[i+1][0] 
	
	print "distribution", distribution
	print "distribution mean", round(float(distribution)/len(snp_hap_seed_dict_sorted_list), 3)
	print "distribution mean over chr", round(float(distribution)/len(snp_hap_seed_dict_sorted_list)/float(chr_length_dict[chr_name]), 10)
	#file_name = "refHaplos.txt"
	#for pos, snp in same_to_B_dict.iteritems():
	#	pass
		#print snp_hap_std_dict[pos]
		#calculate_maf(file_name, pos)
	#output_seed("haplotype_A.txt", seed_title_info, snp_hap_seed_dict, same_to_A_dict)
	return (same_to_A_dict, same_to_B_dict)
	
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

def seperate_homo_hetero(same_to_AB_dict):
    homo_dict = {}
    hetero_dict = {}
    for position, snp in same_to_AB_dict.iteritems():
    	#if snp[2] != 'X' and snp[3] != 'X' and snp[2] != 'N' and snp[3] != 'N':
    	if snp[2] == snp[3]:
    		homo_dict[position] = snp
    	else:
    		hetero_dict[position] = snp
    return (homo_dict, hetero_dict)

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
