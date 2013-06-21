#!/usr/bin/python

# location /home/guoxing/tool/morehouse

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

file_path = "/home/guoxing/disk2/solid/common_files/"
data_record_path = "/home/guoxing/disk2/solid/common_files/data_record/"
currentPath = os.getcwd() + '/'

snp_hap_ref_dict = {}
snp_hap_hifi_dict = {}
snp_hap_seed_dict = {}
snp_hap_ref_total_number = 0
snp_hap_ref_total_number_withoutXN = 0
snp_hap_hifi_total_number = 0
snp_hap_seed_total_number = 0

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
parser.add_option("-i", "--seed", type="string", dest="hifiSeed",help = "Input seed file Name", default="null")

(options, args) = parser.parse_args()
chr_name = options.chrName
seed_file = options.hifiSeed

hap_ref_file_name = "ASW_"+chr_name+"_child_hap_refed.txt"	

inputFile_hap_ref = open(file_path + hap_ref_file_name, "r")
inputFile_hap_seed = open(currentPath + seed_file, "r")
seed_file_name = seed_file[:seed_file.find('.')].strip()

for line in inputFile_hap_ref:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		#if elements[2].strip() != "N" and elements[2].strip() != "X" and elements[3].strip() != "N" and elements[3].strip() != "X":
		if elements[2].strip() != "N" and elements[3].strip() != "N":
			position = elements[1].strip()	
			try:
				position = int(position)
				snp_hap_ref_dict[position] = line.strip()
			except ValueError:
				print file_name, position	
				print line
			snp_hap_ref_total_number += 1
		
for line in inputFile_hap_seed:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		if True:
			position = elements[1].strip()	
			try:
				position = int(position)
				snp_hap_seed_dict[position] = line.strip()
			except ValueError:
				print position	
				print line
			snp_hap_seed_total_number += 1	

print "snp_hap_ref_total_number", snp_hap_ref_total_number
print "snp_hap_seed_total_number", snp_hap_seed_total_number

seed_same_to_A = 0
seed_same_to_B = 0
seed_same_to_AB = 0
seed_X = 0
seed_N = 0
seed_not_AB = 0

error_distribution_output_file_name = seed_file_name + "_std_compare.txt"
error_distribution_output_file = open(currentPath + error_distribution_output_file_name, "w")
print >> error_distribution_output_file, "seed file: ", seed_file
print >> error_distribution_output_file, "std hap file: ", hap_ref_file_name

snp_hap_seed_dict_sorted_list = [x for x in snp_hap_seed_dict.iteritems()] 
snp_hap_seed_dict_sorted_list.sort(key=lambda x: x[0]) # sort by key

for snp in snp_hap_seed_dict_sorted_list:
	
	position = snp[0]
		
	# check hifi seeds, these position need to be in ref too
	if position in snp_hap_ref_dict:
		line_ref = snp_hap_ref_dict[position]
		elements_ref = line_ref.strip().split()
		ref_A = elements_ref[2].strip()
		ref_B = elements_ref[3].strip()
		line_seed = snp_hap_seed_dict[position]
		elements_seed = line_seed.strip().split()
		seed_A = elements_seed[2].strip()
		#seed_A = elements_seed[3].strip() # for solid data 4 and 6, chr from mother
		if seed_A == ref_A:
			if seed_A == ref_B:
				seed_same_to_AB += 1
			else:
				seed_same_to_A += 1
			
		elif seed_A == ref_B:
			seed_same_to_B += 1
		elif ref_A == "X" or ref_B == "X":
			seed_X += 1
		elif ref_A == "N" or ref_B == "N":
			seed_N += 1
		else:
			seed_not_AB += 1

print "seed_same_to_A", seed_same_to_A
print "seed_same_to_B", seed_same_to_B
print "homo seed", seed_same_to_AB
print "seed_not_AB", seed_not_AB

print >> error_distribution_output_file, "seed_same_to_A", seed_same_to_A
print >> error_distribution_output_file, "seed_same_to_B", seed_same_to_B
print >> error_distribution_output_file, "homo seed", seed_same_to_AB
print >> error_distribution_output_file, "seed_not_AB", seed_not_AB

inputFile_hap_ref.close()
inputFile_hap_seed.close()
error_distribution_output_file.close()

