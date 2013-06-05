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

#hap_ref_file_name = "NA12878_hap_new_refed.txt"	# simulation data
hap_ref_file_name = "ASW_"+chr_name+"_child_hap_refed.txt"	

hap_hifi_file_name = "imputedhaplotype.txt"

inputFile_hap_ref = open(file_path + hap_ref_file_name, "r")
inputFile_hap_hifi = open(currentPath + hap_hifi_file_name, "r")
inputFile_hap_seed = open(currentPath + seed_file, "r")

x_number = 0
n_number = 0
a_number = 0
t_number = 0
c_number = 0
g_number = 0

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
			"""
		if elements[2].strip() == "N":
			n_number += 1
		if elements[2].strip() == "X":
			x_number += 1
		if elements[2].strip() == "A":
			a_number += 1
		if elements[2].strip() == "T":
			t_number += 1
		if elements[2].strip() == "C":
			c_number += 1
		if elements[2].strip() == "G":
			g_number += 1
	
print "n_number", n_number
print "x_number", x_number	
print "a_number", a_number
print "t_number", t_number
print "c_number", c_number
print "g_number", g_number
"""			
for line in inputFile_hap_hifi:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		if True:
			position = elements[1].strip()	
			try:
				position = int(position)
				snp_hap_hifi_dict[position] = line.strip()
			except ValueError:
				print position	
				print line
			snp_hap_hifi_total_number += 1	

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
print "snp_hap_hifi_total_number", snp_hap_hifi_total_number
print "snp_hap_seed_total_number", snp_hap_seed_total_number

same_position_total_number = 0
different_position_total_number = 0
same_AB_total_number = 0
same_A_total_number = 0
same_B_total_number = 0
not_same_AB_total_number = 0

error_distribution_output_file_name = "hifi_error_distribution.txt"
error_distribution_output_file = open(currentPath + error_distribution_output_file_name, "w")
#error_distribution_output_file.write("position")
#print >> error_distribution_output_file, "position seed_AB seed_A seed_B seed_X seed_N seed_other hifi_AB hifi_A hifi_B hifi_X hifi_N hifi_other"
print >> error_distribution_output_file, "position \t seed_AB \t seed_A \t seed_B \t seed_X \t seed_N \t seed_other \t hifi_AB \t hifi_A \t hifi_B \t hifi_X \t hifi_N \t hifi_other"
#print >> error_distribution_output_file, position, seed_AB_pos, seed_A_pos, seed_B_pos, seed_X_pos, seed_N_pos, seed_other_pos, hifi_AB_pos, hifi_A_pos, hifi_B_pos, hifi_X_pos, hifi_N_pos, hifi_other_pos

y_axis_value = "80 \t"
seed_correct_axis_value = "100 \t"
seed_error_axis_value = "120 \t"
hifi_error_axis_value = "140 \t"

snp_hap_hifi_dict_sorted_list = [x for x in snp_hap_hifi_dict.iteritems()] 
snp_hap_hifi_dict_sorted_list.sort(key=lambda x: x[0]) # sort by key

#for position, line_hifi in snp_hap_hifi_dict.iteritems():
for snp in snp_hap_hifi_dict_sorted_list:
	
	position = snp[0]
	seed_AB_pos = "\t"	
	seed_A_pos = "\t"
	seed_B_pos = "\t"
	seed_X_pos = "\t"
	seed_N_pos = "\t"
	seed_other_pos = "\t"

	hifi_AB_pos = "\t"
	hifi_A_pos = "\t"
	hifi_B_pos = "\t"
	hifi_X_pos = "\t"
	hifi_N_pos = "\t"
	hifi_other_pos = "\t"
	"""
	seed_A_pos = ""
	seed_B_pos = "NA"
	seed_X_pos = "NA"
	seed_N_pos = "NA"
	seed_other_pos = "NA"

	hifi_AB_pos = "NA"
	hifi_A_pos = "NA"
	hifi_B_pos = "NA"
	hifi_X_pos = "NA"
	hifi_N_pos = "NA"
	hifi_other_pos = "NA"
	"""
	# check hifi results
	if position in snp_hap_ref_dict:
		line_ref = snp_hap_ref_dict[position]
		elements_ref = line_ref.strip().split()
		ref_A = elements_ref[2].strip()
		ref_B = elements_ref[3].strip()
		line_hifi = snp_hap_hifi_dict[position]
		elements_hifi = line_hifi.strip().split()
		hifi_A = elements_hifi[2].strip()
		hifi_B = elements_hifi[3].strip()
	
		if hifi_A == ref_A:
			if hifi_B == ref_B:
				hifi_AB_pos = y_axis_value
			else:
				hifi_A_pos = y_axis_value
				same_A_total_number += 1
		elif hifi_B == ref_B:	# error = B + other	assume A is the selected chr. need to update this for solid data song_4 song_6
			hifi_B_pos = y_axis_value
			#hifi_other_pos = y_axis_value
		elif ref_A == "X" or ref_B == "X":
			hifi_X_pos = y_axis_value
		elif ref_A == "N" or ref_B == "N":
			hifi_N_pos = y_axis_value
		else:
			hifi_other_pos = hifi_error_axis_value
			
		# check hifi seeds, these position need to be in ref too
		if position in snp_hap_seed_dict:
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
					seed_AB_pos = seed_correct_axis_value
				else:
					seed_A_pos = seed_correct_axis_value
			elif seed_A == ref_B:
				seed_B_pos = seed_error_axis_value
			elif ref_A == "X" or ref_B == "X":
				seed_X_pos = y_axis_value
			elif ref_A == "N" or ref_B == "N":
				seed_N_pos = y_axis_value
			else:
				seed_other_pos = y_axis_value		
		if hifi_A != hifi_B:
			print >> error_distribution_output_file, str(position)+"\t", seed_AB_pos, seed_A_pos, seed_B_pos, seed_X_pos, seed_N_pos, seed_other_pos, hifi_AB_pos, hifi_A_pos, hifi_B_pos, hifi_X_pos, hifi_N_pos, hifi_other_pos

inputFile_hap_ref.close()
inputFile_hap_hifi.close()
inputFile_hap_seed.close()
error_distribution_output_file.close()
