#!/usr/bin/python

# location /home/guoxing/tool/morehouse

import os, glob, subprocess, random, operator, time
from optparse import OptionParser
from tools import *
#from tools import file_path, program_path, data_record_path, currentPath
#from tools import sort_dict_by_key, load_raw_data, wccount, keywithmaxval


data_record_path = "/home/guoxing/disk2/solid/common_files/data_record/"

snp_hap_ori_dict = {}
hifi_result_dict = {}
hap_std_total_number = 0
hap_std_total_number_withoutXN = 0
hifi_result_total_number = 0

def get_args():
	desc="Compare seed and std hap, to check purity of seed"
	usage = "seed_std_compare -i seed_file -c chr#" 
	parser = OptionParser(usage = usage) 
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="chr11")
	parser.add_option("-i", "--imputed", type="string", dest="hifiResult",help = "Input hifiResult file Name", default="null")
	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.hifiResult == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

def allele_similarity(hifi_result_dict, hifi_std_dict):
	same_to_A = 0
	same_to_B = 0
	for position, elements_hifi in hifi_result_dict.iteritems():
		if position in hifi_std_dict:
			hifi_A = elements_hifi[2].strip()
			hifi_B = elements_hifi[3].strip()
			elements_std = hifi_std_dict[position]
			std_A = elements_std[2].strip()
			std_B = elements_std[3].strip()
			if hifi_A != hifi_B:
				if hifi_A == std_A:
					same_to_A += 1
				if hifi_A == std_B:
					same_to_B += 1
	#print "similarity, same_to_A", same_to_A
	#print "similarity, same_to_B", same_to_B
	if same_to_A >= same_to_B:
		return "similar_to_A"
	else:
		return "similar_to_B" 

def compare_std_result(hifi_result_dict, hifi_std_dict):
	same_to_A_dict= {}
	same_to_B_dict = {}
	same_to_X_dict = {}
	same_to_N_dict = {}
	same_to_AB_dict = {}
	not_same_to_AB_dict = {}
	different_position_dict = {}
	same_position_total_number = 0
	different_position_total_number = 0
	similarity = allele_similarity(hifi_result_dict, hifi_std_dict)
	
	for position, elements_hifi in hifi_result_dict.iteritems():
		if position in hifi_std_dict:
			hifi_A = elements_hifi[2].strip()
			hifi_B = elements_hifi[3].strip()
			elements_std = hifi_std_dict[position]
			std_A = elements_std[2].strip()
			std_B = elements_std[3].strip()
			if similarity == "similar_to_B":		# for solid data 4 and 6, the hifi seed is from mother (B)
				hifi_A, hifi_B = hifi_B, hifi_A
				# the hifi seed is from father, A				
			if hifi_A == std_A:	#A is A
				if hifi_B == std_B:
					same_to_AB_dict[position] = elements_hifi
				else:
					same_to_A_dict[position] = elements_hifi
			elif hifi_B == std_B:
				same_to_B_dict[position] = elements_hifi
			elif std_A == "X" or std_B == "X" or std_A == "N" or std_B == "N":
				same_to_AB_dict[position] = elements_hifi
			else:
				not_same_to_AB_dict[position] = elements_hifi
			same_position_total_number += 1			
		else:
			different_position_total_number += 1
			different_position_dict[position] = elements_hifi
			
	return (same_to_A_dict, same_to_B_dict, same_to_AB_dict, not_same_to_AB_dict, different_position_dict, same_position_total_number, different_position_total_number)


def output_different_snps(different_position_dict):
	difference_output_file_name = "hifi_difference_position.txt"
	difference_output_file = open(currentPath + difference_output_file_name, "w")
	for position, elements in different_position_dict.iteritems():
		"""
		line = ""
		for a in elements:
			line += a + "\t"
		print >> difference_output_file, line
		"""
		print >> difference_output_file, list_to_line(elements)
	difference_output_file.close()


def hifiAccuCheck (hifi_result_file, chr_name):
	#hap_std_file_name = "NA12878_hap_new_refed.txt"	# simulation data
	hap_std_file_name = file_path + "ASW_"+chr_name+"_child_hap_refed.txt"	
	hifi_std_dict = load_raw_data(hap_std_file_name, raw_data_format)[1]
	
	hifi_std_dict = removeN(hifi_std_dict)
	hap_std_total_number = len(hifi_std_dict)
	
	hifi_result_dict = load_raw_data(hifi_result_file, raw_data_format)[1]
	hifi_result_total_number = len(hifi_result_dict)
	
	print "hap_std_total_number", hap_std_total_number
	print "hifi_result_total_number", hifi_result_total_number
	
	compare_tuple = compare_std_result(hifi_result_dict, hifi_std_dict)
	same_to_A_dict= compare_tuple[0]
	#print len(same_to_A_dict)
	same_to_B_dict = compare_tuple[1]
	same_to_AB_dict = compare_tuple[2]
	not_same_to_AB_dict = compare_tuple[3]
	different_position_dict = compare_tuple[4]
	same_position_total_number = compare_tuple[5]
	different_position_total_number = compare_tuple[6]
	
	output_different_snps(different_position_dict)
	
	same_A_total_number = len(same_to_A_dict)
	same_B_total_number = len(same_to_B_dict)
	same_AB_total_number = len(same_to_AB_dict)
	not_same_AB_total_number = len(not_same_to_AB_dict)
	
	pencentage_in_common = format(float(same_position_total_number)/float(hifi_result_total_number)*100, "0.3f")
	accuracy = format(float(same_A_total_number + same_B_total_number + same_AB_total_number)/float(same_position_total_number)*100, "0.3f")	
	
	print "same_position_total_number", same_position_total_number
	print "different_position_total_number", different_position_total_number
	print "same_AB_total_number", same_AB_total_number
	print "same_A_total_number", same_A_total_number
	print "same_B_total_number", same_B_total_number
	print "not_same_AB_total_number", not_same_AB_total_number
	print "pencentage in common", pencentage_in_common
	print "accuracy", accuracy
	
	accuracy_output_file_name = "hifi_accuracy.txt"
	accuracy_output_file = open(currentPath + accuracy_output_file_name, "w")	
	print >> accuracy_output_file, "hifi result file: ", currentPath+hifi_result_file
	print >> accuracy_output_file, "hifi std file: ", hap_std_file_name
	print >> accuracy_output_file, "hap_std_total_number: ", hap_std_total_number
	print >> accuracy_output_file, "hifi_result_total_number: ", hifi_result_total_number
	print >> accuracy_output_file, "same_position_total_number: ", same_position_total_number
	print >> accuracy_output_file, "different_position_total_number: ", different_position_total_number
	print >> accuracy_output_file, "same_A_total_number: ", same_A_total_number
	print >> accuracy_output_file, "same_B_total_number: ", same_B_total_number
	print >> accuracy_output_file, "same_AB_total_number: ", same_AB_total_number
	print >> accuracy_output_file, "not_same_AB_total_number: ", not_same_AB_total_number
	print >> accuracy_output_file, "pencentage in common: ", pencentage_in_common
	print >> accuracy_output_file, "accuracy: ", accuracy
	
	accuracy_output_file.close()
	
	# record data
	data_record_file_name = "solid_process_4.txt"
	data_record_file = open(data_record_path + data_record_file_name, "a")
	print >> data_record_file, currentPath, same_position_total_number, (same_A_total_number+same_B_total_number+not_same_AB_total_number), same_AB_total_number, hifi_result_total_number, accuracy
	data_record_file.close()
	cmd = "grep hifi_data hifi_accuracy.txt >> " + data_record_path + data_record_file_name
	#print cmd
	#os.system(cmd)

if __name__=='__main__':
	options = get_args()
	chr_name = options.chrName
	hifi_result_file = options.hifiResult
	hifiAccuCheck(hifi_result_file, chr_name)
