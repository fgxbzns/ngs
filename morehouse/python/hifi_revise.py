#!/usr/bin/python
#######################################################################################
# hifi v2
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *
from hifiAccuCheck_v2 import hifiAccuCheck
from seed_std_compare import seed_std_compare

def hifi_test(seed_input_file):
	#hifi = program_path + "hifi_fu_revise " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " " + str(maf_step) + ""	# make sure the other hifi processes are finished
	hifi = program_path + "hifi_fu_revise " + seed_input_file
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()

def get_window_info():
	data = {}
	file_name = "pos.out"
	
	fp = open(file_name, "r")
	for line in fp:
		if line.startswith("Imputing"):
			elements = line.strip().split()
			try:			
				data[int(elements[2])] = elements
			except ValueError:
				#print "error in ", line
				pass
	fp.close()
	return data


def load_hap_std(file_name):
    hap_std_dict = {}
    temp_data_dict = load_raw_data(file_name, raw_data_format)[1]
    for position, elements in temp_data_dict.iteritems():
        """ ??? N X """
        if elements[2].strip() != "N" and elements[3].strip() != "N": 
            try:
                position = elements[1].strip()
                position = int(position)
                hap_std_dict[position] = elements
            except ValueError:
                print "error in ", file_name, position
    return hap_std_dict  

def analyze_data():
	
	hap_std_file = file_path + "ASW_"+chr_name+"_child_hap_refed.txt" 
	hap_std_dict = load_hap_std(hap_std_file)
	
	window_info_dict = get_window_info()
	window_info_sorted_list = sort_dict_by_key(window_info_dict)
	print window_info_sorted_list[0][0], window_info_sorted_list[0][1]
	window_size_total = 0
	
	correct_snp = 0
	error_snp = 0
	window_size_total_correct = 0
	window_size_total_error = 0
	
	maf_total_c = 0.
	maf_total = 0.
	for data in window_info_sorted_list:
		pos = data[0]
		info = data[1]
		window_size_total += int(info[7])
		if pos in hap_std_dict and hap_std_dict[pos][2] != 'X':
			if info[5][0] == hap_std_dict[pos][2] and info[5][1] == hap_std_dict[pos][3]:
				correct_snp += 1
				window_size_total_correct += int(info[7])
				maf_total_c += float(info[9])
			#elif float(info[9]) < 0.2:
			else:
				error_snp += 1
				window_size_total_error += int(info[7])
				print info[5], hap_std_dict[pos][2], hap_std_dict[pos][3], info[7], info[9]
				maf_total += float(info[9])
			
	
	print "window_size_total", window_size_total
	print "window_size_average", window_size_total/len(window_info_dict)
	print "correct_snp", correct_snp
	print "correct_snp", window_size_total_correct/correct_snp
	print "error_snp", error_snp
	print "error_snp", window_size_total_error/error_snp
	
	print "maf_total_c", maf_total_c
	print "maf_total_c", maf_total_c/correct_snp	
	
	print "maf_total", maf_total
	print "maf_total", maf_total/error_snp	
	
def get_args():
	desc = "Compare seed and std hap, to check purity of seed"

	usage = "seed_std_compare -i seed_file -c chr#" 
	parser = OptionParser(usage=usage, description=desc) 
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
	parser.add_option("-i", "--seed", type="string", dest="hifiSeed", help="Input seed file Name", default="null")
	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.hifiSeed == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	global chr_name
	
	options = get_args()
	chr_name = options.chrName
	seed_input_file = options.hifiSeed	
	
	start_time = time.time()
	#hifi_test(seed_input_file)
	#hifiAccuCheck("imputed_" + seed_input_file, chr_name)
	#hifiAccuCheck("imputedhaplotype_1.txt", chr_name)
	
	analyze_data()
	
	
	
	elapse_time = time.time() - start_time
	print "***********************to_impute_window time is: " + str(format(elapse_time, "0.3f")) + "s"












