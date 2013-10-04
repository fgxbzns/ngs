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


def analyze_data():
	window_info_dict = get_window_info()
	window_info_sorted_list = sort_dict_by_key(window_info_dict)
	print window_info_sorted_list[0][0], window_info_sorted_list[0][1]
	window_size_total = 0
	for data in window_info_sorted_list:
		pos = data[0]
		info = data[1]
		window_size_total += int(info[7])
	
	print "window_size_total", window_size_total
	print "window_size_average", window_size_total/len(window_info_dict)
		
	


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
	
	options = get_args()
	chr_name = options.chrName
	seed_input_file = options.hifiSeed	
	
	start_time = time.time()
	#hifi_test(seed_input_file)
	#hifiAccuCheck("imputed_" + seed_input_file, chr_name)
	
	analyze_data()
	
	
	
	elapse_time = time.time() - start_time
	print "***********************to_impute_window time is: " + str(format(elapse_time, "0.3f")) + "s"












