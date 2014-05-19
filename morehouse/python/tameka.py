#!/usr/bin/python
#######################################################################################
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *

def process_file(file_name):
	with open(file_name + "_processed", "w") as output_file:
		with open(file_name, "r") as input_file:
			for line in input_file:
				elements = line.strip().split()
				temp_list = []
				temp_list.append(elements[0])
				for data in elements[1:]:
					temp_data = float(data)
					if temp_data != 0 and temp_data > 1:
						temp_data = 1/temp_data
						temp_list.append(temp_data)
					else:
						temp_list.append(temp_data)
				print >> output_file, list_to_line(temp_list)

def get_args():
	desc="calculate_maf"
	usage = "" 
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--ref", type="string", dest="ref_name", help="Input ref file name", default="null")
	(options, args) = parser.parse_args()
	if options.ref_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	options = get_args()
	file_name = options.ref_name
	process_file(file_name)















			
