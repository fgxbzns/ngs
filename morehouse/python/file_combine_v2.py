#!/usr/bin/python

import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser

""" this works well for small files, saves data in dict"""

currentPath = os.getcwd() + '/'

def sort_dict_by_key(input_dict):
	sorted_list = []
	sorted_list = [x for x in input_dict.iteritems()] 
	sorted_list.sort(key=lambda x: x[0]) # sort by key
	return sorted_list

def list_to_line(list):
	line = ""
	for a in list:
		line += a + "\t"
	return line

def get_file_name_list():
	file_name_list = []
	for infile in glob.glob(os.path.join(currentPath,'*.vcf')):
		file_name = infile[(infile.find('NA')):].strip()
		if file_name not in file_name_list:
			fp = open(file_name, "r")
			line = fp.readline()
			while not line.startswith("chr"):
				line = fp.readline()
			elements = line.strip().split()
			elements_length = len(elements[2:])
			file_name_list.append((file_name, elements_length))
			fp.close()
	for file_info in file_name_list:
		print file_info[0], "column length: ", file_info[1]
	return file_name_list

def load_raw_data(file_name, data_dict):			
	fp = open(file_name, "r")
	for line in fp:
		if line.startswith("chr"):
			elements = line.strip().split()
			try:
				chr_name = elements[0]
				position = elements[1]
				if chr_name not in data_dict:
					data_dict[chr_name] = {}
				if position not in data_dict[chr_name]:
					data_dict[chr_name][position] = []
				data_dict[chr_name][position].append((file_name, elements))
			except ValueError:
				#print "error in ", line
				pass
	fp.close()
	return data_dict

def ouput_data(output_file_name, data_dict):
	output_file = open(currentPath + output_file_name, "w")
	title_info = "chr \t pos \t"
	for file_info in file_name_list:
		file_name = file_info[0]
		elements_length = file_info[1]
		for i in range (0, elements_length):
			title_info += file_name + "\t"
	print >> output_file, title_info
	
	for chr, chr_dict in data_dict.iteritems():
		chr_sorted_list = sort_dict_by_key(chr_dict)
		for data in chr_sorted_list:
			line = ""
			position = data[0]
			line += chr + "\t"
			line += position + "\t"
			info_list = data[1]
			for file_info in file_name_list:
				file_name = file_info[0]
				elements_length = file_info[1]
				for info in info_list:
					input_file_name = info[0]
					elements = info[1][2:]	# remove chr and position
					if file_name == input_file_name:
						line += list_to_line(elements) + "\t"
					else:
						for i in range (0, elements_length):
							line += "\t"
			print >> output_file, line
	output_file.close()

def file_combine():
	global data_dict
	global file_name_list
	data_dict = {}
	file_name_list = []

	file_name_list = get_file_name_list()


	#list = ["NA12877_S1_XY_5000.txt", "NA12878_S1_XY_5000.txt"]
	for file_info in file_name_list:
		file_name = file_info[0]
		data_dict = load_raw_data(file_name, data_dict)
	
	output_file_name = "combined_vxf.txt"
	ouput_data(output_file_name, data_dict)
	
def get_args():
	desc="combine vcf file"
	usage = "file_combine" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
	parser.add_option("-i", "--file", type="string", dest="file_name",help = "Input file Name", default="null")
	(options, args) = parser.parse_args()
	"""
	if options.chrName == "null" or options.file_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	"""
	return options

if __name__=='__main__':
	options = get_args()
	chr_name = options.chrName
	file_name = options.file_name	
	
	start = time.time()	
	file_combine()
	end = time.time()
	run_time = str(format((end - start), "0.3f"))
	print "run time is: " + run_time + "s"
	
	

