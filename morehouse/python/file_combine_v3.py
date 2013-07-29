#!/usr/bin/python

import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser

""" this works well for small files """

currentPath = os.getcwd() + '/'

def list_to_line(list):
	line = ""
	for a in list:
		line += a + "\t"
	return line

def output_title():
	output_file_name = "combined_vxf.txt"
	output_file = ouput_data(output_file_name, "w")
	title_info = "@chr \t pos \t"
	for file_info in file_name_list:
		file_name = file_info[0]
		elements_length = file_info[1]
		for i in range (0, elements_length):
			title_info += file_name + "\t"
	print >> output_file, title_info
	output_file.close()

def get_file_name_list(file_name_1, file_name_2, chr_name):
	file_name_list = []
	for file_name in (file_name_1, file_name_2):
		if file_name not in file_name_list:
			fp = open(file_name, "r")
			line = fp.readline()
			while not line.startswith(chr_name):
				line = fp.readline()
			elements = line.strip().split()
			elements_length = len(elements[2:])
			file_name_list.append((file_name, elements_length))
			fp.close()
	for file_info in file_name_list:
		print file_info[0], "column length: ", file_info[1]
	return file_name_list

def combine_data(file_name_1, file_name_2, chr_name):
	output_file_name = "combined_vcf.txt"
	output_file = ouput_data(output_file_name, "w")
	elements_length_1 = file_name_list[0][1]
	elements_length_2 = file_name_list[1][1]

	for file_info in file_name_list:
		file_name = file_info[0]
		elements_length = file_info[1]
		for i in range (0, elements_length):
			title_info += file_name + "\t"
	print >> output_file, title_info
	
	fp_1 = open(file_name_1, "r")
	line_1 = fp_1.readline()
	while not line_1.startswith(chr_name):
		line_1 = fp_1.readline()
	
	fp_2 = open(file_name_2, "r")
	line_2 = fp_2.readline()
	while not line_2.startswith(chr_name):
		line_2 = fp_2.readline()
		
	while line_1 != "" or line_2 != "":
		list_1 = line_1.strip().split()
		chr_name_1 = list_1[0]
		position_1 = int(list_1[1])
		
		list_2 = line_1.strip().split()
		chr_name_2 = list_2[0]
		position_2 = int(list_2[1])
		
		if chr_name_1 == chr_name and chr_name_2 == chr_name:
			if position_1 < position_2:
				line = chr_name + "\t" + position_1 + "\t"
				line += list_to_line(list_1[2:])
				for i in range (0, len(list_2[2:])):
					line += "\t"
				print >> output_file, line
				line_1 = fp_1.readline()
			elif position_1 == position_2:
				line = chr_name + "\t" + position_1 + "\t"
				line += list_to_line(list_1[2:])
				line += list_to_line(list_2[2:])
				print >> output_file, line
				line_1 = fp_1.readline()
				line_2 = fp_2.readline()
			elif position_1 > position_2:
				line = chr_name + "\t" + position_2 + "\t"
				for i in range (0, len(list_1[2:])):
					line += "\t"
				line += list_to_line(list_2[2:])
				print >> output_file, line
				line_2 = fp_2.readline()
		elif chr_name_1 == chr_name:
			line = chr_name + "\t" + position_1 + "\t"
				line += list_to_line(list_1[2:])
				for i in range (0, elements_length_2)):
					line += "\t"
				print >> output_file, line
				line_1 = fp_1.readline()
			
		
		

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
	
	

