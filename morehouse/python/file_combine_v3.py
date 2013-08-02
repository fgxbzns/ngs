#!/usr/bin/python

import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser

""" this works well for large, multiple files """

currentPath = os.getcwd() + '/'

def list_to_line(list):
	line = ""
	for a in list:
		line += a + "\t"
	return line

def output_title():
	output_file_name = "combined_vcf.txt"
	output_file = ouput_data(output_file_name, "w")
	title_info = "@chr \t pos \t"
	for file_info in file_name_list:
		file_name = file_info[0]
		elements_number = file_info[1]
		for i in range (0, elements_number):
			title_info += file_name + "\t"
	print >> output_file, title_info
	output_file.close()

def get_file_name_list():
	title_info = "chr \t pos \t"
	file_name_list = []
	for infile in glob.glob(os.path.join(currentPath,'*.vcf')):
		#file_name = infile[(infile.find('NA')):].strip()
		file_name = os.path.basename(infile)
		if file_name not in file_name_list:
			fp = open(file_name, "r")
			line = fp.readline()
			while not line.startswith(chr_name):
				line = fp.readline()
			elements = line.strip().split()
			elements_number = len(elements[2:])
			file_name_list.append((file_name, elements_number))
			fp.close()
	for file_info in file_name_list:
		file_name = file_info[0]
		elements_number = file_info[1]
		print file_name, "column number: ", elements_number
		for i in range (0, elements_number):
			title_info += file_name + "\t"		
	return (file_name_list, title_info)

def get_file_name_list_2(file_name_1, file_name_2):
	title_info = "chr \t pos \t"
	file_name_list = []
	for file_name in (file_name_1, file_name_2):
		if file_name not in file_name_list:
			fp = open(file_name, "r")
			line = fp.readline()
			while not line.startswith(chr_name):
				line = fp.readline()
			elements = line.strip().split()
			elements_number = len(elements[2:])
			file_name_list.append((file_name, elements_number))
			fp.close()
	for file_info in file_name_list:
		file_name = file_info[0]
		elements_number = file_info[1]
		print file_name, "column number: ", elements_number
		for i in range (0, elements_number):
			title_info += file_name + "\t"		
	return (file_name_list, title_info)

def combine_data(file_1, file_2, output_file_name):
	file_name_1 = file_1[0]
	elements_number_1 = file_1[1]
	
	file_name_2 = file_2[0]
	elements_number_2 = file_2[1]
	
	#output_file_name = file_name_1[:6] + "_" + file_name_2[:6] + ".txt"
	output_file = open(output_file_name, "w")
	print >> output_file, title_info
	
	fp_1 = open(file_name_1, "r")
	line_1 = fp_1.readline()
	while not line_1.startswith(chr_name):
		line_1 = fp_1.readline()
	
	fp_2 = open(file_name_2, "r")
	line_2 = fp_2.readline()
	while not line_2.startswith(chr_name):
		line_2 = fp_2.readline()
		
	while (line_1 != "" and line_1.startswith(chr_name)) or (line_2 != "" and line_2.startswith(chr_name)):
		if line_1 != "":
			list_1 = line_1.strip().split()
			try:
				chr_name_1 = list_1[0]
				position_1 = int(list_1[1])
			except:
				#print "line_1", line_1, file_name_1
				pass
		else:
			chr_name_1 = ""
		if line_2 != "":
			list_2 = line_2.strip().split()
			try:
				chr_name_2 = list_2[0]
				position_2 = int(list_2[1])
			except:
				#print "line_2", line_2, file_name_2
				pass
		else:
			chr_name_2 = ""
		
		if chr_name_1 == chr_name and chr_name_2 != chr_name:
			line = chr_name + "\t" + str(position_1) + "\t"
			line += list_to_line(list_1[2:])
			for i in range (0, elements_number_2):
				line += "\t"
			print >> output_file, line
			line_1 = fp_1.readline()
		elif chr_name_1 == chr_name and chr_name_2 == chr_name:
			if position_1 < position_2:
				line = chr_name + "\t" + str(position_1) + "\t"
				line += list_to_line(list_1[2:])
				for i in range (0, elements_number_2):
					line += "\t"
				print >> output_file, line
				line_1 = fp_1.readline()
			elif position_1 == position_2:
				line = chr_name + "\t" + str(position_1) + "\t"
				line += list_to_line(list_1[2:])
				line += list_to_line(list_2[2:])
				print >> output_file, line
				line_1 = fp_1.readline()
				line_2 = fp_2.readline()
			elif position_1 > position_2:
				line = chr_name + "\t" + str(position_2) + "\t"
				for i in range (0, elements_number_1):
					line += "\t"
				line += list_to_line(list_2[2:])
				print >> output_file, line
				line_2 = fp_2.readline()
		elif chr_name_1 != chr_name and chr_name_2 == chr_name:
			line = chr_name + "\t" + str(position_2) + "\t"
			for i in range (0, elements_number_1):
				line += "\t"
			line += list_to_line(list_2[2:])
			print >> output_file, line
			line_2 = fp_2.readline()
	fp_1.close()
	fp_2.close()
	output_file.close()
			
def file_combine():
	global file_name_list
	global title_info
	global chr_name
	
	options = get_args()
	chr_name = options.chrName
	a_file_name = options.a_file_name
	b_file_name = options.b_file_name
	if a_file_name == "null" and b_file_name == "null":
		file_tuple = get_file_name_list()
	else:
		file_tuple = get_file_name_list_2(a_file_name, b_file_name)
	
	file_name_list = file_tuple[0]
	title_info = file_tuple[1]
	
	#chr_name = "chrX"

	i = 1
	output_file_name = ""
	while i < len(file_name_list):
		if i == 1:
			output_file_name = "temp_" + chr_name + "_" + str(i) + ".txt"
			elements_number = file_name_list[0][1] + file_name_list[1][1]
			print "combining: ", file_name_list[0][0], file_name_list[1][0], "to: ", output_file_name
			combine_data(file_name_list[0], file_name_list[1], output_file_name)
		elif i > 1:
			temp_file = (output_file_name, elements_number)
			if i < len(file_name_list) - 1:
				output_file_name = "temp_" + chr_name + "_" + str(i) + ".txt"
			else:
				output_file_name = "combined_" + chr_name + "_" + str(len(file_name_list)) + "_files.txt"
			elements_number += file_name_list[i][1]
			print "combining: ", temp_file[0], file_name_list[i][0], "to: ", output_file_name
			combine_data(temp_file, file_name_list[i], output_file_name)
		i+=1

def get_args():
	desc="combine vcf file"
	usage = "file_combine -c chr_name" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
	parser.add_option("-a", "--afile", type="string", dest="a_file_name",help = "Input file Name", default="null")
	parser.add_option("-b", "--bfile", type="string", dest="b_file_name",help = "Input file Name", default="null")
	(options, args) = parser.parse_args()
	
	if options.chrName == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	
	return options

if __name__=='__main__':
	start = time.time()	
	file_combine()
	end = time.time()
	run_time = str(format((end - start), "0.3f"))
	print "run time is: " + run_time + "s"
	
	

