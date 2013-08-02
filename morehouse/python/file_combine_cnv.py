#!/usr/bin/python

import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser

""" remove snps falls in to cnv """

currentPath = os.getcwd() + '/'

# file 1, snp file, file_2, cnv file
def remove_snp_in_cnv(file_name_1, file_name_2):
	
	fp_1 = open(file_name_1, "r")
	line_1 = fp_1.readline()
	title_info = line_1
	while not line_1.startswith(chr_name):
		line_1 = fp_1.readline()
	
	output_file_name = file_name_1 + "_cnv_removed"
	output_file = open(output_file_name, "w")
	print >> output_file, title_info
	
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
				start_position = int(list_2[1])
				end_position = int(list_2[2])
			except:
				#print "line_2", line_2, file_name_2
				pass
		else:
			chr_name_2 = ""
		
		if chr_name_1 == chr_name and chr_name_2 != chr_name:
			print >> output_file, line_1.strip()
			line_1 = fp_1.readline()
			
		elif chr_name_1 == chr_name and chr_name_2 == chr_name:
			if position_1 < start_position:
				print >> output_file, line_1.strip()
				line_1 = fp_1.readline()
			elif position_1 >= start_position and position_1 <= end_position:
				print "remove:", position_1, "in", start_position, end_position
				line_1 = fp_1.readline()
			else:
				line_2 = fp_2.readline()
			
		elif chr_name_1 != chr_name and chr_name_2 == chr_name:
			line_2 = fp_2.readline()
		else:
			print "error", line_1, line_2
	fp_1.close()
	fp_2.close()
	output_file.close()
		
def file_combine():
	global chr_name
	
	options = get_args()
	chr_name = options.chrName
	a_file_name = options.a_file_name
	b_file_name = options.b_file_name
	
	remove_snp_in_cnv(a_file_name, b_file_name)

def get_args():
	desc="combine vcf file"
	usage = "file_combine -c chr_name" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
	parser.add_option("-a", "--afile", type="string", dest="a_file_name",help = "Input file Name", default="null")
	parser.add_option("-b", "--bfile", type="string", dest="b_file_name",help = "Input file Name", default="null")
	(options, args) = parser.parse_args()
	
	if options.chrName == "null" or options.a_file_name == "null" or options.b_file_name == "null":
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
	
	

