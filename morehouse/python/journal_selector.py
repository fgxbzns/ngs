#!/usr/bin/python
#######################################################################################
# Guoxing Fu Nov 14, 2013
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

def load_data(input_filename):
	with open(input_filename, "r") as input_file:
		line = input_file.readline()
		while not line.startswith("EFETCH RESULT"):
			line = input_file.readline()
		while line != "":
			# journal line
			line = input_file.readline().strip()
			print line
			line = input_file.readline()
			while line != "\n" and line != "":
				line = input_file.readline()
			#print line
			
			# title line
			line = input_file.readline().strip()
			#print line
			line = input_file.readline()
			while line != "\n" and line != "":
				line = input_file.readline().strip()

			# author line
			line = input_file.readline().strip()
			line = input_file.readline()			
			while line != "\n" and line != "":
				line = input_file.readline().strip()
			
			# author info or abstract
			abstract = ""
			line = input_file.readline().strip()
			if line.startswith("Author information"):
				#line = input_file.readline()
				while line != "\n" and line != "":
					line = input_file.readline().strip()
				line = input_file.readline().strip()
												
			#line = input_file.readline()
			if line.startswith("PMID"):						# no author info, no abstract
				line = input_file.readline()
				line = input_file.readline()
			else:											# no author info, with abstract
				#abstract += line
				while line != "\n" and line != "":
					abstract += line + " "
					line = input_file.readline().strip()
					print line
				#print abstract
				
				line = input_file.readline()	
				#print line
				if line.startswith("PMID"):
					line = input_file.readline()
					line = input_file.readline()
				else:
					#print input_file.tell()
					print "error", line
				
			
			"""
			black_line = input_file.readline()
			print author_line
			black_line = input_file.readline()
			
			line = input_file.readline()
			print line
			"""

			
			
			"""
			journal_elements = line.split()
			year_index = journal_elements.index("2014")
			print year_index
			
			#journal_name = " ".join(journal_elements[2, year_index])
			#print journal_name
			while not line.startswith("PMID"):
				line = input_file.readline()
			line = input_file.readline()
			line = input_file.readline()
			"""
				
def get_args():
	desc=""
	usage = "" 
	parser = OptionParser(usage = usage, description=desc)
	parser.add_option("-i", "--ifile", type="string", dest="input_filename",help = "Input file Name", default="null")
	(options, args) = parser.parse_args()
 	"""
	parser.add_option("-e", "--eight", type="string", dest="hg18_name",help = "Input file name", default="null")
	parser.add_option("-n", "--nine", type="string", dest="hg19_name",help = "Input file name", default="null")
	parser.add_option("-d", "--del", type="string", dest="del_name",help = "Input file name", default="null")
	
	
	if options.hg18_name == "null" or options.hg19_name == "null" or options.del_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	"""
	return options

if __name__=='__main__':
	options = get_args()
	input_filename = options.input_filename
	
	start_time = time.time()
	load_data(input_filename)
	print "run time is: ", round((time.time() - start_time), 3), "s"
	
	
	
	
	
	
	
	