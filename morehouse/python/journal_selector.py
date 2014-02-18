#!/usr/bin/python
# -*- coding: iso-8859-15 -*-
#######################################################################################
# Guoxing Fu Feb 17, 2014
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *


def update_keywords_dict(keyword, keywords_dict):
	if keyword in keywords_dict:
		keywords_dict[keyword] = keywords_dict[keyword] + 1
	else:
		keywords_dict[keyword] = 1

def load_data(input_filename):
	
	keywords_dict = {}
	exclusive_list = {}
	article_total_number = 0
	
	with open(input_filename, "r") as input_file:
		line = input_file.readline()
		while not line.startswith("EFETCH RESULT"):
			line = input_file.readline()
		while line != "":
			# journal line
			
			journal_line = ""
			line = input_file.readline().strip()
			while line != "\n" and line != "":
					journal_line += line + " "
					line = input_file.readline().strip()
			#print journal_line
			"""
			old method
			#print line
			journal_elements = journal_line.split()
			try:
				year_index = journal_elements.index("2014")
			except:
				print "year_index error", journal_line
			#	break
			
			#print year_index
			
			journal_name = " ".join(journal_elements[1:year_index])
			#print journal_name
			"""
			dot_index = journal_line.find(".")
			year_index = journal_line.find("201")
			journal_name = journal_line[dot_index+2:year_index-1]
			article_total_number += 1
			#print journal_name
						
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
			
			# Collaborators
			abstract = ""
			line = input_file.readline().strip()
			if line.startswith("Collaborators"):
				while line != "\n" and line != "":
					line = input_file.readline().strip()
				line = input_file.readline().strip()
			
			# author info or abstract
			if line.startswith("Author information"):
				#line = input_file.readline()
				while line != "\n" and line != "":
					line = input_file.readline().strip()
				line = input_file.readline().strip()
												
			#line = input_file.readline()
			
			if line.startswith("Copyright"):
				line = input_file.readline()
				line = input_file.readline()
			if line.startswith("PMID"):						# no author info, no abstract
				line = input_file.readline()
				line = input_file.readline()
				if line.startswith("]"):
					line = input_file.readline()
					line = input_file.readline()
					line = input_file.readline()	
					line = input_file.readline()	
			else:											# no author info, with abstract
				#abstract += line
				while line != "\n" and line != "":
					abstract += line + " "
					line = input_file.readline().strip()
					#print line
				#print abstract
				
				line = input_file.readline()	
				#print line 
				if line.startswith("Copyright"):
					line = input_file.readline()
					line = input_file.readline()
				elif not line.startswith("PMID"):
					while not line.startswith("PMID") and line != "":
						line = input_file.readline()
				
				if line.startswith("PMID"):
					line = input_file.readline()
					line = input_file.readline()
					if line.startswith("]"):
						line = input_file.readline()
						line = input_file.readline()
						line = input_file.readline()
						line = input_file.readline()
				else:
					
					#print input_file.tell()
					#print "error", line
					pass
					
				# update dicts
				if True:
				#if journal_name == "Ann Lab Med.":
					for keyword in abstract.strip().split():
						update_keywords_dict(keyword, keywords_dict)
					
	keywords_sorted_list = sort_dict_by_value(keywords_dict)
	for data in keywords_sorted_list:
		print data[0], data[1]
	#for i in range(500):
	#	print keywords_sorted_list[i][0], keywords_sorted_list[i][1]
			
	print "article_total_number", article_total_number
				
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
	
	
	
	
	
	
	
	