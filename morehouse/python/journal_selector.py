#!/usr/bin/python
#######################################################################################
# Guoxing Fu Feb 17, 2014
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

class data_class():
	def __init__(self):
		self.keywords_dict = {}
		self.pre_keywords_dict = {}
		self.non_keywords_dict = {}
		self.symbol_dict = {}

		self.article_file_name = ""
		self.pre_keywords_file_name = "pre_keywords_list.txt"
		self.non_keywords_file_name = "non_keywords_list.txt"
		self.symbol_file_name = "symbol.txt"

def update_keywords_dict(keyword, keywords_dict):
	if keyword not in data.non_keywords_dict and \
					keyword not in data.pre_keywords_dict and keyword not in data.symbol_dict:
		if keyword in keywords_dict:
			keywords_dict[keyword] = keywords_dict[keyword] + 1
		else:
			keywords_dict[keyword] = 1

def load_keywords(input_filename):
	dict = {}
	with open(input_filename, "r") as input_file:
		for line in input_file:
			dict[line.strip().split()[0]] = ""
	return dict

def output_keywords(file_name, keywords_dict):
	with open(file_name, "w") as output_file:
		for keyword, number in keywords_dict.iteritems():
			print >> output_file, keyword, number

def keyword_process(keyword):
	keyword = keyword.lower()
	for symbol in data.symbol_dict.keys():
		if symbol in keyword:
			return "remove_keyword"
			#keyword = keyword.replace(symbol, "")
			#print "removed", symbol, "from", keyword
	# remove number
	try:
		keyword = int(keyword)
		return "remove_keyword"
	except:
		try:
			first_letter = int(keyword[0])
			return "remove_keyword"
		except:
			return keyword
	# first letter non digit


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
			old methodnew_keywords_dict.txt
			#print line
			journal_elements = journal_line.split()
			try:
				year_index = journal_elements.index("2014")
			except:
				print "year_index error", journal_line
			#	break
			
			#print year_index
			
			journal_name = " "with open(input_filename, "r") as input_file:.join(journal_elements[1:year_index])
			#print journal_name
			"""
			dot_index = journal_line.find(".")
			year_index = journal_line.find("201")
			journal_name = journal_line[dot_index + 2:year_index - 1]
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
			if line.startswith("PMID"):  # no author info, no abstract
				line = input_file.readline()
				line = input_file.readline()
				if line.startswith("]"):
					line = input_file.readline()
					line = input_file.readline()
					line = input_file.readline()
					line = input_file.readline()
			else:  # no author info, with abstract
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
						keyword = keyword_process(keyword)
						if len(keyword) > 1 and keyword != "remove_keyword":
							update_keywords_dict(keyword, keywords_dict)

	keywords_sorted_list = sort_dict_by_value(keywords_dict)
	for data in keywords_sorted_list:
		if data[1] > 1:
			print data[0], data[1]
	#for i in range(50):
	#	print keywords_sorted_list[i][0], keywords_sorted_list[i][1]

	print "article_total_number", article_total_number


def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--ifile", type="string", dest="input_filename", help="Input file Name", default="null")
	parser.add_option("-n", "--nfile", type="string", dest="input_filename", help="Input file Name", default="null")
	parser.add_option("-p", "--pfile", type="string", dest="input_filename", help="Input file Name", default="null")

	(options, args) = parser.parse_args()
	"""
	if options.hg18_name == "null" or options.hg19_name == "null" or options.del_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	"""
	return options


if __name__ == '__main__':
	options = get_args()
	input_filename = options.input_filename

	start_time = time.time()
	global data
	data = data_class()
	data.pre_keywords_dict = load_keywords(data.pre_keywords_file_name)
	data.non_keywords_dict = load_keywords(data.non_keywords_file_name)
	data.symbol_dict = load_keywords(data.symbol_file_name)
	#print data.symbol_dict
	#print len(data.symbol_file_name)
	"""
	new_dict = dict_substract(data.non_keywords_dict, data.pre_keywords_dict)
	print len(new_dict)
	output_keywords("new_keywords_dict.txt", new_dict)
	"""
	load_data(input_filename)
	print "run time is: ", round((time.time() - start_time), 3), "s"
	
	
	
	
