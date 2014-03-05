#!/usr/bin/python
# -*- coding: utf-8 -*-
#######################################################################################
# Guoxing Fu Feb 17, 2014
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy, string
from optparse import OptionParser

import unicodedata as ud

from tools import *
#print type('β')

#print ud.normalize('NFC', u'\xce\xA9')

greek_alphabet = {
u'\u0391': 'Alpha',
u'\u0392': 'Beta',
u'\u0393': 'Gamma',
u'\u0394': 'Delta',
u'\u0395': 'Epsilon',
u'\u0396': 'Zeta',
u'\u0397': 'Eta',
u'\u0398': 'Theta',
u'\u0399': 'Iota',
u'\u039A': 'Kappa',
u'\u039B': 'Lamda',
u'\u039C': 'Mu',
u'\u039D': 'Nu',
u'\u039E': 'Xi',
u'\u039F': 'Omicron',
u'\u03A0': 'Pi',
u'\u03A1': 'Rho',
u'\u03A3': 'Sigma',
u'\u03A4': 'Tau',
u'\u03A5': 'Upsilon',
u'\u03A6': 'Phi',
u'\u03A7': 'Chi',
u'\u03A8': 'Psi',
u'\u03A9': 'Omega',
u'\u03B1': 'alpha',
u'\u03B2': 'beta',
u'\u03B3': 'gamma',
u'\u03B4': 'delta',
u'\u03B5': 'epsilon',
u'\u03B6': 'zeta',
u'\u03B7': 'eta',
u'\u03B8': 'theta',
u'\u03B9': 'iota',
u'\u03BA': 'kappa',
u'\u03BB': 'lamda',
u'\u03BC': 'mu',
u'\u03BD': 'nu',
u'\u03BE': 'xi',
u'\u03BF': 'omicron',
u'\u03C0': 'pi',
u'\u03C1': 'rho',
u'\u03C3': 'sigma',
u'\u03C4': 'tau',
u'\u03C5': 'upsilon',
u'\u03C6': 'phi',
u'\u03C7': 'chi',
u'\u03C8': 'psi',
u'\u03C9': 'omega',
}

class keywords_class():
	def __init__(self):
		self.number = 0
		self.journal_name_dict = {}

class data_class():
	def __init__(self):
		self.keywords_dict = {}
		self.pre_keywords_dict = {}
		self.non_keywords_dict = {}
		self.symbol_dict = {}
		self.new_symbol_dict = {}
		self.new_symbol_list = []

		self.article_file_name = ""
		self.pre_keywords_file_name = "pre_keywords_list.txt"
		self.non_keywords_file_name = "non_keywords_list.txt"
		self.symbol_file_name = "symbol.txt"
		self.new_symbol_file_name = "new_symbol.txt"

		self.new_symbol_file = open(self.new_symbol_file_name, 'w')

		self.letter_digit = dict.fromkeys(string.ascii_lowercase, 0)
		for i in range(10):
			self.letter_digit[str(i)] = 0
		#print self.letter_digit

		#self.english_dict = enchant.Dict("en_US")


def update_keywords_dict(keyword, keywords_dict):
	# update keywords from all papers.
	if keyword not in data.non_keywords_dict and \
					keyword not in data.pre_keywords_dict and keyword not in data.symbol_dict:
		if keyword in keywords_dict:
			keywords_dict[keyword] += 1
		else:
			keywords_dict[keyword] = 1

def update_journal_keywords_dict(keyword, journal_name, keywords_dict):
	# update keywords in all papers and in different journals.
	if keyword not in data.non_keywords_dict and \
					keyword not in data.pre_keywords_dict and keyword not in data.symbol_dict:
		if keyword in keywords_dict:
			keywords_dict[keyword][0] += 1
			if journal_name in keywords_dict[keyword][1]:
				keywords_dict[keyword][1][journal_name] += 1
			else:
				keywords_dict[keyword][1][journal_name] = 1
		else:
			journal_name_dict = {}
			keywords_dict[keyword] = [1, journal_name_dict]
			keywords_dict[keyword][1][journal_name] = 1

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
	'''
	# replace alpha with a
	for greek_letter in greek_alphabet:
		if greek_letter in keyword:
			keyword = keyword.replace(greek_letter, greek_alphabet[greek_letter])
			print "replaced", greek_letter, "from", keyword
	'''
	# remove keyword with symbol
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
			# first letter non digit
			first_letter = int(keyword[0])
			return "remove_keyword"
		except:
			for letter in keyword:
				if letter not in data.letter_digit:

						#data.new_symbol_list.append(keyword)
						'''
						if keyword not in data.new_symbol_dict:
							data.non_keywords_dict[keyword] = 1
						else:
							data.non_keywords_dict[keyword] += 1
						'''
						#print keyword
						#print keyword, "is not english"
						print >> data.new_symbol_file, keyword
						return "remove_keyword"
			return keyword



def load_data(input_filename):
	keywords_dict = {}
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
			'''
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
			'''
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
							#update_keywords_dict(keyword, keywords_dict)
							update_journal_keywords_dict(keyword, journal_name, keywords_dict)

	keywords_sorted_list = sort_dict_by_value(keywords_dict)

	for data in keywords_sorted_list:
		if data[1][0] > 2:
			print data[0], data[1][0], data[1][1]
	"""
	#for i in range(50):
	#	print keywords_sorted_list[i][0], keywords_sorted_list[i][1]


	#print data.new_symbol_dict
	for item in data.new_symbol_list:
		#print item.encode('utf-8')
		print item
	print data.new_symbol_list
	"""
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
	#path = "/home/guoxing/node1/disk2/ncbi_download/"
	#os.chdir(path)

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

	#input_filename = "2014[year]_ncbi_output_5000.txt"
	load_data(input_filename)
	data.new_symbol_file.close()
	print "run time is: ", round((time.time() - start_time), 3), "s"
	
	
	
	
