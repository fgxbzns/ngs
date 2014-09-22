#!/usr/bin/python

# calculate the seed percentage and hifi error in a centain range.

import re
import os, glob, subprocess, random, operator, time
from optparse import OptionParser
from tools import *

class data_class():
	def __init__(self):
		self.chr_name = ""

		self.seed_file_name = ""
		self.seed_file_prefix = ""
		self.seed_title_info = ""
		self.seed_dict = {}

def removeCT_ref(ref_file_name):
	with open(ref_file_name[:-3]+"_bis.fa", "w") as bis_file:
		with open(ref_file_name, "r") as ref_file:
			for line in ref_file:
				line = line.strip() if line.startswith(">") else re.sub('[ctCT]', '', line).strip()
				print >> bis_file, line

def removeCT_fastq(fastq_file_name):
	fastq_file = open(fastq_file_name, "r")
	fastq_bis_file = open(fastq_file_name[:-3] + "_bis.fastq", "w")
	line = fastq_file.readline()
	while line != "":
		line = fastq_file.readline()


def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
	parser.add_option("-i", "--fastq", type="string", dest="fastq", help="Input fastq file Name", default="null")
	parser.add_option("-r", "--ref", type="string", dest="ref", help="Input ref file Name", default="null")
	(options, args) = parser.parse_args()
	"""
	if options.chrName == "null" or options.hifiSeed == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	"""
	return options


if __name__ == '__main__':
	options = get_args()

	start_time = time.time()
	global data
	data = data_class()
	"""
	data.chr_name = options.chrName
	data.seed_file_name = options.hifiSeed
	data.hifi_result_file = options.hifiResult
	"""
	ref_file_name = options.ref
	removeCT_ref(ref_file_name)

	print "run time is: ", round((time.time() - start_time), 3), "s"
