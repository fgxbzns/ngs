#!/usr/bin/python
#######################################################################################
# check rough similarity of two chr
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser

import re
from tools import *


def load_chr(file_name):
	chr_seq = ""
	fp = open(file_name, "r")
	for line in fp:
		chr_seq += line.strip()
	fp.close()
	return chr_seq

def select_seq(chr_seq, point_num):
	seq_dict = {}
	total_num = len(chr_seq)
	
	number_ceilling = int(math.ceil(float(total_num)/100)*100)
	num_in_each_part = number_ceilling/point_num
	
	for i in range(point_num):
		seq_dict[i*num_in_each_part] = chr_seq[i*num_in_each_part:(i*num_in_each_part+200)]
	return seq_dict

def search_seq(chr_seq, seq_dict):
	pass

def get_args():
	desc="calculate_maf"
	usage = "" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-e", "--en", type="string", dest="en_name",help = "Input enzyme file name", default="null")
	parser.add_option("-s", "--seq", type="string", dest="seq_name",help = "Input seq file name", default="null")
	(options, args) = parser.parse_args()
	if options.en_name == "null" or options.seq_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__=='__main__':
	#options = get_args()
	#file_name = options.ref_name
	enzyme_seq = "[CG][CA]"
	DNA_seq = "CCCGGTCCGACAAAAATGA"
	enzyme_search(enzyme_seq, DNA_seq)
	enzyme_seq = "GACNNNNRTGA"
	enzyme_seq = seq_convert(enzyme_seq)
	print enzyme_seq
	enzyme_search(enzyme_seq, DNA_seq)

