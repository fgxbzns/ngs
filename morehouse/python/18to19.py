#!/usr/bin/python
#######################################################################################
# Guoxing Fu Aug 29, 2013
# find enzyme cut seq in dna seq
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser

import re
from tools import *

def load_delete_pos(del_name):
	del_pos = {}
	with open(del_name, "r") as fp:
		for line in fp:
			if line.startswith('chr'):
				pos = int(line[(line.strip().index('-')+1):])
				print pos
				
def get_args():
	desc="./18to19.py -e hg18 -n hg19 -d delete"
	usage = "" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-e", "--e", type="string", dest="hg18_name",help = "Input file name", default="null")
	parser.add_option("-n", "--nine", type="string", dest="hg19_name",help = "Input file name", default="null")
	parser.add_option("-d", "--del", type="string", dest="del_name",help = "Input file name", default="null")
	(options, args) = parser.parse_args()
	"""
	if options.hg18_name == "null" or options.hg19_name == "null" or options.del_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	"""
	return options

if __name__=='__main__':
	options = get_args()
	hg18_name = options.hg18_name
	hg19_name = options.hg19_name
	del_name = options.del_name


	load_delete_pos(del_name)