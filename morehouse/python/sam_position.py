#!/usr/bin/python

#######################################################################################
# Find the reads that cover a position
# Author: Guoxing Fu
# Sep 26, 2006
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser

def get_args():
	desc="variation call"
	usage = "sam_position -s sam_file -p position" 
	parser = OptionParser(usage = usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chr",help = "chr", default="null")
	parser.add_option("-p", "--pos", type="string", dest="pos",help = "", default="null")
	(options, args) = parser.parse_args()
	if options.samFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

def sam_process():
	options = get_args()
	sam_file_name = options.samFile
	chr = options.chr
	pos = int(options.pos)
	
	sam_file = open(sam_file_name, "r")
	print "chr: ", chr
	print "searching position: ", pos
	for lines in sam_file:
		if not lines.startswith("@"):
			elements = lines.strip().split()
			if elements[2] == chr:
				try:
					read_ID = elements[0]
					start_postion = int(elements[3])
					read_seq = elements[9]
					read_length = len(read_seq)
					if pos >= start_postion and pos <= (start_postion + read_length):
						print read_ID, start_postion, read_seq, elements[10]
				except ValueError:
					print "error in ", lines
				#pass

if __name__=='__main__':
	sam_process()
		

