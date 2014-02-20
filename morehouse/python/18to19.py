#!/usr/bin/python
#######################################################################################
# Guoxing Fu Nov 14, 2013
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

def load_delete_pos(del_name):
	del_pos_dict = {}
	with open(del_name, "r") as fp:
		for line in fp:
			if line.startswith('chr'):
				del_chr = line[3:line.index(':')]
				pos = int(line[(line.strip().index('-')+1):])
				del_pos_dict[pos] = del_chr
	return del_pos_dict

def combine_files(hg18_name, hg19_name, del_name):
	del_pos_dict = load_delete_pos(del_name)
	hg18_file = open(hg18_name, "r")
	hg19_file = open(hg19_name, "r")
	output_file = open(hg18_name + "_combined", "w")
	
	hg18_line_ori = hg18_file.readline().strip().split()	# skip first line
	print >> output_file, hg18_line_ori[0] + "\t POS_18 \t POS_19 \t" + list_to_line(hg18_line_ori[2:])
	
	hg18_line_ori = hg18_file.readline()
	hg19_line = hg19_file.readline().strip()
	
	k = 0
	while hg18_line_ori != "":
		k+=1
		if k%100000==0:
			print 'Processing line ', str(k)
		hg18_line = hg18_line_ori.strip().split()
		try:
			hg18_chr = hg18_line[0]
			hg18_pos = int(hg18_line[1])
			while hg18_pos in del_pos_dict and hg18_chr == del_pos_dict[hg18_pos]:
				hg18_line_ori = hg18_file.readline()
				hg18_line = hg18_line_ori.strip().split()
				hg18_chr = hg18_line[0]
				hg18_pos = int(hg18_line[1])
					
			hg19_chr = hg19_line[3:hg19_line.index(':')]
			hg19_pos = int(hg19_line[(hg19_line.index('-')+1):])
			if hg18_chr == hg19_chr:
				print >> output_file, hg18_line[0] + "\t" + str(hg18_pos) + "\t" + str(hg19_pos) + "\t" + list_to_line(hg18_line[2:])
			else:
				print "chr not match", hg18_line_ori, hg19_line
				sys.exit(1)
		except:
			print "error in", hg18_line_ori
			sys.exit(1)
		hg18_line_ori = hg18_file.readline()
		hg19_line = hg19_file.readline().strip()
	
	hg18_file.close()	
	output_file.close()	
	hg19_file.close()	
	
def get_args():
	desc="./18to19.py -e hg18 -n hg19 -d delete"
	usage = ""
	parser = OptionParser(usage = usage, description=desc)
	parser.add_option("-e", "--eight", type="string", dest="hg18_name",help = "Input file name", default="null")
	parser.add_option("-n", "--nine", type="string", dest="hg19_name",help = "Input file name", default="null")
	parser.add_option("-d", "--del", type="string", dest="del_name",help = "Input file name", default="null")
	(options, args) = parser.parse_args()
	
	if options.hg18_name == "null" or options.hg19_name == "null" or options.del_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__=='__main__':
	options = get_args()
	hg18_name = options.hg18_name
	hg19_name = options.hg19_name
	del_name = options.del_name
	
	start_time = time.time()
	combine_files(hg18_name, hg19_name, del_name)
	elapsed_time = time.time() - start_time
	print "elapsed_time is: " + str(format(elapsed_time, "0.3f")) + "s"		
	
	
	
	
	
	
	