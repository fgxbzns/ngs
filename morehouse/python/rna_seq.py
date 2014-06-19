#!/usr/bin/python
#######################################################################################
# Guoxing Fu Nov 14, 2013
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *


def get_size(start_path):
	print "getting size of", start_path
	total_size = 0
	MB=1024*1024.0
	GB=1024*MB

	for dirpath, dirnames, filenames in os.walk(start_path):
		for f in filenames:
			fp = os.path.join(dirpath, f)
			#print fp
			total_size += os.path.getsize(fp)
	unit="GB" if total_size>=GB else "MB"
	total_size=round(total_size/GB, 3) if total_size>=GB else round(total_size/MB, 3)
	return total_size, unit

def get_args():
	desc = "./18to19.py -e hg18 -n hg19 -d delete"
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-f", "--folder", type="string", dest="folder_name", help="Input folder name", default="null")
	parser.add_option("-n", "--ninety", type="string", dest="hg19_name", help="Input file name", default="null")
	parser.add_option("-d", "--del", type="string", dest="del_name", help="Input file name", default="null")
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
	folder_name = options.folder_name

	start_time = time.time()
	print get_size(folder_name)
	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"
	
	
	
	
	
	
	