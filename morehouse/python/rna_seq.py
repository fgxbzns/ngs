#!/usr/bin/python
#######################################################################################
# Guoxing Fu Nov 14, 2013
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

def convert_unit(size):
	MB=1024*1024.0
	GB=1024*MB
	unit = "GB" if size >= GB else "MB"
	new_size = round(size/GB, 2) if size >= GB else round(size/MB, 2)
	return new_size, unit

def get_size(folder_name):
	print "getting size of", folder_name
	total_size = 0
	MB=1024*1024.0
	GB=1024*MB
	file_size_list=[]

	for dirpath, dirnames, filenames in os.walk(folder_name):
		#print dirpath
		#print dirnames
		#print filenames
		for f in filenames:
			fp = os.path.join(dirpath, f)
			#print fp
			file_size = os.path.getsize(fp)
			total_size += file_size

			#print f, convert_unit(file_size)
			file_size_list.append((f, convert_unit(file_size)))

	total_size = convert_unit(total_size)
	#return total_size, unit
	return file_size_list

def output_size(folder_name):
	#folder_list = os.walk(folder_name)

	with open(folder_name + "_size.txt", "w") as output_file:
		print >> output_file, folder_name
		for subfolder in range(25):
			file_size_list = get_size(folder_name+"/"+str(subfolder))
			print >> output_file, subfolder
			for data in file_size_list:
				#print >> output_file, data[0], data[1][0], data[1][1]
				print >> output_file, data[1][0],
			print >> output_file, ""

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
	output_size(folder_name)
	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"
	
	
	
	
	
	
	