#!/usr/bin/python

# location /home/guoxing/tool/morehouse

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

file_path = "/home/guoxing/disk2/solid/common_files/"
data_record_path = "/home/guoxing/disk2/solid/common_files/data_record/"
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-a", "--all", type="string", dest="alldelFile",help = "Input all error deleted file Name", default="null")
parser.add_option("-i", "--one", type="string", dest="errorFile",help = "Input error file Name", default="null")

(options, args) = parser.parse_args()
alldel_name = options.alldelFile
error_file_name = options.errorFile

alldel_file = open(currentPath + alldel_name, "r")
error_file = open(currentPath + error_file_name, "r")
output_file = open(currentPath + "error_compare.txt", "w")
print >> output_file, "position \t alldel \t single_error_only"

# skip the title line
alldel_line = alldel_file.readline()
error_line = error_file.readline()

alldel_line = alldel_file.readline()
error_line = error_file.readline()

same_position_number = 0

while alldel_line != "":
	alldel_elements = alldel_line.strip().split()
	position = alldel_elements[0].strip()
	alldel_hifi_error = ""
	error_hifi_error = ""
	#if len(alldel_elements) > 2:
	if True:
		try:
			if alldel_elements[-1].strip() == "140":
				alldel_hifi_error = alldel_elements[-1].strip()
		except:
			pass
	error_elements = error_line.strip().split()
	#if len(error_elements) > 2:
	if True:
		try:
			if error_elements[-1].strip() == "140":
				error_hifi_error = error_elements[-1].strip()
		except:
			pass
#	if alldel_hifi_error != "" and error_hifi_error != "" and alldel_hifi_error == error_hifi_error:
	if alldel_hifi_error == "140" and error_hifi_error == "140":
		same_position_number += 1
		error_hifi_error = ""
	new_line = position + "\t" +alldel_hifi_error +"\t" + error_hifi_error
	print >> output_file, new_line
	alldel_line = alldel_file.readline()
	error_line = error_file.readline()

print same_position_number
	
alldel_file.close()
error_file.close()
output_file.close()
