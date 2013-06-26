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
output_record_file = open(currentPath + "error_record.txt", "w")
print >> output_file, "position \t alldel \t single_error_only"

# skip the title line
alldel_line = alldel_file.readline()
error_line = error_file.readline()

alldel_line = alldel_file.readline()
error_line = error_file.readline()

same_position_number = 0
error_position_list = []

while alldel_line != "":
	alldel_elements = alldel_line.strip().split()
	position = alldel_elements[0].strip()
	alldel_hifi_error = ""
	error_hifi_error = ""
	#if len(alldel_elements) > 2:
	if True:
		try:
			if alldel_elements[-1].strip() == "80":
				alldel_hifi_error = alldel_elements[-1].strip()
		except:
			pass
	error_elements = error_line.strip().split()
	#if len(error_elements) > 2:
	if True:
		try:
			if error_elements[-1].strip() == "80":
				error_hifi_error = error_elements[-1].strip()
		except:
			pass
#	if alldel_hifi_error != "" and error_hifi_error != "" and alldel_hifi_error == error_hifi_error:
	if alldel_hifi_error == "80" and error_hifi_error == "80":
		same_position_number += 1
		error_hifi_error = ""	
	new_line = position + "\t" +alldel_hifi_error +"\t" + error_hifi_error
	print >> output_file, new_line
	if alldel_hifi_error != "80" and error_hifi_error == "80":
		error_position_list.append(int(position))
	alldel_line = alldel_file.readline()
	error_line = error_file.readline()

print "same_error_position_number ", same_position_number
print "Total error caused by this position", len(error_position_list)
print "smallest error position", min(error_position_list)
print "largest error position", max(error_position_list)

print >> output_record_file, "error caused by hifi ", same_position_number
print >> output_record_file, "Total error caused by this position", len(error_position_list)
print >> output_record_file, "smallest error position", min(error_position_list)
print >> output_record_file, "largest error position", max(error_position_list)
	
alldel_file.close()
error_file.close()
output_file.close()
output_record_file.close()
