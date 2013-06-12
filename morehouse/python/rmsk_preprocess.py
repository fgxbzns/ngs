#!/usr/bin/python

# remove entries in rmsk reference other than "Alu" and "Simple_repeat"

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

file_path = "/home/guoxing/disk2/solid/common_files/"
currentPath = os.getcwd() + '/'

# start time
start = time.time()		
		
# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--rmsk", type="string", dest="rmskFile",help = "Input rmsk File Name", default="null")
(options, args) = parser.parse_args()

rmsk_file = options.rmskFile

if rmsk_file == "null":
	print "Please input the rmsk file name"

rmsk_file_name = rmsk_file[:(len(rmsk_file)-4)]
print "rmsk file: ", rmsk_file_name

inputFile_rmsk = open(currentPath + rmsk_file, "r")
outputFile_rmsk = open(currentPath + rmsk_file_name + "_processed.txt", "w")

ori_rmsk_total_number = 0
processed_rmsk_total_number = 0

for line in inputFile_rmsk:
	ori_rmsk_total_number += 1	
	if "Alu" in line or "Simple_repeat" in line:
		processed_rmsk_total_number += 1	
		print >> outputFile_rmsk, line.strip()	
		
	

print "ori_rmsk_total_number: ", ori_rmsk_total_number
print "processed_rmsk_total_number: ", processed_rmsk_total_number

end = time.time()
run_time = str(format((end - start), "0.3f"))
print "run time is: " + run_time + "s"

inputFile_rmsk.close()
outputFile_rmsk.close()








