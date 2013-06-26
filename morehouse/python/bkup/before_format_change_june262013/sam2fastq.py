#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# convert bwa sam file to fastq

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

file_path = "/home/guoxing/disk2/solid/common_files/"
currentPath = os.getcwd() + '/'

# start time
start = time.time()		
		
# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input sam File Name", default="null")
(options, args) = parser.parse_args()

sam_file = options.samFile

if sam_file == "null":
	print "Please input the sam file name"

sam_file_name = sam_file[:(len(sam_file)-4)]

print "sam file: ", sam_file_name

inputFile_sam = open(currentPath + sam_file, "r")

outputFile_fastq = open(currentPath + sam_file_name + ".fastq", "w")
total_reads_num = 0

for read in inputFile_sam:
	if not read.startswith("@"):
		total_reads_num += 1	
		elements_first = read.strip().split()
		try:
			ID = elements_first[0].strip()
			read_seq = elements_first[9].strip()
			qual_line = elements_first[10].strip()
			read = "@"+ID+"\n"+read_seq+"\n+\n"+qual_line
			print >> outputFile_fastq, read.strip()		
		except:
			#print "error in line: ", line
			pass

print "total_reads_num: ", total_reads_num

end = time.time()
run_time = str(format((end - start), "0.3f"))
print "run time is: " + run_time + "s"

inputFile_sam.close()
outputFile_fastq.close()








