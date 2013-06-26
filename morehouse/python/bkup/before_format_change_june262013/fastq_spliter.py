#!/usr/bin/python

# run bwa

import os, glob, subprocess, random, operator, time, math
from optparse import OptionParser

program_path = "/home/guoxing/disk2/ngs/morehouse/python/"
other_path = "/home/guoxing/disk2/ngs/morehouse/other/"
samtools_path = other_path + "samtools-0.1.18/"
ref_path = "/home/guoxing/disk2/UCSC_hg18_index_lm/"
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 

parser.add_option("-i", "--aFile", type="string", dest="aFile",help = "Input file Name", default="null")

(options, args) = parser.parse_args()
fastq_file = options.aFile
fastq_file_name = fastq_file[:fastq_file.find('.')].strip()

# start time
start = time.time()	

def wccount(filename):
    out = subprocess.Popen(['wc', '-l', filename],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

print "fastq_file is :", fastq_file

total_line_number = int(wccount(fastq_file))
print "total_line_number: ", total_line_number

number_of_subfile = 10
total_line_number_ceil = int(math.ceil(float(total_line_number)/100)*100)
print "total_line_number_ceil: ", total_line_number_ceil
# first n-1 files will contain line_in_each_subfile
line_in_each_subfile = total_line_number_ceil/4/number_of_subfile*4
print "line_in_each_subfile: ", line_in_each_subfile
# last file contains
line_in_last_subfile = int(math.fmod(total_line_number, line_in_each_subfile))
print "line_in_last_subfile: ", line_in_last_subfile

fastq_input_file = open(currentPath + fastq_file, "r")
line = fastq_input_file.readline()

file_number = 0
while file_number < number_of_subfile:
	line_number = 0
	output_subfile_name = fastq_file_name + "_" + str(file_number) + ".fastq"
	print "output: ", output_subfile_name
	output_subfile = open(currentPath + output_subfile_name, "w")
	while line_number < line_in_each_subfile and line != "":
		output_subfile.write(line)
		line = fastq_input_file.readline()
		line_number += 1
	file_number += 1
	output_subfile.close()
		
fastq_input_file.close()	


end = time.time()
run_time = str(format((end - start), "0.3f"))
print "run time is: " + run_time + "s"
