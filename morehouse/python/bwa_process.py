#!/usr/bin/python

# run bwa

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

program_path = "/home/guoxing/disk2/ngs/morehouse/python/"
other_path = "/home/guoxing/disk2/ngs/morehouse/other/"
samtools_path = other_path + "samtools-0.1.18/"
ref_path = "/home/guoxing/disk2/UCSC_hg18_index_lm/"
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 

parser.add_option("-i", "--inputfile", type="string", dest="inputfile",help = "Input File Name", default="null")
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="all")
(options, args) = parser.parse_args()

# start time
start = time.time()	

fastq_file = options.inputfile
fastq_file = "song_5_prem_0.fastq"
fastq_file_name = fastq_file[:(fastq_file.find('fastq')-2)].strip()
chr_name = options.chrName
chr_name = "chr11"

number_of_subfile = 10
"""
# bwa
file_number = 5
while file_number < number_of_subfile:
	subfile_name = fastq_file_name + str(file_number) + ".fastq"
	print "current file: ", subfile_name
	bwa = program_path + "bwaRun.py -i " + subfile_name + " -c " + chr_name + " &"
	#print bwa
	#os.system(bwa)
	#subprocess.call(bwa)
	bwa_process = subprocess.Popen(bwa, shell=True)
	bwa_process.wait()
	file_number += 1

#os.wait()
"""
# grep chr
file_number = 0
while file_number < number_of_subfile:
	subfile_sam_name = fastq_file_name + str(file_number) + "_all"  + ".sam"
	print "current file: ", subfile_sam_name
	if file_number == 0:
		#grep = "grep chr " + " " + subfile_sam_name + " > " + fastq_file_name + chr_name + ".sam &"
		grep = "grep " + chr_name + " " + subfile_sam_name + " > " + fastq_file_name + chr_name + ".sam &"
	else:
		#grep = "grep chr " + " " + subfile_sam_name + " >> " + fastq_file_name + chr_name + ".sam &"
		grep = "grep " + chr_name + " " + subfile_sam_name + " >> " + fastq_file_name + chr_name + ".sam &"
	print grep
	#os.system(grep)
	#os.wait()
	grep_Process = subprocess.Popen(grep, shell=True)
	grep_Process.wait()
	grep_Process.poll()
	file_number += 1
	
#os.wait()

end = time.time()
run_time = str(format((end - start), "0.3f"))
print "run time is: " + run_time + "s"

data_record_file_name = fastq_file_name + "bwaTotalTime.txt"
data_record_file = open(currentPath + data_record_file_name, "a")
print >> data_record_file, "run time is: " + run_time + "s"
data_record_file.close()
