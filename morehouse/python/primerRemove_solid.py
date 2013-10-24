#!/usr/bin/python

# location /home/guoxing/tool/morehouse/python

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--aFile", type="string", dest="aFile",help = "Input File Name", default="null")
parser.add_option("-p", "--ps", type="string", dest="primerSequence",help = "Input primer sequence", default="null")

(options, args) = parser.parse_args()

a_file_name = options.aFile
#primer_length = int(options.primerLength)

currentPath = os.getcwd() + '/'

primer_dict = {}

a_file = open(currentPath + a_file_name,'r')
b_file = open(currentPath + a_file_name[:a_file_name.find('.')] + "_prem.fastq",'w')
#c_file = open(currentPath + a_file_name + "_noPri",'w')

primer_seq_1 = "TGTGTTGGGTGTGTTTGG"	# keep
primer_length_1 = len(primer_seq_1)
primer_seq_2 = "CGCCTTGGCCGTACAGCA"	# remove
primer_length_2 = len(primer_seq_2)

total_reads_number = 0
reads_with_N_number = 0
reads_with_correct_primer_number = 0
reads_with_wrong_primer_number = 0
reads_without_primer_number = 0

line = a_file.readline()

while line != "":
	if line.startswith("@song"):
		keep_this_line = True
		total_reads_number += 1
		title = line
		read_seq = a_file.readline().strip()
		plus_symbol = a_file.readline().strip()
		qual_line = a_file.readline().strip()
		if "N" in read_seq:
			keep_this_line = False
			reads_with_N_number += 1
		elif read_seq.startswith(primer_seq_2):
			keep_this_line = False
			reads_with_wrong_primer_number += 1
		elif read_seq.startswith(primer_seq_1):
			read_seq = read_seq[len(primer_seq_1):]
			qual_line = qual_line[len(primer_seq_1):]
			reads_with_correct_primer_number += 1
			b_file.write(title + "\n")
			b_file.write(read_seq + "\n")
			b_file.write("+" + "\n")
			b_file.write(qual_line + "\n")
		else:
			reads_without_primer_number += 1
			b_file.write(title + "\n")
			b_file.write(read_seq + "\n")
			b_file.write("+" + "\n")
			b_file.write(qual_line + "\n")
	line = a_file.readline().strip()


print "total_reads_number", total_reads_number
print "reads_with_N_number", reads_with_N_number
print "reads_with_correct_primer_number", reads_with_correct_primer_number
print "reads_with_wrong_primer_number", reads_with_wrong_primer_number
print "reads_without_primer_number", reads_without_primer_number

a_file.close()
b_file.close()
#c_file.close()
