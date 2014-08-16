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
parser = OptionParser(usage=usage)

parser.add_option("-i", "--aFile", type="string", dest="aFile", help="Input file Name", default="null")

(options, args) = parser.parse_args()
chr_file = options.aFile

# chr_file = "hg18chr.fa"
chr_file_name = chr_file[:chr_file.find('.')].strip()

print chr_file
#chr_input_file = open(ref_path + chr_file, "r")
chr_input_file = open(currentPath + chr_file, "r")

output_in_process = False

for line in chr_input_file:
	if line.startswith(">") and output_in_process:
		chr_output_file.close()
		output_in_process = False
	if line.startswith(">") and not output_in_process:
		output_in_process = True
		chr_name = line[(line.find('>') + 1):].strip()
		print chr_name
		chr_output_file = open(currentPath + chr_file_name + "_" + chr_name + ".fa", "w")
		chr_output_file.write(line.strip() + "\n")
	else:
		chr_output_file.write(line.strip() + "\n")

chr_input_file.close()
