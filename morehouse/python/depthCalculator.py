#!/usr/bin/python

# location /home/guoxing/tool/morehouse/python

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

currentPath = os.getcwd() + '/'
file_path = "/home/guoxing/disk2/solid/common_files/"

# Reading options
usage = "usage: %prog [options] arg1"
parser = OptionParser(usage=usage)
parser.add_option("-s", "--aFile", type="string", dest="samFile", help="Input File Name", default="null")
parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="chr11")

chr_length_dict = {"chr1": 247249719, "chr2": 242951149, "chr3": 199501827, "chr4": 191273063, "chr5": 180857866,
                   "chr6": 170899992, "chr7": 158821424, "chr8": 146274826, "chr9": 140273252, "chr10": 135374737,
                   "chr11": 134452384, "chr12": 132349534, "chr13": 114142980, "chr14": 106368585, "chr15": 100338915,
                   "chr16": 88827254, "chr17": 78774742, "chr18": 76117153, "chr19": 63811651, "chr20": 62435964,
                   "chr21": 46944323, "chr22": 49691432, "chrX": 154913754, "chrY": 57772988}

(options, args) = parser.parse_args()

sam_file = options.samFile
chr_name = options.chrName

total_reads_number = 0
total_base_number = 0

# a_file = open(currentPath + a_file_name,'r')
inputFile_sam = open(currentPath + sam_file, "r")
b_file = open(file_path + "depth_calculator.txt", 'a')

# for mapped, unique, matched chr depth
sam_line_first = inputFile_sam.readline()

while sam_line_first != '':
	if not sam_line_first.startswith("@"):
		total_reads_number += 1
		elements_first = sam_line_first.strip().split()
		try:
			rName_first = elements_first[2].strip()
		except:
			print "no chr name in read:", sam_line_first
		if (rName_first == chr_name):
			read_sequence_first = elements_first[9].strip()
			read_length_first = len(read_sequence_first)
			total_base_number += read_length_first
	sam_line_first = inputFile_sam.readline()

inputFile_sam.close()

# for raw data depth
"""
line = a_file.readline()
while line != "":
	#if line.startswith("@>"):	# 454 data
	if line.startswith("@song"):	# solid
		total_reads_number += 1
		read_seq = a_file.readline().strip()
		total_base_number += len(read_seq)	
	line = a_file.readline().strip()
"""

print "chr_name: ", sam_file
print "chr_name: ", chr_name
print "chr_base_number: ", chr_length_dict[chr_name]
print "total_reads_number: ", total_reads_number
print "total_base_number: ", total_base_number
depth = float(total_base_number) / float(chr_length_dict[chr_name])
print "depth: ", format(depth, "0.4f")

print >> b_file, "file name: ", sam_file
print >> b_file, "chr_name: ", chr_name
print >> b_file, "chr_base_number: ", chr_length_dict[chr_name]
print >> b_file, "total_reads_number: ", total_reads_number
print >> b_file, "total_base_number: ", total_base_number
print >> b_file, "depth: ", format(depth, "0.4f")
print >> b_file, ""

b_file.close()

