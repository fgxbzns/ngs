#!/usr/bin/python

# location /home/guoxing/tool/morehouse/python

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--gFile", type="string", dest="gFile",help = "Input File Name", default="null")
parser.add_option("-n", "--ps", type="string", dest="name",help = "Input name", default="null")

(options, args) = parser.parse_args()

a_file_name = options.aFile
name = options.name

currentPath = os.getcwd() + '/'


genotype_input_file_name = a_file_name

print genotype_input_file_name
genotype_input_file = open(file_path + genotype_input_file_name, "r")
#genotype_input_file = open(currentPath + genotype_input_file_name, "r")

genotype_output_file_name = "genotype.txt"
genotype_output_file = open(currentPath + genotype_output_file_name, "w")
#genotype_output_file.write("rsID \t phys_position \t NA12878 \n")

title_genotype = ""

genotype_dict = {}

for line in genotype_input_file:
	if line.startswith("rsID"):
		title_genotype = line.strip()
	if not line.startswith("rsID"):
		elements = line.strip().split()
		position = elements[1].strip()							
		try:
			position = int(position)
			genotype_dict[position] = line.strip()
		except ValueError:
			print "error in ", line

genotype_output_file.write(title_genotype + "\n")
