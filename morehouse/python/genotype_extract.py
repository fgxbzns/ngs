#!/usr/bin/python

# location /home/guoxing/tool/morehouse/python

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

file_path = "/home/guoxing/disk2/solid/common_files/"
NA10847_F12146_M12239_path = file_path + "NA10847_F12146_M12239/"
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--gFile", type="string", dest="gFile",help = "Input File Name", default="null")
parser.add_option("-n", "--ps", type="string", dest="name",help = "Input name", default="null")
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")

(options, args) = parser.parse_args()

a_file_name = options.gFile
person_ID = options.name
#person_ID = "NA12146"
chr_name = options.chrName

genotype_input_file_name = a_file_name

print genotype_input_file_name
genotype_input_file = open(currentPath + genotype_input_file_name, "r")
#genotype_input_file = open(currentPath + genotype_input_file_name, "r")

genotype_output_file_name = "genotype_"+person_ID+"_"+chr_name+".txt"
genotype_output_file = open(NA10847_F12146_M12239_path + genotype_output_file_name, "w")
print >> genotype_output_file, "rsID" + "\t" + "phys_position" + "\t" + person_ID

ID_index = 0

genotype_dict = {}

for line in genotype_input_file:
	if line.startswith("rs#"):
		elements = line.strip().split()
		try:
			ID_index = elements.index(person_ID)
			print ID_index
		except:
			print person_ID, "not found in file"
			break
	else:
		elements = line.strip().split()
		try:
			rsID = elements[0].strip()							
			position = elements[3].strip()
			genotype = elements[ID_index].strip()
			print >> genotype_output_file, rsID + "\t" + position + "\t" + genotype
		except ValueError:
			print "error in ", line
genotype_output_file.close()
