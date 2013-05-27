#!/usr/bin/python

# location /home/guoxing/tool/morehouse/python

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--aFile", type="string", dest="aFile",help = "Input File Name", default="null")

(options, args) = parser.parse_args()

a_file_name = options.aFile


#file_path = "/home/guoxing/morehouse_new_hap/fils/"
currentPath = os.getcwd() + '/'

a_file = open(currentPath + a_file_name,'r')
b_file = open(currentPath + a_file_name + "_nonLetter",'w')

#base_dict = {'A':'', 'T':'', 'C':'', 'G':'', 'N':'', '*':''}

base_dict = {}

total_base_number = 0
total_non_base_number = 0
total_character_number = 0

for line in a_file:
	if not line.startswith(">"):
		line = line.strip()
		for base in line:
			if base not in base_dict:
				base_dict[base] = ""
				print base
				total_non_base_number += 1
				b_file.write(base + "\n")
				b_file.write(line + "\n")
			total_character_number += 1	

"""
for line in a_file:
	if not line.startswith(">"):
		line = line.strip()
		for base in line:
			if base.upper() in base_dict:
				total_base_number += 1
			else:
				total_non_base_number += 1
				b_file.write(line + "\n")
				break
			total_character_number += 1
"""	
print "dict size: ", len(base_dict)	
print "total_base_number", total_base_number
print "total_non_base_number", total_non_base_number
print "total_character_number", total_character_number

a_file.close()
b_file.close()
