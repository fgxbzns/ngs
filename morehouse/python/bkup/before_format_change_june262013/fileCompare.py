#!/usr/bin/python

# location /home/guoxing/tool/morehouse

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-a", "--aFile", type="string", dest="aFile",help = "Input File Name", default="null")
parser.add_option("-b", "--bFile", type="string", dest="bFile",help = "Input File Name", default="null")

(options, args) = parser.parse_args()

a_file_name = options.aFile
b_file_name = options.bFile


file_path = "/home/guoxing/morehouse_new_hap/fils/"
currentPath = os.getcwd() + '/'

a_file = open(currentPath + a_file_name,'r')
b_file = open(currentPath + b_file_name,'r')

a_line = a_file.readline()
b_line = b_file.readline()

same_line_number = 0
different_line_number = 0
total_element_number = 0

while a_line!='':
	if not a_line.startswith('rsID'):
		if a_line.strip() == b_line.strip():
			a_line = a_file.readline()
			b_line = b_file.readline()
			same_line_number += 1
		else:
			print a_line
			print b_line
			a_line = a_file.readline()
			b_line = b_file.readline()
			different_line_number += 1		
	else:
		print a_line
		print b_line
		a_line = a_file.readline()
		b_line = b_file.readline()
	total_element_number += 1
		
print "same_line_number", same_line_number
print "different_line_number", different_line_number
print "total_element_number", total_element_number
print "per", float(different_line_number)/float(total_element_number)

a_file.close()
b_file.close()
