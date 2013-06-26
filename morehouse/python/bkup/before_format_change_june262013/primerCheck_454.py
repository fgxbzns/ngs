#!/usr/bin/python

# location /home/guoxing/tool/morehouse/python

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--aFile", type="string", dest="aFile",help = "Input File Name", default="null")
parser.add_option("-l", "--pl", type="string", dest="primerLength",help = "Input primer length", default="8")


(options, args) = parser.parse_args()

a_file_name = options.aFile
primer_length = int(options.primerLength)

currentPath = os.getcwd() + '/'

primer_dict = {}
primer_length = 6

b_file = open(currentPath + a_file_name + "_" + str(primer_length),'w')

base_list = ['A', 'T', 'C', 'G', 'N']

#primer_dict[primer_seq] = 0
#print "dict size: ", len(primer_dict)	
#print primer_dict


while primer_length <= 30:
	a_file = open(currentPath + a_file_name,'r')
	total_reads_number = 0

	line = a_file.readline()

	while line != "":
		if line.startswith(">"):
			total_reads_number += 1
			read_seq = ""
			line = a_file.readline().strip()			
			while not line.startswith(">") and line != "":
				read_seq += line
				line = a_file.readline().strip()
			#print read_seq
			primer_seq_begining = read_seq[:primer_length]		
			#print primer_seq_begining
			if primer_seq_begining not in primer_dict:
				primer_dict[primer_seq_begining] = 1
			else:
				primer_dict[primer_seq_begining] += 1
	
			primer_seq_end = read_seq[(len(read_seq)-primer_length):]	
			#print primer_seq_end
			if primer_seq_end not in primer_dict:
				primer_dict[primer_seq_end] = 1
			else:
				primer_dict[primer_seq_end] += 1
		#line = a_file.readline().strip()

	print "primer length: ", primer_length
	print "total_reads_number", total_reads_number

	print "dict size: ", len(primer_dict)
	b_file.write("primer length: " + str(primer_length) + "\n")	
	b_file.write("total reads number is: " + str(total_reads_number) + "\n")
	b_file.write("total primer number is: " + str(len(primer_dict)) + "\n")

	primer_list = [x for x in primer_dict.iteritems()] 
	primer_list.sort(key=lambda x: x[1], reverse=True) # sort by value in reverse order. Max first

	i = 0
	while i < 3:
		print primer_list[i], float(primer_list[i][1])/float(total_reads_number)
		b_file.write(str(primer_list[i]) + "\t" + str(float(primer_list[i][1])/float(total_reads_number)) + "\n")
		i += 1
	primer_length += 1
	primer_dict.clear()	# empty the dict after each round
	a_file.close()



b_file.close()
