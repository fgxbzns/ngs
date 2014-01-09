#!/usr/bin/python

# location /home/guoxing/tool/morehouse/python

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--aFile", type="string", dest="aFile",help = "Input File Name", default="null")
parser.add_option("-l", "--pl", type="string", dest="primerLength",help = "Input primer length", default="16")
parser.add_option("-u", "--up", type="string", dest="upbound",help = "Input upper bound of primer length", default="22")

(options, args) = parser.parse_args()

a_file_name = options.aFile
primer_length = int(options.primerLength)

currentPath = os.getcwd() + '/'

primer_dict = {}
primer_length = 16

b_file = open(currentPath + a_file_name + "_" + str(primer_length),'w')

base_list = ['A', 'T', 'C', 'G', 'N']

def generate_permutations(chars = 6) :
     
	allowed_chars= ['A', 'T', 'C', 'G', 'N']
   
	status = []
	for tmp in range(chars) :
		status.append(0)
		last_char = len(allowed_chars)
		rows = []
	for x in xrange(last_char ** chars) :
		rows.append("")
		for y in range(chars - 1 , -1, -1) :
			key = status[y]
			rows[x] = allowed_chars[key] + rows[x]
				
		for pos in range(chars - 1, -1, -1) :
			if(status[pos] == last_char - 1) :
				status[pos] = 0
			else :
				status[pos] += 1
				break;
           
	return rows
  
#print len(generate_permutations())


#primer_dict[primer_seq] = 0
#print "dict size: ", len(primer_dict)	
#print primer_dict


while primer_length <= 24:
	a_file = open(currentPath + a_file_name,'r')
	total_reads_number = 0

	line = a_file.readline()

	while line != "":
		if line.startswith("@ILLUMINA"):	# for quake data
			total_reads_number += 1
			read_seq = a_file.readline().strip()
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
		line = a_file.readline().strip()

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
