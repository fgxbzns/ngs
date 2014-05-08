#!/usr/bin/python

# location /home/guoxing/tool/morehouse/python

import os
from optparse import OptionParser

currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--aFile", type="string", dest="aFile", help="Input File Name", default="null")

(options, args) = parser.parse_args()

a_file_name = options.aFile



chr_dict = {}

a_file = open(currentPath + a_file_name, 'r')
b_file = open(currentPath + a_file_name + "_percentage", 'w')

total_reads_number = 0

line = a_file.readline()

for line in a_file:
	if not line.startswith("@"):
		total_reads_number += 1
		elements = line.strip().split()
		chr_name = elements[2].strip()
		if chr_name not in chr_dict:
			chr_dict[chr_name] = 1
		else:
			chr_dict[chr_name] += 1

chr_list = [x for x in chr_dict.iteritems()] 
chr_list.sort(key=lambda x: x[1], reverse=True)  # sort by value in reverse order. Max first

b_file.write("total_reads_number" + "\t" + str(total_reads_number) + "\n")

for chr in chr_list:
	b_file.write(chr[0] + "\t" + str(chr[1]) + "\t" + str(float(chr[1])*100/float(total_reads_number))+ "\n")

"""
for name, number in chr_dict.iteritems():
	b_file.write(name + "\t" + str(number) + "\t" + str(float(number)*100/float(total_reads_number))+ "\n")
"""
print "total_reads_number", total_reads_number

a_file.close()
b_file.close()
