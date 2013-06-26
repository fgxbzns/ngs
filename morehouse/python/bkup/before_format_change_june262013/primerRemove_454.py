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

read_seq = "TGTGTTGGGTGTGTTTGGTGTGTTGGGTGTGTTTGGGTGTTTTTTGGGTGTGTTGGGTGTGTTTGGAATGCGCTCGTCATTCATGAGCTTTC"
qual_line = "TGTGTTGGGTGTGTTTGGTGTGTTGGGTGTGTTTGGGTGTTTTTTGGGTGTGTTGGGTGTGTTTGGAATGCGCTCGTCATTCATGAGCTTTC"
primer_seq_1 = "TGTGTTGGGTGTGTTTGG"

read_seq = "TTCCAAATTACTCCTTAGTACCAAACCCAACCAAACACACCCAACACAACAACACCCCAAACACACCCAACACAACAA"
qual_line = "TTCCAAATTACTCCTTAGTACCAAACCCAACCAAACACACCCAACACAACAACACCCCAAACACACCCAACACAACAA"
primer_seq_1 = "CCAAACACACCCAACACA"

"""
def removePrimer(read_seq, qual_line, primer_seq):
	if primer_seq in read_seq:
		new_start_position = read_seq.find(primer_seq) + len(primer_seq)
		read_seq = read_seq[new_start_position:]
		qual_line = qual_line[new_start_position:]
		#print read_seq
		#print qual_line
		removePrimer(read_seq, qual_line, primer_seq)
	else:
		print read_seq
		return (read_seq, qual_line)
"""
def removePrimer(read_seq, qual_line, primer_seq):
	while primer_seq in read_seq:
		new_start_position = read_seq.find(primer_seq) + len(primer_seq)
		read_seq = read_seq[new_start_position:]
		qual_line = qual_line[new_start_position:]
	return (read_seq, qual_line)
"""	
a = removePrimer(read_seq, qual_line, primer_seq_1)
print a[0]
print a [1]
print a

print read_seq
read_seq = read_seq[::-1]
print read_seq
qual_line = qual_line[::-1]
primer_seq_1 = primer_seq_1[::-1]
"""

a = removePrimer(read_seq[::-1], qual_line[::-1], primer_seq_1[::-1])
print a[0][::-1]
print a [1][::-1]
#print a


a_file = open(currentPath + a_file_name,'r')
b_file = open(currentPath + a_file_name[:a_file_name.find('.')] + "_prem.fastq",'w')
c_file = open(currentPath + a_file_name[:a_file_name.find('.')] + "_record.txt",'w')

primer_seq_1 = "TGTGTTGGGTGTGTTTGG"	# keep  in the begining
primer_length_1 = len(primer_seq_1)
primer_seq_2 = "CCAAACACACCCAACACA"	# keep 	this one is in the end
primer_length_2 = len(primer_seq_2)

total_reads_number = 0
reads_with_N_number = 0
reads_with_primer_1 = 0
reads_with_primer_2 = 0
reads_without_primer_number = 0

line = a_file.readline()

while line != "":
	if line.startswith("@>"):
		keep_this_line = True
		total_reads_number += 1
		title = line
		read_seq = a_file.readline().strip()
		plus_symbol = a_file.readline().strip()
		qual_line = a_file.readline().strip()
		if "N" in read_seq:
			keep_this_line = False
			reads_with_N_number += 1
		else:
			if primer_seq_1 in read_seq:
				seq_qual = removePrimer(read_seq, qual_line, primer_seq_1)
				read_seq = seq_qual[0]
				qual_line = seq_qual[1]
				reads_with_primer_1 += 1
			if primer_seq_2 in read_seq:	# reverse the seq, remove primer, then reverse it back 
				seq_qual = removePrimer(read_seq[::-1], qual_line[::-1], primer_seq_2[::-1])
				read_seq = seq_qual[0][::-1]
				qual_line = seq_qual[1][::-1]
				reads_with_primer_2 += 1
			else:
				reads_without_primer_number += 1
			b_file.write(title + "\n")
			b_file.write(read_seq + "\n")
			b_file.write("+" + "\n")
			b_file.write(qual_line + "\n")
	line = a_file.readline().strip()


print "total_reads_number", total_reads_number
print "reads_with_N_number", reads_with_N_number
print "reads_with_primer_1", reads_with_primer_1
print "reads_with_primer_2", reads_with_primer_2
print "reads_without_primer_number", reads_without_primer_number

print >> c_file, "total_reads_number", total_reads_number
print >> c_file, "reads_with_N_number", reads_with_N_number
print >> c_file, "reads_with_primer_1", reads_with_primer_1
print >> c_file, "reads_with_primer_2", reads_with_primer_2
print >> c_file, "reads_without_primer_number", reads_without_primer_number

a_file.close()
b_file.close()
c_file.close()

