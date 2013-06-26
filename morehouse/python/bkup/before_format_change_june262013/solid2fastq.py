#!/usr/bin/python

# location /home/guoxing/tool/morehouse/python

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-c", "--cFile", type="string", dest="csfastaFile",help = "Input csfasta Name", default="null")
parser.add_option("-q", "--qFile", type="string", dest="qualFile",help = "Input qual Name", default="null")
parser.add_option("-o", "--output", type="string", dest="fastqFile",help = "output fastq file ", default="null")


(options, args) = parser.parse_args()

csfasta_file_name = options.csfastaFile
qual_file_name = options.qualFile
fastq_file_name = options.fastqFile

currentPath = os.getcwd() + '/'


code_list = [[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]]
base_list = ['A', 'C', 'G', 'T']

decode_dict = {}

for i in range(0, 4):
	base_dict = {}
	for j in range(0, 4):
		base_dict[base_list[j]] = base_list[code_list[i][j]]
	decode_dict[i]=base_dict

#print decode_dict[2]['T']

def convert_reads(solid_seq):
	fastaq_seq = ""
	anchor = ""
	solid_code_list = list(solid_seq.strip())
	anchor = solid_code_list.pop(0)		# assign first anchor tag and remove it from the list
	for solid_code in solid_code_list:
		if solid_code == ".":
			fastaq_seq += "N"
			anchor = "N"
		elif anchor == "N":
			fastaq_seq += "N"
		else:
			base_code = decode_dict[int(solid_code)][anchor]
			fastaq_seq += base_code
			anchor = base_code
	return fastaq_seq

#print convert_reads("T20130311211003031210002211022222212121122210023222")

def convert_qual(solid_qual):
	fastaq_qual = ""
	solid_qual_list = solid_qual.strip().split()
	for solid_qual_code in solid_qual_list:
		if solid_qual_code == "-1":
			fastaq_qual += chr(0+33)
		else:
			try:
				solid_qual_code = int(solid_qual_code)
			except ValueError:
				print "error in ", solid_qual
				solid_qual_code = 0	
			fastaq_qual += chr(solid_qual_code+33)
	return fastaq_qual

#print convert_qual("6 6 32 31 8 12 4 4 20 17 12 5 4 4 4 6 10 10 4 8 6 12 14 8 22 19 4 11 15 30 5 4 7 12 6 14 4 9 -1 22 22 19 5 8 22 20 4 5 5 11")
"""
csfasta_file_name = "song_1_5_run1_bcSample1_F3_song_5.csfasta"
qual_file_name = "song_1_5_run1_bcSample1_F3_QV_song_5.qual"
#csfasta_file_name = "song_1_5_run1_bcSample1_F3_song_5_500.txt"
#qual_file_name = "song_1_5_run1_bcSample1_F3_QV_song_5_500.txt"
fastq_file_name = "song_5_2_withN"
"""
csfasta_file = open(currentPath + csfasta_file_name,'r')
qual_file = open(currentPath + qual_file_name,'r')
fastq_file = open(currentPath + fastq_file_name + ".fastq",'w')

csfasta_line = ""
while not csfasta_line.startswith(">"):
	csfasta_line = csfasta_file.readline().strip()
print csfasta_line

qual_line = ""
while not qual_line.startswith(">"):
	qual_line = qual_file.readline().strip()
print qual_line	

total_reads_number = 0
reads_with_N_number = 0
	
while csfasta_line != "":
	if csfasta_line.startswith(">") and qual_line.startswith(">"):
		if csfasta_line == qual_line:
			total_reads_number += 1
			reads_complete = True
			title = "@" + fastq_file_name + ":" + csfasta_line + "/1"
			csfasta_line = csfasta_file.readline().strip()
			fastaq_seq = convert_reads(csfasta_line)
			qual_line = qual_file.readline().strip()
			fastaq_qual = convert_qual(qual_line)
			"""
			if "N" in fastaq_seq:
				reads_with_N_number +=1
				reads_complete = False
				index_of_first_N = fastaq_seq.find('N')
				fastaq_seq = fastaq_seq[:index_of_first_N]
				fastaq_qual = fastaq_qual[:index_of_first_N]
			"""
			if reads_complete:
				fastq_file.write(title + "\n")
				fastq_file.write(fastaq_seq + "\n")
				fastq_file.write("+" + "\n")
				fastq_file.write(fastaq_qual + "\n")
		else:
			#csfasta_line = csfasta_file.readline().strip()
			#qual_line = qual_file.readline().strip()
			print "read name does not match"
			print "csfasta_line", csfasta_line
			print "qual_line", qual_line
	csfasta_line = csfasta_file.readline().strip()
	qual_line = qual_file.readline().strip()


print "total_reads_number", total_reads_number
print "reads_with_N_number", reads_with_N_number


csfasta_file.close()
qual_file.close()
fastq_file.close()
