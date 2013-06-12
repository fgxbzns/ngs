#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for pre-processing the sam file before snpPick

import os, glob, subprocess, random, operator, time
from optparse import OptionParser


file_path = "/home/guoxing/disk2/solid/common_files/"
currentPath = os.getcwd() + '/'


# start time
start = time.time()		
		
# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input sam File Name", default="null")

(options, args) = parser.parse_args()

sam_file = options.samFile


indel_info = "39M62S"
read_seq = "GTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAGAGGCCCTTCCCCCCCCCAAAACAAAAAAAAAAAAATGCTTTCCCCGGGGGAAAAAAAAACATCC"

qual_line = "C@BFFFFFHHHHHJJJJJJJHFDDDDDDDDDD#####################################################################"


# keep mismatched seq >= 25

mismatch_threshold = 25

def removePrimer(read_seq, qual_line, indel_info):
	indel_length = ""
	indel_type = ""
	read_seq_final = ""
	qual_line_final = ""
	for info in indel_info:
		try:
			int(info)
			indel_length += info
		except:
			indel_type = info
			indel_length = int(indel_length)
			if indel_type == "M":	# matched, remove this part 
				read_seq = read_seq[indel_length:]
				qual_line = qual_line[indel_length:]
			elif indel_type == "S" or indel_type == "H":	# mismatch, keep this part
				if indel_length >= mismatch_threshold:
					read_seq_final += read_seq[:indel_length]
					qual_line_final += qual_line[:indel_length]
					read_seq = read_seq[indel_length:]
					qual_line = qual_line[indel_length:]
				else:
					read_seq = read_seq[indel_length:]
					qual_line = qual_line[indel_length:]
			elif indel_type == "I":	# insertion, remove
				read_seq = read_seq[indel_length:]
				qual_line = qual_line[indel_length:]
			elif indel_type == "D":	# deletion, insert "N" in read_seq, insert "!" (quality score = 0) in qual_line
				# do not keep deletion
				read_seq = read_seq[indel_length:]
				qual_line = qual_line[indel_length:]
				"""
				i = 0
				while i < indel_length:
					read_seq_final += "N"
					qual_line_final += "!"
					i += 1
				"""
			else:
				print "error in ", read_seq, " info type ", indel_type
			indel_length = ""
			indel_type = ""	
	return (read_seq_final, qual_line_final)


read_qual = removePrimer(read_seq, qual_line, indel_info)
print read_qual[0]
print read_qual[1]
read = ">G16U2IC01CIWF1 0       chr11   223477  25      39M62S *       0       0       \
		GTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAGAGGCCCTTCCCCCCCCCAAAACAAAAAAAAAAAAATGCTTTCCCCGGGGGAAAAAAAAACATCC   \
		C@BFFFFFHHHHHJJJJJJJHFDDDDDDDDDD#####################################################################   AS:i:59 XS:i:0  XF:i:3  XE:i:1  NM:i:5"
print "ori ", read

read = read.replace(read_seq, read_qual[0])
read = read.replace(qual_line, read_qual[1])
#read = read.replace("TCCAACAACACTTGAAGCACCCACAGAAAGATATGGACCTTCTTCCACCTTTAATACTGACAGAAAACACAGCAGCAAAACAAAA", read_qual[0])
print "new ", read




if sam_file == "null":
	print "Please input the sam file name"

sam_file_name = sam_file[:(len(sam_file)-4)]

print "sam file: ", sam_file_name

inputFile_sam = open(currentPath + sam_file, "r")

outputFile_sam = open(currentPath + sam_file_name + "_hanging.sam", "w")
total_reads_num = 0
hanging_reads_num = 0

for read in inputFile_sam:
	if not read.startswith("@"):
		total_reads_num += 1	
		elements_first = read.strip().split()
		try:
			indel_info = elements_first[5].strip()
			read_seq = elements_first[9].strip()
			qual_line = elements_first[10].strip()
			read_qual = removePrimer(read_seq, qual_line, indel_info)
			if len(read_qual[0]) > 0:
				hanging_reads_num += 1
				read = read.replace(read_seq, read_qual[0])
				read = read.replace(qual_line, read_qual[1])
				print >> outputFile_sam, read.strip()		
		except:
			#print "error in line: ", line
			pass

print "total_reads_num: ", total_reads_num
print "hanging_reads_num: ", hanging_reads_num

inputFile_sam.close()
outputFile_sam.close()







