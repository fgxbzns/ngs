#!/usr/bin/python
#######################################################################################
# wli 12878 data
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

def primer_remove():
	output_file = open(input_file_name[:input_file_name.find('.')] + "_priRem.fastq",'w')

	data_tag = "@ILLUMINA"
	#data_tag = "@song"
	total_reads_number = 0
	
	with open(input_file_name,'r') as input_file:
		line = input_file.readline()
		while line != "":
			if line.startswith(data_tag):
				total_reads_number += 1
				title = line.strip()
				read_seq = input_file.readline().strip()
				plus_symbol = input_file.readline().strip()
				qual_line = input_file.readline().strip()
				
				read_seq, qual_line = process_seq(read_seq, qual_line)
				if read_seq != "" and qual_line != "" and len(read_seq) >= 25:
					output_file.write(title + "\n")
					output_file.write(read_seq + "\n")
					output_file.write("+" + "\n")
					output_file.write(qual_line + "\n")
			line = input_file.readline().strip()
	
	print "total_reads_number", total_reads_number
	output_file.close()

def process_seq(read_seq, qual_line):
	for primer_seq in primer_list:
		while primer_seq in read_seq:
			seq_length = len(read_seq)
			primer_position = read_seq.find(primer_seq)
			if primer_position <= (seq_length/2):
				new_start_position = primer_position + len(primer_seq) + 8
				read_seq = read_seq[new_start_position:]
				qual_line = qual_line[new_start_position:]
			else:
				read_seq = read_seq[:primer_position]
				qual_line = qual_line[:primer_position]
	
	while 'N' in read_seq:
		seq_length = len(read_seq)
		N_position = read_seq.find('N')
		if N_position <= (seq_length/2):
			read_seq = read_seq[(N_position + 1):]
			qual_line = qual_line[(N_position + 1):]
		else:
			read_seq = read_seq[:N_position]
			qual_line = qual_line[:N_position]
	return (read_seq, qual_line)

def get_args():
	desc = "remove primer"
	usage = "primerRemove -i input_fastq" 
	parser = OptionParser(usage=usage, description=desc) 
	parser.add_option("-i", "--inputFile", type="string", dest="inputFile",help = "Input File Name", default="null")
	parser.add_option("-p", "--ps", type="string", dest="primerSequence",help = "Input primer sequence", default="null")
	(options, args) = parser.parse_args()
	if options.inputFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	global input_file_name
	global primer_list
	
	primer_list = ["TGTGTTGGGTGTGTTTGG", "TGTNTTGGGTGTGTTTGG", "TGTNTTGGGGTGTTTGG", "TGTTGGGTGTGTTTGG"]
	#primer_list = ["TGTGTTGGGTGTGTTTGG", "TGTGTTGGGTGTGTTTGG", "CGCCTTGGCCGTACAGCA"]
		
	options = get_args()
	input_file_name = options.inputFile
	
	start_time = time.time()
	primer_remove()
	elapse_time = time.time() - start_time
	print "run time is: ", round((time.time() - start_time), 3), "s"

