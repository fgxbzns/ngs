#!/usr/bin/python
#######################################################################################
# wli 12878 data
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

def primer_remove_single_end():
	output_file = open(input_file_name[:input_file_name.find('.')] + "_priRem.fastq", 'w')

	data_tag = "@ILLUMINA"          #wli data
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
				#if read_seq != "" and qual_line != "" and len(read_seq) >= 25:
				#if read_seq != "" and qual_line != "":      # for wli
				if True:
					output_file.write(title + "\n")
					output_file.write(read_seq + "\n")
					output_file.write("+" + "\n")
					output_file.write(qual_line + "\n")
			line = input_file.readline().strip()
	
	print "total_reads_number", total_reads_number
	output_file.close()

def primer_remove_pair_end():
	output_file = open(input_file_name[:input_file_name.find('.')] + "_priRem.fastq", 'w')

	data_tag = "@ILLUMINA"          #wli data
	total_reads_number = 0

	with open(input_file_name,'r') as input_file:
		line = input_file.readline()
		while line != "":
			if line.startswith(data_tag):
				total_reads_number += 1
				first_title = line.strip()
				first_read_seq = input_file.readline().strip()
				plus_symbol = input_file.readline().strip()
				first_qual_line = input_file.readline().strip()

				first_read_seq, first_qual_line = process_seq(first_read_seq, first_qual_line)
				if first_read_seq != "" and first_qual_line != "" and len(first_read_seq) >= 25:
					second_line = input_file.readline().strip()


					output_file.write(first_title + "\n")
					output_file.write(first_read_seq + "\n")
					output_file.write("+" + "\n")
					output_file.write(first_qual_line + "\n")
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

def remove_N(read_seq, qual_line):
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


def sequence_remove_single_end():
	# for lima's fastq file to remove certain number of bases from fastq reads
	# single chromosome
	# pair end illumina GAII sample 12878

	# also used for mimi yang data

	removed_base_number = 35    # for pairend read 2
	#removed_base_number = 28    # for pairend read 2
	output_file = open(input_file_name[:input_file_name.find('.')] + "_seqRem.fastq", 'w')
	#data_tag = "@ILLUMINA"
	data_tag = "@SRR"
	total_reads_number = 0
	print "primer remove:", input_file_name

	with open(input_file_name, 'r') as input_file:
		line = input_file.readline()
		while line != "":
			if line.startswith(data_tag):
				total_reads_number += 1
				title = line.strip()
				read_seq = input_file.readline().strip()
				plus_symbol = input_file.readline().strip()
				qual_line = input_file.readline().strip()

				output_file.write(title + "\n")

				read_seq = read_seq[removed_base_number:]
				qual_line = qual_line[removed_base_number:]
				read_seq, qual_line = remove_N(read_seq, qual_line)
				output_file.write(read_seq + "\n")
				#output_file.write(read_seq[35:] + "\n")
				#output_file.write(read_seq[::-1] + "\n")
				output_file.write("+" + "\n")
				output_file.write(qual_line + "\n")
				#output_file.write(qual_line[35:] + "\n")
				#output_file.write(qual_line[::-1] + "\n")
			line = input_file.readline().strip()

	print "total_reads_number", total_reads_number
	output_file.close()

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
	#primer_remove_single_end()
	sequence_remove_single_end()
	print "run time is: ", round((time.time() - start_time), 3), "s"

