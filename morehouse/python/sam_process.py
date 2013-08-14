#!/usr/bin/python

# this is for pre-processing the sam file before snpPick

import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser
from tools import *

def is_multiple_maping(elements):
	multiple_maping = False
	try:
		XA = elements_first[21].strip()
		multiple_maping_first = True
	except:
		pass
	return multiple_maping

def pair_end_filter(sam_file):
	"""filter the reads by pair_end info"""
	sam_file_name = sam_file[:(len(sam_file)-4)]
	print "pair_end_filter: ", sam_file_name
	output_file = open(sam_file_name + "_pairend.sam", "w")
	inputfile_sam = open(currentPath + sam_file, "r")
	sam_line_first = inputfile_sam.readline() # the first read line in a pair
	total_reads_num = 0
	covered_snp_total_number = 0
	
	insert_size_lower_bond = 100
	insert_size_upper_bond = 1000

	while sam_line_first!='':
		keep_this_pair = False	
		if not sam_line_first.startswith("@"):
			elements_first = sam_line_first.strip().split()
			try:
				read_ID_first = elements_first[0].strip()
				chrName_first = elements_first[2].strip()
				insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for second read is negative
			except:
				print "error in first read:", sam_line_first
			
			if chrName_first.startswith("chr") and (insert_size_first >= insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):			
				# only keep the reads mapped to chr 
				# if the first read is within insert size limit, check the second read
				# the insert_size for a pair is the same. If the first read is passed, the second will be passed, too.
				sam_line_second = inputfile_sam.readline()
				
				elements_second = sam_line_second.strip().split()
				try:
					read_ID_second = elements_second[0].strip()
					chrName_second = elements_second[2].strip()
					insert_size_second = abs(int(elements_second[8].strip()))			#  insert_size for second read is negative
				except:
					print "error in second read:", sam_line_second
				if read_ID_first == read_ID_second:		# check if the two reads belong to the same pair
					""" keep the pair as long as one read is not multiple mapping"""
					""" check if the reads from the same pair are mapped to the same chr	"""
					if (chrName_first == chrName_second) and ((not is_multiple_maping(elements_first)) or (not is_multiple_maping(elements_second))): 
						total_reads_num += 2
						print >> output_file, sam_line_first.strip()
						print >> output_file, sam_line_second.strip()
				else:
					print "first and second read ID do not match", read_ID_first, read_ID_second					
		sam_line_first = inputfile_sam.readline()
	inputfile_sam.close()
	output_file.close()
	print "total_reads_num: ", total_reads_num
	#return total_reads_num

def pair_end_indel(sam_file):
	""" process the reads by both pair_end and indel"""
	sam_file_name = sam_file[:(len(sam_file)-4)]
	print "pair_end_indel: ", sam_file_name
	output_file = open(sam_file_name + "_pairend_indel.sam", "w")
	inputfile_sam = open(currentPath + sam_file, "r")
	sam_line_first = inputfile_sam.readline() # the first read line in a pair
	total_reads_num = 0
	reads_after_process_total_number = 0
	
	insert_size_lower_bond = 100
	insert_size_upper_bond = 1000

	while sam_line_first!='':
		keep_this_pair = False	
		if not sam_line_first.startswith("@"):
			total_reads_num += 1
			elements_first = sam_line_first.strip().split()
			try:
				read_ID_first = elements_first[0].strip()
				chrName_first = elements_first[2].strip()
				insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for second read is negative
				indel_info_first = elements_first[5].strip()
				read_seq_first = elements_first[9].strip()
				qual_line_first = elements_first[10].strip()
			except:
				print "error in first read:", sam_line_first
			
			if chrName_first.startswith("chr") and (insert_size_first >= insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):			
				# only keep the reads mapped to chr 
				# if the first read is within insert size limit, check the second read
				# the insert_size for a pair is the same. If the first read is passed, the second will be passed, too.
				sam_line_second = inputfile_sam.readline()
				total_reads_num += 1
				
				elements_second = sam_line_second.strip().split()
				try:
					read_ID_second = elements_second[0].strip()
					chrName_second = elements_second[2].strip()
					insert_size_second = abs(int(elements_second[8].strip()))			#  insert_size for second read is negative
					indel_info_second = elements_second[5].strip()
					read_seq_second = elements_second[9].strip()
					qual_line_second = elements_second[10].strip()
				except:
					print "error in second read:", sam_line_second
				if read_ID_first == read_ID_second:		# check if the two reads belong to the same pair
					""" keep the pair as long as one read is not multiple mapping"""
					""" check if the reads from the same pair are mapped to the same chr	"""
					if (chrName_first == chrName_second) and ((not is_multiple_maping(elements_first)) or (not is_multiple_maping(elements_second))): 
						reads_after_process_total_number += 2
						
						read_qual_first = indel_correction(read_seq_first, qual_line_first, indel_info_first)
						sam_line_first = sam_line_first.replace(read_seq_first, read_qual_first[0])
						sam_line_first = sam_line_first.replace(qual_line_first, read_qual_first[1])
						print >> output_file, sam_line_first.strip()
						
						read_qual_second = indel_correction(read_seq_second, qual_line_second, indel_info_second)
						sam_line_second = sam_line_second.replace(read_seq_second, read_qual_second[0])
						sam_line_second= sam_line_second.replace(qual_line_second, read_qual_second[1])
						print >> output_file, sam_line_second.strip()
				else:
					print "first and second read ID do not match", read_ID_first, read_ID_second					
		sam_line_first = inputfile_sam.readline()
	inputfile_sam.close()
	output_file.close()
	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number: ", reads_after_process_total_number
	#return total_reads_num

def is_indel(cigar):
        n_indel = cigar.count("I") + cigar.count("D") + cigar.count("N") + cigar.count("S")
	if n_indel > 0 :
		return True
	else :
		return False	

def indel_correction(read_seq, qual_line, indel_info):
	indel_length = ""
	indel_type = ""
	read_seq_final = ""
	qual_line_final = ""
	if is_indel(indel_info):
		for info in indel_info:
			try:
				int(info)
				indel_length += info
			except:
				indel_type = info
				indel_length = int(indel_length)
				if indel_type == "S":	# mismatch, remove this part
					read_seq = read_seq[indel_length:]
					qual_line = qual_line[indel_length:]
				elif indel_type == "M":	# matched, keep this part
					read_seq_final += read_seq[:indel_length]
					qual_line_final += qual_line[:indel_length]
					read_seq = read_seq[indel_length:]
					qual_line = qual_line[indel_length:]
				elif indel_type == "I":	# insertion, remove
					read_seq = read_seq[indel_length:]
					qual_line = qual_line[indel_length:]
				elif indel_type == "D":	# deletion, insert "N" in read_seq, insert "!" (quality score = 0) in qual_line
					i = 0
					while i < indel_length:
						read_seq_final += "N"
						qual_line_final += "!"
						i += 1
				else:
					print "error in ", read_seq, " info type ", indel_type
				indel_length = ""
				indel_type = ""
	else:
		read_seq_final = read_seq
		qual_line_final = qual_line
	return (read_seq_final, qual_line_final)

def indel_process(sam_file):
	sam_file_name = sam_file[:(len(sam_file)-4)]
	print "indel_process: ", sam_file_name
	inputFile_sam = open(currentPath + sam_file, "r")
	outputFile_sam = open(currentPath + sam_file_name + "_indel.sam", "w")
	total_reads_num = 0
	
	for read in inputFile_sam:
		if not read.startswith("@"):
			total_reads_num += 1	
			elements_first = read.strip().split()
			try:
				indel_info = elements_first[5].strip()
				read_seq = elements_first[9].strip()
				qual_line = elements_first[10].strip()
				read_qual = indel_correction(read_seq, qual_line, indel_info)
				read = read.replace(read_seq, read_qual[0])
				read = read.replace(qual_line, read_qual[1])
				print >> outputFile_sam, read.strip()		
			except:
				#print "error in line: ", line
				pass
	print "total_reads_num: ", total_reads_num
	inputFile_sam.close()
	outputFile_sam.close()

def pair_end_indel_multiple():
	"""problem cannot call itself"""
	for infile in glob.glob(os.path.join(currentPath,'*.sam')):
		file_name = os.path.basename(infile)
		print "processing: ", file_name
		cmd = program_path + "sam_process.py -m sp -s " + file_name + " &"
		print cmd
		os.system(cmd)

def combine_files():
	cmd = "cat "
	combined_file = ""
	file_number = 0
	for infile in glob.glob(os.path.join(currentPath,'*.sam')):
		file_name = os.path.basename(infile)
		combined_file = file_name[:9]
		cmd += file_name + " "
		file_number += 1
	cmd = cmd + "> " + combined_file + "_" + str(file_number) + "_combined.sam"
	print cmd
	os.system(cmd)

def seperate_by_chr(sam_file):
	cmd = bash_path + "seperate_sam_by_chr.sh " + sam_file
	print cmd
	os.system(cmd)

def get_args():
	desc="variation call"
	usage = "snpPick_fish -s sam_file" 
	parser = OptionParser(usage = usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode",help = "", default="null")
	(options, args) = parser.parse_args()
	if options.samFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

def sam_process():
	options = get_args()
	sam_file = options.samFile
	mode = options.mode
	sam_file_name = sam_file[:(len(sam_file)-4)]
	
	if mode == "sp":
		pair_end_indel(sam_file)
	elif mode == "mp":
		pair_end_indel_multiple()
	elif mode == "sep":
		seperate_by_chr(sam_file)
	
	
if __name__=='__main__':
	
	#start = time.time()
	sam_process()
		
	
	#combine_files()
	#seperate_by_chr(sam_file)
	#pair_end_indel(sam_file)
	
	
	
	"""
	pair_end_filter(sam_file)
	sam_file = sam_file_name + "_pairend.sam"
	indel_process(sam_file)
	
	

	end = time.time()
	run_time = str(format((end - start), "0.3f"))
	print "run time is: " + run_time + "s"
	"""

