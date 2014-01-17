#!/usr/bin/python

# this is for pre-processing the sam file before snpPick
# It can process single-end file, pair-end file

import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser
from tools import *

class parameters:
	def __init__(self):
		self.chr_name = ""
		self.sam_file_name = ""
		self.insert_size_lower_bond = 0
		self.insert_size_upper_bond = 1000

def is_multiple_maping(elements):
	multiple_maping = False
	XA = ""
	try:
		XA = elements_first[21].strip()
		multiple_maping = True
	except:
		pass
	return multiple_maping, XA

def pair_end_filter(sam_file, chr_name):
	"""filter the reads by pair_end info"""
	print "pair_end_filter: ", parameter.sam_file_name
	output_file = open(parameter.sam_file_name + "_pairend.sam", "w")
	inputfile_sam = open(currentPath + sam_file, "r")
	sam_line_first = inputfile_sam.readline() # the first read line in a pair
	total_reads_num = 0
	covered_snp_total_number = 0

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
			# process all chr or one particular chr
			check_chr_name = chrName_first.startswith("chr") if (chr_name == "chr") else (chr_name == chrName_first)
			
			if check_chr_name and (insert_size_first > parameter.insert_size_lower_bond) and (insert_size_first <= parameter.insert_size_upper_bond):			
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
					if (chrName_first == chrName_second) and ((not is_multiple_maping(elements_first)[0]) or (not is_multiple_maping(elements_second)[0])): 
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

def pair_end_indel(sam_file, chr_name):
	""" process the reads by both pair_end and indel, does not process XA info here"""
	print "pair_end_indel: ", parameter.sam_file_name
	output_file = open(parameter.sam_file_name + "_pairend_indel.sam", "w")
	inputfile_sam = open(currentPath + sam_file, "r")
	sam_line_first = inputfile_sam.readline() # the first read line in a pair
	total_reads_num = 0
	reads_after_process_total_number = 0

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
			
			# process all chr or one particular chr
			check_chr_name = chrName_first.startswith("chr") if (chr_name == "chr") else (chr_name == chrName_first)
			if check_chr_name and (insert_size_first >= parameter.insert_size_lower_bond) and (insert_size_first <= parameter.insert_size_upper_bond):					# only keep the reads mapped to chr 
				# if the first read is within insert size limit, check the second read
				# the insert_size for a pair is the same. If the first read is passed, the second will be passed, too.
				sam_line_second = inputfile_sam.readline()
				total_reads_num += 1
				
				elements_second = sam_line_second.strip().split()
				try:
					read_ID_second = elements_second[0].strip()
					chrName_second = elements_second[2].strip()
					insert_size_second = abs(int(elements_second[8].strip()))		#  insert_size for second read is negative
					indel_info_second = elements_second[5].strip()
					read_seq_second = elements_second[9].strip()
					qual_line_second = elements_second[10].strip()
				except:
					print "error in second read:", sam_line_second
				if read_ID_first == read_ID_second:		# check if the two reads belong to the same pair
					""" keep the pair as long as one read is not multiple mapping"""
					""" check if the reads from the same pair are mapped to the same chr	"""
					if (chrName_first == chrName_second) and ((not is_multiple_maping(elements_first)[0]) or (not is_multiple_maping(elements_second)[0])): 
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

def single_end_indel(sam_file, chr_name):
	""" process the reads by both single_end and indel"""
	print "single_end_indel: ", parameter.sam_file_name
	output_file = open(parameter.sam_file_name + "_sginle_indel.sam", "w")
	inputfile_sam = open(currentPath + sam_file, "r")
	sam_line_first = inputfile_sam.readline() # the first read line in a pair
	total_reads_num = 0
	reads_after_process_total_number = 0

	while sam_line_first != '':
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
			
			# process all chr or one particular chr
			check_chr_name = chrName_first.startswith("chr") if (chr_name == "chr") else (chr_name == chrName_first)
			if check_chr_name and (insert_size_first > parameter.insert_size_lower_bond) and (insert_size_first <= parameter.insert_size_upper_bond):					# only keep the reads mapped to chr 
				reads_after_process_total_number += 1
				read_qual_first = indel_correction(read_seq_first, qual_line_first, indel_info_first)
				sam_line_first = sam_line_first.replace(read_seq_first, read_qual_first[0])
				sam_line_first = sam_line_first.replace(qual_line_first, read_qual_first[1])
				print >> output_file, sam_line_first.strip()
									
		sam_line_first = inputfile_sam.readline()
	inputfile_sam.close()
	output_file.close()
	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number: ", reads_after_process_total_number
	#return total_reads_num

def filter_by_chr():
	""" filter the reads by chr name and insert size """
	print "filter_by_chr: ", parameter.sam_file_name
	output_file = open(parameter.sam_file_name + "_prefiltered.sam", "w")
	inputfile_sam = open(currentPath + parameter.sam_file, "r")
	sam_line_first = inputfile_sam.readline() # the first read line in a pair
	total_reads_num = 0
	reads_after_process_total_number = 0

	while sam_line_first != '':
		if not sam_line_first.startswith("@"):
			total_reads_num += 1
			elements_first = sam_line_first.strip().split()
			try:
				chrName_first = elements_first[2].strip()
				insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for one of the read is negative
			except:
				print "error in first read:", sam_line_first
			
			# process all chr or one particular chr
			check_chr_name = chrName_first.startswith("chr") if (parameter.chr_name == "chr") else (parameter.chr_name == chrName_first)
			if check_chr_name and (insert_size_first > parameter.insert_size_lower_bond) and (insert_size_first <= parameter.insert_size_upper_bond):					# only keep the reads mapped to chr 
				reads_after_process_total_number += 1
				print >> output_file, sam_line_first.strip()
									
		sam_line_first = inputfile_sam.readline()
	inputfile_sam.close()
	output_file.close()
	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number: ", reads_after_process_total_number

def match_pairend():
	""" find the matched read pair in a disordered sam file """
	print "match_pairend: ", parameter.sam_file_name
	
	reads_dict = {}
	previous_size = 0
	total_reads_num = 0
	reads_after_process_total_number = 0
	
	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_pairend.sam", "w") as output_file:
			sam_line_first = inputfile_sam.readline() # the first read line in a pair
			while sam_line_first != '':
				if not sam_line_first.startswith("@"):
					total_reads_num += 1
					elements_first = sam_line_first.strip().split()
					try:
						read_ID_first = elements_first[0].strip()
					except:
						print "error in first read:", sam_line_first
					if read_ID_first not in reads_dict:
						reads_dict[read_ID_first] = sam_line_first.strip()
					else:
						reads_after_process_total_number += 2
						print >> output_file, sam_line_first.strip()
						print >> output_file, reads_dict[read_ID_first]
						del reads_dict[read_ID_first]
					if len(reads_dict)%100000 == 0:
						if len(reads_dict) != previous_size:
							print "current reads_dict size: ", len(reads_dict)
						previous_size = len(reads_dict)
					if 	total_reads_num%1000000 == 0:
						print "current line number : ", total_reads_num				
				sam_line_first = inputfile_sam.readline()
	
	with open(parameter.sam_file_name + "_pairend_leftover.sam", "w") as output_file_lefover:
		for key in reads_dict.keys():
			print >> output_file_lefover, reads_dict[key]
	
	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number paired: ", reads_after_process_total_number
	print "reads_after_process_total_number left: ", len(reads_dict)

def filter_match_pairend():
	""" find the matched read pair in a disordered sam file, combined filter_by_chr and match_pairend"""
	print "filter_match_pairend: ", parameter.sam_file_name
	
	reads_dict = {}
	previous_size = 0
	total_reads_num = 0
	reads_after_process_total_number = 0
	
	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_pairend.sam", "w") as output_file:
			sam_line_first = inputfile_sam.readline() # the first read line in a pair
			while sam_line_first != '':
				if not sam_line_first.startswith("@"):
					total_reads_num += 1
					elements_first = sam_line_first.strip().split()
					try:
						read_ID_first = elements_first[0].strip()
						chrName_first = elements_first[2].strip()
						insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for one of the read is negative
					except:
						print "error in first read:", sam_line_first
						
					# process all chr or one particular chr
					check_chr_name = chrName_first.startswith("chr") if (parameter.chr_name == "chr") else (parameter.chr_name == chrName_first)
					if check_chr_name and (insert_size_first > parameter.insert_size_lower_bond) \
						and (insert_size_first <= parameter.insert_size_upper_bond):					# only keep the reads mapped to chr 
						if read_ID_first not in reads_dict:
							reads_dict[read_ID_first] = sam_line_first.strip()
						else:
							reads_after_process_total_number += 2
							print >> output_file, sam_line_first.strip()
							print >> output_file, reads_dict[read_ID_first]
							del reads_dict[read_ID_first]
					if len(reads_dict)%100000 == 0:
						if len(reads_dict) != previous_size:
							print "current reads_dict size: ", len(reads_dict)
						previous_size = len(reads_dict)
					if 	total_reads_num%1000000 == 0:
						print "current line number : ", total_reads_num				
				sam_line_first = inputfile_sam.readline()
	
	with open(parameter.sam_file_name + "_pairend_leftover.sam", "w") as output_file_lefover:
		for key in reads_dict.keys():
			print >> output_file_lefover, reads_dict[key]
	
	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number paired: ", reads_after_process_total_number
	print "reads_after_process_total_number left: ", len(reads_dict)

def filter_by_XA():
	"""
	The file is already processed by filter_match_pairend, if both the read in a pair have XA
	ane they belong to the same chr in XA, check if the start position distance are within
	the insert size range >= 0 and <= 1000 at this moment
	"""
	print "filter_match_pairend: ", parameter.sam_file_name
	
	total_reads_num = 0
	reads_after_process_total_number = 0
	
	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_pairend.sam", "w") as output_file:
			sam_line_first = inputfile_sam.readline() # the first read line in a pair
			while sam_line_first != '':
				if not sam_line_first.startswith("@"):
					total_reads_num += 1
					elements_first = sam_line_first.strip().split()
					try:
						read_ID_first = elements_first[0].strip()
						chrName_first = elements_first[2].strip()
						insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for one of the read is negative
						indel_info_first = elements_first[5].strip()
						read_seq_first = elements_first[9].strip()
						qual_line_first = elements_first[10].strip()
						first_is_XA, first_XA_info = is_multiple_maping(elements_first)
					except:
						print "error in first read:", sam_line_first
						
					# process all chr or one particular chr, keep these steps for other files that do not need pair match
					check_chr_name = chrName_first.startswith("chr") if (chr_name == "chr") else (chr_name == chrName_first)
					if check_chr_name and (insert_size_first >= parameter.insert_size_lower_bond) and (insert_size_first <= parameter.insert_size_upper_bond):					# only keep the reads mapped to chr 
						# if the first read is within insert size limit, check the second read
						# the insert_size for a pair is the same. If the first read is passed, the second will be passed, too.
						sam_line_second = inputfile_sam.readline()
						total_reads_num += 1
						
						elements_second = sam_line_second.strip().split()
						try:
							read_ID_second = elements_second[0].strip()
							chrName_second = elements_second[2].strip()
							insert_size_second = abs(int(elements_second[8].strip()))		#  insert_size for one of the read is negative
							indel_info_second = elements_second[5].strip()
							read_seq_second = elements_second[9].strip()
							qual_line_second = elements_second[10].strip()
							second_is_XA, second_XA_info = is_multiple_maping(elements_second)
						except:
							print "error in second read:", sam_line_second
						if read_ID_first == read_ID_second:		# check if the two reads belong to the same pair
							""" check if the reads from the same pair are mapped to the same chr	"""
							if (chrName_first == chrName_second):
								keep_this_pair = True
								"""only check XA when both of the reads have XA"""							
								if first_is_XA and second_is_XA:
									pass

								if keep_this_pair:
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
								print "first and second read mapped to different chr", read_ID_first, chrName_first, read_ID_second, chrName_second
						else:
							print "first and second read ID do not match", read_ID_first, read_ID_second					
				sam_line_first = inputfile_sam.readline()
	
	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number: ", reads_after_process_total_number
	#return total_reads_num
	
def check_XA(first_XA_info, second_XA_info):
	# eg: XA:Z:chrY,+10113,101M,1; first to remove "XA:Z:"
	#XA:Z:chrY,+10113,101M,1;
	#XA:Z:chrY,-10264,101M,2;
	keep_this_pair = True
	first_XA_info = first_XA_info[5:]
	second_XA_info = second_XA_info[5:]
	first_XA_element = first_XA_info.split(";")
	second_XA_element = second_XA_info.split(";")
	first_XA_dict = {}
	second_XA_dict = {}
	for info in first_XA_element:
		data = info.split(",")
		chr_name = data[0]
		pos = data[1][1:]
		first_XA_dict[chr_name] = pos
	for info in second_XA_element:
		data = info.split(",")
		chr_name = data[0]
		pos = data[1][1:]
		second_XA_dict[chr_name] = pos
	
	common_chr = [chr for chr in first_XA_dict.keys() if chr in second_XA_dict.keys()]
	for chr in common_chr:
		distance = abs(abs(int(first_XA_dict[chr])) - abs(int(second_XA_dict[chr])))
		if distance >= parameter.insert_size_lower_bond and distance <= parameter.insert_size_upper_bond:
			keep_this_pair = False
	return keep_this_pair

def is_indel(cigar):
	n_indel = cigar.count("I") + cigar.count("D") + cigar.count("N") + cigar.count("S")
	is_indel = True if n_indel > 0 else False
	return is_indel

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

def seperate_by_chr():
	cmd = bash_path + "seperate_sam_by_chr.sh " + parameter.sam_file
	print cmd
	os.system(cmd)

def sam_process(sam_file, chr_name, mode):
	if mode == "single":
		single_end_indel(sam_file, chr_name)
	elif mode == "pair_single":
		pair_end_indel(sam_file, chr_name)
	elif mode == "pair_mutiple":
		pair_end_indel_multiple()
	elif mode == "sep_chr":
		seperate_by_chr()
	elif mode == "combine":
		combine_files()
	elif mode == "prefilter":	
		filter_by_chr()
	elif mode == "match":		
		match_pairend()
	elif mode == "fm":		
		filter_match_pairend()

def get_args():
	desc="variation call"
	usage = "snpPick_fish -s sam_file -c chr -m mode" 
	parser = OptionParser(usage = usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="chr")
	parser.add_option("-m", "--mode", type="string", dest="mode",help = "", default="null")
	(options, args) = parser.parse_args()
	if options.mode == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__=='__main__':
	
	global parameter
	parameter = parameters()
	
	options = get_args()
	parameter.sam_file = options.samFile
	parameter.chr_name = options.chrName
	parameter.sam_file_name = parameter.sam_file[:(len(parameter.sam_file)-4)]
	parameter.mode = options.mode
	
	start = time.time()
	sam_process(parameter.sam_file, parameter.chr_name, parameter.mode)
	#run_time = str(format((end - start), "0.3f"))
	print "run time is: ", round((time.time() - start), 3), "s"
	
	"""
	#combine_files()
	#seperate_by_chr(sam_file)
	#pair_end_indel(sam_file)
	pair_end_filter(sam_file)
	sam_file = sam_file_name + "_pairend.sam"
	indel_process(sam_file)
	"""

