#!/usr/bin/python

# this is for pre-processing the sam file before snpPick
# It can process single-end file, pair-end file

import os, glob, subprocess, random, operator, time, sys, math, string
from optparse import OptionParser
#from tools import *
from repeatRemove_sorted import *

class parameters:
	def __init__(self):
		self.chr_name = ""
		self.sam_file_name = ""
		self.insert_size_lower_bond = 0
		self.insert_size_upper_bond = 1000

def is_multiple_maping(elements):
	multiple_maping = False
	XA = elements[-1].strip()
	if XA.startswith('XA'):
		multiple_maping = True
	else:
		XA = ""
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
	output_file = open(parameter.sam_file_name + "_indel.sam", "w")
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
	"""
	:param sam_file:
	:param chr_name:
	:return:
	process the reads by both single_end and indel
	remove reads with "N" in cigar and correct "ID"
	"""
	print "single_end_indel: ", parameter.sam_file_name
	output_file = open(parameter.sam_file_name + "_indel.sam", "w")
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
			check_chr_name = chrName_first.startswith("Chr") if (chr_name == "Chr") else (chr_name == chrName_first)    #jiang
			#check_chr_name = chrName_first.startswith("chr") if (chr_name == "chr") else (chr_name == chrName_first)
			if check_chr_name:					# only keep the reads mapped to chr

				if is_indel(indel_info_first):
					if indel_info_first.count("N") == 0:
						reads_after_process_total_number += 1
						read_qual_first = indel_correction(read_seq_first, qual_line_first, indel_info_first)
						sam_line_first = sam_line_first.replace(read_seq_first, read_qual_first[0])
						sam_line_first = sam_line_first.replace(qual_line_first, read_qual_first[1])
						print >> output_file, sam_line_first.strip()
				else:
					reads_after_process_total_number += 1
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
	output_file = open(parameter.sam_file_name + parameter.chr_name[3:] + "_chr.sam", "w")
	output_removed_file = open(parameter.sam_file_name + "_" + parameter.chr_name + "_nonchr.sam", "w")
	inputfile_sam = open(currentPath + parameter.sam_file, "r")
	sam_line_first = inputfile_sam.readline() # the first read line in a pair
	total_reads_num = 0
	reads_after_process_total_number = 0

	while sam_line_first != '':
		if sam_line_first.startswith("@"):
			print >> output_file, sam_line_first.strip()
		else:
			total_reads_num += 1
			elements_first = sam_line_first.strip().split()
			try:
				chrName_first = elements_first[2].strip()
				insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for one of the read is negative
			except:
				print "error in first read:", sam_line_first
			
			# process all chr or one particular chr
			check_chr_name = chrName_first.startswith("chr") if (parameter.chr_name == "chr") else (parameter.chr_name == chrName_first)
			# to remove reads mapped to non chrs
			try:
				chr_number = int(chrName_first[3:])
			except:
				#print "chr_number", chr_number
				check_chr_name = False
				print >> output_removed_file, sam_line_first.strip()

			if check_chr_name and (insert_size_first >= parameter.insert_size_lower_bond) and (insert_size_first <= parameter.insert_size_upper_bond):					# only keep the reads mapped to chr
				reads_after_process_total_number += 1
				print >> output_file, sam_line_first.strip()
			#else:
			#	print >> output_removed_file, sam_line_first.strip()

		sam_line_first = inputfile_sam.readline()
	inputfile_sam.close()
	output_file.close()
	output_removed_file.close()
	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number: ", reads_after_process_total_number
	print "reads removed : ", total_reads_num - reads_after_process_total_number
	print "reads kept % : ", round(float(reads_after_process_total_number)*100/total_reads_num, 2)

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
				#if not sam_line_first.startswith("@"):  # mali data
				if sam_line_first.startswith("@ILLUMINA"):  # wli data
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
					if total_reads_num%1000000 == 0:
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
					if check_chr_name and (insert_size_first >= parameter.insert_size_lower_bond) \
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
					if total_reads_num%1000000 == 0:
						print "current line number : ", total_reads_num				
				sam_line_first = inputfile_sam.readline()
	
	with open(parameter.sam_file_name + "_pairend_removed.sam", "w") as output_file_lefover:
		for key in reads_dict.keys():
			print >> output_file_lefover, reads_dict[key]
	
	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number paired: ", reads_after_process_total_number
	print "reads_after_process_total_number left: ", len(reads_dict)

def single_end_xa():
	""" remove mutiple mapped reads in single end file.
	Assume the file is already processed by chr, insert size"""

	print "filter_match_pairend: ", parameter.sam_file_name

	total_reads_num = 0
	reads_after_process_total_number = 0
	xa_total = 0

	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_XA.sam", "w") as output_file:
			sam_line_first = inputfile_sam.readline() # the first read line in a pair
			while sam_line_first != '':
				if sam_line_first.startswith("@"):
					# the @header is needed for samtoos sorting
					print >> output_file, sam_line_first.strip()
				else:
					total_reads_num += 1
					elements_first = sam_line_first.strip().split()
					try:
						first_is_XA, first_XA_info = is_multiple_maping(elements_first)
					except:
						print "error in first read:", sam_line_first

					if not first_is_XA:
						print >> output_file, sam_line_first.strip()
						reads_after_process_total_number += 1
					else:
						xa_total += 1
				sam_line_first = inputfile_sam.readline()

	print "total_reads_num: ", total_reads_num
	print "xa_total: ", xa_total
	print "reads_after_process_total_number: ", reads_after_process_total_number
	print "kept percentage", round(float(reads_after_process_total_number)*100/total_reads_num, 2)
	#return total_reads_num

def filter_by_XA():
	"""
	The file is already processed by filter_match_pairend, if both the read in a pair have XA
	ane they belong to the same chr in XA, check if the start position distance are within
	the insert size range > 0 and <= 1000 at this moment
	"""
	print "filter_match_pairend: ", parameter.sam_file_name
	
	total_reads_num = 0
	reads_after_process_total_number = 0
	
	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_XA.sam", "w") as output_file:
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
					check_chr_name = chrName_first.startswith("Chr") if (parameter.chr_name == "Chr") else (parameter.chr_name == chrName_first)    #jiang
					#check_chr_name = chrName_first.startswith("chr") if (parameter.chr_name == "chr") else (parameter.chr_name == chrName_first)

					if check_chr_name and (insert_size_first > parameter.insert_size_lower_bond) and (insert_size_first <= parameter.insert_size_upper_bond):					# only keep the reads mapped to chr
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
									keep_this_pair = check_XA(first_XA_info, second_XA_info)
									
								if keep_this_pair:
									reads_after_process_total_number += 2
									
									read_qual_first = indel_correction(read_seq_first, qual_line_first, indel_info_first)
									sam_line_first = sam_line_first.replace(read_seq_first, read_qual_first[0])
									sam_line_first = sam_line_first.replace(qual_line_first, read_qual_first[1])
									#print sam_line_first
									print >> output_file, sam_line_first.strip()
									
									read_qual_second = indel_correction(read_seq_second, qual_line_second, indel_info_second)
									sam_line_second = sam_line_second.replace(read_seq_second, read_qual_second[0])
									sam_line_second= sam_line_second.replace(qual_line_second, read_qual_second[1])
									#print sam_line_second
									print >> output_file, sam_line_second.strip()
							else:
								print "first and second read mapped to different chr", read_ID_first, chrName_first, read_ID_second, chrName_second
						else:
							print "first and second read ID do not match", read_ID_first, read_ID_second					
				sam_line_first = inputfile_sam.readline()
	
	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number: ", reads_after_process_total_number
	#return total_reads_num

def filter_by_XA_mimi():
	"""
	The file is already processed by filter_match_pairend, if both the read in a pair have XA
	ane they belong to the same chr in XA, check if the start position distance are within
	the insert size range >= 0 and <= 1000 at this moment
	This is for mimi project. the indel is not processed for sorting purpose
	"""
	print "filter_match_pairend: ", parameter.sam_file_name

	total_reads_num = 0
	reads_after_process_total_number = 0

	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_XA.sam", "w") as output_file:
			sam_line_first = inputfile_sam.readline() # the first read line in a pair
			while sam_line_first != '':
				if sam_line_first.startswith("@"):
					# the @header is needed for samtoos sorting
					print >> output_file, sam_line_first.strip()
					#sam_line_first = inputfile_sam.readline()
				else:
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
					check_chr_name = chrName_first.startswith("chr") if (parameter.chr_name == "chr") else (parameter.chr_name == chrName_first)
					if check_chr_name and (insert_size_first > parameter.insert_size_lower_bond) and (insert_size_first <= parameter.insert_size_upper_bond):					# only keep the reads mapped to chr
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
									keep_this_pair = check_XA(first_XA_info, second_XA_info)

								if keep_this_pair:
									reads_after_process_total_number += 2

									#read_qual_first = indel_correction(read_seq_first, qual_line_first, indel_info_first)
									#sam_line_first = sam_line_first.replace(read_seq_first, read_qual_first[0])
									#sam_line_first = sam_line_first.replace(qual_line_first, read_qual_first[1])
									#print sam_line_first
									print >> output_file, sam_line_first.strip()

									#read_qual_second = indel_correction(read_seq_second, qual_line_second, indel_info_second)
									#sam_line_second = sam_line_second.replace(read_seq_second, read_qual_second[0])
									#sam_line_second= sam_line_second.replace(qual_line_second, read_qual_second[1])
									#print sam_line_second
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
		try:
			data = info.split(",")
			chr_name = data[0]
			pos = data[1][1:]
			first_XA_dict[chr_name] = pos
		except:
			#print info
			pass
	for info in second_XA_element:
		try:
			data = info.split(",")
			chr_name = data[0]
			pos = data[1][1:]
			second_XA_dict[chr_name] = pos
		except:
			#print info
			pass
	
	#print first_XA_dict, second_XA_dict
	common_chr = [chr for chr in first_XA_dict.keys() if chr in second_XA_dict.keys()]
	
	for chr in common_chr:
		distance = abs(abs(int(first_XA_dict[chr])) - abs(int(second_XA_dict[chr])))
		#print distance
		if distance > parameter.insert_size_lower_bond and distance <= parameter.insert_size_upper_bond:
			keep_this_pair = False
			break
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

def indel_process_rnaseq(sam_file):
	# for rna seq data, to process indel information

	sam_file_name = sam_file[:(len(sam_file)-4)]
	print "indel_process: ", sam_file_name
	total_reads_num = 0
	with open(currentPath + sam_file_name + "_indel.sam", "w") as outputFile_sam:
		with open(currentPath + sam_file, "r") as inputFile_sam:
			for read in inputFile_sam:
				if not read.startswith("@"):
					total_reads_num += 1	
					elements_first = read.strip().split()
					try:
						chr_name = elements_first[2]
						start_pos = elements_first[3]

						cigar = elements_first[5].strip()
						read_seq = elements_first[9].strip()
						qual_line = elements_first[10].strip()
						read_qual = indel_correction(read_seq, qual_line, cigar)
						read = read.replace(read_seq, read_qual[0])
						read = read.replace(qual_line, read_qual[1])
						print >> outputFile_sam, read.strip()		
					except:
						#print "error in line: ", line
						pass
	print "total_reads_num: ", total_reads_num

def indel_process(sam_file):
	sam_file_name = sam_file[:(len(sam_file)-4)]
	print "indel_process: ", sam_file_name
	total_reads_num = 0
	with open(currentPath + sam_file_name + "_indel.sam", "w") as outputFile_sam:
		with open(currentPath + sam_file, "r") as inputFile_sam:
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


def pair_end_indel_multiple():
	"""problem cannot call itself"""
	for infile in glob.glob(os.path.join(currentPath,'*.sam')):
		file_name = os.path.basename(infile)
		print "processing: ", file_name
		cmd = program_path + "sam_process.py -m sp -s " + file_name + " &"
		print cmd
		os.system(cmd)

def combine_files():
	# combine sam files
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

def load_repeat(repeat_cnv_file, chr_index, start_pos_index, end_pos_index):
	# This is used load repeat, cnv from files for combine_cnv_repeat
	print "combine_cnv_repeat: ", repeat_cnv_file
	repeat_list = []
	#read_length = 101
	repeat_length_cutoff = 15       # a cutoff for the repeat, the overlap range for repeat and read
	repeat_cnv_file_name = repeat_cnv_file[:(len(repeat_cnv_file)-4)]
		
	with open(currentPath + repeat_cnv_file, "r") as inputfile_repeat:
		for line in inputfile_repeat:
			if not line.startswith("#"):
				elements = line.strip().split()
				start_pos = 0
				end_pos = 0
				chr = ""
				try:
					start_pos = int(elements[start_pos_index].strip())
					end_pos = int(elements[end_pos_index].strip())
					chr = elements[chr_index].strip()
				except:
					print "error in line:", line
				#if end_pos - start_pos >= read_length:
				if end_pos - start_pos >= repeat_length_cutoff:
					repeat_list.append((start_pos, end_pos, chr, repeat_cnv_file_name))
	print "list size:", repeat_cnv_file, len(repeat_list)
	return repeat_list

def combine_repeat_list(repeat_list_1, repeat_list_2):
	repeat_list = []
	i = 0
	j = 0
	start_pos_1 = 0
	end_pos_1 = 0
	start_pos_2 = 0
	end_pos_2 = 0

	while i < len(repeat_list_1) or j < len(repeat_list_2):
		
		if i == len(repeat_list_1) - 1 and j < len(repeat_list_2) - 1:
			repeat_list.append(repeat_list_2[j])
			j += 1
		elif j == len(repeat_list_2) - 1 and i < len(repeat_list_1) - 1:
			#print "here"
			repeat_list.append(repeat_list_1[i])
			i += 1
		elif j == len(repeat_list_2) - 1 and i == len(repeat_list_1) - 1:
			if j >= i:
				repeat_list.append(repeat_list_2[j])
			else:
				repeat_list.append(repeat_list_1[i])
			i += 1
			j += 1
		else:
			if i < len(repeat_list_1):
				start_pos_1 = repeat_list_1[i][0]
				end_pos_1 = repeat_list_1[i][1]
	
			if j < len(repeat_list_2):
				start_pos_2 = repeat_list_2[j][0]
				end_pos_2 = repeat_list_2[j][1]
		
			if end_pos_1 <= start_pos_2:
				repeat_list.append(repeat_list_1[i])
				i += 1
			elif start_pos_1 < start_pos_2:
				if end_pos_1 <= end_pos_2:
					repeat_list.append(repeat_list_1[i])
					i += 1
				elif end_pos_1 > end_pos_2:	# repeat 1 contains repeat 2. 2 can be skipped
					j += 1
					#print len(repeat_list), start_pos_1, end_pos_1, start_pos_2, end_pos_2
			elif start_pos_1 >= start_pos_2:
				if start_pos_1 < end_pos_2:
					if end_pos_1 <= end_pos_2:	# repeat 2 contains repeat 1. 1 can be skipped
						i += 1
						#print len(repeat_list), start_pos_1, end_pos_1, start_pos_2, end_pos_2
					elif end_pos_1 > end_pos_2:
						repeat_list.append(repeat_list_2[j])
						j += 1
				elif start_pos_1 >= end_pos_2:
					try:
						repeat_list.append(repeat_list_2[j])
						j += 1
					except:
						print start_pos_1, end_pos_2, i, j-1
						sys.exit(1)
					
			else:
				print "combine_repeat_list error in", start_pos_1, end_pos_1, start_pos_2, end_pos_2
			#print len(repeat_list), start_pos_1, end_pos_1, start_pos_2, end_pos_2
	return repeat_list		

def combine_cnv_repeat(file_1, file_2):
	
	print "combine_cnv_repeat: ", file_1, file_2
	file_1_name = file_1[:(len(file_1)-4)]
	file_2_name = file_2[:(len(file_2)-4)]
	"""
	#MultSNPs_chrX (same pos with DGV_chrX)
	chr_index = 1
	start_pos_index = 2
	end_pos_index = 3
	repeat_list_1 = load_repeat(file_1, chr_index, start_pos_index, end_pos_index)
	"""
	#rmsk_chrX_hg19
	chr_index = 5
	start_pos_index = 6
	end_pos_index = 7
	repeat_list_1 = load_repeat(file_1, chr_index, start_pos_index, end_pos_index)
	"""
	#SegDups_chrX
	chr_index = 1
	start_pos_index = 2
	end_pos_index = 3
	repeat_list_2 = load_repeat(file_2, chr_index, start_pos_index, end_pos_index)
	"""
	#combined
	chr_index = 2
	start_pos_index = 0
	end_pos_index = 1
	repeat_list_2 = load_repeat(file_2, chr_index, start_pos_index, end_pos_index)

	
	if len(repeat_list_1) > 0 and len(repeat_list_2) > 0:
		repeat_list = combine_repeat_list(repeat_list_1, repeat_list_2)
		print "repeat_list combined size before cleanup", len(repeat_list)
		
		repeat_list = cleanup_repeat_list(repeat_list)
		print "repeat_list combined size", len(repeat_list)
		# output combined file
		title_line = "#start_pos	end_pos chr	ori_file"
		with open(currentPath + file_1_name + "_" + file_2_name + ".txt", "w") as outputfile_combine:
			print >> outputfile_combine, title_line
			for repeat in repeat_list:
				print >> outputfile_combine, repeat[0], repeat[1], repeat[2], repeat[3]
	else:
		print "list size is zero.", len(repeat_list_1), len(repeat_list_2)

def cleanup_repeat_list(ori_repeat_list):
	# compare the repeat inside the list, remove the repeats that are contained in other repeats
	new_repeat_list = []
	current_repeat = [0, 0]
	for repeat in ori_repeat_list:
		if repeat[0] >= current_repeat[0] and repeat[1] <= current_repeat[1]:
			# if the next repeat is covered by the previous repeat, skip it.
			pass
		elif repeat[0] <= current_repeat[0] and repeat[1] >= current_repeat[1]:
			# if the next repeat covers the previous repeat, remove the previous repeat, keep the current.
			del new_repeat_list[-1]
			new_repeat_list.append(repeat)
			current_repeat = list(repeat)
		else:
			new_repeat_list.append(repeat)
			current_repeat = list(repeat)
	return new_repeat_list
		
def filter_by_repeat(repeat_cnv_file):	
	"""
	The file is already processed by insert size, XA, chr
	"""	
	print "filter_by_repeat: ", repeat_cnv_file
	repeat_cnv_file_name = repeat_cnv_file[:(len(repeat_cnv_file)-4)]
	
	#MultSNPs_chrX
	chr_index = 1
	start_pos_index = 2
	end_pos_index = 3
	
	repeat_list = load_repeat(repeat_cnv_file, chr_index, start_pos_index, end_pos_index)
	
	total_reads_num = 0
	reads_after_process_total_number = 0
	
	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_" + repeat_cnv_file_name + ".sam", "w") as output_file:
			sam_line_first = inputfile_sam.readline() # the first read line in a pair
			while sam_line_first != '':
				if not sam_line_first.startswith("@"):
					total_reads_num += 1
					elements_first = sam_line_first.strip().split()
					chrName_first = elements_first[2].strip()
					insert_size_first = abs(int(elements_first[8].strip()))
					try:
						read_ID_first = elements_first[0].strip()
						#chrName_first = elements_first[2].strip()
						start_pos_first = int(elements_first[3].strip())
						#insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for one of the read is negative
						#indel_info_first = elements_first[5].strip()
						#read_seq_first = elements_first[9].strip()
						#qual_line_first = elements_first[10].strip()
						#first_is_XA, first_XA_info = is_multiple_maping(elements_first)
					except:
						print "error in first read:", sam_line_first
						
					# process all chr or one particular chr, keep these steps for other files that do not need pair match
					check_chr_name = chrName_first.startswith("chr") if (parameter.chr_name == "chr") else (parameter.chr_name == chrName_first)
					if check_chr_name and (insert_size_first > parameter.insert_size_lower_bond) and (insert_size_first <= parameter.insert_size_upper_bond):					# only keep the reads mapped to chr 
						# if the first read is within insert size limit, check the second read
						# the insert_size for a pair is the same. If the first read is passed, the second will be passed, too.
						sam_line_second = inputfile_sam.readline()
						total_reads_num += 1

def add_header(sam_file):
	# for samtools sort
	head_path = "/home/guoxing/disk2/lima/hg_header/"
	hg18_header = "headerhg18.txt"
	hg19_header = "headerhg19.txt"

	cat = "cat " + head_path + hg19_header + " " + sam_file + " > " + sam_file + "_wheader"
	grep_Process = subprocess.Popen(cat, shell=True)
	grep_Process.wait()

	os.system("mv " + sam_file + "_wheader " + sam_file)

def samtools_sort(sam_file):
	other_path = "/home/guoxing/disk2/ngs/morehouse/other/"
	samtools_path = other_path + "samtools-0.1.18/"

	# sort the sam file
	print "sorting running"

	sort_input = sam_file[:len(sam_file)-4]
	sam2bam = samtools_path + "samtools view -bS " + sam_file + " > " + sort_input + ".bam"
	print sam2bam
	sam2bam_Process = subprocess.Popen(sam2bam, shell=True)
	sam2bam_Process.wait()

	sortbam = samtools_path + "samtools sort " + sort_input + ".bam " + sort_input + "_sorted"
	print sortbam
	sortbam_Process = subprocess.Popen(sortbam, shell=True)
	sortbam_Process.wait()

	bam2sam = samtools_path + "samtools view " + sort_input + "_sorted.bam > " + sort_input + "_sorted.sam"
	print bam2sam
	bam2sam_Process = subprocess.Popen(bam2sam, shell=True)
	bam2sam_Process.wait()

	#os.system("rm " + sort_input + ".bam")
	#os.system("rm " + sort_input + "_sorted.bam")

def os_sort(sam_file):
	print "sorting", sam_file
	cmd = "sort -k 3,3 -k 4,4n " + parameter.sam_file + " > " + parameter.sam_file_name + "_sorted.sam"
	os.system(cmd)
	print "done"

def mimi_solid_base_remove(sam_file):
	# for solid data, to remove 10 bases after the primer, which are not accurate.

	read_length_cutoff = 32
	total_reads_num = 0
	reads_after_process_total_number = 0

	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_base_cleaned.sam", "w") as output_file:
			sam_line_first = inputfile_sam.readline() # the first read line in a pair
			while sam_line_first != '':
				if not sam_line_first.startswith("@"):
					total_reads_num += 1
					elements_first = sam_line_first.strip().split()
					try:
						read_seq_first = elements_first[9].strip()
						if len(read_seq_first) > read_length_cutoff:
							print >> output_file, sam_line_first.strip()
							reads_after_process_total_number += 1
					except:
						print "error in first read:", sam_line_first

				sam_line_first = inputfile_sam.readline()

	print "total_reads_num: ", total_reads_num
	print "reads_after_process_total_number: ", reads_after_process_total_number
	print "removed reads: ", total_reads_num-reads_after_process_total_number
	print "removed percentage: ", round(float(total_reads_num-reads_after_process_total_number)/total_reads_num, 2)


def sep_pairend():
	"""
	seperate pariend from simulation data
	:return:
	"""
	n = 0
	with open(parameter.sam_file_name + "_1.txt", "w") as sam_output_first:
		with open(parameter.sam_file_name + "_2.txt", "w") as sam_output_second:
			with open(parameter.sam_file, "r") as sam_input:
				for line in sam_input:
					n += 1
					if line.startswith("@"):
						print >> sam_output_first, line.strip()
						print >> sam_output_second, line.strip()
					else:
						if n % 2 == 0:
							print >> sam_output_first, line.strip()
						elif n % 2 == 1:
							print >> sam_output_second, line.strip()

def add_length():
	target_length = 101
	with open(parameter.sam_file_name + "_" + str(target_length) + ".txt", "w") as sam_output:
		with open(parameter.sam_file, "r") as sam_input:
			for line in sam_input:
				if line.startswith("@"):
					print >> sam_output, line.strip()
				else:
					line = line.strip()
					elements = line.split()
					cigar = elements[5]
					ori_seq = elements[9]
					ori_qs = elements[10]
					seq = elements[9]
					qs = elements[10]

					seq_length = len(seq)
					for i in range(target_length - seq_length):
						seq += "N"
						qs += "R"
					line = string.replace(line, cigar, str(target_length)+"M", 1)
					line = string.replace(line, ori_seq, seq, 1)
					line = string.replace(line, ori_qs, qs, 1)

					print >> sam_output, line

def meth_pos_process():
	# adjust position in crick read
	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_pos_adjusted.sam", "w") as output_file:
			for line in inputfile_sam:
				line = line.strip()
				elements = line.split()
				start_pos = int(elements[3].strip())
				read_seq = elements[9].strip()
				new_start_pos = start_pos - len(read_seq) + 1
				line = line.replace(str(start_pos), str(new_start_pos))
				print >> output_file, line

def meth_crick_process():
	# adjust position in crick read and reverse complement the read seq, reverse the qual
	with open(currentPath + parameter.sam_file, "r") as inputfile_sam:
		with open(parameter.sam_file_name + "_crick_processed.sam", "w") as output_file:
			for line in inputfile_sam:
				line = line.strip()
				elements = line.split()
				start_pos = int(elements[3].strip())
				read_seq = elements[9].strip()
				qual_line = elements[10].strip()
				new_start_pos = start_pos - len(read_seq) + 1
				new_read_seq = reverse_complementary(read_seq)
				new_qual_line = qual_line[::-1]
				line = line.replace(str(start_pos), str(new_start_pos))
				line = line.replace(read_seq, new_read_seq)
				line = line.replace(qual_line, new_qual_line)
				print >> output_file, line

def sam_process(sam_file, chr_name, mode):
	if mode == "single":
		single_end_indel(sam_file, chr_name)
	elif mode == "pair_indel":
		pair_end_indel(sam_file, chr_name)
	elif mode == "pair_mutiple":
		pair_end_indel_multiple()
	elif mode == "sep_chr":
		seperate_by_chr()
	elif mode == "combine":
		combine_files()
	elif mode == "chr":
		filter_by_chr()
	elif mode == "match":	
		match_pairend()
	elif mode == "fm":
		filter_match_pairend()
	elif mode == "xa":		
		filter_by_XA()
	elif mode == "xamimi":
		filter_by_XA_mimi()
	elif mode == "indel":		
		indel_process(sam_file)
	elif mode == "single_xa":
		single_end_xa()
	elif mode == "sort":
		#add_header(sam_file)
		samtools_sort(sam_file)
	elif mode == "ossort":
		os_sort(sam_file)
	elif mode == "base_remove":
		# to remove bases following primer in solid data
		mimi_solid_base_remove(sam_file)
	elif mode == "ext_single":
		extract_single_overlapped_read(parameter.sam_file)
	elif mode == "add_length":
		add_length()
	elif mode == "sep_pairend":
		sep_pairend()
	elif mode == "meth_crick":
		meth_crick_process()
	elif mode == "mimi":
		ori_sam_file_name = parameter.sam_file_name
		"""
		print "1. filter chr and find pairend",  parameter.sam_file
		filter_match_pairend()

		parameter.sam_file_name = parameter.sam_file_name + "_pairend"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "2. filter by XA",  parameter.sam_file
		filter_by_XA_mimi()

		parameter.sam_file_name = parameter.sam_file_name + "_XA"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "3. sorting",  parameter.sam_file
		add_header(parameter.sam_file)
		samtools_sort(parameter.sam_file)
		"""
		parameter.sam_file_name = parameter.sam_file_name + "_pairend"
		parameter.sam_file_name = parameter.sam_file_name + "_XA"

		parameter.sam_file_name = parameter.sam_file_name + "_sorted"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "4. repeat remove",  parameter.sam_file
		rmsk_file = "/home/guoxing/disk2/lima/rmsk_chrX_hg19_MultSNPs_chrX_SegDups_chrX.txt"
		repeat_remove_mimi(rmsk_file, parameter.sam_file)

		parameter.sam_file_name = parameter.sam_file_name + "_rmsk"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "5. find matched pairend after repeat remove",  parameter.sam_file
		filter_match_pairend()

		#parameter.sam_file_name = parameter.sam_file_name[:len(parameter.sam_file_name)-5]
		#parameter.sam_file = parameter.sam_file_name + ".sam"
		print "6. extract_single_overlapped_read",  parameter.sam_file
		extract_single_overlapped_read(parameter.sam_file)

		parameter.sam_file_name = parameter.sam_file_name + "_combined"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "7. find matched pairend from combined file and process indel",  parameter.sam_file
		pair_end_indel(parameter.sam_file, parameter.chr_name)

		#snpPick_mimi -s NA12893_S1_ChrXnew_pairend_XA_sorted_rmsk_combined_indel.sam -c chrX -m update -d NA12893_S1_chrX

		print "8. clean up"
		# keep the pairend_XA_sorted.sam and pairend_XA_sorted_rmsk.sam
		os.system("rm " + ori_sam_file_name + ".sam")
		os.system("rm " + ori_sam_file_name + "_pairend.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_removed.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA_sorted_record.txt")
		sorted_rmsk_name = ori_sam_file_name + "_pairend_XA_sorted_rmsk"
		os.system("rm " + sorted_rmsk_name + "_pairend*.sam")
		os.system("rm " + sorted_rmsk_name + "_recovered*.sam")
		os.system("rm " + sorted_rmsk_name + "_removed.sam")
		os.system("rm " + sorted_rmsk_name + "_combined.sam")

	elif mode == "exon":
		ori_sam_file_name = parameter.sam_file_name

		print "1. filter chr and find pairend",  parameter.sam_file
		filter_match_pairend()

		parameter.sam_file_name = parameter.sam_file_name + "_pairend"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "2. filter by XA",  parameter.sam_file
		filter_by_XA_mimi()

		parameter.sam_file_name = parameter.sam_file_name + "_XA"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "3. sorting",  parameter.sam_file

		cmd = "sort -k 3,3 -k 4,4n " + parameter.sam_file + " > " + parameter.sam_file_name + "_sorted.sam"
		os.system(cmd)

		parameter.sam_file_name = parameter.sam_file_name + "_sorted"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "4. repeat remove",  parameter.sam_file
		rmsk_file = "/home/guoxing/disk2/solid/common_files/hg18_rmsk_chr/rmsk_" + parameter.chr_name + ".txt"
		print rmsk_file
		repeat_remove_fish_wli(rmsk_file, parameter.sam_file)

		parameter.sam_file_name = parameter.sam_file_name + "_rmsk"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "5. find matched pairend after repeat remove",  parameter.sam_file
		filter_match_pairend()

		print "6. extract_single_overlapped_read",  parameter.sam_file
		extract_single_overlapped_read(parameter.sam_file)

		parameter.sam_file_name = parameter.sam_file_name + "_combined"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "7. find matched pairend from combined file and process indel",  parameter.sam_file
		pair_end_indel(parameter.sam_file, parameter.chr_name)

		#snpPick_mimi -s NA12893_S1_ChrXnew_pairend_XA_sorted_rmsk_combined_indel.sam -c chrX -m update -d NA12893_S1_chrX

		print "8. clean up"
		# keep the pairend_XA_sorted.sam and pairend_XA_sorted_rmsk.sam
		#os.system("rm " + ori_sam_file_name + ".sam")
		"""
		os.system("rm " + ori_sam_file_name + "_pairend.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_removed.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA_sorted_record.txt")
		sorted_rmsk_name = ori_sam_file_name + "_pairend_XA_sorted_rmsk"
		os.system("rm " + sorted_rmsk_name + "_pairend*.sam")
		os.system("rm " + sorted_rmsk_name + "_recovered*.sam")
		os.system("rm " + sorted_rmsk_name + "_removed.sam")
		os.system("rm " + sorted_rmsk_name + "_combined.sam")
		"""

	elif mode == "fish_wli":
		ori_sam_file_name = parameter.sam_file_name

		print "1. filter chr and find pairend",  parameter.sam_file
		filter_match_pairend()

		parameter.sam_file_name = parameter.sam_file_name + "_pairend"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "2. filter by XA",  parameter.sam_file
		filter_by_XA_mimi()

		parameter.sam_file_name = parameter.sam_file_name + "_XA"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "3. sorting",  parameter.sam_file

		cmd = "sort -k 3,3 -k 4,4n " + parameter.sam_file + " > " + parameter.sam_file_name + "_sorted.sam"
		os.system(cmd)

		parameter.sam_file_name = parameter.sam_file_name + "_sorted"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "4. repeat remove",  parameter.sam_file
		rmsk_file = "/home/guoxing/disk2/wli/rmsk/rmsk_" + parameter.chr_name + "_sorted.txt"
		print rmsk_file
		repeat_remove_fish_wli(rmsk_file, parameter.sam_file)

		parameter.sam_file_name = parameter.sam_file_name + "_rmsk"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "5. find matched pairend after repeat remove",  parameter.sam_file
		filter_match_pairend()

		print "6. extract_single_overlapped_read",  parameter.sam_file
		extract_single_overlapped_read(parameter.sam_file)

		parameter.sam_file_name = parameter.sam_file_name + "_combined"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "7. find matched pairend from combined file and process indel",  parameter.sam_file
		pair_end_indel(parameter.sam_file, parameter.chr_name)


		#snpPick_mimi -s NA12893_S1_ChrXnew_pairend_XA_sorted_rmsk_combined_indel.sam -c chrX -m update -d NA12893_S1_chrX

		print "8. clean up"
		"""
		# keep the pairend_XA_sorted.sam and pairend_XA_sorted_rmsk.sam
		#os.system("rm " + ori_sam_file_name + ".sam")
		os.system("rm " + ori_sam_file_name + "_pairend.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_removed.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA_sorted_record.txt")
		sorted_rmsk_name = ori_sam_file_name + "_pairend_XA_sorted_rmsk"
		os.system("rm " + sorted_rmsk_name + "_pairend*.sam")
		os.system("rm " + sorted_rmsk_name + "_recovered*.sam")
		os.system("rm " + sorted_rmsk_name + "_removed.sam")
		os.system("rm " + sorted_rmsk_name + "_combined.sam")
		"""

	elif mode == "solid_lima":
		# solid mimi process, single end, no XA_filter needed, hg18, already sorted.
		ori_sam_file_name = parameter.sam_file_name
		"""
		#parameter.sam_file_name = parameter.sam_file_name + "_sorted"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "1. repeat remove",  parameter.sam_file
		rmsk_file = "hg18_rmsk.txt_original"
		#repeat_remove(rmsk_file, parameter.sam_file)

		parameter.sam_file_name = parameter.sam_file_name + "_rmsk"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "2. process indel",  parameter.sam_file
		single_end_indel(parameter.sam_file, chr_name)
		"""
		parameter.sam_file_name = parameter.sam_file_name + "_rmsk_single_indel"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "3. remove base after primer",  parameter.sam_file
		mimi_solid_base_remove(parameter.sam_file)

		#snpPick_mimi -s NA12893_S1_ChrXnew_pairend_XA_sorted_rmsk_combined_indel.sam -c chrX -m update -d NA12893_S1_chrX
		"""
		print "8. clean up"
		# keep the pairend_XA_sorted.sam and pairend_XA_sorted_rmsk.sam
		os.system("rm " + ori_sam_file_name + ".sam")
		os.system("rm " + ori_sam_file_name + "_pairend.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_removed.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA_sorted_record.txt")
		sorted_rmsk_name = ori_sam_file_name + "_pairend_XA_sorted_rmsk"
		os.system("rm " + sorted_rmsk_name + "_pairend*.sam")
		os.system("rm " + sorted_rmsk_name + "_recovered*.sam")
		os.system("rm " + sorted_rmsk_name + "_removed.sam")
		os.system("rm " + sorted_rmsk_name + "_combined.sam")
		"""
	elif mode == "yang_mimi":

		ori_sam_file_name = parameter.sam_file_name
		print "1. filter by chr",  parameter.sam_file
		filter_by_chr()

		parameter.sam_file_name = parameter.sam_file_name + parameter.chr_name[3:]
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "2. filter by XA",  parameter.sam_file
		# singend data
		single_end_xa()

		parameter.sam_file_name = parameter.sam_file_name + "_XA"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "3. sorting",  parameter.sam_file
		#add_header(parameter.sam_file)
		samtools_sort(parameter.sam_file)

		parameter.sam_file_name = parameter.sam_file_name + "_sorted"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "4. repeat remove",  parameter.sam_file
		rmsk_file = "hg18_rmsk.txt_original"
		repeat_remove(rmsk_file, parameter.sam_file)

		parameter.sam_file_name = parameter.sam_file_name + "_rmsk"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "5. process indel",  parameter.sam_file
		single_end_indel(parameter.sam_file, parameter.chr_name)

		parameter.sam_file_name = parameter.sam_file_name + "_single_indel"
		parameter.sam_file = parameter.sam_file_name + ".sam"
		print "6. mv sam file",  parameter.sam_file
		os.system("cp " + parameter.sam_file + " ../mimi_yang_sam/ &")

		"""

		print "7. find matched pairend from combined file and process indel",  parameter.sam_file

		#snpPick_mimi -s NA12893_S1_ChrXnew_pairend_XA_sorted_rmsk_combined_indel.sam -c chrX -m update -d NA12893_S1_chrX

		print "8. clean up"
		# keep the pairend_XA_sorted.sam and pairend_XA_sorted_rmsk.sam
		os.system("rm " + ori_sam_file_name + ".sam")
		os.system("rm " + ori_sam_file_name + "_pairend.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_removed.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA.sam")
		os.system("rm " + ori_sam_file_name + "_pairend_XA_sorted_record.txt")
		sorted_rmsk_name = ori_sam_file_name + "_pairend_XA_sorted_rmsk"
		os.system("rm " + sorted_rmsk_name + "_pairend*.sam")
		os.system("rm " + sorted_rmsk_name + "_recovered*.sam")
		os.system("rm " + sorted_rmsk_name + "_removed.sam")
		os.system("rm " + sorted_rmsk_name + "_combined.sam")
		"""

def get_args():
	desc="variation call"
	usage = "sam_process"
	parser = OptionParser(usage = usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile", help="Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="chr")
	parser.add_option("-m", "--mode", type="string", dest="mode", help="", default="null")
	
	parser.add_option("-a", "--repeat", type="string", dest="repeatFile", help="", default="null")
	parser.add_option("-b", "--combined", type="string", dest="combinedFile", help="", default="null")
	(options, args) = parser.parse_args()
	if options.mode == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	
	start_time = time.time()
	global parameter
	parameter = parameters()
	
	options = get_args()
	parameter.mode = options.mode
	if parameter.mode == "com_repeat":
		# combine repeat files.
		# sam_process -a repeat_file -b already_combined_file
		repeat_cnv_file = options.repeatFile
		combined_file = options.combinedFile
		combine_cnv_repeat(repeat_cnv_file, combined_file)
	else:
		parameter.sam_file = options.samFile
		parameter.chr_name = options.chrName
		parameter.sam_file_name = parameter.sam_file[:(len(parameter.sam_file)-4)]
		sam_process(parameter.sam_file, parameter.chr_name, parameter.mode)
	print "run time is: ", round((time.time() - start_time), 3), "s"
	
	"""
	#combine_files()
	#seperate_by_chr(sam_file)
	#pair_end_indel(sam_file)
	pair_end_filter(sam_file)
	sam_file = sam_file_name + "_pairend.sam"
	indel_process(sam_file)
	"""

