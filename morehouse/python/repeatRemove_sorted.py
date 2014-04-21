#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for solid data, to remove reads covered by repeats

import os,sys, glob, subprocess, random, operator, time
from tools import *
from optparse import OptionParser

file_path = "/home/guoxing/disk2/solid/common_files/"

currentPath = os.getcwd() + '/'

# store repeat info from rmsk file
class repeat:
	def __init__(me, chr_name, repeat_start, repeat_end, matched_repeat, repeat_class, repeat_family):
		me.chr_name = chr_name
		me.repeat_start = repeat_start
		me.repeat_end = repeat_end
		me.matched_repeat = matched_repeat
		me.repeat_class = repeat_class
		me.repeat_family = repeat_family

def is_multiple_maping(elements):
	multiple_maping = False
	#XA = ""
	XA = elements[-1].strip()
	if XA.startswith('XA'):
		multiple_maping = True
	else:
		XA = ""
	return multiple_maping

def output_files(file_name, data_list):
	with open(file_name, "w") as output_file:
		for line in data_list:
			print >> output_file, line

def repeat_remove(rmsk_file, sam_file):
	# For removing repeat for solid data

	start = time.time()
	rmsk_list = []
	multmap_list = []
	removed_list = []
	record_list = []

	repeat_name = rmsk_file[10:(len(rmsk_file) - 4)]
	sam_file_name = sam_file[:(len(sam_file) - 4)]

	print "rmsk file: ", rmsk_file
	print "sam file: ", sam_file_name

	inputFile_rmsk = open(file_path + rmsk_file, "r")
	inputFile_sam = open(currentPath + sam_file, "r")

	total_reads_number = 0
	multiple_mapping_number = 0
	kept_reads_nmuber = 0
	removed_reads_number = 0

	#overlap = 15
	overlap = 25

	repeat_line = inputFile_rmsk.readline()

	read_line = inputFile_sam.readline()

	while repeat_line != '' and read_line != '':
		repeat_elements = repeat_line.strip().split()
		repeat_chr = repeat_elements[5].strip()

		read_line = read_line.strip()
		read_elements = read_line.split()
		read_chr = read_elements[2].strip()

		# remove reads with mulitple mapping
		multiple_maping = is_multiple_maping(read_elements)

		if not multiple_maping:
			if repeat_chr == "chrX":
				repeat_chr = "chr23"
			if repeat_chr == "chrY":
				repeat_chr = "chr24"
			if read_chr == "chrX":
				read_chr = "chr23"
			if read_chr == "chrY":
				read_chr = "chr24"

			if repeat_chr == read_chr:

				repeat_start = int(repeat_elements[6].strip())
				repeat_end = int(repeat_elements[7].strip())
				matched_repeat = repeat_elements[10].strip()
				repeat_class = repeat_elements[11].strip()
				repeat_family = repeat_elements[12].strip()

				read_start = int(read_elements[3].strip())
				read_seq = read_elements[9].strip()
				read_length = len(read_seq)
				read_end = read_start + read_length

				if read_end <= (repeat_start+overlap):
					#outputFile_sam.write(read_line.strip() + "\n")
					rmsk_list.append(read_line)
					read_line = inputFile_sam.readline()
					total_reads_number += 1
					kept_reads_nmuber += 1
					#print "kept", repeat_start
				elif read_end > (repeat_start+overlap) and read_end <= (repeat_end+read_length-overlap) and (repeat_end-repeat_start) >= overlap: # need to consider length of repeat
					#print "removed", repeat_start
					#outputFile_removed.write(read_chr + "\t" + str(repeat_start) + "\t" + str(repeat_end) + "\t" + repeat_class + "\t" + repeat_family + "\t" + read_line.strip() + "\n")
					removed_line = read_chr + "\t" + str(repeat_start) + "\t" + str(repeat_end) + "\t" + repeat_class + "\t" + repeat_family + "\t" + read_line.strip()
					removed_list.append(removed_line)
					read_line = inputFile_sam.readline()
					total_reads_number += 1
					removed_reads_number += 1
				else:
					repeat_line = inputFile_rmsk.readline()
			elif int(repeat_chr[3:]) < int(read_chr[3:]): # to match the repeat chr and read chr
				#print repeat_chr[3:], int(read_chr[3:])
				repeat_line = inputFile_rmsk.readline()
			else: # repeat is finished, more reads left. Keep them all
				#outputFile_sam.write(read_line.strip() + "\n")
				rmsk_list.append(read_line)
				read_line = inputFile_sam.readline()
				total_reads_number += 1
				kept_reads_nmuber += 1
		else:
			multmap_list.append(read_line)
			multiple_mapping_number += 1
			read_line = inputFile_sam.readline()
			total_reads_number += 1


	print "total_reads_number", total_reads_number
	print "multiple_mapping_number", multiple_mapping_number
	print "kept_reads_nmuber", kept_reads_nmuber
	print "removed_reads_number", removed_reads_number
	print "removed percentage", round(float(removed_reads_number+multiple_mapping_number)/total_reads_number, 3)
	print "kept reads percentage", round(float(kept_reads_nmuber)/total_reads_number, 3)



	end = time.time()
	run_time = str(end - start)
	run_time = run_time[:(run_time.find('.') + 3)]
	print "run time is: " + run_time + "s"

	record_list.append("total_reads_number: " + str(total_reads_number))
	record_list.append("multiple_mapping_number: " + str(multiple_mapping_number))
	record_list.append("kept_reads_nmuber: " + str(kept_reads_nmuber))
	record_list.append("removed_reads_number: " + str(removed_reads_number))
	record_list.append("run time is: " + str(run_time))

	output_files(sam_file_name + "_rmsk.sam", rmsk_list)
	output_files(sam_file_name + "_mulmap.sam", multmap_list)
	output_files(sam_file_name + "_rmsk_removed.sam", removed_list)
	output_files(sam_file_name + "_rmsk_record.txt", record_list)

def repeat_remove_mimi(rmsk_file, sam_file):
	# For removing repeat for mimi data
	# there might be a glitch here. If the last chr in repeat is finished first, will the
	# remaining reads be kept?

	start = time.time()
	rmsk_list = []
	multmap_list = []
	removed_list = []
	record_list = []

	repeat_name = rmsk_file[10:(len(rmsk_file) - 4)]
	sam_file_name = sam_file[:(len(sam_file) - 4)]

	print "rmsk file: ", rmsk_file
	print "sam file: ", sam_file_name

	#outputFile_mulmap = open(currentPath + sam_file_name + "_mulmap.sam", "w")
	outputFile_sam = open(currentPath + sam_file_name + "_rmsk.sam", "w")
	outputFile_removed = open(currentPath + sam_file_name + "_rmsk_removed.sam", "w")

	#inputFile_rmsk = open(file_path + rmsk_file, "r")  # solid data
	inputFile_rmsk = open(rmsk_file, "r")  # mimi data lima
	inputFile_sam = open(currentPath + sam_file, "r")

	total_reads_number = 0
	multiple_mapping_number = 0
	kept_reads_nmuber = 0
	removed_reads_number = 0

	overlap = 15

	repeat_line = inputFile_rmsk.readline() # skip title line
	repeat_line = inputFile_rmsk.readline()

	read_line = inputFile_sam.readline()

	while repeat_line != '' and read_line != '':
		repeat_elements = repeat_line.strip().split()
		repeat_chr = repeat_elements[2].strip()

		read_line = read_line.strip()
		read_elements = read_line.split()
		read_chr = read_elements[2].strip()

		# remove reads with multiple mapping
		multiple_maping = is_multiple_maping(read_elements)
		# multiple mapping is already checked with sam_process XA
		if True:
			if repeat_chr == "chrX":
				repeat_chr = "chr23"
			if repeat_chr == "chrY":
				repeat_chr = "chr24"
			if repeat_chr == "chrM":
				repeat_chr = "chr25"
			if read_chr == "chrX":
				read_chr = "chr23"
			if read_chr == "chrY":
				read_chr = "chr24"
			if read_chr == "chrM":
				read_chr = "chr25"
			if read_chr == "*":         # chrM is converted to * by samtools sort
				read_chr = "chr25"

			try:
				repeat_chr_int = int(repeat_chr[3:])
				read_chr_int = int(read_chr[3:])
			except:
				print repeat_chr, repeat_line
				print read_line, read_chr

			if repeat_chr == read_chr:

				repeat_start = int(repeat_elements[0].strip())
				repeat_end = int(repeat_elements[1].strip())
				#matched_repeat = repeat_elements[10].strip()
				repeat_class = repeat_elements[3].strip()
				#repeat_family = repeat_elements[12].strip()

				read_start = int(read_elements[3].strip())
				read_seq = read_elements[9].strip()
				read_length = len(read_seq)
				read_end = read_start + read_length

				if read_end <= (repeat_start+overlap):
					outputFile_sam.write(read_line.strip() + "\n")
					#rmsk_list.append(read_line)
					read_line = inputFile_sam.readline()
					total_reads_number += 1
					kept_reads_nmuber += 1
					#print "kept", repeat_start
				elif read_end > (repeat_start+overlap) and read_end <= (repeat_end+read_length-overlap) and (repeat_end-repeat_start) >= overlap: # need to consider length of repeat
					#print "removed", repeat_start
					# old output, include repeat info
					#removed_line = read_chr + "\t" + str(repeat_start) + "\t" + str(repeat_end) + "\t" + repeat_class + "\t" + read_line.strip()
					#outputFile_removed.write(removed_line + "\n")
					#removed_list.append(removed_line)

					# to recover read pair with only one ready overlaps with repeat.
					outputFile_removed.write(read_line.strip() + "\n")

					read_line = inputFile_sam.readline()
					total_reads_number += 1
					removed_reads_number += 1
				else:
					repeat_line = inputFile_rmsk.readline()
			elif int(repeat_chr[3:]) < int(read_chr[3:]): # to match the repeat chr and read chr
				#print repeat_chr[3:], int(read_chr[3:])
				repeat_line = inputFile_rmsk.readline()
			else: # repeat is finished, more reads left. Keep them all
				outputFile_sam.write(read_line.strip() + "\n")
				#rmsk_list.append(read_line)
				read_line = inputFile_sam.readline()
				total_reads_number += 1
				kept_reads_nmuber += 1
		else:
			#multmap_list.append(read_line)
			#outputFile_mulmap.write(read_line.strip() + "\n")
			multiple_mapping_number += 1
			read_line = inputFile_sam.readline()
			total_reads_number += 1


	print "total_reads_number", total_reads_number
	print "multiple_mapping_number", multiple_mapping_number
	print "removed_reads_number", removed_reads_number
	print "kept_reads_nmuber", kept_reads_nmuber
	print "percentage", round(float(kept_reads_nmuber)/total_reads_number, 3)

	end = time.time()
	run_time = str(end - start)
	run_time = run_time[:(run_time.find('.') + 3)]
	print "run time is: " + run_time + "s"

	record_list.append("total_reads_number: " + str(total_reads_number))
	record_list.append("multiple_mapping_number: " + str(multiple_mapping_number))
	record_list.append("kept_reads_nmuber: " + str(kept_reads_nmuber))
	record_list.append("removed_reads_number: " + str(removed_reads_number))
	record_list.append("run time is: " + str(run_time))

	# takes time if the file is big
	#output_files(sam_file_name + "_rmsk.sam", rmsk_list)
	#output_files(sam_file_name + "_mulmap.sam", multmap_list)
	#output_files(sam_file_name + "_removed.sam", removed_list)
	output_files(sam_file_name + "_record.txt", record_list)

	outputFile_sam.close()
	outputFile_removed.close()
	#outputFile_mulmap.close()

def repeat_remove_mimi_solid(rmsk_file, sam_file):
	# For removing repeat for solid mimi data
	# for solid data, remove rmsk_hg18 only, do not combine with others
	# non-repeat region, >= 7, overlap > 25

	start = time.time()
	rmsk_list = []
	multmap_list = []
	removed_list = []
	record_list = []

	repeat_name = rmsk_file[10:(len(rmsk_file) - 4)]
	sam_file_name = sam_file[:(len(sam_file) - 4)]

	print "rmsk file: ", rmsk_file
	print "sam file: ", sam_file_name

	#outputFile_mulmap = open(currentPath + sam_file_name + "_mulmap.sam", "w")
	outputFile_sam = open(currentPath + sam_file_name + "_rmsk.sam", "w")
	outputFile_removed = open(currentPath + sam_file_name + "_rmsk_removed.sam", "w")

	#inputFile_rmsk = open(file_path + rmsk_file, "r")  # solid data
	inputFile_rmsk = open(rmsk_file, "r")  # mimi data lima
	inputFile_sam = open(currentPath + sam_file, "r")

	total_reads_number = 0
	multiple_mapping_number = 0
	kept_reads_nmuber = 0
	removed_reads_number = 0

	overlap = 25

	repeat_line = inputFile_rmsk.readline() # skip title line
	repeat_line = inputFile_rmsk.readline()

	read_line = inputFile_sam.readline()

	while repeat_line != '' and read_line != '':
		repeat_elements = repeat_line.strip().split()
		repeat_chr = repeat_elements[2].strip()

		read_line = read_line.strip()
		read_elements = read_line.split()
		read_chr = read_elements[2].strip()

		# remove reads with multiple mapping
		multiple_maping = is_multiple_maping(read_elements)
		# multiple mapping is already checked with sam_process XA
		if True:
			if repeat_chr == "chrX":
				repeat_chr = "chr23"
			if repeat_chr == "chrY":
				repeat_chr = "chr24"
			if repeat_chr == "chrM":
				repeat_chr = "chr25"
			if read_chr == "chrX":
				read_chr = "chr23"
			if read_chr == "chrY":
				read_chr = "chr24"
			if read_chr == "chrM":
				read_chr = "chr25"
			if read_chr == "*":         # chrM is converted to * by samtools sort
				read_chr = "chr25"

			try:
				repeat_chr_int = int(repeat_chr[3:])
				read_chr_int = int(read_chr[3:])
			except:
				print repeat_chr, repeat_line
				print read_line, read_chr

			if repeat_chr == read_chr:

				repeat_start = int(repeat_elements[0].strip())
				repeat_end = int(repeat_elements[1].strip())
				#matched_repeat = repeat_elements[10].strip()
				repeat_class = repeat_elements[3].strip()
				#repeat_family = repeat_elements[12].strip()

				read_start = int(read_elements[3].strip())
				read_seq = read_elements[9].strip()
				read_length = len(read_seq)
				read_end = read_start + read_length

				if read_end <= (repeat_start + overlap):
					outputFile_sam.write(read_line.strip() + "\n")
					#rmsk_list.append(read_line)
					read_line = inputFile_sam.readline()
					total_reads_number += 1
					kept_reads_nmuber += 1
					#print "kept", repeat_start
				elif read_end > (repeat_start + overlap) and read_end <= (repeat_end + read_length - overlap) \
						and (repeat_end - repeat_start) >= overlap: # need to consider length of repeat
					#print "removed", repeat_start
					# old output, include repeat info
					#removed_line = read_chr + "\t" + str(repeat_start) + "\t" + str(repeat_end) + "\t" + repeat_class + "\t" + read_line.strip()
					#outputFile_removed.write(removed_line + "\n")
					#removed_list.append(removed_line)

					# to recover read pair with only one ready overlaps with repeat.
					outputFile_removed.write(read_line.strip() + "\n")

					read_line = inputFile_sam.readline()
					total_reads_number += 1
					removed_reads_number += 1
				else:
					repeat_line = inputFile_rmsk.readline()
			elif int(repeat_chr[3:]) < int(read_chr[3:]): # to match the repeat chr and read chr
				#print repeat_chr[3:], int(read_chr[3:])
				repeat_line = inputFile_rmsk.readline()
			else: # repeat is finished, more reads left. Keep them all
				outputFile_sam.write(read_line.strip() + "\n")
				#rmsk_list.append(read_line)
				read_line = inputFile_sam.readline()
				total_reads_number += 1
				kept_reads_nmuber += 1
		else:
			#multmap_list.append(read_line)
			#outputFile_mulmap.write(read_line.strip() + "\n")
			multiple_mapping_number += 1
			read_line = inputFile_sam.readline()
			total_reads_number += 1


	print "total_reads_number", total_reads_number
	print "multiple_mapping_number", multiple_mapping_number
	print "removed_reads_number", removed_reads_number
	print "kept_reads_nmuber", kept_reads_nmuber
	print "kept reads percentage", round(float(kept_reads_nmuber)/total_reads_number, 3)

	end = time.time()
	run_time = str(end - start)
	run_time = run_time[:(run_time.find('.') + 3)]
	print "run time is: " + run_time + "s"

	record_list.append("total_reads_number: " + str(total_reads_number))
	record_list.append("multiple_mapping_number: " + str(multiple_mapping_number))
	record_list.append("kept_reads_nmuber: " + str(kept_reads_nmuber))
	record_list.append("removed_reads_number: " + str(removed_reads_number))
	record_list.append("run time is: " + str(run_time))

	# takes time if the file is big
	#output_files(sam_file_name + "_rmsk.sam", rmsk_list)
	#output_files(sam_file_name + "_mulmap.sam", multmap_list)
	#output_files(sam_file_name + "_removed.sam", removed_list)
	output_files(sam_file_name + "_record.txt", record_list)

	outputFile_sam.close()
	outputFile_removed.close()
	#outputFile_mulmap.close()

def load_mimi_data(file_name):
	# used in mimi_combine_files
	data = {}
	with open(file_name, "r") as fp:
		for line in fp:
			elements = line.strip().split()
			try:
				data[int(elements[0])] = elements[:7]
			except:
				#print "error in ", line, file_name
				pass
	return data

def mimi_combine_files():
	# combine genotype from variation call
	file_number_dict = {}
	for infile in glob.glob(os.path.join("NA128??_S1_chrX_0_158375978_filtered.txt")):
		mimi_file_number = int(infile[2:7].strip())
		#print mimi_file_number
		file_number_dict[mimi_file_number] = {}

	title_line = "pos \t"
	second_line = "pos \t"
	all_snp_dict = {}
	for mimi_file_number in file_number_dict.keys():
		file_name = "NA" + str(mimi_file_number) + "_S1_chrX_0_158375978_filtered.txt"
		file_number_dict[mimi_file_number] = load_mimi_data(file_name)
		print "NA" + str(mimi_file_number), len(file_number_dict[mimi_file_number])
		title_line += "NA" + str(mimi_file_number) + "\t\t\t\t"
		second_line += "A \t T \t C \t G \t"
		for pos in file_number_dict[mimi_file_number]:
			if pos not in all_snp_dict:
				all_snp_dict[pos] = ""
	print "total pos number: ", len(all_snp_dict)
	#print title_line

	for pos in all_snp_dict.keys():
		for pos_data in file_number_dict.values():
			if pos in pos_data:
				#print pos, pos_data[pos]
				all_snp_dict[pos] = all_snp_dict[pos] + list_to_line(pos_data[pos][3:]) + "\t"
			else:
				all_snp_dict[pos] += "\t\t\t\t"
				#pass
	#print all_snp_dict[3343118]

	with open("combined_mimi.txt", "w") as fp:
		print >> fp, title_line
		print >> fp, second_line
		for pos, data in all_snp_dict.iteritems():
			print >> fp, str(pos) + "\t" + data

def map_position_repeat_combined(rmsk_file, pos_file):
	# use -s samfile to input pos_file
	# this is used to map the pos to the combined repeat file.
	pos_list = []
	with open(pos_file, "r") as pos_fp:
		for line in pos_fp:
			try:
				pos = int(line.strip().split()[0])
				pos_list.append(pos)
			except:
				print "error in map_position_repeat", pos
	pos_list_ori_prder = list(pos_list)
	pos_list.sort()
	pos_dict = {pos: "" for pos in pos_list}

	inputFile_rmsk = open(rmsk_file, "r")  # mimi data lima

	mapped_to_repeat_nmuber = 0
	not_mapped_to_repeat_nmuber = 0

	repeat_line = inputFile_rmsk.readline() # skip title line
	repeat_line = inputFile_rmsk.readline()

	i = 0
	while repeat_line != '' and i < len(pos_list):
		# only if the snp pos is mapped to a repeat, the repeat info will be added
		# if the repeat finish earlier, the remaining pos will be output as unmapped.
		repeat_elements = repeat_line.strip().split()
		repeat_chr = repeat_elements[2].strip()

		snp_pos = pos_list[i]

		# chr is not checked here. Each time check only one chr
		repeat_start = int(repeat_elements[0].strip())
		repeat_end = int(repeat_elements[1].strip())
		repeat_class = repeat_elements[3].strip()

		if snp_pos < repeat_start:
			i += 1
			not_mapped_to_repeat_nmuber += 1
		elif snp_pos >= repeat_start and snp_pos <= repeat_end:
			mapped_to_repeat_nmuber += 1
			pos_dict[snp_pos] = str(repeat_start) + "\t" + str(repeat_end) + "\t" + repeat_class
			i += 1
			#print snp_pos, pos_dict[snp_pos]
		elif snp_pos > repeat_end:
			repeat_line = inputFile_rmsk.readline()

	print "len(pos_list)", len(pos_list)
	print "mapped_to_repeat_nmuber", mapped_to_repeat_nmuber
	print "not_mapped_to_repeat_nmuber", not_mapped_to_repeat_nmuber
	print "percentage", round(float(mapped_to_repeat_nmuber)/len(pos_list), 3)
	"""
	with open("mimi_pos_mapped.txt", "w") as output_file:
		for pos in pos_list_ori_prder:
			print >> output_file, pos + "\t" + pos_dict[pos]
	"""
	with open(pos_file[:len(pos_file)-4]+"_mapped_combined.txt", "w") as output_file:
		for pos in pos_list:
			print >> output_file, pos + "\t" + pos_dict[pos]

def map_position_repeat_rmsk(rmsk_file, pos_file):
	# use -s samfile to input pos_file
	# this is used to map the pos to the rmsk_chrX_hg19.txt. After snpPick_mimi
	pos_list = []
	with open(pos_file, "r") as pos_fp:
		for line in pos_fp:
			try:
				pos = int(line.strip().split()[0])
				pos_list.append(pos)
			except:
				print "error in map_position_repeat", pos
	pos_list_ori_prder = list(pos_list)
	pos_list.sort()
	pos_dict = {pos: "" for pos in pos_list}

	inputFile_rmsk = open(rmsk_file, "r")  # mimi data lima

	mapped_to_repeat_nmuber = 0
	not_mapped_to_repeat_nmuber = 0

	repeat_line = inputFile_rmsk.readline() # skip title line
	repeat_line = inputFile_rmsk.readline()

	i = 0
	while repeat_line != '' and i < len(pos_list):
		# only if the snp pos is mapped to a repeat, the repeat info will be added
		# if the repeat finish earlier, the remaining pos will be output as unmapped.
		repeat_elements = repeat_line.strip().split()
		repeat_chr = repeat_elements[2].strip()

		snp_pos = pos_list[i]

		# chr is not checked here. Each time check only one chr
		repeat_start = int(repeat_elements[6].strip())
		repeat_end = int(repeat_elements[7].strip())
		repeat_class = repeat_elements[10].strip() + "\t" + repeat_elements[11].strip() + "\t" + repeat_elements[12].strip()

		if snp_pos < repeat_start:
			i += 1
			not_mapped_to_repeat_nmuber += 1
		elif snp_pos >= repeat_start and snp_pos <= repeat_end:
			mapped_to_repeat_nmuber += 1
			pos_dict[snp_pos] = str(repeat_start) + "\t" + str(repeat_end) + "\t" + repeat_class
			i += 1
			#print snp_pos, pos_dict[snp_pos]
		elif snp_pos > repeat_end:
			repeat_line = inputFile_rmsk.readline()

	print "len(pos_list)", len(pos_list)
	print "mapped_to_repeat_nmuber", mapped_to_repeat_nmuber
	print "not_mapped_to_repeat_nmuber", not_mapped_to_repeat_nmuber
	print "percentage", round(float(mapped_to_repeat_nmuber)/len(pos_list), 3)

	with open(pos_file[:len(pos_file)-4]+"_mapped_rmsk.txt", "w") as output_file:
		for pos in pos_list_ori_prder:
			print >> output_file, pos + "\t" + pos_dict[pos]

	with open(pos_file[:len(pos_file)-4]+"_mapped_ordered_rmsk.txt", "w") as output_file:
		for pos in pos_list:
			print >> output_file, pos + "\t" + pos_dict[pos]

def extract_single_overlapped_read(sam_file):
	# to extract_single_overlapped_read from pairend_removed.sam file

	sam_file_name = sam_file[:(len(sam_file) - 4)]
	rmsk_file_name = sam_file
	rmsk_pairend_file_name = sam_file_name + "_pairend.sam"
	rmsk_removed = sam_file_name + "_removed.sam"
	rmsk_pairend_removed = sam_file_name + "_pairend_removed.sam"
	recoverd_overl_read_name = sam_file_name + "_recovered_overlap_read.sam"

	rmsk_pairend_rm_dict = {}
	with open(rmsk_pairend_removed, "r") as rmsk_pairend_rm_file:
		for read in rmsk_pairend_rm_file:
			rmsk_pairend_rm_dict[read.strip().split()[0]] = read.strip()
	print "rmsk_removed_pairend_removed size", len(rmsk_pairend_rm_dict)

	with open(recoverd_overl_read_name, "w") as recovered_file:
		with open(rmsk_removed, "r") as rmsk_removed_file:
			for read in rmsk_removed_file:
				read_id = read.strip().split()[0]
				if len(rmsk_pairend_rm_dict) > 0:
					if read_id in rmsk_pairend_rm_dict:
						print >> recovered_file, read.strip()
						print >> recovered_file, rmsk_pairend_rm_dict[read_id]
						del rmsk_pairend_rm_dict[read_id]
				else:
					break

	combined_file_name = sam_file_name + "_combined.sam"
	cmd = "cat " + recoverd_overl_read_name + " " + rmsk_pairend_file_name + " > " + combined_file_name
	os.system(cmd)

def get_args():
	desc = "Compare seed and std hap, to check purity of seed"

	usage = "repeatRemove_sorted -s sam_file -r repeat_file"
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-r", "--rmsk", type="string", dest="rmskFile", help="Input rmsk File Name", default="null")
	parser.add_option("-s", "--sam", type="string", dest="samFile", help="Input sam File Name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode", help="", default="null")

	(options, args) = parser.parse_args()
	if options.samFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	start_time = time.time()

	options = get_args()
	rmsk_file = options.rmskFile
	sam_file = options.samFile
	#rmsk_file = "hg18_rmsk.txt"
	mode = options.mode

	if mode == "solid":
		repeat_remove(rmsk_file, sam_file)
	elif mode == "mimi_remove":
		rmsk_file = "/home/guoxing/disk2/lima/rmsk_chrX_hg19_MultSNPs_chrX_SegDups_chrX.txt"
		repeat_remove_mimi(rmsk_file, sam_file)
		extract_single_overlapped_read(sam_file)
	elif mode == "mimi_map_rmsk":
		# map to rmsk_chrX_hg19.txt
		rmsk_file = "/home/guoxing/disk2/lima/repeat_chrx/rmsk_chrX_hg19.txt"
		map_position_repeat_rmsk(rmsk_file, sam_file)
	elif mode == "mimi_map_combined":
		# map to combined repeat file
		rmsk_file = "/home/guoxing/disk2/lima/rmsk_chrX_hg19_MultSNPs_chrX_SegDups_chrX.txt"
		map_position_repeat_combined(rmsk_file, sam_file)
	elif mode == "mimi_snp_combine":
		mimi_combine_files()
	elif mode == "mimi_solid":
		#rmsk_file = "hg18_rmsk.txt"
		rmsk_file = "hg18_rmsk.txt_original"
		repeat_remove(rmsk_file, sam_file)

	elapse_time = time.time() - start_time
	print "run time: ", round(elapse_time, 3), "s"







