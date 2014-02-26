#!/usr/bin/python

# remove entries in rmsk reference other than "Alu" and "Simple_repeat"

import os, glob, subprocess, random, operator, time
from optparse import OptionParser
from tools import *


class parameters:
	def __init__(self):
		self.ori_rmsk_file = ""
		self.repeat_rmsk_file = ""

		self.sam_file = ""

		self.record_file_name = ""

	def update(self):
		self.repeat_rmsk_file_name = self.repeat_rmsk_file[:(len(self.repeat_rmsk_file) - 4)]
		self.sam_file_name = self.sam_file[:(len(self.sam_file) - 4)]


	#outputFile_record = open(currentPath + self.sam_file + "_record.txt", "w")


def get_repeat_type(ori_rmsk_file):
	repeat_type_dict = {}
	with open(file_path + ori_rmsk_file, "r") as fp:
		for line in fp:
			try:
				elements = line.strip().split()
				type = elements[12]
				if type not in repeat_type_dict:
					repeat_type_dict[type] = ""
			except:
				print "error in", line
	print "repeat type:", len(repeat_type_dict)
	return repeat_type_dict


def output_repeat_file(ori_rmsk_file, output_repeat):
	output_repeat_name = "hg18_rmsk_" + output_repeat + ".txt"
	print "rmsk file: ", output_repeat_name
	ori_rmsk_total_number = 0
	processed_rmsk_total_number = 0
	with open(output_repeat_name, "w") as output_rmsk:
		with open(file_path + ori_rmsk_file, "r") as ori_rmsk:
			for line in ori_rmsk:
				ori_rmsk_total_number += 1
				if output_repeat in line:
					processed_rmsk_total_number += 1
					print >> output_rmsk, line.strip()
	print "ori_rmsk_total_number: ", ori_rmsk_total_number
	print "processed_rmsk_total_number: ", processed_rmsk_total_number


def repeatRemove_sorted():
	sam_file = parameter.sam_file
	rmsk_file = parameter.repeat_rmsk_file
	sam_file_name = parameter.sam_file_name

	print "rmsk file: ", parameter.repeat_rmsk_file_name
	print "sam file: ", sam_file_name

	inputFile_rmsk = open(currentPath + rmsk_file, "r")
	inputFile_sam = open(currentPath + sam_file, "r")

	outputFile_mulmap = open(currentPath + sam_file_name + "_mulmap.sam", "w")
	outputFile_sam = open(currentPath + sam_file_name + "_rmsk.sam", "w")
	#outputFile_removed = open(currentPath + sam_file_name + "_removed.sam", "w")
	outputFile_record = open(currentPath + sam_file_name + "_record.txt", "w")

	total_reads_number = 0
	multiple_mapping_number = 0
	kept_reads_nmuber = 0
	removed_reads_number = 0

	overlap = 15

	repeat_line = inputFile_rmsk.readline()

	read_line = inputFile_sam.readline()
	multiple_maping = False

	while repeat_line != '' and read_line != '':
		repeat_elements = repeat_line.strip().split()
		repeat_chr = repeat_elements[5].strip()

		read_elements = read_line.strip().split()
		read_chr = read_elements[2].strip()

		# remove reads with mulitple mapping
		try:
			XA = read_elements[19].strip()
			#print XA
			multiple_mapping_number += 1
			multiple_maping = True
		except:
			multiple_maping = False

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
				keep_this_read = True

				repeat_start = int(repeat_elements[6].strip())
				repeat_end = int(repeat_elements[7].strip())
				matched_repeat = repeat_elements[10].strip()
				repeat_class = repeat_elements[11].strip()
				repeat_family = repeat_elements[12].strip()

				read_start = int(read_elements[3].strip())
				read_seq = read_elements[9].strip()
				read_length = len(read_seq)
				read_end = read_start + read_length

				if read_end <= (repeat_start + overlap):
					outputFile_sam.write(read_line.strip() + "\n")
					read_line = inputFile_sam.readline()
					total_reads_number += 1
					kept_reads_nmuber += 1
				#print "kept", repeat_start
				elif read_end > (repeat_start + overlap) and read_end <= (repeat_end + read_length - overlap) and (
					repeat_end - repeat_start) >= overlap:  # need to consider length of repeat
					keep_this_read = False
					#print "removed", repeat_start
					#outputFile_removed.write(read_chr + "\t" + str(repeat_start) + "\t" + str(repeat_end) + "\t" + repeat_class + "\t" + repeat_family + "\t" + read_line.strip() + "\n")
					read_line = inputFile_sam.readline()
					total_reads_number += 1
					removed_reads_number += 1
				else:
					repeat_line = inputFile_rmsk.readline()
			elif int(repeat_chr[3:]) < int(read_chr[3:]):  # to match the repeat chr and read chr
				#print repeat_chr[3:], int(read_chr[3:])
				repeat_line = inputFile_rmsk.readline()
			else:  # repeat is finished, more reads left. Keep them all
				outputFile_sam.write(read_line.strip() + "\n")
				read_line = inputFile_sam.readline()
				total_reads_number += 1
				kept_reads_nmuber += 1
		else:
			outputFile_mulmap.write(read_line.strip() + "\n")
			multiple_maping = False
			read_line = inputFile_sam.readline()
			total_reads_number += 1

	print "total_reads_number", total_reads_number
	print "multiple_mapping_number", multiple_mapping_number
	print "kept_reads_nmuber", kept_reads_nmuber
	print "removed_reads_number", removed_reads_number

	outputFile_sam.close()
	#outputFile_removed.close()
	outputFile_record.close()
	outputFile_mulmap.close()

	pass


def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--rmsk", type="string", dest="rmskFile", help="Input rmsk File Name", default="null")
	parser.add_option("-s", "--sam", type="string", dest="samFile", help="Input sam File Name", default="null")
	(options, args) = parser.parse_args()
	if options.rmskFile == "rmskFile":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options


if __name__ == '__main__':
	start_time = time.time()
	global parameter
	parameter = parameters()
	"""
	# for split repeat types
	ori_rmsk_file = "hg18_rmsk.txt_original"
	repeat_type_dict = get_repeat_type(ori_rmsk_file)
	for repeat in repeat_type_dict.keys():
		output_repeat_file(ori_rmsk_file, repeat)
	"""
	options = get_args()

	parameter.repeat_rmsk_file = options.rmskFile
	parameter.sam_file = options.samFile
	parameter.update()
	print parameter.repeat_rmsk_file_name, parameter.sam_file_name
	print parameter.sam_file

	a = parameter.sam_file
	a = 0
	print parameter.sam_file

	print "run time is: ", round((time.time() - start_time), 3), "s"




