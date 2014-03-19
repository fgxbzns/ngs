#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for solid data, to remove reads covered by repeats

import os,sys, glob, subprocess, random, operator, time
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

	overlap = 15

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

			if repeat_chr == read_chr :

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
	output_files(sam_file_name + "_removed.sam", removed_list)
	output_files(sam_file_name + "_record.txt", record_list)

def repeat_remove_mimi(rmsk_file, sam_file):
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

	#inputFile_rmsk = open(file_path + rmsk_file, "r")  # solid data
	inputFile_rmsk = open(currentPath + rmsk_file, "r")  # mimi data lima
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
					#outputFile_sam.write(read_line.strip() + "\n")
					rmsk_list.append(read_line)
					read_line = inputFile_sam.readline()
					total_reads_number += 1
					kept_reads_nmuber += 1
					#print "kept", repeat_start
				elif read_end > (repeat_start+overlap) and read_end <= (repeat_end+read_length-overlap) and (repeat_end-repeat_start) >= overlap: # need to consider length of repeat
					#print "removed", repeat_start
					removed_line = read_chr + "\t" + str(repeat_start) + "\t" + str(repeat_end) + "\t" + repeat_class + "\t" + read_line.strip()
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
	output_files(sam_file_name + "_removed.sam", removed_list)
	output_files(sam_file_name + "_record.txt", record_list)

def get_args():
	desc = "Compare seed and std hap, to check purity of seed"

	usage = "seed_std_compare -i seed_file -c chr#"
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-r", "--rmsk", type="string", dest="rmskFile", help="Input rmsk File Name", default="null")
	parser.add_option("-s", "--sam", type="string", dest="samFile", help="Input sam File Name", default="null")
	(options, args) = parser.parse_args()
	if options.samFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options


if __name__ == '__main__':
	start = time.time()

	options = get_args()
	rmsk_file = options.rmskFile
	sam_file = options.samFile
	#rmsk_file = "hg18_rmsk.txt"

	#repeat_remove(rmsk_file, sam_file)
	repeat_remove_mimi(rmsk_file, sam_file)







