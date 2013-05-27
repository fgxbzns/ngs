#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for simulation NGS data

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

rmsk_path = "/home/guoxing/disk2/repeatMasker/"
file_path = "/home/guoxing/disk2/UCSC_hg18_index_lm/"

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

# start time
start = time.time()		
		
# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-r", "--rmsk", type="string", dest="rmskFile",help = "Input rmsk File Name", default="null")
parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input sam File Name", default="null")

(options, args) = parser.parse_args()

rmsk_file = options.rmskFile
sam_file = options.samFile

#if haplotype_file == "null":
	#print "Please input the haplotype file name"

rmsk_file = "hg18_rmsk.txt"
sam_file = "song_5_2_noPri_chr.sam"
sam_file = "song_5_2_removePri_chr.sam"


if sam_file == "null":
	print "Please input the sam file name"
	#sam_file = ".sam"

sam_file_name = sam_file[:(len(sam_file)-4)]

print "rmsk file: ", rmsk_file
print "sam file: ", sam_file_name

inputFile_rmsk = open(rmsk_path + rmsk_file, "r")
inputFile_sam = open(file_path + sam_file, "r")

outputFile_sam = open(currentPath + sam_file_name + "_kept.sam", "w")
outputFile_removed = open(currentPath + sam_file_name + "_removed.sam", "w")


# get rmsk repeat info
chr_dict={}

"""
for line in inputFile_rmsk:
	if True:
		elements = line.strip().split()
		chr_name = elements[5].strip()
		start = int(elements[6].strip())
		end = int(elements[7].strip())
		matched_repeat = elements[10].strip()
		repeat_class = elements[11].strip()
		repeat_family = elements[12].strip()
		if not chr_name in chr_dict:
			chr_dict[chr_name] = {}
		if not start in chr_dict[chr_name]:
			chr_dict[chr_name][start] = repeat(chr_name, start, end, matched_repeat, repeat_class, repeat_family)
		else:
			print start, " already exists in ", chr_name, " the longer repeat is kept."
			if (end - start) > (chr_dict[chr_name][start].end - chr_dict[chr_name][start].start):
				chr_dict[chr_name][start] = repeat(chr_name, start, end, matched_repeat, repeat_class, repeat_family)
"""
chr_start_dict = {}
chr_end_dict = {}
total_reads_number = 0
kept_reads_nmuber = 0
removed_reads_number = 0

for line in inputFile_rmsk:
	if True:
		elements = line.strip().split()
		chr_name = elements[5].strip()
		repeat_start = int(elements[6].strip())
		repeat_end = int(elements[7].strip())
		matched_repeat = elements[10].strip()
		repeat_class = elements[11].strip()
		repeat_family = elements[12].strip()
		# start postion dict
		if not chr_name in chr_start_dict:
			chr_start_dict[chr_name] = {}
		if not repeat_start in chr_start_dict[chr_name]:
			chr_start_dict[chr_name][repeat_start] = repeat(chr_name, repeat_start, repeat_end, matched_repeat, repeat_class, repeat_family)
		else:
			#print "repeat start ", repeat_start, " already exists in ", chr_name, " the longer repeat is kept."
			if (repeat_end - repeat_start) > (chr_start_dict[chr_name][repeat_start].repeat_end - chr_start_dict[chr_name][repeat_start].repeat_start):
				chr_start_dict[chr_name][repeat_start] = repeat(chr_name, repeat_start, repeat_end, matched_repeat, repeat_class, repeat_family)
		# end position dict
		if not chr_name in chr_end_dict:
			chr_end_dict[chr_name] = {}
		if not repeat_end in chr_end_dict[chr_name]:
			chr_end_dict[chr_name][repeat_end] = chr_start_dict[chr_name][repeat_start]
		else:
			#print "repeat end ", repeat_end, " already exists in ", chr_name, " the longer repeat is kept."
			if (repeat_end - repeat_start) > (chr_end_dict[chr_name][repeat_end].repeat_end - chr_end_dict[chr_name][repeat_end].repeat_start):
				chr_end_dict[chr_name][repeat_end] = chr_start_dict[chr_name][repeat_start]

# find max length in each chr
for chr_name, repeat_dict in chr_start_dict.iteritems():
	max_length = 0
	for repeat_start, repeat in chr_start_dict[chr_name].iteritems():
		if (chr_start_dict[chr_name][repeat_start].repeat_end - chr_start_dict[chr_name][repeat_start].repeat_start) > max_length:
			max_length = (chr_start_dict[chr_name][repeat_start].repeat_end - chr_start_dict[chr_name][repeat_start].repeat_start)
	chr_start_dict[chr_name]["max_length"] = max_length
	print chr_name, chr_start_dict[chr_name]["max_length"]
		
				
inputFile_rmsk.close()		

print len(chr_start_dict)
print (chr_start_dict['chr1'][168198715].matched_repeat)
print (chr_start_dict['chr2'][225229171].matched_repeat)

print len(chr_end_dict)

overlap = 15
for line in inputFile_sam:
	if not line.startswith("@"):
		total_reads_number += 1
		keep_this_read = True
		elements = line.strip().split()
		try:
			read_chr = elements[2].strip()
			read_start = int(elements[3].strip())
			read_seq = elements[9].strip()
			read_length = len(read_seq)
			read_end = read_start + read_length
		except:
			print "error in read:", line
		# check reads start
		i = 0
		while i < (read_length - overlap):
			if (read_start+i) in chr_start_dict[read_chr] and chr_start_dict[read_chr][read_start+i].repeat_end >= read_end:
				keep_this_read = False
				removed_reads_number += 1
				current_repeat = chr_start_dict[read_chr][read_start+i]
				outputFile_removed.write("start: \t" + read_chr + "\t" + str(current_repeat.repeat_start) + "\t" + str(current_repeat.repeat_end) + "\t" + current_repeat.repeat_class + "\t" + current_repeat.repeat_family + "\t" + line.strip() + "\n")
				break
			i += 1
		# check reads end
		if keep_this_read:
			i = 0
			while i < (read_length - overlap):
				if (read_end-i) in chr_end_dict[read_chr] and chr_end_dict[read_chr][read_end-i].repeat_start <= read_start:
					keep_this_read = False
					removed_reads_number += 1
					current_repeat = chr_end_dict[read_chr][read_end-i]
					outputFile_removed.write("end: \t" + read_chr + "\t" + str(current_repeat.repeat_start) + "\t" + str(current_repeat.repeat_end) + "\t" + current_repeat.repeat_class + "\t" + current_repeat.repeat_family + "\t" + line.strip() + "\n")
					break
				i += 1
		# check repeat that covers the whole read
		if keep_this_read:
			i = 0
			while i < (int(chr_start_dict[read_chr]["max_length"]) - read_length):
				if (read_start-i) > 0 and (read_start-i) in chr_start_dict[read_chr] and chr_start_dict[read_chr][read_start-i].repeat_end >= read_end:
					keep_this_read = False
					removed_reads_number += 1
					current_repeat = chr_start_dict[read_chr][read_start-i]
					outputFile_removed.write("cover all: \t" + read_chr + "\t" + str(current_repeat.repeat_start) + "\t" + str(current_repeat.repeat_end) + "\t" + current_repeat.repeat_class + "\t" + current_repeat.repeat_family + "\t" + line.strip() + "\n")
					break
				i += 1		
		if keep_this_read:
			kept_reads_nmuber += 1
			outputFile_sam.write(line.strip() + "\n")
	else:
		outputFile_sam.write(line.strip() + "\n")
		
print "total_reads_number", total_reads_number
print "kept_reads_nmuber", kept_reads_nmuber
print "removed_reads_number", removed_reads_number

outputFile_sam.close()	
outputFile_removed.close()	

















