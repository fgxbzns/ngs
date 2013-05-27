#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for simulation NGS data

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

rmsk_path = "/home/guoxing/disk2/dgv_struct_val/"
file_path = "/home/guoxing/disk2/UCSC_hg18_index_lm/" # for repeatRemove, after bwa
file_path = "/home/guoxing/disk2/repeatMasker/remove_rep_v2/" # for dgv

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

rmsk_file = "dgv_hg18.txt"
sam_file = "song_5_2_noPri_chr.sam"
sam_file = "song_5_2_removePri_chr.sam"
sam_file = "song_5_2_removePri_chr_sorted_kept.sam"
#sam_file = "song_5_2_noPri_chr_sorted_kept.sam"



if sam_file == "null":
	print "Please input the sam file name"
	#sam_file = ".sam"

sam_file_name = sam_file[:(len(sam_file)-4)]

print "rmsk file: ", rmsk_file
print "sam file: ", sam_file_name

inputFile_rmsk = open(rmsk_path + rmsk_file, "r")
inputFile_sam = open(file_path + sam_file, "r")


#outputFile_mulmap = open(currentPath + sam_file_name + "_mulmap.sam", "w")
outputFile_sam = open(currentPath + sam_file_name + "_dgv.sam", "w")
outputFile_removed = open(currentPath + sam_file_name + "_removed.sam", "w")
outputFile_record = open(currentPath + sam_file_name + "_record.txt", "w")


total_reads_number = 0
multiple_mapping_number = 0
kept_reads_nmuber = 0
removed_reads_number = 0

overlap = 15

repeat_line = inputFile_rmsk.readline()

read_line = inputFile_sam.readline()
multiple_maping = False
#read_elements = read_line.strip().split()
#read_chr = read_elements[2].strip()
#print read_chr


while repeat_line != '' and read_line != '':
	repeat_elements = repeat_line.strip().split()
	repeat_chr = repeat_elements[1].strip()
	
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
		"""
		# for dgv chr5	
		if 	len(repeat_chr) > 5:
			repeat_chr = repeat_chr[:repeat_chr.find('_')]
		"""	
		
		if repeat_chr == read_chr :
			keep_this_read = True
			
			repeat_start = int(repeat_elements[2].strip())
			repeat_end = int(repeat_elements[3].strip())
			matched_repeat = repeat_elements[4].strip()
			repeat_class = repeat_elements[5].strip()
			repeat_family = repeat_elements[10].strip()
			
			read_start = int(read_elements[3].strip())
			read_seq = read_elements[9].strip()
			read_length = len(read_seq)
			read_end = read_start + read_length
			
			if read_end <= (repeat_start+overlap):
				outputFile_sam.write(read_line.strip() + "\n")
				read_line = inputFile_sam.readline()
				total_reads_number += 1
				kept_reads_nmuber += 1
				#print "kept", repeat_start			
			elif read_end > (repeat_start+overlap) and read_end <= (repeat_end+read_length-overlap) and (repeat_end-repeat_start) >= overlap: # need to consider length of repeat
				keep_this_read = False
				#print "removed", repeat_start
				outputFile_removed.write(read_chr + "\t" + str(repeat_start) + "\t" + str(repeat_end) + "\t" + repeat_class + "\t" + repeat_family + "\t" + read_line.strip() + "\n")
				read_line = inputFile_sam.readline()
				total_reads_number += 1
				removed_reads_number += 1	
			else:
				repeat_line = inputFile_rmsk.readline()
		elif len(repeat_chr) > 5:	# skip repeat such as chr5_h2_hap1, chr1_random . etc
			repeat_line = inputFile_rmsk.readline()
		elif int(repeat_chr[3:]) < int(read_chr[3:]): # to match the repeat chr and read chr
			#print repeat_chr[3:], int(read_chr[3:])
			repeat_line = inputFile_rmsk.readline()		
		else: # repeat is finished, more reads left. Keep them all
			outputFile_sam.write(read_line.strip() + "\n")
			read_line = inputFile_sam.readline()
			total_reads_number += 1	
			kept_reads_nmuber += 1
	else:
		#outputFile_mulmap.write(read_line.strip() + "\n")
		multiple_maping = False
		read_line = inputFile_sam.readline()
		total_reads_number += 1	
		
print "total_reads_number", total_reads_number
print "multiple_mapping_number", multiple_mapping_number
print "kept_reads_nmuber", kept_reads_nmuber
print "removed_reads_number", removed_reads_number

end = time.time()
run_time = str(end - start)
run_time = run_time[:(run_time.find('.')+3)]
print "run time is: " + run_time + "s"

outputFile_record.write("total_reads_number: " + str(total_reads_number)+"\n")
outputFile_record.write("multiple_mapping_number: " + str(multiple_mapping_number)+"\n")
outputFile_record.write("kept_reads_nmuber: " + str(kept_reads_nmuber)+"\n")
outputFile_record.write("removed_reads_number: " + str(removed_reads_number)+"\n")
outputFile_record.write("run time is: " + str(run_time)+"\n")

outputFile_sam.close()	
outputFile_removed.close()	
outputFile_record.close()
#outputFile_mulmap.close()
















