#!/usr/bin/python

# for zebra fish data

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser
from tools import *

# A from Father, B from Mother
class snps:
	def __init__(self):
		self.rsID = ""
		self.position = 0
		self.A = ""
		self.B = ""
		self.covered_reads_list = []
		self.allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
		self.consistence = True
		self.max_allele = ""
		self.max_allele_number = 0
		self.max_allele_percentage = 0

# class to store reads from sam file
class read:
	def __init__(self, qName, flag, chrName, start_position, read_sequence, quality_score_sequence, read_length, covered_snp):
		self.qName = qName
		self.flag = flag
		self.chrName = chrName
		self.start_position = start_position
		self.read_sequence = read_sequence
		self.quality_score_sequence = quality_score_sequence
		self.read_length = read_length
		self.covered_snp = covered_snp
		
def get_args():
	desc="variation call"
	usage = "snpPick_fish -s sam_file" 
	parser = OptionParser(usage = usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
	(options, args) = parser.parse_args()
	if options.samFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

def is_multiple_maping(elements):
	multiple_maping = False
	try:
		XA = elements_first[21].strip()
		multiple_maping_first = True
	except:
		pass
	return multiple_maping

def variant_call_pair_end(sam_file, chr_dict):
	
	inputfile_sam = open(currentPath + sam_file, "r")
	sam_line_first = inputfile_sam.readline() # the first read line in a pair
	total_reads_num = 0
	covered_snp_total_number = 0
	
	insert_size_lower_bond = 100
	insert_size_upper_bond = 1000

	while sam_line_first!='':
		if not sam_line_first.startswith("@"):
			total_reads_num += 1	
			elements_first = sam_line_first.strip().split()
			try:
				read_ID_first = elements_first[0].strip()
				chrName_first = elements_first[2].strip()
				insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for second read is negative
			except:
				print "error in first read:", sam_line_first
				
			if (insert_size_first >= insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):
				
				# if the first read is within insert size limit, check the second read
				sam_line_second = inputfile_sam.readline()
				total_reads_num += 1
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
					if chrName_first.startswith("chr") and (chrName_first == chrName_second) and ((not is_multiple_maping(elements_first)) or (not is_multiple_maping(elements_second))): 
						
						if chrName_first not in chr_dict:
							chr_dict[chrName_first] = {}
						
						# first read
						qName_first = elements_first[0].strip()
						flag_first = elements_first[1].strip()
						start_position_first = int(elements_first[3].strip())
						read_sequence_first = elements_first[9].strip()
						read_length_first = len(read_sequence_first)
						quality_score_sequence_first = elements_first[10].strip()
						i = 0
						while i < read_length_first:
							current_base_position = start_position_first+i
							if (current_base_position) not in chr_dict[chrName_first]:
								chr_dict[chrName_first][current_base_position] = snps()
							covered_snp = read_sequence_first[i]			# ith position is the covered snp
							quality_score_symbol = quality_score_sequence_first[i]
							if (not covered_snp == 'N') and ((ord(quality_score_symbol)-33) > quality_score_threshold):	# check quality_score
								covered_snp_total_number += 1
								chr_dict[chrName_first][current_base_position].covered_reads_list.append(read(qName_first, flag_first, chrName_first, start_position_first, read_sequence_first, quality_score_sequence_first, read_length_first, covered_snp))					
								for base, value in chr_dict[chrName_second][current_base_position].allele_dict.iteritems():
									if base == covered_snp:
										chr_dict[chrName_second][current_base_position].allele_dict[base] += 1		
							i += 1
						
						# second read
						qName_second = elements_second[0].strip()
						flag_second = elements_second[1].strip()
						start_position_second = int(elements_second[3].strip())
						read_sequence_second = elements_second[9].strip()
						read_length_second = len(read_sequence_second)
						quality_score_sequence_second = elements_second[10].strip()
						i = 0
						while i < read_length_second:
							current_base_position = start_position_second+i
							if (current_base_position) not in chr_dict[chrName_second]:
								chr_dict[chrName_second][current_base_position] = snps()
							covered_snp = read_sequence_second[i]			# ith position is the covered snp
							quality_score_symbol = quality_score_sequence_second[i]
							if (not covered_snp == 'N') and ((ord(quality_score_symbol)-33) > quality_score_threshold):
								covered_snp_total_number += 1
								chr_dict[chrName_first][current_base_position].covered_reads_list.append(read(qName_second, flag_second, chrName_second, start_position_second, read_sequence_second, quality_score_sequence_second, read_length_second, covered_snp))					
								for base, value in chr_dict[chrName_second][current_base_position].allele_dict.iteritems():
									if base == covered_snp:
										chr_dict[chrName_second][current_base_position].allele_dict[base] += 1	
							i += 1
				else:
					print "first and second read ID do not match", read_ID_first, read_ID_second					
		sam_line_first = inputfile_sam.readline()
	inputfile_sam.close()
	return total_reads_num

def update_max_allele_number(max_allele_number, dict):
	dict[max_allele_number] = 1 if max_allele_number not in dict else (dict[max_allele_number] + 1)

def output_coverage_info(chr_dict):
	output_file = open(currentPath + sam_file_name + "_genotype.txt", "w")
	print >> output_file, "chr", "pos", "coverage","max_allele", "max_allele_number", "A", "T", "C", "G"
	for chr, chr_sub_dict in chr_dict.iteritems():
		chr_sub_sorted_list = sort_dict_by_key(chr_sub_dict)
		for data in chr_sub_sorted_list:
			position = data[0]
			snp = data[1]
			coverage = len(snp.covered_reads_list)
			allele_dict = snp.allele_dict
			max_allele = keywithmaxval(allele_dict)
			max_allele_number = allele_dict[max_allele]
			print >> output_file, chr, position, coverage,max_allele, max_allele_number, allele_dict["A"], allele_dict["T"], allele_dict["C"], allele_dict["G"]
			for reads in snp.covered_reads_list:
				print >> output_file, reads.qName, reads.flag, reads.start_position, reads.covered_snp, reads.read_sequence, reads.quality_score_sequence
	output_file.close()

def snpPick(sam_file):
	# start time
	start = time.time()		
		
	global quality_score_threshold
	global max_allele_percentage_threshold
	global sam_file_name
	global chr_dict

	chr_dict = {}
	quality_score_threshold = 13
		
	sam_file_name = sam_file[:(len(sam_file)-4)]
	
	print "sam file: ", sam_file_name
	"""
	data_record_file_name = sam_file_name + "_data_record.txt"
	data_record_file = open(currentPath + data_record_file_name, "w")
"""

	total_reads_num = variant_call_pair_end(sam_file, chr_dict)
	print "total_reads_num", total_reads_num
	print "chr_dict", len(chr_dict)
	
	for chr, chr_sub_dict in chr_dict.iteritems():
		print chr, len(chr_sub_dict)
	
	print chr_dict['chr12'][23241515].covered_reads_list[0].read_sequence
	print chr_dict['chr12'][23241515].allele_dict
	
	output_coverage_info(chr_dict)

	end = time.time()
	run_time = str(format((end - start), "0.3f"))
	print "run time is: " + run_time + "s"

	#data_record_file.write("run time is: " + run_time + "s \n")
	#data_record_file.close()
	
if __name__=='__main__':
	options = get_args()
	sam_file = options.samFile
	
	"""
	for infile in glob.glob(os.path.join(file_path,"*"+chr_name+"_???.phased")):	# add * for chrX
		ref_file_name = file_path + infile[(infile.find("hapmap3")):].strip()
		print ref_file_name
	"""
	snpPick(sam_file)

