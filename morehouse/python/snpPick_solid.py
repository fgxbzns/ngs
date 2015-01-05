#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for solid data May 09 2013
# June 22, 2013, added purity. >= 0.80 max_allele_number will be kept as seed. A. B, X, homo will be kept as seeds.
# quality score check. ord('!')-33 > 13

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *
from seed_std_compare import seed_std_compare

# A from Father, B from Mother
class snps:
	def __init__(self):
		self.rsID = ""
		self.position = 0
		self.A = ""
		self.B = ""
		self.covered_reads_list = []
		self.allele_dict = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
		self.max_allele = ""
		self.max_allele_number = 0
		self.max_allele_percentage = 0


# class to store reads from sam file
class read:
	def __init__(self, qName, flag, chrName, start_position, read_sequence, quality_score_sequence, read_length,
	             covered_snp):
		self.qName = qName
		self.flag = flag
		self.chrName = chrName
		self.start_position = start_position
		self.read_sequence = read_sequence
		self.quality_score_sequence = quality_score_sequence
		self.read_length = read_length
		self.covered_snp = covered_snp

class parameters:
	def __init__(self):
		self.chr_name = ""
		self.depth_threshold = 1

def load_hap_std(file_name):
	hap_std_dict = {}
	title_haplotype = load_raw_data(file_name, raw_data_format)[0]
	data_dict = load_raw_data(file_name, raw_data_format)[1]
	for position, elements in data_dict.iteritems():
		if True:
			try:
				snp = snps()
				snp.rsID = elements[0]
				snp.position = int(position)
				# for wli nea data, for solid data. remove the comments
				#snp.A = elements[2].strip()
				#snp.B = elements[3].strip()
				hap_std_dict[position] = snp
			except ValueError:
				print "error in ", file_name, position
	return title_haplotype, hap_std_dict

def variant_call_single_end(sam_file, hap_std_dict, chr_name):
	inputfile_sam = open(sam_file, "r")
	sam_line_first = inputfile_sam.readline()  # the first read line in a pair
	total_reads_num = 0
	covered_snp_total_number = 0

	# solid data, single end. no insert size
	while sam_line_first != '':
		if not sam_line_first.startswith("@"):
			total_reads_num += 1
			elements_first = sam_line_first.strip().split()
			try:
				read_ID_first = elements_first[0].strip()
				chrName_first = elements_first[2].strip()
			except:
				print "error in first read:", sam_line_first
			# first read
			qName_first = elements_first[0].strip()
			flag_first = elements_first[1].strip()
			start_position_first = int(elements_first[3].strip())
			read_sequence_first = elements_first[9].strip()
			read_length_first = len(read_sequence_first)
			quality_score_sequence_first = elements_first[10].strip()
			if chrName_first == chr_name:
				for i in range(read_length_first):
					current_base_position = start_position_first + i
					if current_base_position in hap_std_dict:
						covered_snp = read_sequence_first[i]  # ith position is the covered snp
						try:
							quality_score_symbol = quality_score_sequence_first[i]
							if (not covered_snp == 'N') and (
								(ord(quality_score_symbol) - 33) > quality_score_threshold):  # check quality_score
								covered_snp_total_number += 1
								hap_std_dict[current_base_position].covered_reads_list.append(
									read(qName_first, flag_first, chrName_first, \
									     start_position_first, read_sequence_first, quality_score_sequence_first,
									     read_length_first, covered_snp))
								for base, value in hap_std_dict[current_base_position].allele_dict.iteritems():
									if base == covered_snp:
										hap_std_dict[current_base_position].allele_dict[base] += 1
						except:
							print "different length of read and quality score", sam_line_first

		sam_line_first = inputfile_sam.readline()
	inputfile_sam.close()
	hap_std_sorted_list = sort_dict_by_key(hap_std_dict)
	return (total_reads_num, covered_snp_total_number, hap_std_sorted_list)

def variant_call_single_end_bkup(sam_file, hap_std_dict, chr_name):
	# working for solid, backuped March 19 2014
	inputfile_sam = open(sam_file, "r")
	sam_line_first = inputfile_sam.readline()  # the first read line in a pair
	total_reads_num = 0
	covered_snp_total_number = 0

	# solid data, single end. no insert size
	while sam_line_first != '':
		if not sam_line_first.startswith("@"):
			total_reads_num += 1
			elements_first = sam_line_first.strip().split()
			try:
				read_ID_first = elements_first[0].strip()
				chrName_first = elements_first[2].strip()
			except:
				print "error in first read:", sam_line_first
			# first read
			qName_first = elements_first[0].strip()
			flag_first = elements_first[1].strip()
			start_position_first = int(elements_first[3].strip())
			read_sequence_first = elements_first[9].strip()
			read_length_first = len(read_sequence_first)
			quality_score_sequence_first = elements_first[10].strip()
			i = 0
			while i < read_length_first:
				current_base_position = start_position_first + i
				if current_base_position in hap_std_dict:
					covered_snp = read_sequence_first[i]  # ith position is the covered snp
					quality_score_symbol = quality_score_sequence_first[i]
					if (chrName_first == chr_name) and (not covered_snp == 'N') and (
						(ord(quality_score_symbol) - 33) > quality_score_threshold):  # check quality_score
						covered_snp_total_number += 1
						hap_std_dict[current_base_position].covered_reads_list.append(
							read(qName_first, flag_first, chrName_first, \
							     start_position_first, read_sequence_first, quality_score_sequence_first,
							     read_length_first, covered_snp))
						for base, value in hap_std_dict[current_base_position].allele_dict.iteritems():
							if base == covered_snp:
								hap_std_dict[current_base_position].allele_dict[base] += 1
				i += 1
		sam_line_first = inputfile_sam.readline()
	inputfile_sam.close()
	hap_std_sorted_list = sort_dict_by_key(hap_std_dict)
	return (total_reads_num, covered_snp_total_number, hap_std_sorted_list)

"""
def is_multiple_maping(elements_first):
	multiple_maping = False
	try:
		XA = elements_first[21].strip()
		multiple_maping_first = True
	except:
		pass
	return multiple_maping
"""

def is_multiple_maping(elements):
	multiple_maping = False
	XA = elements[-1].strip()
	if XA.startswith('XA'):
		multiple_maping = True
	else:
		XA = ""
	return multiple_maping

def variant_call_pair_end(sam_file, hap_std_dict):
	inputfile_sam = open(currentPath + sam_file, "r")
	sam_line_first = inputfile_sam.readline()  # the first read line in a pair
	total_reads_num = 0
	covered_snp_total_number = 0

	insert_size_lower_bond = 100
	insert_size_upper_bond = 1000

	while sam_line_first != '':
		if not sam_line_first.startswith("@"):
			total_reads_num += 1
			elements_first = sam_line_first.strip().split()
			try:
				read_ID_first = elements_first[0].strip()
				chrName_first = elements_first[2].strip()
				insert_size_first = abs(int(elements_first[8].strip()))  #  insert_size for second read is negative
			except:
				print "error in first read:", sam_line_first
			if (insert_size_first >= insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):
				multiple_maping_first = is_multiple_maping(elements_first)

				# if the first read is within insert size limit, check the second read
				sam_line_second = inputfile_sam.readline()
				total_reads_num += 1
				elements_second = sam_line_second.strip().split()
				try:
					read_ID_second = elements_second[0].strip()
					chrName_second = elements_second[2].strip()
					insert_size_second = abs(
						int(elements_second[8].strip()))  #  insert_size for second read is negative
				except:
					print "error in second read:", sam_line_second
				if read_ID_first == read_ID_second:  # check if the two reads belong to the same pair
					multiple_maping_second = is_multiple_maping(elements_second)

					if (not multiple_maping_first) or (
					not multiple_maping_second):  # keep the pair as long as one read is not multiple mapping
						# first read
						qName_first = elements_first[0].strip()
						flag_first = elements_first[1].strip()
						start_position_first = int(elements_first[3].strip())
						read_sequence_first = elements_first[9].strip()
						read_length_first = len(read_sequence_first)
						quality_score_sequence_first = elements_first[10].strip()
						i = 0
						while i < read_length_first:
							current_base_position = start_position_first + i
							if current_base_position in hap_std_dict:
								covered_snp = read_sequence_first[i]  # ith position is the covered snp
								quality_score_symbol = quality_score_sequence_first[i]
								if (chrName_first == chr_name) and (not covered_snp == 'N') and (
									(ord(quality_score_symbol) - 33) > quality_score_threshold):  # check quality_score
									covered_snp_total_number += 1
									hap_std_dict[current_base_position].covered_reads_list.append(
										read(qName_first, flag_first, chrName_first, start_position_first,
										     read_sequence_first, quality_score_sequence_first, read_length_first,
										     covered_snp))
									for base, value in hap_std_dict[current_base_position].allele_dict.iteritems():
										if base == covered_snp:
											hap_std_dict[current_base_position].allele_dict[base] += 1
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
							current_base_position = start_position_second + i
							if current_base_position in hap_std_dict:
								covered_snp = read_sequence_second[i]  # ith position is the covered snp
								quality_score_symbol = quality_score_sequence_second[i]
								if (chrName_second == chr_name) and (not covered_snp == 'N') and (
									(ord(quality_score_symbol) - 33) > quality_score_threshold):
									covered_snp_total_number += 1
									hap_std_dict[current_base_position].covered_reads_list.append(
										read(qName_second, flag_second, chrName_second, start_position_second,
										     read_sequence_second, quality_score_sequence_second, read_length_second,
										     covered_snp))
									for base, value in hap_std_dict[current_base_position].allele_dict.iteritems():
										if base == covered_snp:
											hap_std_dict[current_base_position].allele_dict[base] += 1
							i += 1

				else:
					print "first and second read ID do not match", read_ID_first, read_ID_second
		sam_line_first = inputfile_sam.readline()
	inputfile_sam.close()
	hap_std_sorted_list = sort_dict_by_key(hap_std_dict)
	return (total_reads_num, covered_snp_total_number, hap_std_sorted_list)

def output_coverage_info(hap_std_sorted_list):
	outputFile_reads = open(currentPath + sam_file_name + "_" + chr_name + "_reads.txt", "w")
	outputFile_reads.write("SNP position \t Depth \n")
	outputfile_allele = open(currentPath + sam_file_name + "_" + chr_name + "_allele.txt", "w")
	outputfile_allele.write("Chromosome \t position \t Total Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")
	for snp_data in hap_std_sorted_list:
		snp = snp_data[1]
		if len(snp.covered_reads_list) > parameter.depth_threshold:
			outputfile_allele.write(
				chr_name + "\t" + str(snp.position) + "\t" + str(len(snp.covered_reads_list)) + "\t" + str(
					snp.allele_dict['A']) \
				+ "\t" + str(snp.allele_dict['T']) + "\t" + str(snp.allele_dict['C']) + "\t" + str(
					snp.allele_dict['G']) + "\n")
			outputFile_reads.write(
				"@_" + snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\t" + str(
					len(snp.covered_reads_list)) + "\n")
			for reads in snp.covered_reads_list:
				outputFile_reads.write(
					reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" \
					+ reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
	outputFile_reads.close()
	outputfile_allele.close()

def update_max_allele_number(max_allele_number, dict):
	dict[max_allele_number] = 1 if max_allele_number not in dict else (dict[max_allele_number] + 1)

def output_ABX_files(file_name, title_info, data_dict):
	file = open(file_name, "w")
	print >> file, title_info
	for position, snp in data_dict.iteritems():
		print >> file, snp.rsID, snp.position, snp.A, snp.B, snp.max_allele, snp.max_allele_number, snp.max_allele_percentage
		for reads in snp.covered_reads_list:
			print >> file, reads.qName, reads.flag, reads.chrName, reads.start_position, reads.covered_snp, reads.read_sequence, reads.quality_score_sequence
	file.close()

def output_seed_file(file_name, title_info, data_dict):
	file = open(file_name, "w")
	print >> file, title_info
	for position, snp in data_dict.iteritems():
		print >> file, snp.rsID, snp.position, snp.max_allele
	file.close()

def get_called_seed_dict():
	# for solid data
	prefiltered_seed_dict = {}
	called_seed_dict = {}
	for snp_data in hap_std_sorted_list:
		snp = snp_data[1]
		position = snp.position
		if len(snp.covered_reads_list) > parameter.depth_threshold:
			max_allele = keywithmaxval(snp.allele_dict)
			max_allele_number = snp.allele_dict[max_allele]
			snp.max_allele = max_allele
			snp.max_allele_number = max_allele_number
			snp.max_allele_percentage = float(max_allele_number) / float(len(snp.covered_reads_list))

			if snp.max_allele_percentage >= max_allele_percentage_threshold and max_allele_number > parameter.depth_threshold:
				prefiltered_seed_dict[position] = snp
				# check genotype to remove called base that does not in genotype
				if position in geno_dict and max_allele in geno_dict[position][2]:
					called_seed_dict[position] = snp
	return called_seed_dict, prefiltered_seed_dict

def get_called_seed_dict_neandertal():
	# for neandertal data wli
	prefiltered_seed_dict = {}
	called_seed_dict = {}
	for snp_data in hap_std_sorted_list:
		snp = snp_data[1]
		position = snp.position
		if len(snp.covered_reads_list) > parameter.depth_threshold:
			max_allele = keywithmaxval(snp.allele_dict)
			max_allele_number = snp.allele_dict[max_allele]
			snp.max_allele = max_allele
			snp.max_allele_number = max_allele_number
			snp.max_allele_percentage = float(max_allele_number) / float(len(snp.covered_reads_list))

			if snp.max_allele_percentage >= max_allele_percentage_threshold and max_allele_number > parameter.depth_threshold:
				called_seed_dict[position] = snp
	return called_seed_dict

def compare_with_std_hap():
	distribution_file = open(sam_file_name + "_distribution.txt", "w")
	distribution_file.write("rsID \t phys_position \t snp.A	\t snp.B \t distribution \n")

	pure_total = 0
	not_in_genotype = 0
	A_dict = {}
	B_dict = {}
	X_dict = {}
	not_ABX_dict = {}
	not_in_genotype_dict_num = {}

	same_to_A_dict = {}
	same_to_B_dict = {}
	same_to_AB_dict = {}
	same_to_X_dict = {}
	same_to_N_dict = {}
	not_in_genotype_dict = {}
	for position, snp in prefiltered_seed_dict.iteritems():
		if position in called_seed_dict:
			max_allele = snp.max_allele
			max_allele_number = snp.max_allele_number
			if max_allele == snp.A or max_allele == snp.B:  #check genotype
				if max_allele == snp.A and max_allele != snp.B:
					same_to_A_dict[position] = snp
					distribution_file.write(
						snp.rsID + "\t" + str(position) + "\t" + snp.A + "\t" + snp.B + "\t" + "A" + "\n")
					update_max_allele_number(max_allele_number, A_dict)
				if max_allele == snp.B and max_allele != snp.A:
					same_to_B_dict[position] = snp
					distribution_file.write(
						snp.rsID + "\t" + str(position) + "\t" + snp.A + "\t" + snp.B + "\t" + "B" + "\n")
					update_max_allele_number(max_allele_number, B_dict)
				if max_allele == snp.B and max_allele == snp.A:
					same_to_AB_dict[position] = snp
					distribution_file.write(
						snp.rsID + "\t" + str(position) + "\t" + snp.A + "\t" + snp.B + "\t" + "" + "\n")
			elif snp.A == "X" or snp.B == "X":  # keep these reads too
				same_to_X_dict[position] = snp
				update_max_allele_number(max_allele_number, X_dict)
			else:
				same_to_N_dict[position] = snp
				update_max_allele_number(max_allele_number, not_ABX_dict)
		else:  # not in genotype, sequencing error???
			not_in_genotype_dict[position] = snp
			update_max_allele_number(max_allele_number, not_in_genotype_dict_num)

	same_to_A = len(same_to_A_dict)
	same_to_B = len(same_to_B_dict)
	same_to_AB = len(same_to_AB_dict)
	X_AB = len(same_to_X_dict)
	not_ABX = len(same_to_N_dict)
	not_in_genotype = len(not_in_genotype_dict)

	A_in_hetero = format((float(same_to_A) / float(same_to_A + same_to_B)) * 100, "0.2f")
	B_in_hetero = format((float(same_to_B) / float(same_to_A + same_to_B)) * 100, "0.2f")
	correct_rate_in_hetero = A_in_hetero if A_in_hetero >= B_in_hetero else B_in_hetero

	A_in_all = format((float(same_to_A) / float(called_seed_total)) * 100, "0.2f")
	B_in_all = format((float(same_to_B) / float(called_seed_total)) * 100, "0.2f")
	correct_rate_in_all_seed = A_in_all if A_in_all >= B_in_all else B_in_all

	hete_seed_rate = format(float(same_to_A + same_to_B) / float(called_seed_total) * 100, "0.2f")

	print "same_to_A", same_to_A
	print "same_to_B", same_to_B
	print "same_to_AB: ", same_to_AB
	print "X_AB", X_AB
	print "not_ABX", not_ABX
	print "A_in_hetero", A_in_hetero, "%"
	print "B_in_hetero", B_in_hetero, "%"
	print "correct_rate_in_hetero: ", correct_rate_in_hetero
	print "correct_rate_in_all_seed: ", correct_rate_in_all_seed
	print "(A+B)/seed: ", hete_seed_rate
	print "not_in_genotype", not_in_genotype
	"""
	print >> data_record_file, "same_to_A: ", same_to_A
	print >> data_record_file, "same_to_B: ", same_to_B
	print >> data_record_file, "same_to_AB: ", same_to_AB
	print >> data_record_file, "X_AB: ", X_AB
	print >> data_record_file, "not_ABX: ", not_ABX
	print >> data_record_file, "A_in_hetero: ", A_in_hetero, "%"
	print >> data_record_file, "B_in_hetero: ", B_in_hetero, "%"
	print >> data_record_file, "correct_rate_in_hetero: ", correct_rate_in_hetero
	print >> data_record_file, "correct_rate_in_all_seed: ", correct_rate_in_all_seed
	print >> data_record_file, "hete_seed_rate: ", hete_seed_rate
	print >> data_record_file, "not_in_genotype: ", not_in_genotype
	
	print >> data_record_file, "all", "sam_file", "chr", "depth_threshold", "same_to_A", "same_to_B", "correct_rate_in_hetero", "same_to_AB", "X_AB", "not_ABX", "called_seed_total", "correct_rate_in_all_seed", "hete_seed_rate", "pure_total", "not_in_genotype"
	print >> data_record_file, "data", sam_file, chr_name, depth_threshold, same_to_A, same_to_B, correct_rate_in_hetero, same_to_AB, X_AB, not_ABX, called_seed_total, correct_rate_in_all_seed, hete_seed_rate, pure_total, not_in_genotype

	title_info = "rsID \t phys_position \t std_A \t std_B \t max_allele \t max_allele_number \t max_allele_percentage"
	file_name = sam_file_name + "_A.txt"
	output_ABX_files(file_name, title_info, same_to_A_dict)
	file_name = sam_file_name + "_B.txt"
	output_ABX_files(file_name, title_info, same_to_B_dict)
	file_name = sam_file_name + "_AB.txt"
	output_ABX_files(file_name, title_info, same_to_AB_dict)
	file_name = sam_file_name + "_X.txt"
	output_ABX_files(file_name, title_info, same_to_X_dict)
	file_name = sam_file_name + "_notABX.txt"
	output_ABX_files(file_name, title_info, same_to_N_dict)
	file_name = sam_file_name + "_notInGenotype.txt"
	output_ABX_files(file_name, title_info, not_in_genotype_dict)
	"""


def coverage_distribution():
	A_dict = {}
	B_dict = {}
	X_dict = {}
	not_ABX_dict = {}

	""" caulculate the coverage distribution for each snp """
	coverage_distribution_file = open(currentPath + sam_file_name + "_coverage_distribution.txt", "w")
	coverage_distribution_file.write("coverage depth \t A \t B \t X \t not_ABX \t percentage \n")

	for key, value in A_dict.iteritems():
		B_value = 0
		X_value = 0
		notABX_value = 0
		coverage_distribution_file.write(str(key) + "\t" + str(value) + "\t")
		if key in B_dict:
			coverage_distribution_file.write(str(B_dict[key]) + "\t")
			B_value = B_dict[key]
		else:
			coverage_distribution_file.write("\t")
		if key in X_dict:
			coverage_distribution_file.write(str(X_dict[key]) + "\t")
			X_value = X_dict[key]
		else:
			coverage_distribution_file.write("\t")
		if key in not_ABX_dict:
			coverage_distribution_file.write(str(not_ABX_dict[key]) + "\t")
			notABX_value = not_ABX_dict[key]
		else:
			coverage_distribution_file.write("\t")
		coverage_distribution_file.write(
			format(float(value) / float(value + B_value + X_value + notABX_value), "0.3f") + "\n")
	coverage_distribution_file.close()


def combine_called_seed_geno(called_seed_dict, geno_homo_dict):
	for position, geno_data in geno_homo_dict.iteritems():
		if position not in called_seed_dict and geno_data[2][0] != 'N':
			snp = snps()
			snp.rsID = geno_data[0]
			snp.position = int(position)
			snp.max_allele = geno_data[2][0]
			called_seed_dict[position] = snp


def snpPick_solid(sam_file, depth_threshold, chr_name):
	# for solid data
	#global chr_name
	global quality_score_threshold
	global max_allele_percentage_threshold

	global geno_dict
	global hap_std_dict
	global prefiltered_seed_dict
	global called_seed_dict
	global geno_homo_dict
	global called_seed_total
	global hap_ref_dict

	global sam_file_name
	global title_haplotype
	global hap_std_sorted_list

	global parameter
	parameter = parameters()
	parameter.chr_name = chr_name
	parameter.depth_threshold = depth_threshold

	quality_score_threshold = 13    # solid data
	#quality_score_threshold = 30    # neand data
	max_allele_percentage_threshold = 0.8

	#haplotype_file = "ASW_" + chr_name + "_child_hap_refed.txt"  # for solid and 454 NA10847
	#genotype_file = "genotype_NA10847_" + chr_name + ".txt"  # for solid and 454 NA10847

	#haplotype_file = "NA12878_chr4_haplotype_std_hg18.txt"	# for illumina hg18 NA12878 chr4
	#genotype_file = "genotype_NA12878_"+chr_name+".txt"	# for illumina hg18 NA12878 chr4


	haplotype_file = "NA12878_hap_new_refed.txt"	# for simulation hg18 NA12878 chr6
	genotype_file = "genotype_NA12878_"+chr_name+".txt"	# for simulation
	#haplotype_file = "NA12878_hap_new_refed.txt" # for quake data
	#genotype_file = "genotype_NA12878_chr6.txt"	# for quake data


	sam_file_name = sam_file[:(len(sam_file) - 4)] + "_" + str(parameter.depth_threshold)

	print "haplotype file: ", haplotype_file
	print "genotype_file : ", genotype_file
	print "sam file: ", sam_file_name

	data_record_file_name = sam_file_name + "_data_record.txt"
	data_record_file = open(data_record_file_name, "w")

	title_haplotype, hap_std_dict = load_hap_std(file_path + haplotype_file)

	geno_dict = load_raw_data(file_path + genotype_file, raw_data_format)[1]

	total_reads_num, covered_snp_total_number, hap_std_sorted_list = variant_call_single_end(sam_file, hap_std_dict,
	                                                                                         chr_name)

	print "haplotye list size is: ", len(hap_std_dict)
	print "total_reads_num", total_reads_num

	seed_tuple = get_called_seed_dict()
	called_seed_dict = seed_tuple[0]
	prefiltered_seed_dict = seed_tuple[1]

	called_seed_file_name = sam_file_name + "_called_seed.txt"
	title_info = title_haplotype
	output_seed_file(called_seed_file_name, title_info, called_seed_dict)
	called_seed_total = len(called_seed_dict)
	pure_total = len(prefiltered_seed_dict)

	print "called_seed_dict", called_seed_total
	print "pure_total (include alleles pure but may not in genotype)", pure_total
	print >> data_record_file, "haplotype file: ", haplotype_file
	print >> data_record_file, "haplotye list size is:", len(hap_std_dict)
	print >> data_record_file, "genotype_file file: ", genotype_file
	print >> data_record_file, "sam file: ", sam_file
	print >> data_record_file, "reads list size is: ", total_reads_num
	print >> data_record_file, "called_seed_dict: ", called_seed_total
	print >> data_record_file, "pure_total (include alleles pure but may not in genotype): ", pure_total

	# correct errors in seed, add homo from genotype to seed. for testing purpose
	geno_sorted_list = sort_dict_by_key(geno_dict)

	group_tuple = group_seed(geno_dict, geno_dict)
	geno_homo_dict = group_tuple[0]
	#geno_hetero_dict = group_tuple[1]
	print len(geno_homo_dict)
	#print len(geno_hetero_dict)
	"""
	group_tuple = group_seed(called_seed_dict, geno_dict)
	called_seed_homo_dict = group_tuple[0]
	called_seed_hetero_dict = group_tuple[1]
	print len(called_seed_homo_dict)
	print len(called_seed_hetero_dict)
	"""

	print "called_seed_dict", len(called_seed_dict)
	combine_called_seed_geno(called_seed_dict, geno_homo_dict)
	print "called_seed_dict", len(called_seed_dict)
	combined_seed_file_name = sam_file_name + "_combined_seed.txt"
	title_info = title_haplotype
	output_seed_file(combined_seed_file_name, title_info, called_seed_dict)

	data_record_file.close()

	seed_std_compare(called_seed_file_name, chr_name)
	#seed_std_compare(combined_seed_file_name, chr_name)
	"""
	refMerger(combined_seed_file_name, chr_name)
	file_name = "haplotype.txt"
	hifi_run(file_name, chr_name)
	hifiAccuCheck("imputed_"+file_name, chr_name)
	"""

def snpPick_neandertal(sam_file, chr_name):
	# for neandertal data wli
	global quality_score_threshold
	global max_allele_percentage_threshold

	global geno_dict
	global hap_std_dict
	global prefiltered_seed_dict
	global called_seed_dict
	global geno_homo_dict
	global called_seed_total
	global hap_ref_dict

	global sam_file_name
	global title_haplotype
	global hap_std_sorted_list

	global parameter
	parameter = parameters()
	parameter.chr_name = chr_name
	parameter.depth_threshold = 0

	#quality_score_threshold = 13  # solid data
	quality_score_threshold = 30   # neand wli data
	max_allele_percentage_threshold = 0.8

	haplotype_file = "ASW_" + chr_name + "_child_hap_refed.txt"  # for solid and 454 NA10847

	sam_file_name = sam_file[:(len(sam_file) - 4)] + "_" + str(parameter.depth_threshold)

	print "haplotype file: ", haplotype_file
	print "sam file: ", sam_file_name

	data_record_file_name = sam_file_name + "_data_record.txt"
	data_record_file = open(data_record_file_name, "w")
	hap_std_path = "/home/wli/nfs1_node2/ASW_hap/"
	title_haplotype, hap_std_dict = load_hap_std(hap_std_path + haplotype_file)


	total_reads_num, covered_snp_total_number, hap_std_sorted_list = variant_call_single_end(sam_file, hap_std_dict,
	                                                                                         chr_name)

	print "haplotye list size is: ", len(hap_std_dict)
	print "total_reads_num", total_reads_num

	called_seed_dict = get_called_seed_dict_neandertal()

	called_seed_file_name = sam_file_name + "_called_seed_" + chr_name + ".txt"
	title_info = title_haplotype
	output_seed_file(called_seed_file_name, title_info, called_seed_dict)
	called_seed_total = len(called_seed_dict)

	print "called_seed_dict", called_seed_total
	print >> data_record_file, "haplotype file: ", haplotype_file
	print >> data_record_file, "haplotye list size is:", len(hap_std_dict)
	print >> data_record_file, "sam file: ", sam_file
	print >> data_record_file, "reads list size is: ", total_reads_num
	print >> data_record_file, "called_seed_dict: ", called_seed_total

	data_record_file.close()

	output_coverage_info(hap_std_sorted_list)

def get_args():
	desc = "variation call"
	usage = "snpPick_solid -s sam -c chr# -d depth_threshold"
	parser = OptionParser(usage=usage)
	parser.add_option("-i", "--haplotype", type="string", dest="haplotypeFile", help="Input File Name", default="null")
	parser.add_option("-s", "--sam", type="string", dest="samFile", help="Input File Name", default="null")
	parser.add_option("-d", "--threshold", type="string", dest="threshold", help="Input the depth threshold", default="1")
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode", help="solid or nea", default="null")
	#parser.add_option("-h", "--hstd", type="string", dest="hstd", help="haplotype std file", default="null")

	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.samFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	options = get_args()
	sam_file = options.samFile
	depth_threshold = int(options.threshold)
	global chr_name
	chr_name = options.chrName
	mode = options.mode

	start_time = time.time()
	if mode == "solid":

		snpPick_solid(sam_file, depth_threshold, chr_name)
		seed_std_compare(sam_file_name + "_combined_seed.txt", chr_name)
	elif mode == "illumina":

		snpPick_solid(sam_file, depth_threshold, chr_name)
		seed_std_compare(sam_file_name + "_combined_seed.txt", chr_name)
	elif mode == "nea":
		#haplotype_file = options.hstd
		snpPick_neandertal(sam_file, chr_name)

	elapsed_time = time.time() - start_time
	print "Elapsed time is: ", round(elapsed_time, 3), "s"
	#compare_with_std_hap()




















# compare with ori hap file			
"""
unchanged_snp_file_name = sam_file_name + "_unchanged_snp.txt"
unchanged_snp_file = open(currentPath + unchanged_snp_file_name, "w")		
unchanged_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

changed_snp_file_name = sam_file_name + "_changed_snp.txt"
changed_snp_file = open(currentPath + changed_snp_file_name, "w")		
changed_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

# need to update for A or B 	
for snp_data in hap_std_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) > depth_threshold:	
		consisitent = True
		for reads in snp.covered_reads_list:
			if snp.A != reads.covered_snp:
			#if snp.B != reads.covered_snp:
				consisitent = False
		if consisitent:
			unchanged_snp_file.write(str(snp.position) + "\t" + snp.A + "\t" + str(len(snp.covered_reads_list)) + "\t" +"\n")
			#unchanged_snp_file.write(str(snp.position) + "\t" + snp.B + "\t" + str(len(snp.covered_reads_list)) + "\t" +"\n")
		if not consisitent:
			changed_snp_file.write(str(snp.position) + "\t" + snp.A + "\t" + str(len(snp.covered_reads_list)) +"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+ "\n")
			#changed_snp_file.write(str(snp.position) + "\t" + snp.B + "\t" + str(len(snp.covered_reads_list)) +"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+ "\n")

unchanged_snp_file.close()
changed_snp_file.close()
"""


# check symmetric and asymmetric
"""
same_A_file = open(currentPath + sam_file_name + "_same_A.txt", "w")
same_A_file.write("rsID \t phys_position \t NA12878_A \t selected_base \t reads	\n")

pair_A_file = open(currentPath + sam_file_name + "_pair_A.txt", "w")
pair_A_file.write("rsID \t phys_position \t NA12878_A \t selected_base \t reads	\n")

different_A_file = open(currentPath + sam_file_name + "_different_A.txt", "w")
different_A_file.write("rsID \t phys_position \t NA12878_A \t selected_base \t reads \n")

same_A = 0
pair_A = 0
different_A = 0


for snp_data in hap_std_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) > depth_threshold:	
		max_allele = keywithmaxval(snp.allele_dict)
		max_allele_number = snp.allele_dict[max_allele]

		pure = True
		for base in base_list:
			if base != max_allele:
				if snp.allele_dict[base] != 0:
					pure = False
		if pure:
			pure_total += 1
			hap_base = snp.A
			if hap_base == "A":
				if max_allele == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_allele == 'T':
					pair_A += 1
					pair_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B +"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						pair_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				else:
					different_A += 1
					different_A_file.positionID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						different_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")		
			elif hap_base == "T":
				if max_allele == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_allele == 'A':
					pair_A += 1
					pair_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list)) + "\n")
					for reads in snp.covered_reads_list:
						pair_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				else:
					different_A += 1
					different_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						different_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
			elif hap_base == "C":
				if max_allele == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_allele == 'G':
					pair_A += 1
					pair_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list)) + "\n")
					for reads in snp.covered_reads_list:
						pair_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				else:
					different_A += 1
					different_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B +"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						different_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
			elif hap_base == "G":
				if max_allele == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A  +"\t"+ snp.B+"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_allele == 'C':
					pair_A += 1
					pair_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B +"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						pair_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				else:
					different_A += 1
					different_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B + "\t" + str(len(snp.covered_reads_list)) + "\n")
					for reads in snp.covered_reads_list:
						different_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.chrName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")	
			#else:
				#print hap_base
			


print "same A", same_A
print "symmetric to A is", pair_A
print "asymmetric to A is", different_A

data_record_file.write("same A is: " + str(same_A) + "\n")
data_record_file.write("symmetric to A is: " + str(pair_A) + "\n")
data_record_file.write("asymmetric to A is: " + str(different_A) + "\n")
"""
