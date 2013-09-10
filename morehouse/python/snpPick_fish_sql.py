#!/usr/bin/python

# for zebra fish data

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser
from tools import *
import sqlite as sql


# A from Father, B from Mother
class snps:
	def __init__(self):
		self.rsID = ""
		self.position = 0
		self.A = ""
		self.B = ""
		self.depth = 0
		self.covered_reads_list = []
		self.allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
		self.consistence = True
		self.max_allele = ""
		self.max_allele_number = 0
		self.max_allele_percentage = 0

# class to store reads from sam file
class reads:
	def __init__(self):
		self.qName = ""
		self.flag = ""
		self.chrName = ""
		self.start_position = 0
		self.read_sequence = ""
		self.quality_score_sequence = ""
		self.read_length = 0
		self.covered_snp = ""
		
def variant_call_pair_end(sam_file, chr_dict):
	"""the sequence pair has already been processed
	now treat the read as single end """
	
	global table_name
	
	
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
	
				if True:
					if True:
					#if chrName_first.startswith("chr8"): 
						
						if chrName_first not in chr_dict:
							chr_dict[chrName_first] = {}
						
						# first read
						qName_first = elements_first[0].strip()
						flag_first = elements_first[1].strip()
						start_position_first = int(elements_first[3].strip())
						read_sequence_first = elements_first[9].strip()
						read_length_first = len(read_sequence_first)
						quality_score_sequence_first = elements_first[10].strip()
						
						for i in range(read_length_first):
							current_base_position = start_position_first+i
							#geno_allele = ""
							#total_depth = 0
							A_depth = 0
							T_depth = 0
							C_depth = 0
							G_depth = 0
							#max_allele = ""
							#max_allele_number = 0
							#max_allele_percentage = 0.0
							
							
							
	
	
							covered_snp = read_sequence_first[i]			# ith position is the covered snp
							quality_score_symbol = quality_score_sequence_first[i]
							if (not covered_snp == 'N') and ((ord(quality_score_symbol)-33) > quality_score_threshold):	# check quality_score
								if covered_snp == "A":
									A_depth += 1
								elif covered_snp == "T":
									T_depth += 1
								elif covered_snp == "C":
									C_depth += 1
								elif covered_snp == "G":
									G_depth += 1
								
								
								querry = "INSERT INTO " + table_name + " (position, read_ID, chr, A_depth,\
								T_depth, C_depth, G_depth ) VALUES (" + str(current_base_position) + ",'" + \
								read_ID_first + "','" + chrName_first + "'," + str(A_depth) + "," + str(T_depth) \
								 + "," + str(C_depth) + "," + str(G_depth) + ")"
								print querry
								sql.execute_querry(db_name, querry)
								
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
			#coverage = len(snp.covered_reads_list)
			allele_dict = snp.allele_dict
			depth = snp.depth
			max_allele = keywithmaxval(allele_dict)
			max_allele_number = allele_dict[max_allele]
			print >> output_file, chr, position, depth, max_allele, max_allele_number, allele_dict["A"], allele_dict["T"], allele_dict["C"], allele_dict["G"]
			#for reads in snp.covered_reads_list:
			#	print >> output_file, reads.qName, reads.flag, reads.start_position, reads.covered_snp, reads.read_sequence, reads.quality_score_sequence
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
	
	#print chr_dict['chr12'][23241515].covered_reads_list[0].read_sequence
	#print chr_dict['chr12'][23241515].allele_dict
	
	#output_coverage_info(chr_dict)

	end = time.time()
	run_time = str(format((end - start), "0.3f"))
	print "run time is: " + run_time + "s"

	#data_record_file.write("run time is: " + run_time + "s \n")
	#data_record_file.close()

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
	
if __name__=='__main__':
	options = get_args()
	sam_file = options.samFile
	
	"""
	for infile in glob.glob(os.path.join(file_path,"*"+chr_name+"_???.phased")):	# add * for chrX
		ref_file_name = file_path + infile[(infile.find("hapmap3")):].strip()
		print ref_file_name
	"""
	
	global table_name
	db_name = "/home/guoxing/disk2/ngs.db"
	table_name = "zebra_fish"
	attribute = "position INT PRIMARY KEY, read_ID TEXT, chr TEXT, geno_allele TEXT, total_depth INT, 	\
	A_depth INT, T_depth INT, C_depth INT, G_depth INT, max_allele TEXT, max_allele_number INT, max_allele_percentage FLOAT"
	sql.creat_table(db_name, table_name, attribute)
	
	
	snpPick(sam_file)
	sql.retrive_data(db_name, table_name)

