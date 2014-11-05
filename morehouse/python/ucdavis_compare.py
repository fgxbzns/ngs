#!/usr/bin/python
#######################################################################################
# Guoxing Fu July 7, 2013
# for rna-seq data
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

def read_nature_data(file):
	nature_pos_dict = {}
	nature_geneid_dict = {}
	nature_geneid_line_dict = {}

	nature_genename_dict = {}

	with open(file, "r") as input_file:
		for line in input_file:
			if not line.startswith("#"):
				elements = line.strip().split(',')
				try:
					gene_id = elements[10][:-3]
					gene_name = elements[11]

					start = elements[-2]
					end = elements[-1]

					strand = "+" if int(start) <= int(end) else "-"
					if strand == "-":
						start = elements[-1]
						end = elements[-2]

					if gene_id not in nature_geneid_dict:
						nature_geneid_dict[gene_id] = gene_name
						nature_geneid_line_dict[gene_id] = line.strip()
					else:
						#print gene_id
						pass
					if start not in nature_pos_dict:
						nature_pos_dict[start] = gene_id
					else:
						#print start
						pass
					"""
					if gene_name not in nature_genename_dict:
						nature_genename_dict[gene_name] = gene_id
					else:
						#print start
						pass
					"""

				except:
					print "error in ", line
					sys.exit(1)
	print len(nature_pos_dict)
	print len(nature_geneid_dict)
	return nature_pos_dict, nature_geneid_dict, nature_geneid_line_dict

def read_cons_data(file):
	cons_dict = {}
	with open(file, "r") as input_file:
			for line in input_file:
				elements = line.strip().split()
				chr = elements[3]
				start_position = chr[chr.index(":")+1:chr.index("-")]
				gene_name = elements[2]
				#print start_position
				if gene_name not in cons_dict:
					cons_dict[gene_name] = line.strip()
				else:
					#print start_position
					pass
	return cons_dict

def read_davis_data(file):
	davis_dict = {}
	with open(file, "r") as input_file:
			for line in input_file:
				elements = line.strip().split()
				gene_id = elements[0]
				#print gene_id
				if gene_id not in davis_dict:
					davis_dict[gene_id] = line.strip()
				else:
					print gene_id
					pass
	return davis_dict

def compare_data(cons_dict, davis_dict, nature_pos_dict, nature_geneid_dict):
	geneid_not_in_nature = []
	pos_indavis_not_incons = []
	num = 0
	for gene_id in davis_dict.keys():
		if gene_id in nature_geneid_dict:
			gene_name = nature_geneid_dict[gene_id]
			#print gene_name

			if gene_name in cons_dict.keys():
				#print gene_id
				#print cons_dict[gene_name]
				#print davis_dict[gene_id]
				#print nature_geneid_line_dict[gene_id]
				num += 1
			else:
				pos_indavis_not_incons.append(gene_name)


		else:
			geneid_not_in_nature.append(gene_id)
	print len(geneid_not_in_nature)
	print geneid_not_in_nature

	print len(pos_indavis_not_incons)
	#print pos_indavis_not_incons

	print "num", num

# following functions are for combining report

def read_nature_data(file):
	# to output combined report
	nature_geneid_dict = {}

	with open(file, "r") as input_file:
		for line in input_file:
			if not line.startswith("#"):
				elements = line.strip().split(',')
				try:
					#gene_id = elements[10][:-3]
					gene_trans_id = elements[10]

					if gene_trans_id not in nature_geneid_dict:
						nature_geneid_dict[gene_trans_id] = line.strip()
					else:
						print gene_trans_id
				except:
					print "error in", line

	print len(nature_geneid_dict)
	return nature_geneid_dict

def read_conseidon_data(file):
	cons_dict = {}
	with open(file, "r") as input_file:
			for line in input_file:
				elements = line.strip().split()
				gene_trans_id = elements[2]

				if gene_trans_id not in cons_dict:
					cons_dict[gene_trans_id] = line.strip()
				else:
					#print gene_trans_id
					pass
	return cons_dict

def combine(file1, file2):
	nature_geneid_dict = read_nature_data(file1)
	cons_dict = read_conseidon_data(file2)

	with open("combined_gene.txt", "w") as output_file:
		print >> output_file, "gene_id gene_name transcript_id chr_number start_pos end_pos sample_1 sample_2", \
		                      "status FPKM_1 FPKM_2 log2(fold_change) test_stat p_value q_value significant", \
							  "go_number"
		for gene_trans_id in cons_dict.keys():
			gene_id = gene_trans_id[:-3]
			trans_id = gene_trans_id[-2:]

			cons_elements = cons_dict[gene_trans_id].split()

			chr_pos = cons_elements[3]
			chr_name = chr_pos[:chr_pos.index(":")]
			start_position = chr_pos[chr_pos.index(":")+1:chr_pos.index("-")]
			end_position = chr_pos[chr_pos.index("-")+1:]

			sample_1 = cons_elements[4]
			sample_2 = cons_elements[5]
			status = cons_elements[6]
			value_1 = cons_elements[7]
			value_2 = cons_elements[8]
			logFC = cons_elements[9]
			test_stat = cons_elements[10]

			if gene_trans_id in nature_geneid_dict:

				try:
					elements = nature_geneid_dict[gene_trans_id].split(",")

					gene_name = elements[11]
					gene_product = elements[12]
					go_number = elements[8]
					gene_function = elements[9]

					print >> output_file, gene_id, gene_name, trans_id, chr_name, start_position, end_position,
					for i in range(4, 14):
						print >> output_file, cons_elements[i],

					#print >> output_file, gene_product, go_number, gene_function
					print >> output_file, go_number
				except:
					print nature_geneid_dict[gene_trans_id]

def get_args():
	desc = "./18to19.py -e hg18 -n hg19 -d delete"
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-f", "--folder", type="string", dest="folder_name", help="Input folder name", default="null")
	parser.add_option("-n", "--ninety", type="string", dest="hg19_name", help="Input file name", default="null")
	parser.add_option("-d", "--del", type="string", dest="del_name", help="Input file name", default="null")
	(options, args) = parser.parse_args()
	"""
	if options.hg18_name == "null" or options.hg19_name == "null" or options.del_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	"""
	return options

if __name__ == '__main__':
	options = get_args()
	#file = options.folder_name
	path = "/home/guoxing/disk3/vsUCDavis/"
	start_time = time.time()
	mode = "combine"
	if mode == "compare":
		cons_dict = read_cons_data(path + "cons.txt")

		davis_dict = read_davis_data(path + "davis.txt")

		global nature_geneid_line_dict
		nature_pos_dict, nature_geneid_dict, nature_geneid_line_dict = read_nature_data(path + "nature12817-s3_sorted_5.csv")

		compare_data(cons_dict, davis_dict, nature_pos_dict, nature_geneid_dict)

	if mode == "combine":
		#combine("nature12817-s3_sorted_5.csv", "gene_exp.diff_significant")
		combine("nature12817-s3_sorted_5.csv", "NIBM.diff")

	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"
	
	
	
	
	

