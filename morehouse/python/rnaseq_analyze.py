#!/usr/bin/python
#######################################################################################
# Guoxing Fu July 7, 2013
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

chr_dict = {}

class Chr:
	def __init__(self):
		self.chr_name = ""
		self.length = 0
		self.gene_dict = {}
		self.gene_number = 0
		self.transcript_dict = {}
		self.transcript_number = 0

class Transcript:
	def __init__(self):
		self.transcript_id = ""
		self.transcript_start = 0
		self.transcript_end = 0
		self.transcript_length = 0
		self.exon_dict = {}
		self.exon_number = 0
		self.FPKM = 0

class Exon:
	def __init__(self):
		self.exon_id = ""
		self.exon_start = 0
		self.exon_end = 0
		self.exon_length = 0

def convert_file(file):
	#with open(file + "_output.txt", "w") as output_file:
	with open(file, "r") as input_file:
		line = input_file.readline()
		while line != "":
			elements = line.strip().split()
			chr_name = elements[0]
			if chr_name not in chr_dict:
				chr_dict[chr_name] = Chr()
			element_type = elements[2]
			start = elements[3]
			end = elements[4]
			length = int(end) - int(start)
			gene_id = elements[9][1:-2]
			if gene_id not in chr_dict[chr_name].gene_dict:
				chr_dict[chr_name].gene_dict[gene_id] = ""
			transcript_id = elements[11][1:-2]
			FPKM = elements[11][1:-2]

			if element_type == "transcript":
				temp_trans = Transcript()
				temp_trans.transcript_id = transcript_id
				temp_trans.transcript_start = start
				temp_trans.transcript_end = end
				temp_trans.transcript_length = length
				temp_trans.FPKM = FPKM

				line = input_file.readline()



def process_line(line):

	elements = line.strip().split()
	chr_name = elements[0]


def covert_line(line):
	#position_key = 0
	gtf_line = ""
	elements = line.strip().split(',')
	try:

		gene_id = elements[10][:-3]
		chr_name = "chr" + gene_id[2]
		#transcript_id = elements[10][-1:]
		transcript_id = elements[10]
		source_name = "ncbi"
		seq_type = "exon"
		start = elements[-2]
		end = elements[-1]

		gene_name = elements[11]
		gene_product = elements[12]
		go_number = elements[8]
		gene_function = elements[9]

		strand = "+" if int(start) <= int(end) else "-"
		#position_key = int(start) if int(start) >= int(end) else int(end)
		if strand == "-":
			start = elements[-1]
			end = elements[-2]

		gtf_line = chr_name + "\t" + source_name + "\t" + seq_type + "\t" \
			+ start + "\t" + end + "\t.\t" + strand + "\t.\t" + "gene_id \"" + gene_id \
			+ "\"; " + "transcript_id \"" + transcript_id \
			+ "\"; " + "gene_name \"" + gene_name \
			+ "\"; " + "gene_product \"" + gene_product  \
			+ "\"; " + "GO_number \"" + go_number + "\"; " \
			+ "gene_function \"" + gene_function + "\"; "

	except:
		print "error in ", line
		sys.exit(1)

	return chr_name, gtf_line

def convert_file(file):
	with open(file + ".gtf", "w") as output_file:
		with open(file, "r") as input_file:
			for line in input_file:
				if not line.startswith("#"):
					chr_name, gtf_line = covert_line(line)
					#gtf_line = covert_line_att(line)
					if not chr_name == "chr_":
						print >> output_file, gtf_line

def convert_ref_transcriptome():
	"""
	convert transcriptome from refbeet-1.2 to refbeet-1.1
	"""

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
	file = options.folder_name

	start_time = time.time()
	convert_file(file)
	os.system("sort -k 1,1 -k 4,4n " + file + ".gtf > sorted_ref.gtf")
	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"
	
	
	
	
	
	
	