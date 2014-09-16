#!/usr/bin/python
#######################################################################################
# Guoxing Fu July 7, 2013
# for rna-seq data
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

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
	
	
	
	
	
	
	