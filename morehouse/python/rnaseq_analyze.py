#!/usr/bin/python
#######################################################################################
# Guoxing Fu July 7, 2014
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *



class Chr:
	def __init__(self):
		self.chr_name = ""
		self.length = 0
		self.gene_dict = {}
		self.gene_number = 0
		self.transcript_dict = {}
		self.transcript_number = 0

		self.transcript_length_list = []
		self.max_transcript_length = 0
		self.min_transcript_length = 0
		self.avg_transcript_length = 0
		self.stdev_transcript_length = 0

		self.exon_number_list = []
		self.max_exon_number = 0
		self.min_exon_number = 0
		self.avg_exon_number = 0
		self.stdev_exon_number = 0

		self.exon_length_list = []
		self.max_exon_length = 0
		self.min_exon_length = 0
		self.avg_exon_length = 0
		self.stdev_exon_length = 0

class Transcript:
	def __init__(self):
		self.transcript_id = ""
		self.transcript_start = 0
		self.transcript_end = 0
		self.transcript_length = 0
		self.exon_list = []
		self.exon_number = 0
		self.FPKM = 0

class Exon:
	def __init__(self):
		self.exon_id = ""
		self.exon_start = 0
		self.exon_end = 0
		self.exon_length = 0

def convert_file(file):
	with open(file, "r") as input_file:
		line = input_file.readline()
		while line != "":
			try:
				elements = line.strip().split()
				chr_name = elements[0]
				if chr_name not in chr_dict:
					chr_dict[chr_name] = Chr()
				element_type = elements[2]
				start = elements[3]
				end = elements[4]
				length = int(end) - int(start) + 1
				gene_id = elements[9][1:-2]
				if gene_id not in chr_dict[chr_name].gene_dict:
					chr_dict[chr_name].gene_dict[gene_id] = ""
				transcript_id = elements[11][1:-2]
				FPKM = elements[13][1:-2]

				if element_type == "transcript":
					temp_transcript = Transcript()
					temp_transcript.transcript_id = transcript_id
					temp_transcript.transcript_start = start
					temp_transcript.transcript_end = end
					temp_transcript.transcript_length = length
					temp_transcript.FPKM = FPKM

					chr_dict[chr_name].transcript_dict[temp_transcript.transcript_id] = temp_transcript

					line = input_file.readline()
					elements = line.strip().split()
					chr_name = elements[0]

					element_type = elements[2]
					while line != "" and element_type == "exon":
						temp_exon = Exon()
						temp_exon.exon_start = elements[3]
						temp_exon.exon_end = elements[4]
						temp_exon.exon_length = int(temp_exon.exon_end) - int(temp_exon.exon_start) + 1
						temp_exon.exon_id = elements[13][1:-2]
						temp_exon.FPKM = elements[15][1:-2]
						chr_dict[chr_name].transcript_dict[temp_transcript.transcript_id].exon_list.append(temp_exon)


						line = input_file.readline()
						if line != "":
							elements = line.strip().split()
							element_type = elements[2]
				else:
					print "exon not within transcript"
					sys.exit(1)
			except:
				print "error in ", line

def process_data():

	for chr_name, chr_data in chr_dict.iteritems():
		chr_data.gene_number = len(chr_data.gene_dict.keys())
		chr_data.transcript_number = len(chr_data.transcript_dict.keys())
		for transcript_id, transcript_data in chr_data.transcript_dict.iteritems():
			transcript_data.exon_number = len(transcript_data.exon_list)
			chr_data.transcript_length_list.append(transcript_data.transcript_length)
			chr_data.exon_number_list.append(transcript_data.exon_number)

			for exon in transcript_data.exon_list:
				chr_data.exon_length_list.append(exon.exon_length)

	for chr_name, chr_data in chr_dict.iteritems():
		chr_data.max_transcript_length = max(chr_data.transcript_length_list)
		chr_data.min_transcript_length = min(chr_data.transcript_length_list)
		chr_data.avg_transcript_length = get_average(chr_data.transcript_length_list)
		chr_data.stdev_transcript_length = get_stdev(chr_data.transcript_length_list)

		chr_data.max_exon_number = max(chr_data.exon_number_list)
		chr_data.min_exon_number = min(chr_data.exon_number_list)
		chr_data.avg_exon_number = get_average(chr_data.exon_number_list)
		chr_data.stdev_exon_number = get_stdev(chr_data.exon_number_list)

		chr_data.max_exon_length = max(chr_data.exon_length_list)
		chr_data.min_exon_length = min(chr_data.exon_length_list)
		chr_data.avg_exon_length = get_average(chr_data.exon_length_list)
		chr_data.stdev_exon_length = get_stdev(chr_data.exon_length_list)

def output_data(file):
	with open(file[:-4] + "_statistics.txt", "w") as output_file:
		print >> output_file, "chr", "genes_number", "transcript_number", "max_transcript_length", "min_transcript_length",   \
							"avg_transcript_length", "stdev_transcript_length", "max_#_exons", "min_#_exons",   \
							"avg_#_exons", "stdev_#_exons", "max_exons_length", "min_exons_length",   \
							"avg_exons_length", "stdev_exons_length"
		for i in range(1, 10):
			chr = "chr" + str(i)
			chr_data = chr_dict[chr]
			print >> output_file, chr, chr_data.gene_number, chr_data.transcript_number, chr_data.max_transcript_length, \
				chr_data.min_transcript_length, chr_data.avg_transcript_length, chr_data.stdev_transcript_length, \
			chr_data.max_exon_number, chr_data.min_exon_number, chr_data.avg_exon_number, \
			chr_data.stdev_exon_number, chr_data.max_exon_length, chr_data.min_exon_length, \
			chr_data.avg_exon_length, chr_data.stdev_exon_length

def get_args():
	desc = ""
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

	global chr_dict
	chr_dict = {}

	start_time = time.time()
	convert_file(file)
	process_data()
	output_data(file)
	print chr_dict.keys()
	"""

	for chr in chr_dict.keys():
		print chr
		print len(chr_dict[chr].gene_dict.keys())
		print len(chr_dict[chr].transcript_dict.keys())
		print chr_dict[chr].max_transcript_length
		print chr_dict[chr].stdev_transcript_length


		print chr_dict[chr].max_exon_number
		print chr_dict[chr].avg_exon_number
		print chr_dict[chr].max_exon_length
		print chr_dict[chr].stdev_exon_length
	"""

	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"

	
	