#!/usr/bin/python

# for zebra fish data

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser
from tools import *
from vcfCall_tools import *

def get_args():
	desc = "variation call"
	usage = "snpPick_fish -s sam_file -c chr -m update -d db_name -q qscore \n" \
	        "snpPick_fish -c chr -m output -b startLine -e endLine -d db_name -q qscore"
	parser = OptionParser(usage=usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile", help="Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode", help="mode", default="null")
	parser.add_option("-b", "--startLine", type="string", dest="startLine", help="start line", default="null")
	parser.add_option("-e", "--endLine", type="string", dest="endLine", help="end line", default="null")
	parser.add_option("-d", "--dbname", type="string", dest="dbname", help="db name", default="null")
	parser.add_option("-q", "--qscore", type="int", dest="qscore", help="qscore", default="40")
	parser.add_option("-p", "--pos", type="string", dest="posList", help="Input position list Name", default="null")
	(options, args) = parser.parse_args()
	if options.mode == "null" or options.chrName == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	start_time = time.time()

	options = get_args()
	global parameters
	parameters = parameter()
	mode = options.mode

	parameters.chr_name = options.chrName

	if (mode == "rat"):
		#parameters.quality_score_threshold = options.qscore
		parameters.quality_score_threshold = 20
		print "quality_score_threshold: ", parameters.quality_score_threshold
		parameters.second_largest_allele_depth_cutoff = 2

		# for jendai mouse data
		ref_path = "/home/guoxing/storage1/Reference_from_wenzhi/Rat/chr/"
		parameters.ref_file = ref_path + "rn5_" + parameters.chr_name + ".fa"
		print "ref file", parameters.ref_file

		if (mode == "rat"):
			parameters.sam_file = options.samFile
			parameters.sam_file_name = parameters.sam_file[:parameters.sam_file.find('.')]
			output_file_name = parameters.sam_file_name + "_qs" + str(parameters.quality_score_threshold) + ".txt"
			parameters.output_file = open(output_file_name, "w")

			snpPick(parameters)
			parameters.output_file.close()

	elif (mode == "exon"):
		parameters.quality_score_threshold = options.qscore
		#parameters.quality_score_threshold = 20
		print "quality_score_threshold: ", parameters.quality_score_threshold
		parameters.second_largest_allele_depth_cutoff = 2

		# for exon data
		ref_path = "/home/guoxing/disk2/UCSC_hg18_index_lm/"
		parameters.ref_file = ref_path + "hg18chr_" + parameters.chr_name + ".fa"
		print "ref file", parameters.ref_file

		if (mode == "exon"):
			parameters.sam_file = options.samFile
			parameters.sam_file_name = parameters.sam_file[:parameters.sam_file.find('.')]
			output_file_name = parameters.sam_file_name + "_qs" + str(parameters.quality_score_threshold) + ".txt"
			parameters.output_file = open(output_file_name, "w")

			snpPick(parameters)
			parameters.output_file.close()

	print "run time: ", round(time.time() - start_time, 3), "s"
	#print "snp_in_mimi", parameters.snp_in_mimi
	#print "total_mimi", parameters.total_mimi



