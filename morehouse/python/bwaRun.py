#!/usr/bin/python

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

currentPath = os.getcwd() + '/'

class parameters:
	def __init__(self):
		self.chr_name = ""
		self.fastq_file_1 = ""
		self.fastq_file_name_1 = ""
		self.sam_file_name = ""

		self.fastq_file_2 = ""
		self.fastq_file_name_2 = ""

		self.ref_path = "/home/guoxing/disk2/UCSC_hg18_index_lm/"
		self.ref_file_name = ""

def bwa_sai(fastq_file):
	print "bwa aln running", parameter.fastq_file
	fastq_file_name = fastq_file[:(fastq_file.find('fastq')-1)].strip()
	sai_file = fastq_file_name + "_" + parameter.chr_name + ".sai"
	sam_file = fastq_file_name + "_" + parameter.chr_name + ".sam"

	bwa_sai = "bwa" + " aln " + parameter.ref_file_name + " " + fastq_file + " > " + sai_file
	print bwa_sai
	bwa_sai_Process = subprocess.Popen(bwa_sai, shell=True)
	bwa_sai_Process.wait()

def bwa_samse(fastq_file):
	print "bwa samse running", parameter.fastq_file
	fastq_file_name = fastq_file[:(fastq_file.find('fastq')-1)].strip()
	sai_file = fastq_file_name + "_" + parameter.chr_name + ".sai"
	sam_file = fastq_file_name + "_" + parameter.chr_name + ".sam"

	bwa_sam = "bwa" + " samse " + parameter.ref_file_name + " " + sai_file + " " + fastq_file + " > " + sam_file
	print bwa_sam
	bwa_sam_Process = subprocess.Popen(bwa_sam, shell=True)
	bwa_sam_Process.wait()

def bwa_bwasw(fastq_file):
	print "bwa samse running", parameter.fastq_file
	fastq_file_name = fastq_file[:(fastq_file.find('fastq')-1)].strip()
	sai_file = fastq_file_name + "_" + parameter.chr_name + ".sai"
	sam_file = fastq_file_name + "_" + parameter.chr_name + ".sam"

	bwa_long = "bwa" + " bwasw " + parameter.ref_file_name + " " + fastq_file + " > " + sam_file
	print bwa_long

	bwa_long_Process = subprocess.Popen(bwa_long, shell=True)
	bwa_long_Process.wait()

def bwa_run():
	options = get_args()
	parameter.fastq_file = options.inputfile1
	parameter.fastq_file_name = parameter.fastq_file[:(parameter.fastq_file.find('fastq')-1)].strip()
	parameter.chr_name = options.chrName
	mode = options.mode

	parameter.ref_file_name = parameter.ref_path + "hg18chr.fa" \
							if parameter.chr_name == "all" else parameter.ref_path + "hg18chr_" + parameter.chr_name + ".fa"
	#if mode == "single":
	if True:
		bwa_sai(parameter.fastq_file)
		bwa_samse(parameter.fastq_file)

def get_args():
	desc = "bwaRun"
	usage = "bwaRun -i input_fastq"
	parser = OptionParser(usage=usage, description=desc) 
	parser.add_option("-i", "--inputfile1", type="string", dest="inputfile1", help="Input File Name", default="null")
	parser.add_option("-n", "--inputfile2", type="string", dest="inputfile2", help="Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName",help="Input chr Name", default="all")
	parser.add_option("-m", "--mode", type="string", dest="mode", help="", default="null")

	(options, args) = parser.parse_args()
	if options.inputfile1 == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	global parameter
	parameter = parameters()

	start_time = time.time()
	bwa_run()
	print "a"

	print "elapse_time is: ", round(time.time() - start_time, 3), "s"


