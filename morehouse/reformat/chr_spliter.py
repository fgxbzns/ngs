#!/usr/bin/python

import os
import sys
import time
from optparse import OptionParser

currentPath = os.getcwd() + '/'

def seprate_ref(chr_file):
	chr_file_name = chr_file[:chr_file.find('.')].strip()
	print "Separate: ", chr_file

	output_in_process = False
	with open(chr_file, "r") as chr_input_file:
		for line in chr_input_file:
			if line.startswith(">") and output_in_process:
				chr_output_file.close()
				print "run time for current chr: ", round(time.time() - chr_start_time, 2), "s"
				output_in_process = False
			if line.startswith(">") and not output_in_process:
				output_in_process = True
				chr_name = line[(line.find('>') + 1):].strip()
				print "****** Start processing: ", chr_name
				chr_output_file = open(chr_file_name + "_" + chr_name + ".fa", "w")
				chr_start_time = time.time()
				chr_output_file.write(line.strip() + "\n")
			else:
				chr_output_file.write(line.strip() + "\n")

def get_args():
	desc = "Separate the reference fasta file into each chromosome"
	usage = "chr_spliter -i reference.fa"

	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--aFile", type="string", dest="aFile", help="Input file Name", default="null")
	(options, args) = parser.parse_args()
	if options.aFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	start_time = time.time()
	options = get_args()

	chr_file = options.aFile
	seprate_ref(chr_file)

	print "Total run time: ", round(time.time() - start_time, 2), "s"
