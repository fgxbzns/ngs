#!/usr/bin/python

# ######################################################################################
# Guoxing Fu Mar 28, 2013
#
#######################################################################################

import os
import time
import sys
from optparse import OptionParser


def chr_percentage(a_file_name):

	chr_dict = {}
	total_reads_number = 0

	with open(a_file_name, 'r') as input_file:
			with open(a_file_name + "_percentage", 'w') as output_file:
				for line in input_file:
					if not line.startswith("@"):
						elements = line.strip().split()
						chr_name = elements[2].strip()
						#insertion_size = elements[8].strip()
						if True:
							total_reads_number += 1
							if chr_name not in chr_dict:
								chr_dict[chr_name] = 1
							else:
								chr_dict[chr_name] += 1

				chr_list = [x for x in chr_dict.iteritems()]
				chr_list.sort(key=lambda x: x[1], reverse=True)  # sort by value in reverse order. Max first
				for chr in chr_list:
					print >> output_file, str(chr[0]) + "\t" + str(chr[1]) + "\t" + str(round(float(chr[1])*100/total_reads_number, 4))
				print >> output_file, "total_reads_number", total_reads_number

	#print "chr6", chr_dict["chr6"]
	print "total_reads_number", total_reads_number

def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--aFile", type="string", dest="aFile", help="Input File Name", default="null")
	(options, args) = parser.parse_args()
	if options.aFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	currentPath = os.getcwd() + '/'

	options = get_args()
	a_file_name = options.aFile
	start_time = time.time()
	chr_percentage(a_file_name)
	print "run time is: ", round((time.time() - start_time), 3), "s"