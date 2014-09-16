#!/usr/bin/python

import sys
from optparse import OptionParser

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

def sam_process():
	options = get_args()
	sam_file = options.samFile
	sam_file_name = sam_file[:(len(sam_file)-4)]

	fop = open(sam_file_name+"_output.csv", "w")
	fp = open(sam_file, "r")
	for line in fp:
		keep_line = False
		if line.startswith("chr"):
			print >> fop, line.strip()
		else:
			elements = line.strip().split(',')
			#print len (elements)
			try:
				data_list = elements[24:]
				#print len (data_list)
				for data in data_list:
					if float(data) < 0.1:
						keep_line = True
				if keep_line:
					print >> fop, line.strip()
			except ValueError:
				pass
	fp.close()
	fop.close()	
	
if __name__=='__main__':
	sam_process()















 
