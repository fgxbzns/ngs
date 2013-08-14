#!/usr/bin/python
#######################################################################################
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *


def ref_convert(file_name):
	hap_ref_dict = load_raw_data(file_name, raw_data_format)[1]
	for pos, ref in hap_ref_dict.iteritems():
		alleles = ref[2:]
		unique_alleles = set(alleles)
		alleles_number = len(unique_alleles)
		if n_alleles == 0 or n_alleles > 2:
			print "error in: ", position
			sys.exit(1)
		else:
			pass
			


def calculate_maf(file_name, position):
	hap_ref_dict = load_raw_data(file_name, raw_data_format)[1]
	alleles = hap_ref_dict[position]
	alleles = alleles[2:]
	#print "total allele: ", len(alleles)
	unique_alleles = set(alleles)
	n_alleles = len(unique_alleles)
	if n_alleles == 0 or n_alleles > 2:
		print "error in: ", position
		sys.exit(1)
	else:
		print position
		for allele in unique_alleles:
			print allele, alleles.count(allele), format((float(alleles.count(allele))/float(len(alleles))), "0.3f")

def get_args():
	desc="calculate_maf"
	usage = "" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-i", "--ref", type="string", dest="ref_name",help = "Input ref file name", default="null")
	parser.add_option("-p", "--pos", type="string", dest="position",help = "Input position", default="null")
	(options, args) = parser.parse_args()
	if options.ref_name == "null" or options.position == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__=='__main__':
	options = get_args()
	file_name = "refHaplos.txt"
	position = int(options.position)
	#hap_ref_dict = load_raw_data(file_name, raw_data_format)[1]
	calculate_maf(file_name, position)
	#calculate_maf(hap_ref_dict, position)















			
