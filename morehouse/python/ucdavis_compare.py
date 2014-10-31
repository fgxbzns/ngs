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

	cons_dict = read_cons_data(path + "cons.txt")

	davis_dict = read_davis_data(path + "davis.txt")

	global nature_geneid_line_dict
	nature_pos_dict, nature_geneid_dict, nature_geneid_line_dict = read_nature_data(path + "nature12817-s3_sorted_5.csv")

	compare_data(cons_dict, davis_dict, nature_pos_dict, nature_geneid_dict)

	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"
	
	
	
	
	

