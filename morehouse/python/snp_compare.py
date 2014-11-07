#!/usr/bin/python
#######################################################################################
# Guoxing Fu July 7, 2013
# for rna-seq data
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

def read_data_list(file, index):
	list = []
	with open(file, "r") as input_file:
			for line in input_file:
				elements = line.strip().split()
				rs = elements[index]
				list.append(rs)
	return list

def read_data_dict(file, index):
	dict = {}
	with open(file, "r") as input_file:
			for line in input_file:
				elements = line.strip().split()
				rs = elements[index]
				dict[rs] = line.strip()
	return dict

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

def find_common(list1, list2, snp_list_dict):
	common_list = [x for x in snp_list_dict[list1] if x in snp_list_dict[list2]]
	print list1 + " in " + list2, len(common_list)
	#print common_list
	for rs in common_list:
		print rs,
	print ""

def output_common(list1, list2, snp_list_dict):
	with open(list1 + "_in_" + list2 + ".txt", "w") as common:
		with open(list1 + "_not_in_" + list2 + ".txt", "w") as not_common:
			with open(list1 + "_" + list2 + "_notRS.txt", "w") as not_rs:
				for rsID in snp_list_dict[list1].keys():
					if rsID.startswith("rs"):
						if rsID in snp_list_dict[list2]:
							print >> common, snp_list_dict[list1][rsID]
						else:
							print >> not_common, snp_list_dict[list1][rsID]
					else:
						print >> not_rs, snp_list_dict[list1][rsID]

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
	#file = options.folder_name
	path = "/home/guoxing/disk2/snp_test/"
	start_time = time.time()

	snp_list_dict = {}

	sun_list = "sun_list"
	igene = "igene"
	cancer_list = "cancer_list"
	k500_list = "k500_list"
	gwas_list = "gwas_list"
	gwas_dcit = "gwas_dict"

	snp_list_dict[sun_list] = read_data_list(path + "sun.txt", -1)
	snp_list_dict[igene] = read_data_list(path + "igene.txt", -1)
	snp_list_dict[cancer_list] = read_data_list(path + "Cancer_SNP_Panel_Annotation.txt", 0)
	snp_list_dict[k500_list] = read_data_dict(path + "OncoArray-500K_B_GeneAnnotation.txt", 0)

	"""
	snp_list_dict[gwas_dcit] = read_data_dict(path + "gwascatalog_all.txt", 0)
	#print "dict", len(snp_list_dict[gwas_list])
	snp_list_dict[gwas_list] = read_data_list(path + "gwascatalog_all.txt", 0)
	#print "list", len(snp_list_dict[gwas_list])
	"""


	with open("gwascatalog_all.txt", "r") as input:
		with open("gwas_in_k500.txt", "w") as common:
			with open("gwas_notin_k500.txt", "w") as not_common:
				with open("gwas_notRS.txt", "w") as not_rs:
					for line in input:
						line = line.strip()
						if line.startswith("rs"):
							element = line.strip().split()
							rsID_list = element[0].split(",")
							for rsID in rsID_list:
								if rsID in snp_list_dict[k500_list]:
									print >> common, line
								else:
									print >> not_common, line
						else:
							print >> not_rs, line

	print "sun_list", len(snp_list_dict[sun_list])
	print "igene", len(snp_list_dict[igene])
	print "cancer_list", len(snp_list_dict[cancer_list])
	print "k500_list", len(snp_list_dict[k500_list])
	#print "gwas_list", len(snp_list_dict[gwas_list])

	#find_common(sun_list, igene, snp_list_dict)
	#find_common(sun_list, igene, snp_list_dict)
	#for list1 in (sun_list, igene, cancer_list, k500_list):
	list1 = sun_list
	for list2 in (sun_list, igene, cancer_list, k500_list):
		if list1 != list2:
			#find_common(list1, list2, snp_list_dict)
			pass
	#output_common(gwas_list, k500_list, snp_list_dict)

	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"
	
	
	
	
	

