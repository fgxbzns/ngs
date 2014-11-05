#!/usr/bin/python
#######################################################################################
# Guoxing Fu July 7, 2013
# for rna-seq data
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

def read_data(file, index):
	list = []
	with open(file, "r") as input_file:
			for line in input_file:
				elements = line.strip().split()
				rs = elements[index]
				list.append(rs)
	return list




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

	snp_list_dict[sun_list] = read_data(path + "sun.txt", -1)
	snp_list_dict[igene] = read_data(path + "igene.txt", -1)
	snp_list_dict[cancer_list] = read_data(path + "Cancer_SNP_Panel_Annotation.txt", 0)
	snp_list_dict[k500_list] = read_data(path + "OncoArray-500K_B_GeneAnnotation.txt", 0)

	print "sun_list", len(snp_list_dict[sun_list])
	print "igene", len(snp_list_dict[igene])
	print "cancer_list", len(snp_list_dict[cancer_list])
	print "k500_list", len(snp_list_dict[k500_list])

	#find_common(sun_list, igene, snp_list_dict)
	#find_common(sun_list, igene, snp_list_dict)
	#for list1 in (sun_list, igene, cancer_list, k500_list):
	list1 = sun_list
	for list2 in (sun_list, igene, cancer_list, k500_list):
		if list1 != list2:
			find_common(list1, list2, snp_list_dict)

	"""
	sun_list = read_data(path + "sun.txt", -1)
	print "sun size", len(sun_list)
	print sun_list

	cancer_list = read_data(path + "Cancer_SNP_Panel_Annotation.txt", 0)
	print "1421_list size", len(cancer_list)
	#print cancer_list

	fk_list = read_data(path + "OncoArray-500K_B_GeneAnnotation.txt", 0)
	print "fk_list size", len(fk_list)

	igene_list = read_data(path + "igene.txt", -1)
	print "igene_list size", len(igene_list)

	sun_in_cancer = [x for x in sun_list if x in cancer_list]
	print "sun_in_cancer size", len(sun_in_cancer), sun_in_cancer

	sun_in_fk = [x for x in sun_list if x in fk_list]
	print "sun_in_fk size", len(sun_in_fk)

	sun_in_igene = [x for x in sun_list if x in igene_list]
	print "sun_in_igene size", len(sun_in_igene)

	cancer_in_fk = [x for x in cancer_list if x in fk_list]
	print "cancer_in_fk size", len(cancer_in_fk)
	"""




	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"
	
	
	
	
	

