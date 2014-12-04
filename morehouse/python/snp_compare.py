#!/usr/bin/python
# ######################################################################################
# Guoxing Fu Nov 24, 2013
# sun pei data
#######################################################################################

import gzip
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


def load_frequency(file_name):
	"""
	this one load all data into dict
	"""
	maf_chr_dict = {}
	with gzip.open(file_name, "rb") as input_file:
		for line in input_file:
			elements = line.strip().split()
			rs = elements[0]
			#chr = elements[1]
			#pos = elements[2]
			#frequence = elements[-7:]
			#if rs not in maf_chr_dict:
			maf_chr_dict[rs] = line.strip()
	return maf_chr_dict


def load_frequency(file_name, rs_list):
	"""
	this one only load the rs in rs_list into dict
	"""
	maf_chr_dict = {}
	with gzip.open(file_name, "rb") as input_file:
		for line in input_file:
			elements = line.strip().split()
			rs = elements[0]
			#chr = elements[1]
			#pos = elements[2]
			#frequence = elements[-7:]
			#if rs not in maf_chr_dict:
			if rs in rs_list:
				maf_chr_dict[rs] = line.strip()
	#print len(maf_chr_dict)
	return maf_chr_dict


def get_maf(rs_file):
	"""
	get af from hapmap data
	"""
	order_list = []  # to keep the order of snps
	chr_rs_dict = {}

	# to load snp data
	with open(rs_file, "r") as input_file:
		for line in input_file:
			if not line.startswith("id"):
				elements = line.strip().split()
				chr = elements[1]
				if chr not in chr_rs_dict:
					chr_rs_dict[chr] = {}
				rs = elements[2]
				order_list.append((chr, rs))
				if rs not in chr_rs_dict[chr]:
					chr_rs_dict[chr][rs] = ""

	maf_path = "/home/guoxing/disk2/sunpei/chb_chs/"
	#file_name = "allele_freqs_chr1_CHB_phase3.2_nr.b36_fwd.txt.gz"

	maf_dict = {}
	for data in order_list:
		chr, rs = data[0], data[1]
		if chr not in maf_dict:
			file_name = maf_path + "allele_freqs_chr" + chr + "_CHB_phase3.2_nr.b36_fwd.txt.gz"
			#print file_name
			maf_dict[chr] = load_frequency(file_name, chr_rs_dict[chr].keys())

		if rs in chr_rs_dict[chr]:
			print chr, rs,
			if rs in maf_dict[chr]:
				elements = maf_dict[chr][rs].split()
				print elements[2], list_to_line(elements[-7:])
			else:
				print ""

def find_index(file_1, file_2):
	"""
	find index of population id from 1000g data
	"""
	id_list = []
	with open(file_2, "r") as input_file:
		for line in input_file:
			id_list = line.strip().split()
	print "len(id_list)", len(id_list)

	population_list = []
	with open(file_1, "r") as input_file:
		for line in input_file:
			population_list.append(line.strip().split()[0])
	print "len(population_list)", len(population_list)

	id_wo_genotype_list = []
	for id in population_list:
		if id in id_list:
			print id_list.index(id),
		else:
			id_wo_genotype_list.append(id)
	print ""
	print id_wo_genotype_list

def get_rsID_order(rsID_file):
	"""
	get a list of rsID in its original order
	"""
	order_list = []  # to keep the order of snps rsID
	with open(rsID_file, "r") as input_file:
		for line in input_file:
			if line.strip() != "":
				order_list.append(line.strip().split()[0])
	return order_list

def get_af(vcf_file, rsID_file):
	"""
    for 1000g chb chs data
    :return:
    """
	vcf_dict = {}

	with open(vcf_file, "r") as input_file:
		for line in input_file:
			if not line.startswith("id"):
				elements = line.strip().split()
				data = list_to_line(elements[9:])
				num_0 = data.count("0")
				num_1 = data.count("1") if data.count("1") != 0 else data.count("2")

				#print elements[0], elements[1], elements[2], elements[3], num_0, round(float(num_0) / (num_0 + num_1),4), \
				#	elements[4], num_1, round(float(num_1) / (num_0 + num_1), 4), num_0 + num_1
				vcf_dict[elements[2]] = [elements[0], elements[1], elements[2], elements[3], num_0, round(float(num_0) / (num_0 + num_1),4), \
					elements[4], num_1, round(float(num_1) / (num_0 + num_1), 4), num_0 + num_1]

	print "chr \t pos \t rsID \t ref_allele \t ref_num \t ref_AF \t other_allele \t other_num \t other_AF \t total_sample"

	rsID_order_list = get_rsID_order(rsID_file)
	for rsID in rsID_order_list:
		if rsID in vcf_dict:
			print list_to_line(vcf_dict[rsID])
		else:
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

	get_af("chr_rsID.txt", "cancer_SNP_20141203_dup_not_removed.txt")

	"""
	#get_maf("/home/guoxing/disk2/sunpei/cancer_panel20141124.txt")

	path = "/home/guoxing/disk2/sunpei/1000g/"
	chb_file = path + "chb.txt"
	chs_file = path + "chs.txt"
	id_file = path + "id_chr1.txt"
	print "chb id index"
	find_index(chb_file, id_file)

	print "chs id index"
	find_index(chs_file, id_file)
	"""

	"""
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


	snp_list_dict[gwas_dcit] = read_data_dict(path + "gwascatalog_all.txt", 0)
	#print "dict", len(snp_list_dict[gwas_list])
	snp_list_dict[gwas_list] = read_data_list(path + "gwascatalog_all.txt", 0)
	#print "list", len(snp_list_dict[gwas_list])
	"""

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
	"""

	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"






