#!/usr/bin/python

#######################################################################################
# Author: Guoxing Fu
# generate files for genome laser
#######################################################################################



import os, glob, subprocess, random, operator, time
from optparse import OptionParser

from tools import *
from seed_std_compare import seed_std_compare
from hifiAccuCheck_v2 import hifiAccuCheck



def load_hap_ref_data_single(ref_file_name):
	raw_data_format_ref = "list"
	ref_title_info, hap_ref_dict = load_raw_data(ref_file_name, raw_data_format_ref)
	return ref_title_info, hap_ref_dict

def compare_geno_ref(geno_dict, hap_ref_dict):

	geno_x_dict = {}
	geno_n_dict = {}
	geno_ref_not_consistent = {}
	ref_homo_dict = {}

	# to remove homo snps in ref
	for position in hap_ref_dict.keys():
		alleles = hap_ref_dict[position]
		alleles = alleles[2:]
		unique_alleles = list(set(alleles))
		n_alleles = len(unique_alleles)
		if n_alleles == 1:
			ref_homo_dict[position] = list_to_line(unique_alleles)

	# to remove snps that are conflict in ref and geno
	for position, snp in geno_dict.iteritems():
		if position in hap_ref_dict:
			exist_in_ref = False
			geno_A = snp[2][0]
			geno_B = snp[2][1]
			alleles = hap_ref_dict[position]
			alleles = alleles[2:]
			unique_alleles = list(set(alleles))

			n_alleles = len(unique_alleles)
			if n_alleles == 0 or n_alleles > 2:
				pass
			else:
				try:
					if geno_A == geno_B:
						if n_alleles == 2:
							if geno_A == unique_alleles[0] or geno_A == unique_alleles[1]:
								exist_in_ref = True
						if n_alleles == 1:
							if geno_A == unique_alleles[0]:
								exist_in_ref = True
					else:
						if n_alleles == 2:
							if (geno_A == unique_alleles[0] and geno_B == unique_alleles[1]) or (geno_A == unique_alleles[1] and geno_B == unique_alleles[0]):
								exist_in_ref = True
							if n_alleles == 1:	# hetero_geno, homo_ref
								pass
					if not exist_in_ref:
						if geno_A == 'X':
							geno_x_dict[position] = list_to_line(snp)
						if geno_A == 'N':
							geno_n_dict[position] = list_to_line(snp)
						else:
							geno_ref_not_consistent[position] = list_to_line(snp)
				except:
					print position, unique_alleles, n_alleles
					pass

	print "geno_x_dict: ", len(geno_x_dict)
	print "geno_n_dict: ", len(geno_n_dict)
	print "geno_ref_not_consistent: ", len(geno_ref_not_consistent)
	print "ref_homo", len(ref_homo_dict)
	return ref_homo_dict, geno_ref_not_consistent, geno_n_dict

def dict_substract(large_dict, small_dict):
	return {index: value for index, value in large_dict.iteritems() if index not in small_dict}

def ref_preprocess(geno_dict, hap_ref_dict):
	ref_homo_dict, geno_ref_not_consistent, geno_n_dict = compare_geno_ref(geno_dict, hap_ref_dict)
	hap_ref_dict = dict_substract(hap_ref_dict, geno_ref_not_consistent)
	hap_ref_dict = dict_substract(hap_ref_dict, geno_n_dict)
	#hap_ref_dict = dict_substract(hap_ref_dict, ref_homo_dict)
	return hap_ref_dict

def group_seed(seed_dict, geno_dict):
	seed_homo_dict = {}
	seed_hetero_dict = {}
	for position, snp in seed_dict.iteritems():
		if position in geno_dict:  # pos in seed may not be in geno
			geno_allele = geno_dict[position][2]
			if geno_allele[0] == geno_allele[1]:
				seed_homo_dict[position] = snp
			else:
				seed_hetero_dict[position] = snp
		else:
			seed_hetero_dict[position] = snp
	return (seed_homo_dict, seed_hetero_dict)

def output_files(file_name, title_info, dict):
	outpuf_file = open(currentPath + file_name, "w")    # for hifi
	print >> outpuf_file, title_info
	sorted_list = sort_dict_by_key(dict)
	for element in sorted_list:
		print >> outpuf_file, element[1]
	outpuf_file.close()

def make_hifi_files(seed_dict, geno_dict, hap_ref_dict):
	common_snp_number = 0
	hifi_seed_dict = {}
	hifi_geno_dict = {}
	hifi_ref_dict = {}


	last_seed_position = sort_dict_by_key(seed_dict)[-1][0]
	#print "last_haplotype_position", last_seed_position
	hap_ref_sorted_list = sort_dict_by_key(hap_ref_dict)

	seed_homo_dict, seed_hetero_dict = group_seed(seed_dict, geno_dict)

	hap_ref_sorted_list = [x for x in hap_ref_sorted_list if x[0] not in seed_hetero_dict]


	hap_ref_size = len(hap_ref_dict)
	for i in range(int(remPercent*hap_ref_size)):
		if len(hap_ref_sorted_list) > 1:
			random_index = random.randrange(0, (len(hap_ref_sorted_list)-1))
			position = hap_ref_sorted_list[random_index][0]
			if position in hap_ref_dict:
				del hap_ref_dict[int(position)]
				del hap_ref_sorted_list[random_index]

	hap_ref_sorted_list = sort_dict_by_key(hap_ref_dict)

	for snp in hap_ref_sorted_list:
		position = snp[0]
		if position <= int(last_seed_position):
			hifi_ref_dict[position] = list_to_line(hap_ref_dict[position])
			if position in seed_dict:
				hifi_seed_dict[position] = list_to_line(seed_dict[position])
			elif position in geno_dict and geno_dict[position][2][0] == geno_dict[position][2][1]:
				hifi_seed_dict[position] = geno_dict[position][0] + " " + geno_dict[position][1] + " " + geno_dict[position][2][0]
			if position in geno_dict:
				hifi_geno_dict[position] = list_to_line(geno_dict[position])


	output_files("haplotype.txt", seed_title_info, hifi_seed_dict)
	output_files("genotype.txt", geno_title_info, hifi_geno_dict)
	output_files("refHaplos.txt", ref_title_info, hifi_ref_dict)

	seed_not_in_ref = dict_substract(seed_dict, hap_ref_dict)
	output_files("seed_not_in_ref.txt", seed_title_info, seed_not_in_ref)


def refMerger(haplotype_file, chr_name, remPercent):
	global seed_title_info
	global seed_dict
	global geno_title_info
	global geno_dict
	global hap_ref_dict
	global ref_title_info
	global raw_data_format_ref

	raw_data_format_ref = "list"

	ref_title_info, hap_ref_dict = load_hap_ref_data_single(chr_name)

	geno_title_info, geno_dict = load_raw_data(genotype_file, raw_data_format_ref)
	#print "genotype_file", genotype_file
	#print "total_geno_number: ", len(geno_dict)

	seed_file = haplotype_file
	seed_title_info, seed_dict = load_raw_data(seed_file, raw_data_format_ref)
	print "seed_file is :", seed_file

	print "hap_ref_dict before process", len(hap_ref_dict)
	hap_ref_dict = ref_preprocess(geno_dict, hap_ref_dict)
	print "hap_ref_dict after process", len(hap_ref_dict)

	make_hifi_files(remPercent)


def refMerge(person):

	geno_dict = person.genotype_dict
	hap_dict = person.haplotype
	#std_hap_dict = person.std_hap_dict
