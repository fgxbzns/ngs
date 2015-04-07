#!/usr/bin/python

#######################################################################################
# Author: Guoxing Fu
# generate files for genome laser
#######################################################################################



import os, glob, subprocess, random, operator, time, copy
from optparse import OptionParser

from tools import *
from seed_std_compare import seed_std_compare
from hifiAccuCheck_v2 import hifiAccuCheck


class refs:
	def __init__(self):
		self.seed_title = ""
		self.seed_dict = {}
		self.geno_title = ""
		self.geno_dict = {}
		self.ref_title = ""
		self.ref_dict = {}


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
		#print position
		alleles = hap_ref_dict[position]
		#print "aaaaaaaaaaa", alleles
		alleles = alleles[2:]
		unique_alleles = list(set(alleles))
		n_alleles = len(unique_alleles)
		if n_alleles == 1:
			ref_homo_dict[position] = list_to_line(unique_alleles)

	# to remove snps that are conflict in ref and geno
	for position, snp in geno_dict.iteritems():
		if position in hap_ref_dict:
			exist_in_ref = False
			geno_A = snp[0]
			geno_B = snp[1]
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

	# print "geno_x_dict: ", len(geno_x_dict)
	# print "geno_n_dict: ", len(geno_n_dict)
	# print "geno_ref_not_consistent: ", len(geno_ref_not_consistent)
	# print "ref_homo", len(ref_homo_dict)
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
			geno_allele = geno_dict[position]
			if geno_allele[0] == geno_allele[1]:
				seed_homo_dict[position] = snp
			else:
				seed_hetero_dict[position] = snp
		else:
			seed_hetero_dict[position] = snp
	return (seed_homo_dict, seed_hetero_dict)


def output_hap(file_name, parameter, seed_dict):
	with open(file_name, "w") as output:
		print >> output, ref.seed_title
		for pos in sorted(seed_dict.keys()):
			print >> output, parameter.rsID_dict[pos], pos, seed_dict[pos][0]

def output_geno(file_name, parameter, geno_dict):
	with open(file_name, "w") as output:
		print >> output, ref.geno_title
		for pos in sorted(geno_dict.keys()):
			print >> output, parameter.rsID_dict[pos], pos, geno_dict[pos][0]+geno_dict[pos][1]

def output_ref(file_name, ref_dict):
	with open(file_name, "w") as output:
		print >> output, list_to_line(ref.ref_title)
		for pos in sorted(ref_dict.keys()):
			print >> output, ref_dict[pos].strip()

def make_hifi_files(ref, parameter):

	seed_dict = ref.seed_dict
	geno_dict = ref.geno_dict
	hap_ref_dict = ref.ref_dict

	hifi_seed_dict = {}
	hifi_geno_dict = {}
	hifi_ref_dict = {}


	last_seed_position = sort_dict_by_key(seed_dict)[-1][0]

	hap_ref_sorted_list = sort_dict_by_key(hap_ref_dict)

	for snp in hap_ref_sorted_list:
		position = snp[0]
		if position <= int(last_seed_position):
			hifi_ref_dict[position] = list_to_line(hap_ref_dict[position])
			if position in seed_dict:
				hifi_seed_dict[position] = seed_dict[position][0]
			elif position in geno_dict and geno_dict[position][0] == geno_dict[position][1]:
				hifi_seed_dict[position] = geno_dict[position][0]
			if position in geno_dict:
				hifi_geno_dict[position] = geno_dict[position]


	output_hap("haplotype.txt", parameter, hifi_seed_dict)
	output_geno("genotype.txt", parameter, hifi_geno_dict)
	output_ref("refHaplos.txt", hifi_ref_dict)


def find_index(a, b):
	index_list = []
	for str_b in b:
		for id in a:
			if id in str_b:
				index_list.append(b.index(str_b))
	return index_list

def remove_element(index_list, list):
	for index in sorted(index_list, reverse=True):
		del list[index]


def remove_hap_in_ref(id, parameter):
	ref_title = copy.deepcopy(parameter.ori_ref_title)
	ref_dict = copy.deepcopy(parameter.ori_ref_dict)

	id_list = []
	id_list.append(id)
	if parameter.person_dict[id].father != "N/A":
		id_list.append(parameter.person_dict[id].father)
	if parameter.person_dict[id].mather != "N/A":
		id_list.append(parameter.person_dict[id].mather)
	print id, id_list

	index_list = find_index(id_list, ref_title)
	remove_element(index_list, ref_title)

	for pos in ref_dict.keys():
		remove_element(index_list, ref_dict[pos])

	return ref_title, ref_dict


def refMerger(id, parameter):

	person = parameter.person_dict[id]

	global ref
	ref = refs()
	#print id
	ref.ref_title, ref.ref_dict = remove_hap_in_ref(id, parameter)

	ref.seed_title = "rsID pos " + id + "_A"

	for pos in person.haplotype.keys():
		hap = person.haplotype[pos]
		if hap[0] != "N" and hap[0] != "X" and hap[0] != "":
			ref.seed_dict[pos] = hap[0]

	ref.geno_title = "rsID pos " + id
	ref.geno_dict = person.genotype_dict

	ref.ref_dict = ref_preprocess(ref.geno_dict, ref.ref_dict)

	make_hifi_files(ref, parameter)


def refMerger_laser2(id, parameter):

	person = parameter.person_dict[id]

	global ref
	ref = refs()
	#print id
	ref.ref_title, ref.ref_dict = load_hap_ref_data_single("l3_ref.txt")
	#print ref.ref_title

	ref.seed_title = "rsID pos " + id + "_A"

	for pos in person.haplotype.keys():
		hap = person.haplotype[pos]
		if hap[0] != "N" and hap[0] != "X" and hap[0] != "":
			ref.seed_dict[pos] = hap[0]

	ref.geno_title = "rsID pos " + id
	ref.geno_dict = person.genotype_dict

	ref.ref_dict = ref_preprocess(ref.geno_dict, ref.ref_dict)

	make_hifi_files(ref, parameter)