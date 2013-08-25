#!/usr/bin/python
#######################################################################################
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser

from tools import *

def ref_convert(file_name, maf_upper_bound, maf_lower_bound):
	global medium_maf_num
	
	#maf_upper_bound = 0.5
	#maf_lower_bound = 0.1
	total_haplotype_number = 0
	minor_allele_dict = {}
	hap_ref_dict = load_raw_data(file_name, raw_data_format)[1]
	"""
	print list_to_line(hap_ref_dict[36918909])
	print list_to_line(hap_ref_dict[36916028])
	print list_to_line(hap_ref_dict[36918576])
	"""
	key = hap_ref_dict.iterkeys().next()
	total_haplotype_number = len(hap_ref_dict[key][2:])
	#print "total_haplotype_number", total_haplotype_number
	maf_upper_num = int(total_haplotype_number * maf_upper_bound)
	maf_lower_num = int(total_haplotype_number * maf_lower_bound)
	medium_maf_num = int(total_haplotype_number * 0.5)
	#print maf_upper_num, maf_lower_num
	
	for pos, ref in hap_ref_dict.iteritems():
		alleles = ref[2:]
		minor_allele_number = len(alleles)
		unique_alleles = set(alleles)
		alleles_number = len(unique_alleles)
		if alleles_number > 2:		
			print "more than two alleles in: ", pos
			sys.exit(1)
		elif alleles_number == 1:
			pass
		else:
			maf_list = []
			for allele in unique_alleles:
				this_allele_number = alleles.count(allele)
				maf_list.append((allele, this_allele_number))
				#print allele, this_allele_number
				#minor_allele_number = this_allele_number if this_allele_number < minor_allele_number else minor_allele_number
			major_allele = maf_list[0][0]
			minor_allele = maf_list[1][0]
			minor_allele_number = maf_list[1][1]
			if maf_list[0][1] < maf_list[1][1]:
				major_allele = maf_list[1][0]
				minor_allele = maf_list[0][0]
				minor_allele_number = maf_list[0][1]
			if minor_allele_number <= maf_upper_num and minor_allele_number >= maf_lower_num:
				if minor_allele_number not in minor_allele_dict:
					minor_allele_dict[minor_allele_number] = []
				"""To deal with maf 0.5 snps. problem is each snp will be duplicated and will form a pair with itself"""
				if minor_allele_number == medium_maf_num:
					temp_ref = copy.deepcopy(ref)
					for i in range(len(temp_ref)):
						if temp_ref[i] == major_allele:
							temp_ref[i] = "1"
						if temp_ref[i] == minor_allele:
							temp_ref[i] = "2"
					minor_allele_dict[minor_allele_number].append(temp_ref)
					
					temp_ref = copy.deepcopy(ref)
					for i in range(len(temp_ref)):
						if temp_ref[i] == major_allele:
							temp_ref[i] = "2"
						if temp_ref[i] == minor_allele:
							temp_ref[i] = "1"
					minor_allele_dict[minor_allele_number].append(temp_ref)			
				else:
					for i in range(len(ref)):
						temp_ref = ref
						if temp_ref[i] == major_allele:
							temp_ref[i] = "1"
						if temp_ref[i] == minor_allele:
							temp_ref[i] = "2"
					minor_allele_dict[minor_allele_number].append(temp_ref)
	print len(minor_allele_dict), " different alleles"
	return minor_allele_dict
		
def get_cluster(file_name, maf_upper_bound, maf_lower_bound):
	minor_allele_dict = ref_convert(file_name, maf_upper_bound, maf_lower_bound)
	ref_cluster_dict = {}
	for num, snp_list in minor_allele_dict.iteritems():
		ref_cluster_dict[num] = []
		while len(snp_list) > 1:
			cluster_dict = {}
			current_snp = snp_list.pop(0)
			cluster_dict[current_snp[1]] = current_snp
			for snp in snp_list:
				"""to remove snps pair with itself"""
				if current_snp[2:] == snp[2:] and snp[1] not in cluster_dict:		  
					cluster_dict[snp[1]] = snp
			if len(cluster_dict) > 1:
				#print "cluster size: ", len(cluster_dict)
				for pos, snps in cluster_dict.iteritems():
					#print snps[0], snps[1]# list_to_line(snps[2:])
					if snps in snp_list:		
						snp_list.remove(snps)
				ref_cluster_dict[num].append(cluster_dict)	
	#print "size old", len(ref_cluster_dict)
	ref_cluster_dict = clean_half_maf_cluster(ref_cluster_dict)
	#print "size new", len(ref_cluster_dict)
	ref_cluster_dict = convert_to_allele(file_name, ref_cluster_dict)
	#print_cluster(ref_cluster_dict)
	return ref_cluster_dict

def convert_to_allele(file_name, ref_cluster_dict):
	""""convert 1,2 back to ATCG"""
	hap_ref_dict = load_raw_data(file_name, raw_data_format)[1]
	for num, snp_list in ref_cluster_dict.iteritems():
		for i in range(len(snp_list)):
			for pos, snps in snp_list[i].iteritems():
				snp_list[i][pos] = hap_ref_dict[int(pos)]
				#print hap_ref_dict[int(pos)]
				#print snp_list[i][pos]
		ref_cluster_dict[num] = snp_list
	return ref_cluster_dict

def clean_half_maf_cluster(ref_cluster_dict):
	
	if medium_maf_num in ref_cluster_dict:
		medium_maf_list = ref_cluster_dict[medium_maf_num]
		for first_cluster_dict in medium_maf_list:
			for second_cluster_dict in medium_maf_list:
				if first_cluster_dict.keys() == second_cluster_dict.keys():
					medium_maf_list.remove(second_cluster_dict)
		ref_cluster_dict[medium_maf_num] = medium_maf_list
	
	return ref_cluster_dict

def print_cluster(ref_cluster_dict):
	total_cluster_number = 0
	for num, snp_list in ref_cluster_dict.iteritems():
		print "minor allele number: ", num, len(snp_list)
		for cluster_dict in snp_list:
			#print "cluster size: ", len(cluster_dict)
			total_cluster_number += len(cluster_dict)
			for pos, snps in cluster_dict.iteritems():
				#print snps[0], snps[1] # list_to_line(snps[2:])
				pass
	print "total_cluster_number", total_cluster_number
		
def get_args():
	desc="calculate_maf"
	usage = "" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-i", "--ref", type="string", dest="ref_name",help = "Input ref file name", default="null")
	parser.add_option("-u", "--upper", type="float", dest="maf_upper_bound",help = "maf_upper_bound", default=0.5)
	parser.add_option("-l", "--lower", type="float", dest="maf_lower_bound",help = "maf_lower_bound", default=0.1)
	(options, args) = parser.parse_args()
	if options.ref_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__=='__main__':
	options = get_args()
	file_name = options.ref_name
	maf_upper_bound = options.maf_upper_bound
	maf_lower_bound = options.maf_lower_bound
	get_cluster(file_name, maf_upper_bound, maf_lower_bound)
