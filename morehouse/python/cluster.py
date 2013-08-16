#!/usr/bin/python
#######################################################################################
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *

def ref_convert(file_name):
	maf_upper_bound = 0.5
	maf_lower_bound = 0.5
	total_person_number = 0
	minor_allele_dict = {}
	hap_ref_dict = load_raw_data(file_name, raw_data_format)[1]
	print list_to_line(hap_ref_dict[36918909])
	print list_to_line(hap_ref_dict[36916028])
	print list_to_line(hap_ref_dict[36918576])

	#total_person_number = len(list(hap_ref_dict)[0][1])
	total_person_number = 372
	maf_upper_num = int(total_person_number * maf_upper_bound)
	maf_lower_num = int(total_person_number * maf_lower_bound)
	print maf_upper_num, maf_lower_num
	
	for pos, ref in hap_ref_dict.iteritems():
		alleles = ref[2:]
		minor_allele_number = len(alleles)
		unique_alleles = set(alleles)
		alleles_number = len(unique_alleles)
		if alleles_number == 1 or alleles_number > 2:
			pass
			#print "error in: ", pos
			#sys.exit(1)
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
				for i in range(len(ref)):
					if ref[i] == major_allele:
						ref[i] = "1"
					if ref[i] == minor_allele:
						ref[i] = "2"
#				print ref
				minor_allele_dict[minor_allele_number].append(ref)
	print len(minor_allele_dict)
	#for num, snp_list in minor_allele_dict.iteritems():
    #           print num, len(snp_list)
	return minor_allele_dict
		
def find_cluster(minor_allele_dict):
	total_cluster_number = 0
	for num, snp_list in minor_allele_dict.iteritems():
		print num, len(snp_list)
		while len(snp_list) > 1:
			cluster_list = []
			current_snp = snp_list.pop(0)
			cluster_list.append(current_snp)
			for snp in snp_list:
				if current_snp[2:] == snp[2:]:
					cluster_list.append(snp)
#					snp_list.remove(snp)
			if len(cluster_list) > 1:
				total_cluster_number += len(cluster_list)
				print "cluster size: ", len(cluster_list)
				for snps in cluster_list:
					print snps[0], snps[1]# list_to_line(snps[2:])
					if snps in snp_list:		
						snp_list.remove(snps)
	print total_cluster_number

def get_args():
	desc="calculate_maf"
	usage = "" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-i", "--ref", type="string", dest="ref_name",help = "Input ref file name", default="null")
	(options, args) = parser.parse_args()
	if options.ref_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__=='__main__':
	options = get_args()
	file_name = "refHaplos.txt"
	minor_allele_dict = ref_convert(file_name)
	find_cluster(minor_allele_dict)
