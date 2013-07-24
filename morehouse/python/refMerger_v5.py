#!/usr/bin/python

#######################################################################################
# 
#######################################################################################
"""
keep only the snps that are available in ref data
assume last_haplotype_position is the smallest in all three. Make sure the last position in three files is the same. Limited by current version of hifi.
#seed_output_file_name = haplotype_file[:haplotype_file.find(".")]+"_refed.txt"	# for making hap_refed file
"""


import os, glob, subprocess, random, operator, time
from optparse import OptionParser

from tools import *

def output_dict(filename, title_info, dict):
	output_file = open(currentPath + filename, "w")
	print >> output_file, title_info
	dict_sorted_list = sort_dict_by_key(dict)
	for snp in dict_sorted_list:
		line = snp[1]
		print >> seed_new_file, line
	output_file.close()

def get_args():
	desc="Compare seed and std hap, to check purity of seed"

	usage = "seed_std_compare -i seed_file -c chr#" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-i", "--haplotype", type="string", dest="haplotypeFile",help = "Input File Name", default="null")
	parser.add_option("-s", "--genotype", type="string", dest="genotypeFile",help = "Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="chr11")
	(options, args) = parser.parse_args()
	if options.haplotypeFile == "null" or options.chrName == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

def compare_geno_ref(geno_dict, hap_ref_dict):

	geno_x_dict = {}
	geno_n_dict = {}
	geno_ref_not_consistent = {}
	ref_homo_dict = {}
	
	for position, snp in geno_dict.iteritems():
		if position in hap_ref_dict:
			exist_in_ref = False
			snp = snp.split()
			geno_A = snp[2][0]
			geno_B = snp[2][1]
			alleles = hap_ref_dict[position].split()
			alleles = alleles[2:]
			unique_alleles = list(set(alleles))
			n_alleles = len(unique_alleles)
			if n_alleles == 0 or n_alleles > 2:
				print "error in: ", position, n_alleles, unique_alleles, alleles
				sys.exit(1)
			else:
				try:
					if n_alleles == 1:
						ref_homo_dict[position] = list_to_line(snp)
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
							#print position, unique_alleles, n_alleles, geno_A, geno_B, geno_dict[position][2]
				except:
					#print position, unique_alleles, n_alleles
					pass
	print "geno_x_dict: ", len(geno_x_dict)
	print "geno_n_dict: ", len(geno_n_dict)
	print "geno_ref_not_consistent: ", len(geno_ref_not_consistent)
	print "ref_homo", len(ref_homo_dict)
	return (ref_homo_dict, geno_ref_not_consistent, geno_n_dict)

def load_hap_ref_data(chr_name):
	ref_title_info = ""
	hap_ref_dict = {}
	for infile in glob.glob(os.path.join(file_path,"*"+chr_name+"_???.phased")):	# add * for chrX
		ref_file_name = file_path + infile[(infile.find("hapmap3")):].strip()
		print ref_file_name
		if ref_title_info == "" and len(hap_ref_dict) == 0:
			data_tuple = load_raw_data(ref_file_name, raw_data_format_ref)
			ref_title_info = data_tuple[0]
			hap_ref_dict = data_tuple[1]
		else:
			temp_data_tuple = load_raw_data(ref_file_name, raw_data_format_ref)
			ref_title_info += "\t" + list_to_line((temp_data_tuple[0]).split()[2:])
			temp_hap_ref_dict = temp_data_tuple[1]
			for position, line in temp_hap_ref_dict.iteritems():
				if position in hap_ref_dict:
					hap_ref_dict[position] += "\t" + list_to_line(line.split()[2:])
	return (ref_title_info, hap_ref_dict)

def ref_preprocess(geno_dict, hap_ref_dict):
	ref_tulpe = compare_geno_ref(geno_dict, hap_ref_dict)
	ref_homo_dict = ref_tulpe[0]
	geno_ref_not_consistent = ref_tulpe[1]
	geno_n_dict = ref_tulpe[2]
	hap_ref_dict = dict_substract(hap_ref_dict, geno_ref_not_consistent)
	hap_ref_dict = dict_substract(hap_ref_dict, geno_n_dict)
	hap_ref_dict = dict_substract(hap_ref_dict, ref_homo_dict)
	return hap_ref_dict

def output_files(file_name, title_info, dict):
	outpuf_file = open(currentPath + file_name, "w")	# for hifi
	print >> outpuf_file, title_info
	sorted_list = sort_dict_by_key(dict)
	for element in sorted_list:
		print >> outpuf_file, element[1]
	outpuf_file.close()

def make_hifi_files():
	common_snp_number = 0
	hifi_seed_dict = {}
	hifi_geno_dict = {}
	hifi_ref_dict = {}
	
	last_seed_position = sort_dict_by_key(seed_dict)[-1][0]
	#print "last_haplotype_position", last_seed_position
	
	hap_ref_sorted_list = sort_dict_by_key(hap_ref_dict)
	for snp in hap_ref_sorted_list:
		position = snp[0]
		if position <= int(last_seed_position):  
			hifi_ref_dict[position] = hap_ref_dict[position]
			#outputFile.write(snp[1] + "\n")
			if position in seed_dict:
				hifi_seed_dict[position] = seed_dict[position]
				#seed_output_file.write(seed_dict[position] + "\n")
				
			if position in geno_dict:
				hifi_geno_dict[position] = geno_dict[position]
				#genotype_output_file.write(geno_dict[position] + "\n")		
	
	output_files("haplotype.txt", seed_title_info, hifi_seed_dict)
	output_files("genotype.txt", geno_title_info, hifi_geno_dict)
	output_files("refHaplos.txt", ref_title_info, hifi_ref_dict)

def refMerger(haplotype_file, chr_name):
	global seed_title_info
	global seed_dict
	global geno_title_info
	global geno_dict
	global hap_ref_dict
	global ref_title_info
	global raw_data_format_ref
	
	raw_data_format_ref = "string"
		
	ref_tulpe = load_hap_ref_data(chr_name)
	ref_title_info = ref_tulpe[0]
	hap_ref_dict = ref_tulpe[1]		
	
	#total_person_number = len(ref_title_info.strip().split())
	#print "total_person_number", total_element_number/2
	
	#genotype_input_file_name = "genotype_NA12878_chr6.txt"	# for simulation data
	genotype_file = file_path + "genotype_NA10847_"+chr_name+".txt"	# for all
	geno_tulpe = load_raw_data(genotype_file, raw_data_format_ref)
	geno_title_info = geno_tulpe[0]
	geno_dict = geno_tulpe[1]
	print "genotype_file", genotype_file 
	print "total_geno_number: ", len(geno_dict)
	
	seed_file = haplotype_file
	seed_tuple = load_raw_data(seed_file, raw_data_format_ref)
	seed_title_info = seed_tuple[0]
	seed_dict = seed_tuple[1]
	print "seed_file is :", seed_file
	
	print "hap_ref_dict before process", len(hap_ref_dict)
	hap_ref_dict = ref_preprocess(geno_dict, hap_ref_dict)
	print "hap_ref_dict after process", len(hap_ref_dict)
	
	make_hifi_files()


if __name__=='__main__':
	options = get_args()
	haplotype_file = options.haplotypeFile
	chr_name = options.chrName
	refMerger(haplotype_file, chr_name)
	
