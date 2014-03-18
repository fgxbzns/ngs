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
from seed_std_compare import seed_std_compare
from hifiAccuCheck_v2 import hifiAccuCheck

def output_dict(filename, title_info, dict):
	output_file = open(currentPath + filename, "w")
	print >> output_file, title_info
	dict_sorted_list = sort_dict_by_key(dict)
	for snp in dict_sorted_list:
		line = snp[1]
		print >> output_file, line
	output_file.close()

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
				#print "error in: ", position, n_alleles, unique_alleles, alleles
				#sys.exit(1)
				pass
			else:
				try:
					#if n_alleles == 1:
					#	ref_homo_dict[position] = list_to_line(snp)
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
					print position, unique_alleles, n_alleles
					pass
	
	print "geno_x_dict: ", len(geno_x_dict)
	print "geno_n_dict: ", len(geno_n_dict)
	print "geno_ref_not_consistent: ", len(geno_ref_not_consistent)
	print "ref_homo", len(ref_homo_dict)
	#print ref_homo_dict[104922938]
	return (ref_homo_dict, geno_ref_not_consistent, geno_n_dict)
	

def load_hap_ref_data(chr_name):
	ref_title_info = ""
	hap_ref_dict = {}
	for infile in glob.glob(os.path.join(file_path,"*"+chr_name+"_???.phased")):	# add * for chrX
	#for infile in glob.glob(os.path.join(file_path,"*"+chr_name+"_???.phased_12878FM_remd")):	# simulation data NA12878 FM removed add * for chrX
		ref_file_name = file_path + infile[(infile.find("hapmap3")):].strip()
		print ref_file_name
		if ref_title_info == "" and len(hap_ref_dict) == 0:
			ref_title_info, hap_ref_dict = load_raw_data(ref_file_name, raw_data_format_ref)
		else:
			temp_data_tuple = load_raw_data(ref_file_name, raw_data_format_ref)
			ref_title_info += "\t" + list_to_line((temp_data_tuple[0]).split()[2:])
			temp_hap_ref_dict = temp_data_tuple[1]
			for position, line in temp_hap_ref_dict.iteritems():
				if position in hap_ref_dict:
					#hap_ref_dict[position] += "\t" + list_to_line(line[2:])
					hap_ref_dict[position].extend(line[2:])
	return (ref_title_info, hap_ref_dict)

def ref_preprocess(geno_dict, hap_ref_dict):
	ref_homo_dict, geno_ref_not_consistent, geno_n_dict = compare_geno_ref(geno_dict, hap_ref_dict)
	hap_ref_dict = dict_substract(hap_ref_dict, geno_ref_not_consistent)
	hap_ref_dict = dict_substract(hap_ref_dict, geno_n_dict)
	#print len(hap_ref_dict)
	#print len(ref_homo_dict)
	hap_ref_dict = dict_substract(hap_ref_dict, ref_homo_dict)
	#print len(hap_ref_dict)
	return hap_ref_dict

def output_files(file_name, title_info, dict):
	outpuf_file = open(currentPath + file_name, "w")	# for hifi
	print >> outpuf_file, title_info
	sorted_list = sort_dict_by_key(dict)
	for element in sorted_list:
		print >> outpuf_file, element[1]
	outpuf_file.close()

def make_hifi_files(remPercent):
	common_snp_number = 0
	hifi_seed_dict = {}
	hifi_geno_dict = {}
	hifi_ref_dict = {}
	
	
	last_seed_position = sort_dict_by_key(seed_dict)[-1][0]
	#print "last_haplotype_position", last_seed_position
	hap_ref_sorted_list = sort_dict_by_key(hap_ref_dict)
	
	if remPercent > 0:
		seed_homo_dict, seed_hetero_dict = group_seed(seed_dict, geno_dict)
		"""
		print "seed_hetero_dict", len(seed_hetero_dict)
		print len(hap_ref_sorted_list)
		print len(hap_ref_dict)
		"""
		hap_ref_sorted_list = [x for x in hap_ref_sorted_list if x[0] not in seed_hetero_dict]

		#print len(hap_ref_sorted_list)
		
		hap_ref_size = len(hap_ref_dict)
		print "remPercent is", remPercent
		for i in range(int(remPercent*hap_ref_size)):
			if len(hap_ref_sorted_list) > 1:	
				random_index = random.randrange(0,(len(hap_ref_sorted_list)-1))
				position = hap_ref_sorted_list[random_index][0]
				if position in hap_ref_dict:
					del hap_ref_dict[int(position)]
					del hap_ref_sorted_list[random_index]			
		
		hap_ref_sorted_list = sort_dict_by_key(hap_ref_dict)

	for snp in hap_ref_sorted_list:
		position = snp[0]
		if position <= int(last_seed_position):  
			#hifi_ref_dict[position] = " ".join(hap_ref_dict[position])
			hifi_ref_dict[position] = list_to_line(hap_ref_dict[position])
			if position in seed_dict:
				#hifi_seed_dict[position] = " ".join(seed_dict[position])
				hifi_seed_dict[position] = list_to_line(seed_dict[position])
			elif position in geno_dict and geno_dict[position][2][0] == geno_dict[position][2][1]:
				hifi_seed_dict[position] = geno_dict[position][0] + " " + geno_dict[position][1] + " " + geno_dict[position][2][0]	
			#else:
			#	seed_not_in_ref[position] = list_to_line(seed_dict[position])
			if position in geno_dict:
				#hifi_geno_dict[position] = " ".join(geno_dict[position])
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
		
	ref_title_info, hap_ref_dict = load_hap_ref_data(chr_name)	
	
	#total_person_number = len(ref_title_info.strip().split())
	#print "total_person_number", total_element_number/2
	
	#genotype_file = file_path + "genotype_NA12878_chr6.txt"	# for simulation data
	genotype_file = file_path + "genotype_NA10847_"+chr_name+".txt"	# for all
	geno_title_info, geno_dict = load_raw_data(genotype_file, raw_data_format_ref)
	print "genotype_file", genotype_file 
	print "total_geno_number: ", len(geno_dict)
	
	seed_file = haplotype_file
	seed_title_info, seed_dict = load_raw_data(seed_file, raw_data_format_ref)
	print "seed_file is :", seed_file
	
	print "hap_ref_dict before process", len(hap_ref_dict)
	hap_ref_dict = ref_preprocess(geno_dict, hap_ref_dict)
	print "hap_ref_dict after process", len(hap_ref_dict)
	
	make_hifi_files(remPercent)

def extract_genotype(genotype_file_name, chr_name, person_ID):
	with open("genotype_"+person_ID+"_"+chr_name+".txt", "w") as genotype_output_file:
		print >> genotype_output_file, "rsID" + "\t" + "phys_position" + "\t" + person_ID
		with open(genotype_file_name, "r") as genotype_input_file:
			for line in genotype_input_file:
				elements = line.strip().split()
				if line.startswith("rs#"):
					try:
						ID_index = elements.index(person_ID)
						print "ID_index", ID_index
					except:
						print person_ID, "not found in file"
						pass
				else:
					try:
						rsID = elements[0].strip()							
						position = elements[3].strip()
						genotype = elements[ID_index].strip()
						print >> genotype_output_file, rsID + "\t" + position + "\t" + genotype
					except ValueError:
						print "error in ", line

def get_args():
	desc="Compare seed and std hap, to check purity of seed"

	usage = "seed_std_compare -i seed_file -c chr#" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-i", "--haplotype", type="string", dest="haplotypeFile",help = "Input File Name", default="null")
	parser.add_option("-g", "--genotype", type="string", dest="genotypeFile",help = "Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
	parser.add_option("-p", "--pct", type="float", dest="remPercent",help = "percent of ref removed", default=0)
	parser.add_option("-d", "--id", type="string", dest="personID",help = "personID", default=0)
	
	(options, args) = parser.parse_args()
	if options.haplotypeFile == "null" or options.chrName == "null":
		print "parameters missing..."
		print usage
		#sys.exit(1)
	return options

if __name__=='__main__':
	options = get_args()
	
	haplotype_file = options.haplotypeFile
	chr_name = options.chrName
	remPercent = options.remPercent
	start_time = time.time()
	
	refMerger(haplotype_file, chr_name, remPercent)
	
	elapsed_time = time.time() - start_time
	print "Elapsed time is: " + str(format(elapsed_time, "0.3f")) + "s"
	
	"""
	# for extracting genotype
	genotype_file_name = options.genotypeFile
	chr_name = options.chrName
	person_ID = options.personID
	extract_genotype(genotype_file_name, chr_name, person_ID)
	"""
	#file_name = "haplotype.txt"
	#seed_std_compare(file_name, chr_name)
	#hifi_run(file_name, chr_name)
	#hifiAccuCheck("imputed_"+file_name, chr_name)
	
