#!/usr/bin/python

#######################################################################################
# To correct and extend the seed
# based on genotype, divide the seed into hetero and homo groups. Keep the homo, divide the hetero into n parts
# randomly remove some hetero seed each time. (note, the hetero seeds here may also contain homo snps. The genotype is not available for them)
#######################################################################################

import os, glob, subprocess, random, operator, time, math, copy
from optparse import OptionParser

from tools import file_path, program_path, data_record_path, currentPath
from tools import sort_dict_by_key, load_raw_data, wccount, keywithmaxval

class seeds:
	def __init__(me):
		me.rsID = ""
		me.position = 0
		me.allele_ori = ""
		me.allele_new = ""
		me.allele_new_percentage = 0
		me.allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}

def load_seed_data(file_name):
	snp_dict = {}
	data_tulpe = load_raw_data(file_name)
	title_info = data_tulpe[0]
	data_dict = data_tulpe[1]
	for position, elements in data_dict.iteritems():
		seed = seeds()
		seed.rsID = elements[0].strip()
		seed.position = int(elements[1].strip())
		seed.allele_ori = elements[2].strip()
		snp_dict[position] = seed
	return (title_info, snp_dict)

def load_hifi_result(file_name, hifi_dict):
	data_dict = load_raw_data(file_name)[1]
	for position, elements in data_dict.iteritems():	
		try:
			if position not in hifi_dict:
				seed = seeds()
				seed.rsID = elements[0].strip()
				seed.position = int(elements[1].strip())
				hifi_dict[position] = seed							
			allele_new = elements[2].strip()	# need to modify this if seed is from second snp column, mother side
			hifi_dict[position].allele_new = allele_new
			# update the number of the allele from hifi result
			for base, value in hifi_dict[position].allele_dict.iteritems():
					if allele_new == base:
						hifi_dict[position].allele_dict[base] += 1
		except:
			#print "error at ", position
			pass
	print "total hifi snp number in : ", file_name, len(data_dict) 
	return hifi_dict

"""

def load_hifi_result(file_name, hifi_dict):
	#hifi_dict = {}
	#for position, seed in seed_dict.iteritems():
	#	hifi_dict[position] = seed
	total_hifi_snp_number = 0
	input_file = open(currentPath + file_name, "r")	
	for line in input_file:
		if line.startswith("rsID"):
			title_haplotype = line.strip()
		else:
			total_hifi_snp_number += 1
			elements = line.strip().split()
			rsID = elements[0].strip()
			try:
				position = int(elements[1].strip())
				if position not in hifi_dict:
					allele_ori = ""	
					allele_new = ""
					allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
					hifi_dict[position] = seed(rsID, position, allele_ori, allele_new, allele_dict)							
				allele_new = elements[2].strip()	# need to modify this if seed is from second snp column, mother side
				hifi_dict[position].allele_new = allele_new
				if allele_new == 'A':
					hifi_dict[position].allele_dict['A'] += 1
				elif allele_new == 'T':
					hifi_dict[position].allele_dict['T'] += 1
				elif allele_new == 'C':
					hifi_dict[position].allele_dict['C'] += 1
				elif allele_new == 'G':
					hifi_dict[position].allele_dict['G'] += 1
			except:
				#print line.strip()
				pass
	input_file.close()
	print "total hifi snp number in : ", file_name, total_hifi_snp_number 
	return hifi_dict
"""

# start time
start = time.time()		
		
# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--seed", type="string", dest="seedFile",help = "Input seed File Name", default="null")


parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
parser.add_option("-d", "--threshold", type="string", dest="threshold",help = "Input the depth threshold", default="3")
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")

(options, args) = parser.parse_args()

seed_file = options.seedFile
sam_file = options.samFile
depth_threshold = int(options.threshold)
chr_name = options.chrName

print "seed_file is :", seed_file

seed_file_name = seed_file[:seed_file.find('.')].strip()

seed_tulpe = load_seed_data(seed_file)
title_haplotype = seed_tulpe[0]
seed_dict = seed_tulpe[1]


print "title_haplotype", title_haplotype
print "total_seed_number: ", len(seed_dict)

# get genotype info
def load_genotype(genotype_file):
	geno_dict = {}
	inputFile_geno = open(file_path + genotype_file, "r")
	for line in inputFile_geno:
		if not line.startswith("rsID"):
			elements = line.strip().split()
			#rsID = elements[0].strip()
			position = int(elements[1].strip())
			alleles = elements[2].strip()
			geno_dict[position] = alleles
	inputFile_geno.close()
	return geno_dict

genotype_file = "genotype_NA10847_"+chr_name+".txt"	# for all
geno_dict = load_genotype(genotype_file)



def load_hap_std(inputFile_hap_std):
	hap_std_dict = {}
	inputFile_hap_std = open(file_path + inputFile_hap_std, "r")
	for line in inputFile_hap_std:
		if not line.startswith("rsID"):
			elements = line.strip().split()
			if elements[2].strip() != "N" and elements[3].strip() != "N":
				position = elements[1].strip()	
				try:
					position = int(position)
					hap_std_dict[position] = line.strip()
				except ValueError:
					print file_name, position	
					print line
	inputFile_hap_std.close()
	return hap_std_dict
	
hap_std_file = "ASW_"+chr_name+"_child_hap_refed.txt"	
hap_std_dict = load_hap_std(hap_std_file)

seed_homo_dict = {}
seed_hetero_dict = {}
for position, snp in seed_dict.iteritems():
	if position in geno_dict and geno_dict[position][0] == geno_dict[position][1]:
		seed_homo_dict[position] = snp
	else:
		seed_hetero_dict[position] = snp


def seed_spliter(seed_dict, number_of_subfile):
	global seed_homo_dict
	global seed_hetero_dict
	seed_homo_dict = {}
	seed_hetero_dict = {}
	for position, snp in seed_dict.iteritems():
		if position in geno_dict and geno_dict[position][0] == geno_dict[position][1]:
			seed_homo_dict[position] = snp
		else:
			seed_hetero_dict[position] = snp 
	seed_hetero_sorted_list = [x for x in seed_hetero_dict.iteritems()] 
	seed_hetero_sorted_list.sort(key=lambda x: x[0]) # sort by key
	seed_number_ceilling = int(math.ceil(float(len(seed_hetero_sorted_list))/100)*100)
	print "seed_hetero_number_ceilling: ", seed_number_ceilling
	seed_removed_in_each_subfile = seed_number_ceilling/number_of_subfile
	print "seed_hetero_in_each_subfile: ", seed_removed_in_each_subfile
	seed_removed_in_last_subfile = int(math.fmod(len(seed_hetero_sorted_list), seed_removed_in_each_subfile))
	print "seed_removed_in_last_subfile: ", seed_removed_in_last_subfile
		
	seed_homo_sorted_list = [x for x in seed_homo_dict.iteritems()] 
	file_number = 1
	while file_number <= number_of_subfile:

		if True:	
			output_subfile_name = seed_file_name + "_" + str(file_number) + ".txt"
			print "output: ", output_subfile_name
			output_subfile = open(currentPath + output_subfile_name, "w")
			print >> output_subfile, title_haplotype

			#sub_seed_list.append(seed_hetero_sorted_list[-1])	
			seed_hetero_sorted_list_bkup = copy.deepcopy(seed_hetero_sorted_list)
			# random remove seeds
			i = 0
			#while i <  seed_removed_in_each_subfile:
			while i <  300:
				seed_to_be_removed = random.randrange(0,(len(seed_hetero_sorted_list)-1))
				#seed_hetero_sorted_list.remove(seed_hetero_sorted_list[seed_to_be_removed])
				del seed_hetero_sorted_list[seed_to_be_removed]
				i += 1
			"""
			i = 0 
			while i < seed_removed_in_each_subfile:
				try:
					del seed_hetero_sorted_list[i*file_number+file_number]
				except:
					pass
				i += 1
			"""	
			sub_seed_list = seed_hetero_sorted_list + seed_homo_sorted_list	# add the homo snps to seed
			sub_seed_list.sort(key=lambda x: x[0]) # sort by key
			print len(sub_seed_list)
			for seed in sub_seed_list:
				line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
				print >> output_subfile, line	
			output_subfile.close()
			seed_hetero_sorted_list = copy.deepcopy(seed_hetero_sorted_list_bkup)
			
			hifi = program_path + "hifi_fu " + output_subfile_name + " &"
			hifi_process = subprocess.Popen(hifi, shell=True)
			hifi_process.wait()
			
			"""
			cp = "cp imputedhaplotype.txt imputedhaplotype_" + str(file_number)  + ".txt"
			cp_process = subprocess.Popen(cp, shell=True)
			cp_process.wait()
			#i += 1
			"""
		file_number += 1
			

number_of_subfile = 10 
#seed_spliter(seed_dict, number_of_subfile)

hifi_dict = {}
for position, seeds in seed_dict.iteritems():
	hifi_dict[position] = seeds
print len(hifi_dict)
	
file_number = 1
while file_number <= number_of_subfile:	
	input_subfile_name = "imputed_haplotype_" + str(file_number) + ".txt"
	hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)
	print file_number," :", len(hifi_dict)
	file_number += 1

seed_new_file = open(currentPath + "haplotype_new.txt", "w")
print >> seed_new_file, title_haplotype
num_100 = 0

seed_hetero_new_file = open(currentPath + "hetero.txt", "w")
#print >> seed_hetero_new_file, title_haplotype
print >> seed_hetero_new_file, "rsID \t pos \t std_F \t std_M \t seed_ori \t seed_new \t seed_new_perc"

hifi_sorted_list = [x for x in hifi_dict.iteritems()] 
hifi_sorted_list.sort(key=lambda x: x[0]) # sort by key
for snp in hifi_sorted_list:
	position = snp[0]
	seed = snp[1]
	max_base = keywithmaxval(seed.allele_dict)
	max_value = seed.allele_dict[max_base]
	seed.allele_new_percentage = float(max_value)/float(number_of_subfile)
	# snp in geno >=90 not in geno >= 100
	if seed.allele_new_percentage*100 >= 100 and position in seed_dict:
		num_100 += 1
		seed.allele_new = max_base
		line = seed.rsID + "\t" + str(seed.position) + "\t" + max_base
		print >> seed_new_file, line
	if position in seed_hetero_dict and seed.allele_ori != max_base and seed.allele_new_percentage*100 >= 80:
		hap_std = hap_std_dict[position].split()
		allele_dict = hifi_dict[position].allele_dict
		line = seed.rsID + "\t" + str(seed.position) + "\t" + hap_std[2] + "\t" + \
				hap_std[3] + "\t" + seed.allele_ori + "\t" + max_base + "\t" + str(seed.allele_new_percentage) + "\t" + str(allele_dict['A']+allele_dict['T']+allele_dict['C']+allele_dict['G'])
		print >> seed_hetero_new_file, line
		
seed_new_file.close()
seed_hetero_new_file.close()
	
print "num_100", num_100




