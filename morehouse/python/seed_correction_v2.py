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
from seed_std_compare import seed_std_compare
from calculate_maf import calculate_maf

seed_dict = {}
geno_dict = {}
hifi_dict = {}


seed_title_info = ""
seed_file_name = ""
number_of_subfile = 10 

def get_args():
	desc="Compare seed and std hap, to check purity of seed"
	usage = "usage: %prog [options] arg1" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-i", "--seed", type="string", dest="seedFile",help = "Input seed File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode",help = "split or combine", default="null")
	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.seedFile == "null":
		print "param 'C': 2, 'T': 8, 'G'eters missing..."
		print usage
		sys.exit(1)
	return options

class seeds:
	def __init__(self):
		self.rsID = ""
		self.position = 0
		self.allele_ori = ""
		self.allele_new = ""
		self.allele_new_percentage = 0
		self.allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}

def load_seed_data(file_name):
	seed_dict = {}
	data_tuple = load_raw_data(file_name)
	seed_title_info = data_tuple[0]
	data_dict = data_tuple[1]
	for position, elements in data_dict.iteritems():
		seed = seeds()
		seed.rsID = elements[0].strip()
		seed.position = int(elements[1].strip())
		seed.allele_ori = elements[2].strip()
		seed_dict[position] = seed
	return (seed_title_info, seed_dict)

def load_hap_std(file_name):
	hap_std_dict = {}
	data_dict = load_raw_data(file_name)[1]
	for position, elements in data_dict.iteritems():
		if elements[2].strip() != "N" and elements[3].strip() != "N":
			try:
				position = elements[1].strip()
				position = int(position)
				hap_std_dict[position] = elements
			except ValueError:
				print "error in ", file_name, position
	return hap_std_dict	

def load_hifi_result(file_name, hifi_dict):
	data_dict = load_raw_data(file_name)[1]
	for position, elements in data_dict.iteritems():		
		try:
			if position not in hifi_dict:
				seed = seeds()
				seed.rsID = elements[0]
				seed.position = int(elements[1])
				hifi_dict[position] = seed			
			# need to modify this if seed is from second snp column, mother side
			hifi_dict[position].allele_new = elements[2].strip()
			# update the number of the allele from hifi result
			for base, value in hifi_dict[position].allele_dict.iteritems():
				if base == hifi_dict[position].allele_new:
					hifi_dict[position].allele_dict[base] += 1
		except:
				#print "error at ", file_name, position, elements
				pass
	print "total snp number in : ", file_name, len(data_dict) 
	return hifi_dict

# group the seed into homo and hetero groups
def group_seed(seed_dict, geno_dict):
	seed_homo_dict = {}
	seed_hetero_dict = {}
	for position, snp in seed_dict.iteritems():
		if position in geno_dict:	# pos in seed may not be in geno
			geno_allele = geno_dict[position][2]
			if geno_allele[0] == geno_allele[1]:
				seed_homo_dict[position] = snp
			else:
				seed_hetero_dict[position] = snp
		else:
			seed_hetero_dict[position] = snp
	return (seed_homo_dict, seed_hetero_dict)

pos_deleted = 545653
def seed_spliter(seed_dict, number_of_subfile):
	#del seed_hetero_dict[pos_deleted]
	print "seed_hetero_dict new", len(seed_hetero_dict)
	seed_hetero_sorted_list = [x for x in seed_hetero_dict.iteritems()] 
	seed_hetero_sorted_list.sort(key=lambda x: x[0]) # sort by key
	seed_number_ceilling = int(math.ceil(float(len(seed_hetero_sorted_list))/100)*100)
	#print "seed_hetero_number_ceilling: ", seed_number_ceilling
	seed_removed_in_each_subfile = seed_number_ceilling/number_of_subfile
	print "hetero_seed_removed_in_each_subfile: ", seed_removed_in_each_subfile
	#seed_removed_in_last_subfile = int(math.fmod(len(seed_hetero_sorted_list), seed_removed_in_each_subfile))
	#print "seed_removed_in_last_subfile: ", seed_removed_in_last_subfile
		
	seed_homo_sorted_list = [x for x in seed_homo_dict.iteritems()]
	"""
	seed_hetero_pos_list = []
	j = 0
	while j < len(seed_hetero_sorted_list):
		seed_hetero_pos_list.append(seed_hetero_sorted_list[j][0])
		j += 1

	print "seed_hetero_sorted_list", len(seed_hetero_sorted_list)
	pos_index = seed_hetero_pos_list.index(pos_deleted)
	print "pos_index", pos_deleted, pos_index, seed_hetero_pos_list[pos_index], seed_hetero_sorted_list[pos_index][0], seed_hetero_sorted_list[pos_index][1]
	start_pos = pos_index - 100
	if start_pos < 0:
		start_pos = 0
	end_pos = start_pos + 500
	print start_pos, end_pos
	del seed_hetero_sorted_list[pos_index]
	print "seed_hetero_sorted_list new", len(seed_hetero_sorted_list)
	"""
	print "seed_hetero_sorted_list", len(seed_hetero_sorted_list)
	file_number = 1
	while file_number <= number_of_subfile:	
		output_subfile_name = seed_file_name + "_" + str(file_number) + ".txt"
		#print "output: ", output_subfile_name
		output_subfile = open(currentPath + output_subfile_name, "w")
		print >> output_subfile, seed_title_info

		#sub_seed_list.append(seed_hetero_sorted_list[-1])	
		seed_hetero_sorted_list_bkup = copy.deepcopy(seed_hetero_sorted_list)
		
		i = 0 
		while i < int(seed_removed_in_each_subfile):
			try:
				del seed_hetero_sorted_list[i*file_number+file_number]
				del seed_hetero_sorted_list[i*(file_number+1)-file_number]
				#randomly del the second seed in whole range
				random_index = random.randrange(0,(len(seed_hetero_sorted_list)-1))
				while random_index == (i*file_number+file_number) or random_index == (i*file_number+file_number):
					random_index = random.randrange(0,(len(seed_hetero_sorted_list)-1))
				
				del seed_hetero_sorted_list[random_index]
				"""
				random_index = random.randrange(0,(len(seed_hetero_sorted_list)-1))
				while random_index == (i*file_number+file_number) or random_index == (i*file_number+file_number):
					random_index = random.randrange(0,(len(seed_hetero_sorted_list)-1))
				
				del seed_hetero_sorted_list[random_index]
				
				random_index = random.randrange(0,(len(seed_hetero_sorted_list)-1))
				while random_index == (i*file_number+file_number) or random_index == (i*file_number+file_number):
					random_index = random.randrange(0,(len(seed_hetero_sorted_list)-1))
				
				del seed_hetero_sorted_list[random_index]
				"""
			except:
				pass
			i += 1

		print "seed_hetero_sorted_list new", len(seed_hetero_sorted_list)
		sub_seed_list = seed_hetero_sorted_list + seed_homo_sorted_list	# add the homo snps to seed
		sub_seed_list.sort(key=lambda x: x[0]) # sort by key
		print len(sub_seed_list)
		for seed in sub_seed_list:
			line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
			print >> output_subfile, line	
		output_subfile.close()
		seed_hetero_sorted_list = copy.deepcopy(seed_hetero_sorted_list_bkup)
		
		if file_number < (number_of_subfile-1):
			hifi = program_path + "hifi_fu " + output_subfile_name + " &"
		else:
			hifi = program_path + "hifi_fu " + output_subfile_name + ""	# make sure the other hifi processes are finished
		hifi_process = subprocess.Popen(hifi, shell=True)
		hifi_process.wait()
		file_number += 1

def seed_extract(seed_dict):
	revised_seed_dict = {}
	hifi_dict = seed_dict.copy()
	print "hifi_dict initial", len(hifi_dict)
		
	file_number = 1
	while file_number <= number_of_subfile:	
		input_subfile_name = "imputed_haplotype_" + str(file_number) + ".txt"
		hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)
		file_number += 1
	seed_hetero_new_file = open(currentPath + "hetero.txt", "w")
	print >> seed_hetero_new_file, "rsID \t pos \t std_F \t std_M \t seed_ori \t seed_new \t seed_new_perc"

	hifi_sorted_list = sort_dict_by_key(hifi_dict) 
	for snp in hifi_sorted_list:
		position = snp[0]
		seed = snp[1]
		max_base = keywithmaxval(seed.allele_dict)
		max_value = seed.allele_dict[max_base]
		seed.allele_new_percentage = float(max_value)/float(number_of_subfile)
		if seed.allele_new_percentage*100 >= 100 and seed.allele_new == seed.allele_ori and position in seed_dict:
			seed.allele_new = max_base
			revised_seed_dict[position] = seed
		if position in seed_hetero_dict and seed.allele_ori != max_base and seed.allele_new_percentage*100 >= 80:
			hap_std = hap_std_dict[position]
			allele_dict = hifi_dict[position].allele_dict
			line = seed.rsID + "\t" + str(seed.position) + "\t" + hap_std[2] + "\t" + \
					hap_std[3] + "\t" + seed.allele_ori + "\t" + max_base + "\t" + str(seed.allele_new_percentage) + "\t" + str(allele_dict['A']+allele_dict['T']+allele_dict['C']+allele_dict['G'])
			print >> seed_hetero_new_file, line		
		if position == pos_deleted:
			print hifi_dict[position].allele_dict
			
	seed_hetero_new_file.close()	
	print "new seed total number", len(revised_seed_dict)
	output_revised_seed(revised_seed_dict)
	return revised_seed_dict

def output_revised_seed(revised_seed_dict):
	seed_new_file = open(currentPath + "haplotype_new.txt", "w")
	print >> seed_new_file, seed_title_info
	for position, seed in revised_seed_dict.iteritems():
		line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
		print >> seed_new_file, line
	seed_new_file.close()
	
def output_revised_seed_without_error(revised_seed_dict, same_to_B_dict):
	seed_new_file = open(currentPath + "haplotype_without_error.txt", "w")
	print >> seed_new_file, seed_title_info
	for position, seed in revised_seed_dict.iteritems():
		if position not in same_to_B_dict:
			line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
			print >> seed_new_file, line
	seed_new_file.close()

def output_revised_seed_with_error(revised_seed_dict, same_to_B_dict):
	seed_new_file = open(currentPath + "haplotype_with_error.txt", "w")
	print >> seed_new_file, seed_title_info
	for position, seed in revised_seed_dict.iteritems():
		if position in same_to_B_dict:
			line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_ori
			print >> seed_new_file, line
	seed_new_file.close()

def error_seed_distance(seed_dict, same_to_B_dict):
	seed_sorted_list = sort_dict_by_key(seed_dict) 
	seed_pos_list = []
	j = 0
	while j < len(seed_sorted_list):
		seed_pos_list.append(seed_sorted_list[j][0])
		j += 1
	for position, seed in same_to_B_dict.iteritems():
		index = seed_pos_list.index(position)
		print position
		print seed_sorted_list[index-1][0], position - seed_sorted_list[index-1][0]	
		print seed_sorted_list[index+1][0], seed_sorted_list[index+1][0]- position


def seed_correction(seed_file, chr_name, mode):

	global geno_dict
	global hap_std_dict
	global seed_dict
	global seed_homo_dict
	global seed_hetero_dict
	
	global seed_title_info
	global seed_file_name
	
	start = time.time()	
	options = get_args()
	seed_file = options.seedFile
	chr_name = options.chrName
	mode = options.mode	

	print "seed_file is :", seed_file
	seed_file_name = seed_file[:seed_file.find('.')].strip()
	seed_tuple = load_seed_data(seed_file)
	seed_title_info = seed_tuple[0]
	seed_dict = seed_tuple[1]

	#print "seed_title_info", seed_title_info
	print "total_seed_number: ", len(seed_dict)

	genotype_file = file_path + "genotype_NA10847_"+chr_name+".txt"	# for all
	geno_dict = load_raw_data(genotype_file)[1]

	hap_std_file = file_path + "ASW_"+chr_name+"_child_hap_refed.txt"	
	hap_std_dict = load_hap_std(hap_std_file)

	group_tuple = group_seed(seed_dict, geno_dict)
	seed_homo_dict = group_tuple[0]
	seed_hetero_dict = group_tuple[1]
	print "seed_homo_dict", len(seed_homo_dict)
	print "seed_hetero_dict", len(seed_hetero_dict)
			
	if mode == "split":
		seed_spliter(seed_dict, number_of_subfile)
		seed_extract(seed_dict)
		seed_std_compare("haplotype_new.txt", chr_name)
	elif mode == None:
		pass
	else:
		revised_seed_dict = seed_extract(seed_dict)
		same_to_B_dict = seed_std_compare("haplotype_new.txt", chr_name)
		print same_to_B_dict
		file_name = "refHaplos.txt"
		
		for position, seed in same_to_B_dict.iteritems():
			calculate_maf(file_name, position)	
		#error_seed_distance(seed_dict, same_to_B_dict)
		#error_seed_distance(seed_hetero_dict, same_to_B_dict)
		"""
		output_revised_seed_without_error(revised_seed_dict, same_to_B_dict)
		seed_std_compare("haplotype_without_error.txt", chr_name)
		output_revised_seed_with_error(revised_seed_dict, same_to_B_dict)
		"""
"""
	os.system("mkdir -p new")
	os.system("cp haplotype_new* imputed_haplotype_* new")
	os.system("mv haplotype_new.txt haplotype.txt")
"""

#def seed_extent():
	



if __name__=='__main__':
	options = get_args()
	seed_file = options.seedFile
	chr_name = options.chrName
	mode = options.mode	
	for i in range (5,9):
		seed_correction(seed_file, chr_name, mode)
		os.system("mkdir -p " + str(i))
		os.system("cp haplotype_new* imputed_haplotype_* " + str(i))
		os.system("mv haplotype_new.txt haplotype.txt")



