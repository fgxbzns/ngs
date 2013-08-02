#!/usr/bin/python

#######################################################################################
# To correct and extend the seed
# based on genotype, divide the seed into hetero and homo groups. Keep the homo, divide the hetero into n parts
# randomly remove some hetero seed each time. (note, the hetero seeds here may also contain homo snps. The genotype is not available for them)
#######################################################################################

import os, glob, subprocess, random, operator, time, math, copy
from optparse import OptionParser

from tools import *
from seed_std_compare import seed_std_compare
from calculate_maf import calculate_maf
from hifiAccuCheck_v2 import hifiAccuCheck

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
	data_tuple = load_raw_data(file_name, raw_data_format)
	seed_title_info = data_tuple[0]
	data_dict = data_tuple[1]
	for position, elements in data_dict.iteritems():
		seed = seeds()
		seed.rsID = elements[0].strip()
		seed.position = int(elements[1].strip())
		seed.allele_ori = elements[2].strip()
		seed.allele_new = elements[2].strip()
		seed_dict[position] = seed
	return (seed_title_info, seed_dict)

def load_seed_data_from_dict(data_dict):
	seed_dict = {}
	for position, elements in data_dict.iteritems():
		seed = seeds()
		seed.rsID = elements[0].strip()
		seed.position = int(elements[1].strip())
		seed.allele_ori = elements[2].strip()
		seed.allele_new = elements[2].strip()
		seed_dict[position] = seed
	return seed_dict


def load_hap_std(file_name):
	hap_std_dict = {}
	data_dict = load_raw_data(file_name, raw_data_format)[1]
	for position, elements in data_dict.iteritems():
		""" ??? N X """
		if elements[2].strip() != "N" and elements[3].strip() != "N": 
			try:
				position = elements[1].strip()
				position = int(position)
				hap_std_dict[position] = elements
			except ValueError:
				print "error in ", file_name, position
	return hap_std_dict	

def load_hifi_result(file_name, hifi_dict):
	data_dict = load_raw_data(file_name, raw_data_format)[1]
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
def seed_error_remove(seed_dict):
	#del seed_hetero_dict[pos_deleted]
	global seed_hetero_dict
	global seed_homo_dict
	print "seed_hetero_dict new", len(seed_hetero_dict)
	
	seed_hetero_sorted_list = sort_dict_by_key(seed_hetero_dict) 
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
	file_number = 0
	while file_number < number_of_subfile:	
		output_subfile_name = seed_file_name + "_" + str(file_number) + ".txt"
		#print "output: ", output_subfile_name
		output_subfile = open(currentPath + output_subfile_name, "w")
		print >> output_subfile, seed_title_info

		#sub_seed_list.append(seed_hetero_sorted_list[-1])	
		#seed_hetero_sorted_list_bkup = copy.deepcopy(seed_hetero_sorted_list)
		seed_hetero_dict_bkup = seed_hetero_dict.copy()
		
		i = 0 
		while i < int(seed_removed_in_each_subfile):
			try:
				forward_index = i*number_of_subfile+file_number
				backward_index = (i+1)*number_of_subfile-file_number
				position = seed_hetero_sorted_list[forward_index][0]
				del seed_hetero_dict[position]
				#del seed_hetero_sorted_list[forward_index]
				
				position = seed_hetero_sorted_list[backward_index][0]
				del seed_hetero_dict[position]
				#del seed_hetero_sorted_list[backward_index]
				
				#randomly del the second seed in whole range
				random_index = random.randrange(0,(len(seed_hetero_sorted_list)-1))
				while random_index == (forward_index) or random_index == (backward_index):
					random_index = random.randrange(0,(len(seed_hetero_sorted_list)-1))
				
				position = seed_hetero_sorted_list[random_index][0]
				del seed_hetero_dict[position]
				""" error, size of the sorted list will change each time after one element is deleted"""

			except:
				pass
			i += 1

		print "seed_hetero_sorted_list new", len(seed_hetero_sorted_list)
		
		#sub_seed_list = seed_hetero_sorted_list + seed_homo_sorted_list	# add the homo snps to seed
		#sub_seed_list.sort(key=lambda x: x[0]) # sort by key
		sub_seed_dict = dict_add(seed_hetero_dict, seed_homo_dict)				

		sub_seed_list = sort_dict_by_key(sub_seed_dict)
		
		print len(sub_seed_list)
		for seed in sub_seed_list:
			line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
			print >> output_subfile, line	
		output_subfile.close()
		#seed_hetero_sorted_list = copy.deepcopy(seed_hetero_sorted_list_bkup)
		seed_hetero_dict = seed_hetero_dict_bkup.copy()
		
		if file_number < (number_of_subfile-2):
			hifi = program_path + "hifi_fu " + output_subfile_name + " &"
		else:
			hifi = program_path + "hifi_fu " + output_subfile_name + ""	# make sure the other hifi processes are finished
		hifi_process = subprocess.Popen(hifi, shell=True)
		hifi_process.wait()
		file_number += 1

def seed_error_remove_extract(seed_dict):
	revised_seed_dict = {}
	hifi_dict = seed_dict.copy()
	print "hifi_dict initial", len(hifi_dict)
		
	file_number = 0
	while file_number < number_of_subfile:	
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
		""" ??? """
		if seed.allele_new_percentage*100 >= 100 and max_base == seed.allele_ori and position in seed_dict:
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
	output_revised_seed("haplotype_new.txt", revised_seed_dict)
	return revised_seed_dict

def seed_recover(seed_dict, revised_seed_dict):
	removed_seed_dict = dict_substract(seed_hetero_dict, revised_seed_dict)
	
	print "removed_seed_dict", len(removed_seed_dict)
	seed_number_ceilling = int(math.ceil(float(len(removed_seed_dict))/100)*100)
	seed_added_in_each_subfile = seed_number_ceilling/number_of_subfile
	print "hetero_seed_added_in_each_subfile: ", seed_added_in_each_subfile
		
	removed_seed_list = sort_dict_by_key(removed_seed_dict)
	
	seed_homo_sorted_list = [x for x in seed_homo_dict.iteritems()]
	
	print "revised_seed_dict", len(revised_seed_dict)
	file_number = 0
	while file_number < number_of_subfile:	
		output_subfile_name = seed_file_name + "_" + str(file_number) + ".txt"
		#print "output: ", output_subfile_name
		output_subfile = open(currentPath + output_subfile_name, "w")
		print >> output_subfile, seed_title_info
		revised_seed_dict_bkup = copy.deepcopy(revised_seed_dict)	
		i = 0 
		while i < int(seed_added_in_each_subfile):	# remove 5% seed
			try:
				index = i*number_of_subfile+file_number
				position = removed_seed_list[index][0]
				if position not in revised_seed_dict:
					revised_seed_dict[position] = removed_seed_list[index][1]
				else:
					print "seed already in new seed dict"
				
				random_pos = random.randrange(0,number_of_subfile)
				while random_pos == position:
					random_pos = random.randrange(0,number_of_subfile)
				if random_pos not in revised_seed_dict:
					revised_seed_dict[random_pos] = removed_seed_dict[random_pos]
				else:
					print "seed already in new seed dict"			
			except:
				pass
			i += 1

		print "revised_seed_dict new", len(revised_seed_dict)
	
		revised_seed_list = sort_dict_by_key(revised_seed_dict)
		
		for seed in revised_seed_list:
			line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
			print >> output_subfile, line	
		output_subfile.close()
		revised_seed_dict = copy.deepcopy(revised_seed_dict_bkup)
		
		if file_number < (number_of_subfile-1):
			hifi = program_path + "hifi_fu " + output_subfile_name + " &"
		else:
			hifi = program_path + "hifi_fu " + output_subfile_name + ""	# make sure the other hifi processes are finished
		hifi_process = subprocess.Popen(hifi, shell=True)
		hifi_process.wait()
		file_number += 1

def seed_recover_extract(seed_hetero_dict, revised_seed_dict):
	#revised_seed_dict = {}
	recovered_seed_dict = {}
	removed_seed_dict = dict_substract(seed_hetero_dict, revised_seed_dict)
	hifi_dict = {}
	#hifi_dict = seed_dict.copy()
	#print "hifi_dict initial", len(hifi_dict)
	
	file_number = 0
	while file_number < number_of_subfile:	
		input_subfile_name = "imputed_haplotype_ori_" + str(file_number) + ".txt"
		hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)
		file_number += 1
	print "hifi_dict size: ", len(hifi_dict)
	for position, snp in removed_seed_dict.iteritems():
		if position in hifi_dict:
			seed = hifi_dict[position]
			max_base = keywithmaxval(seed.allele_dict)	
			max_value = seed.allele_dict[max_base]
			seed.allele_new_percentage = float(max_value)/float(number_of_subfile)
			if seed.allele_new_percentage*100 >= 90:	# ref
			#if seed.allele_new_percentage*100 >= 80 and max_base == seed_hetero_dict[position].allele_ori:	# ori
				seed.allele_new = max_base
				recovered_seed_dict[position] = seed
	
	revised_seed_dict = dict_add(revised_seed_dict, recovered_seed_dict)
		
	"""
	
	seed_sorted_list = sort_dict_by_key(seed_dict)
	hifi_sorted_list = sort_dict_by_key(hifi_dict)
	
	hifi_sorted_pos_list = []
	j = 0
	while j < len(hifi_sorted_list):
		hifi_sorted_pos_list.append(hifi_sorted_list[j][0])
		j += 1
	print
	print "hifi_sorted_pos_list", len(hifi_sorted_pos_list)
	for snp in seed_sorted_list:
		position = snp[0]
		#print position
		pos_index = hifi_sorted_pos_list.index(position)
		for new_pos in [(pos_index-1),(pos_index+1)]:
			if new_pos < len(hifi_sorted_list):
				seed = hifi_sorted_list[new_pos][1]
				max_base = keywithmaxval(seed.allele_dict)
				max_value = seed.allele_dict[max_base]
				seed.allele_new_percentage = float(max_value)/float(number_of_subfile)	geno_homo_dict = group_tuple[0]
	geno_hetero_dict
				if seed.allele_new_percentage*100 == 100 and position in seed_dict:
					seed.allele_new = max_base
					revised_seed_dict[position] = seed
	
	
		
	seed_hetero_new_file = open(currentPath + "hetero.txt", "w")
	print >> seed_hetero_new_file, "rsID \t pos \t std_F \t std_M \t seed_ori \t seed_new \t seed_new_perc"

	hifi_sorted_list = sort_dict_by_key(hifi_dict) 
	for snp in hifi_sorted_list:
		position = snp[0]
		seed = snp[1]
		max_base = keywithmaxval(seed.allele_dict)
		max_value = seed.allele_dict[max_base]
		seed.allele_new_percentage = float(max_value)/float(number_of_subfile)
		if seed.allele_new_percentage*100 == 100 and position in seed_dict:
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
	"""
	print "new seed total number", len(recovered_seed_dict)
	output_revised_seed("haplotype_recoved.txt", recovered_seed_dict)
	output_revised_seed("haplotype_new.txt", revised_seed_dict)
	return recovered_seed_dict



def seed_expand_qs(seed_file):
	expanded_seed_dict = {}
	hifi_dict = {}
	revised_seed_dict = load_seed_data(seed_file)[1]
	"""
	hifi = program_path + "hifi_fu_qs " + seed_file + ""
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()
	"""
	hifi_file = "imputed_" + seed_file
	hifi_dict = load_hifi_result(hifi_file, hifi_dict)
	hifi_dict = dict_substract(hifi_dict, revised_seed_dict)
	print "hifi_dict: ", len(hifi_dict)
	
	qscore_file = "qscore_" + seed_file
	qscore_dict = load_raw_data(qscore_file, raw_data_format)[1]
	qscore_dict = dict_substract(qscore_dict, revised_seed_dict)
	print "qscore_dict: ", len(qscore_dict)
	a = 0
	for position, seed in hifi_dict.iteritems():
		if float(qscore_dict[position][3]) == 1.0 and position in geno_dict and seed.allele_new in geno_dict[position][2]:
			expanded_seed_dict[position] = seed
		
		elif float(qscore_dict[position][3]) >= 0.6 and position in geno_dict and seed.allele_new in geno_dict[position][2]:
				a += 1
				expanded_seed_dict[position] = seed
		
	print "a", a
	revised_seed_dict = dict_add(revised_seed_dict, expanded_seed_dict)
	
	print "new seed total number", len(expanded_seed_dict)
	output_revised_seed("haplotype_expanded.txt", expanded_seed_dict)
	output_revised_seed("haplotype_t.txt", revised_seed_dict)
	return revised_seed_dict

def seed_expand_ref():

	seed_ref_difference_dict = dict_substract(hap_ref_dict, seed_dict)
	
	print "hap_ref_dict", len(hap_ref_dict)
	print "seed_ref_difference", len(seed_ref_difference_dict)
	ref_number_ceilling = int(math.ceil(float(len(seed_ref_difference_dict))/100)*100)
	ref_removed_in_each_subfile = ref_number_ceilling/number_of_subfile
	print "ref_removed_in_each_subfile: ", ref_removed_in_each_subfile
		
	seed_ref_difference_list = sort_dict_by_key(seed_ref_difference_dict)
	
	seed_ref_same_dict = dict_substract(hap_ref_dict, seed_ref_difference_dict)
	#seed_ref_same_list = sort_dict_by_key(seed_ref_same_dict)
	global geno_dict
	print "geno_dict global", len(geno_dict)	
	
	#print seed_ref_difference_list[0]
	
	
	file_number = 0
	while file_number < number_of_subfile:	
		
		pos_del_record = open("pos_deledted_"+str(file_number), "w")
		
		ref_subfile_name = "refHaplos" + "_" + str(file_number) + ".txt"
		ref_subfile = open(currentPath + ref_subfile_name, "w")
		print >> ref_subfile, ref_title_info
		#seed_ref_difference_list_bkup = copy.deepcopy(seed_ref_difference_list)	
		seed_ref_difference_dict_bkup = seed_ref_difference_dict.copy()
		geno_dict_bkup = geno_dict.copy()
		i = 0 
		while i < int(ref_removed_in_each_subfile):			
			try:
						
				forward_index = i*number_of_subfile+file_number
				position = seed_ref_difference_list[forward_index][0]
				if position in geno_dict:
					del geno_dict[position]
				
				print >> pos_del_record, seed_ref_difference_list[forward_index][0]
						
				del seed_ref_difference_dict[position]
					
				random_index = random.randrange(0,(len(seed_ref_difference_dict)-1))
				while random_index == (forward_index):
					random_index = random.randrange(0,(len(seed_ref_difference_dict)-1))
				
				position = seed_ref_difference_list[random_index][0]
				if position in geno_dict:
					del geno_dict[position]
				print >> pos_del_record, seed_ref_difference_list[random_index][0]
					
				del seed_ref_difference_dict[position]
				
			except:
				#print "error", forward_index
				pass
			i += 1

		print "seed_ref_difference_dict new", len(seed_ref_difference_dict)
		pos_del_record.close()
		
		sub_ref_dict = dict_add(seed_ref_difference_dict, seed_ref_same_dict)	
		sub_ref_list = sort_dict_by_key(sub_ref_dict)	
		
		
		for seed in sub_ref_list:
			print >> ref_subfile, list_to_line(seed[1])	
		ref_subfile.close()
		#seed_ref_difference_list = copy.deepcopy(seed_ref_difference_list_bkup)
		seed_ref_difference_dict = seed_ref_difference_dict_bkup.copy()
		print "seed_ref_difference_dict ori", len(seed_ref_difference_dict)
	
		# generate new genotype file
		geno_subfile_name = "genotype" + "_" + str(file_number) + ".txt"
		geno_subfile = open(currentPath + geno_subfile_name, "w")
		print >> geno_subfile, geno_title_info
		print "geno_dict reduced", len(geno_dict)
		geno_sorted_list = sort_dict_by_key(geno_dict)
		for geno in geno_sorted_list:
			print >> geno_subfile, list_to_line(geno[1])	
		geno_subfile.close()
		geno_dict = geno_dict_bkup.copy()
		print "geno_dict", len(geno_dict)
		
		hap_subfile_name = "haplotype" + "_" + str(file_number) + ".txt"
		
		os.system("cp haplotype.txt " + hap_subfile_name)
		
		if file_number < (number_of_subfile-1):
			hifi = program_path + "hifi_fu_ref " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " &"
		else:
			hifi = program_path + "hifi_fu_ref " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + ""	# make sure the other hifi processes are finished
		hifi_process = subprocess.Popen(hifi, shell=True)
		hifi_process.wait()
		
		
		file_number += 1

def load_qscore_result(file_name, qscore_dict):
	data_dict = load_raw_data(file_name, raw_data_format)[1]
	for position, elements in data_dict.iteritems():		
		try:
			if position not in qscore_dict:
				qscore_dict[position] = float(elements[3])
			else:
				qscore_dict[position] = qscore_dict[position] + float(elements[3].strip())
		except:
				#print "error at ", file_name, position, elements
				pass
	return qscore_dict

	
def seed_recover_extract_ref():
	recovered_seed_dict = {}
	hifi_dict = {}
	qscore_dict = {}
	
	seed_ref_difference_dict = dict_substract(hap_ref_dict, seed_dict)
	
	file_number = 0
	while file_number < number_of_subfile:	
		input_subfile_name = "imputed_haplotype_" + str(file_number) + ".txt"
		hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)
		file_number += 1
	print "hifi_dict size: ", len(hifi_dict)
	
	file_number = 0
	while file_number < number_of_subfile:	
		input_subfile_name = "qscore_haplotype_" + str(file_number) + ".txt"
		qscore_dict = load_qscore_result(input_subfile_name, qscore_dict)
		file_number += 1
	print "qscore_dict size: ", len(qscore_dict)
	
	
	"""
	#a_file = open("expanded_seed.txt", "w")
	for position, snp in seed_ref_difference_dict.iteritems():
		if position in hifi_dict:
			seed = hifi_dict[position]
			max_base = keywithmaxval(seed.allele_dict)
			max_value = seed.allele_dict[max_base]
			seed.allele_new_percentage = float(max_value)/float(number_of_subfile)
			#each difference position should only have consistent value in 9 of 10 files. it is deleted at least twice.
			if seed.allele_new_percentage*100 >= 0: #and position in geno_dict and max_base in geno_dict[position][2]:	# ref
				if True:
				#if position in qscore_dict and qscore_dict[position]/float(number_of_subfile-1) >= 0.80:
					seed.allele_new = max_base
					recovered_seed_dict[position] = seed
	#			print >> a_file, position, max_base, max_value, seed.allele_new_percentage, seed.allele_dict
	#a_file.close()
	"""
	hifi_sorted_list = sort_dict_by_key(hifi_dict)
	hifi_sorted_pos_list = []
	j = 0
	while j < len(hifi_sorted_list):
		hifi_sorted_pos_list.append(hifi_sorted_list[j][0])
		j += 1
	print "hifi_sorted_pos_list", len(hifi_sorted_pos_list)
	
	position_distance = 5000
	expand_range = 2
	#while position_distance <= 5000:
	#	recovered_seed_dict = {}
	print "current position_distance", position_distance
	for position, snp in seed_hetero_dict.iteritems():
		pos_index = hifi_sorted_pos_list.index(position)
		start_position = (pos_index-expand_range) if pos_index-expand_range >=0 else 0
		end_position = (pos_index+expand_range) if (pos_index+expand_range) < len(hifi_sorted_list) else len(hifi_sorted_list)-1
		for new_pos in range (start_position, end_position):
			#if hifi_sorted_list[new_pos][0] in seed_ref_difference_dict: #and new_pos in hifi_sorted_list:
			if hifi_sorted_list[new_pos][0] in seed_ref_difference_dict and math.fabs(hifi_sorted_list[new_pos][0] - position) < position_distance:
				seed = hifi_sorted_list[new_pos][1]
				max_base = keywithmaxval(seed.allele_dict)
				max_value = seed.allele_dict[max_base]
				seed.allele_new_percentage = float(max_value)/float(number_of_subfile)
				if seed.allele_new_percentage*100 >= 0.80: #and position in geno_dict and max_base in geno_dict[position][2]:
					#if True:
					if position in qscore_dict and qscore_dict[position]/float(number_of_subfile-2) >= 0.80:
						seed.allele_new = max_base
						recovered_seed_dict[hifi_sorted_list[new_pos][0]] = seed
	
	
	print "seed_dict seed number", len(seed_dict)
	revised_seed_dict = dict_add(seed_dict, recovered_seed_dict)
	
	print "recovered seed number", len(recovered_seed_dict)
	print "new seed total number", len(revised_seed_dict)
	output_revised_seed("haplotype_recoved.txt", recovered_seed_dict)
	file_name = "haplotype_expanded.txt"
	output_revised_seed(file_name, revised_seed_dict)
	seed_std_compare(file_name, chr_name)
		#position_distance += 500
	
	#hifi_run(file_name, chr_name)
	#hifiAccuCheck("imputed_"+file_name, chr_name)

	return revised_seed_dict


def seed_expand_geno():
	global seed_hetero_dict
	global random_geno_seed_dict
	
	seed_geno_difference_dict = dict_substract(geno_hetero_dict, seed_hetero_dict)
	
	print "seed_geno_difference_dict", len(seed_geno_difference_dict)
	
	seed_geno_difference_sorted_list = sort_dict_by_key(seed_geno_difference_dict)
	
	random_geno_seed_dict = {}
	geno_seed_selected_number = 1000
	for i in range (0, geno_seed_selected_number):
		random_index = random.randrange(0,len(seed_geno_difference_dict))
		position = seed_geno_difference_sorted_list[random_index][0]
		while position in random_geno_seed_dict:
			random_index = random.randrange(0,len(seed_geno_difference_dict))
			position = seed_geno_difference_sorted_list[random_index][0]
		random_geno_seed_dict[position] = seed_geno_difference_dict[position]
				
		
	geno_seed_added_each_file = geno_seed_selected_number/number_of_subfile
	print "geno_seed_added_each_file: ", geno_seed_added_each_file
	
	random_geno_seed_sorted_list = sort_dict_by_key(random_geno_seed_dict)
	
	for geno_allele_index in (0, 1): # 0 for A, 1 for B
	
		file_number = 0
		while file_number < number_of_subfile:		
			output_subfile_name = seed_file_name +  "_" + str(geno_allele_index) +  "_" + str(file_number) + ".txt"
			#print "output: ", output_subfile_name
			output_subfile = open(currentPath + output_subfile_name, "w")
			print >> output_subfile, seed_title_info
			seed_hetero_dict_bkup = seed_hetero_dict.copy()
			print "seed_hetero_sorted_list original", len(seed_hetero_dict)
	
			i = 0 
			while i < int(geno_seed_added_each_file):			
				try:		
					forward_index = i*number_of_subfile+file_number
					position = random_geno_seed_sorted_list[forward_index][0]
					elements = random_geno_seed_sorted_list[forward_index][1]
					if position not in seed_hetero_dict:
						seed = seeds()
						seed.rsID = elements[0].strip()
						seed.position = int(elements[1].strip())
						""" take A from geno first """
						seed.allele_ori = elements[2].strip()[geno_allele_index]
						#print position, seed.position, seed.allele_ori, elements[2].strip()
						seed_hetero_dict[position] = seed
				except:
					#print "error", forward_index
					pass
				i += 1
			"""last position in seed may not be the same with others this needs to be checked"""
			print "seed_hetero_sorted_list new", len(seed_hetero_dict)
			
			sub_seed_dict = dict_add(seed_hetero_dict, seed_homo_dict)				
	
			sub_seed_list = sort_dict_by_key(sub_seed_dict)
			
			print "new sub_seed_list size", len(sub_seed_list)
			for seed in sub_seed_list:
				line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
				print >> output_subfile, line	
			output_subfile.close()
			seed_hetero_dict = seed_hetero_dict_bkup.copy()
			geno_subfile_name = "genotype.txt"
			ref_subfile_name = "refHaplos.txt"
			if file_number < (number_of_subfile-2):
				hifi = program_path + "hifi_fu_ref " + output_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " &"
				#hifi = program_path + "hifi_fu " + output_subfile_name + " &"
			else:
				hifi = program_path + "hifi_fu_ref " + output_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + ""
				#hifi = program_path + "hifi_fu " + output_subfile_name + ""	# make sure the other hifi processes are finished
			hifi_process = subprocess.Popen(hifi, shell=True)
			hifi_process.wait()	
			
			file_number += 1
	
def seed_geno_extract():
	revised_seed_dict = {}
	
	for geno_allele_index in (0, 1): # 0 for A, 1 for B
	
		hifi_dict = {}
		qscore_dict = {}
			
		file_number = 0
		while file_number < number_of_subfile:	
			input_subfile_name = "imputed_haplotype_" + str(geno_allele_index) + "_" + str(file_number) + ".txt"
			hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)
			file_number += 1
		
		file_number = 0
		while file_number < number_of_subfile:	
			input_subfile_name = "qscore_haplotype_" + str(geno_allele_index) + "_" + str(file_number) + ".txt"
			qscore_dict = load_qscore_result(input_subfile_name, qscore_dict)
			file_number += 1
		print "qscore_dict size: ", len(qscore_dict)
		
		random_geno_seed_sorted_list = sort_dict_by_key(random_geno_seed_dict)
		
		"""
		for position, snp in hifi_dict.iteritems():
			if position in qscore_dict and qscore_dict[position]/float(number_of_subfile) == 1.0:
				seed = hifi_dict[position]
				max_base = keywithmaxval(seed.allele_dict)
				max_value = seed.allele_dict[max_base]
				seed.allele_new = max_base
				revised_seed_dict[position] = seed
		
		print "added seed by qs = 1", len(revised_seed_dict)
		"""
		for snp in random_geno_seed_sorted_list:
			position = snp[0]
			if position not in revised_seed_dict and position in hifi_dict:
				seed = hifi_dict[position]
				max_base = keywithmaxval(seed.allele_dict)
				max_value = seed.allele_dict[max_base]
				seed.allele_new_percentage = float(max_value)/float(number_of_subfile)
				if seed.allele_new_percentage*100 >= 90 and max_base == geno_dict[position][2][geno_allele_index]:
					if position in qscore_dict and qscore_dict[position]/float(number_of_subfile) >= 0.90:
						seed.allele_new = max_base
						revised_seed_dict[position] = seed
		print "added seed by", geno_allele_index, len(revised_seed_dict)

	print "added seed total number", len(revised_seed_dict)
	revised_seed_dict = dict_add(seed_dict, revised_seed_dict)
	print "new seed total number", len(revised_seed_dict)
	output_revised_seed("haplotype_expanded.txt", revised_seed_dict)
	return revised_seed_dict	
	
def output_revised_seed(filename, revised_seed_dict):
	seed_new_file = open(currentPath + filename, "w")
	print >> seed_new_file, seed_title_info
	revised_seed_sorted_list = sort_dict_by_key(revised_seed_dict) 	# need to sort the snps by position
	for snp in revised_seed_sorted_list:
		seed = snp[1]
		line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
		print >> seed_new_file, line
	seed_new_file.close()
	
def output_revised_seed_dict(filename, revised_seed_dict):
	seed_new_file = open(currentPath + filename, "w")
	print >> seed_new_file, seed_title_info
	revised_seed_sorted_list = sort_dict_by_key(revised_seed_dict) 	# need to sort the snps by position
	for snp in revised_seed_sorted_list:
		seed = snp[1]
		line = seed[0] + "\t" + str(seed[1]) + "\t" + seed[2]
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


def compare_std_hap_ref():
	ref_file_name = "refHaplos.txt"
	hap_ref_dict = load_raw_data(ref_file_name, raw_data_format)[1]
	
	std_x_dict = {}
	std_not_in_ref_dict = {}
	a = 0
	
	for position, snp in hap_std_dict.iteritems():
		if position in hap_ref_dict:
			exist_in_ref = False
			std_A = snp[2]
			std_B = snp[3]
			alleles = hap_ref_dict[position]
			alleles = alleles[2:]
			unique_alleles = list(set(alleles))
			n_alleles = len(unique_alleles)
			if n_alleles == 0 or n_alleles > 2:
				print "error in: ", position
				sys.exit(1)
			else:
				try:
					if n_alleles == 1:
						a += 1
					if std_A == std_B:
						if n_alleles == 2:
							if std_A == unique_alleles[0] or std_A == unique_alleles[1]:
								exist_in_ref = True
						if n_alleles == 1:
							if std_A == unique_alleles[0]:
								exist_in_ref = True
					else:
						if n_alleles == 2:
							if (std_A == unique_alleles[0] and std_B == unique_alleles[1]) or (std_A == unique_alleles[1] and std_B == unique_alleles[0]):
								exist_in_ref = True
							if n_alleles == 1:	# hetero_std, homo_ref
								pass
					if not exist_in_ref:
						if std_A == 'X':
							std_x_dict[position] = list_to_line(snp)
						else:
							std_not_in_ref_dict[position] = list_to_line(snp)
							print position, unique_alleles, n_alleles, std_A, std_B, geno_dict[position][2]
				except:
					#print position, unique_alleles, n_alleles
					pass
	print "std_x_dict: ", len(std_x_dict)
	print "std_not_in_ref_dict: ", len(std_not_in_ref_dict)
	print a


def compare_geno_ref():
	ref_file_name = "refHaplos.txt"
	hap_ref_dict = load_raw_data(ref_file_name, raw_data_format)[1]
	
	std_x_dict = {}
	std_n_dict = {}
	std_not_in_ref_dict = {}
	ref_homo = 0
	
	for position, snp in geno_dict.iteritems():
		if position in hap_ref_dict:
			exist_in_ref = False
			std_A = snp[2][0]
			std_B = snp[2][1]
			alleles = hap_ref_dict[position]
			alleles = alleles[2:]
			unique_alleles = list(set(alleles))
			n_alleles = len(unique_alleles)
			if n_alleles == 0 or n_alleles > 2:
				print "error in: ", position
				sys.exit(1)
			else:
				try:
					if n_alleles == 1:
						ref_homo += 1
					if std_A == std_B:
						if n_alleles == 2:
							if std_A == unique_alleles[0] or std_A == unique_alleles[1]:
								exist_in_ref = True
						if n_alleles == 1:
							if std_A == unique_alleles[0]:
								exist_in_ref = True
					else:
						if n_alleles == 2:
							if (std_A == unique_alleles[0] and std_B == unique_alleles[1]) or (std_A == unique_alleles[1] and std_B == unique_alleles[0]):
								exist_in_ref = True
							if n_alleles == 1:	# hetero_std, homo_ref
								pass
					if not exist_in_ref:
						if std_A == 'X':
							std_x_dict[position] = list_to_line(snp)
						if std_A == 'N':
							std_n_dict[position] = list_to_line(snp)
						else:
							std_not_in_ref_dict[position] = list_to_line(snp)
							print position, unique_alleles, n_alleles, std_A, std_B, geno_dict[position][2]
				except:
					#print position, unique_alleles, n_alleles
					pass
	print "std_x_dict: ", len(std_x_dict)
	print "std_n_dict: ", len(std_n_dict)
	print "std_not_in_ref_dict: ", len(std_not_in_ref_dict)
	print "ref_homo", ref_homo

def generate_std_seed(seed_number):
		revised_seed_dict = load_seed_data(file_path+"ASW_chr9_child_hap_refed.txt")[1]
		revised_seed_list = sort_dict_by_key(revised_seed_dict)
		
		group_tuple = group_seed(revised_seed_dict, geno_dict)
		seed_homo_dict = group_tuple[0]
		seed_hetero_dict = group_tuple[1]
		
		seed_hetero_list = sort_dict_by_key(seed_hetero_dict)
		
		selected_seed_dict = {}
		i = 0
		while i < seed_number-1:
			random_index = random.randrange(0,(len(seed_hetero_list)-1))
			while revised_seed_list[random_index][0] in selected_seed_dict:
				random_index = random.randrange(0,(len(seed_hetero_list)-1))
			selected_seed_dict[seed_hetero_list[random_index][0]] = seed_hetero_list[random_index][1]
			i += 1
		selected_seed_dict[revised_seed_list[-1][0]] = revised_seed_list[-1][1]
				
		file_name = "haplotype_std.txt"
		output_revised_seed(file_name, selected_seed_dict)
		seed_std_compare(file_name, chr_name)
		hifi_run(file_name, chr_name)
		hifiAccuCheck("imputed_"+file_name, chr_name)

def combine_hifi_seed(ori_seed_file):
	
	revised_seed_dict = {}
	hifi_dict = {}
	
	ori_seed_tuple = load_seed_data(ori_seed_file)
	ori_seed_dict = ori_seed_tuple[1]
		
	for i in range (0,ref_cycle_number):
		input_subfile_name = "haplotype_expanded.txt_" + str(i)
		hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)
	
	print "hifi_dict seed total number", len(hifi_dict)
	
	for position, snp in hifi_dict.iteritems():
		if position not in ori_seed_dict:
			seed = hifi_dict[position]
			max_base = keywithmaxval(seed.allele_dict)
			max_value = seed.allele_dict[max_base]
			seed.allele_new_percentage = float(max_value)/float(ref_cycle_number)
			if seed.allele_new_percentage == 1.00:
				seed.allele_new = max_base
				revised_seed_dict[position] = seed
			else:
				#print seed.allele_new_percentage 
				pass

	print "consistent seed total number", len(revised_seed_dict)
	revised_seed_dict = dict_add(ori_seed_dict, revised_seed_dict)
	print "new seed total number", len(revised_seed_dict)
	output_filename = "haplotype.txt"
	output_revised_seed(output_filename, revised_seed_dict)
	seed_std_compare(output_filename, chr_name)

def multple_ref_expand(seed_file, chr_name, mode):
	#global seed_dict
	global ref_cycle_number
	ref_cycle_number = 3
	sub_cycle = 2
	
	os.system("cp haplotype.txt haplotype_ori.txt")
	
	for i in range (0,1):
		"""
		print "########### overall expand cycle #########", i
		for j in range (0,ref_cycle_number):
			os.system("cp haplotype_ori.txt haplotype.txt")
			mode = "ref"
			print "########### ref expand cycle #########", j
			seed_correction(seed_file, chr_name, mode)
			os.system("cp haplotype_expanded.txt haplotype.txt")
			
			for k in range (0,sub_cycle):
				print "*************** split *************",  k
				mode = "split"
				seed_correction(seed_file, chr_name, mode)
				os.system("cp haplotype_new.txt haplotype.txt")	
			
			os.system("cp haplotype_new.txt haplotype_expanded.txt")
			file_name = "haplotype_expanded.txt"
			output_file = file_name + "_" + str(j)
			os.system("cp haplotype_expanded.txt " + output_file)
		"""
		print "*************** combine expanded hifi results *************"
		combine_hifi_seed("haplotype_ori.txt")
		"""
		for j in range (0,sub_cycle):
			print "*************** split *************",  j
			mode = "split"
			seed_correction(seed_file, chr_name, mode)
			os.system("cp haplotype_new.txt haplotype.txt")	
		"""	
		for j in range (0,sub_cycle):
			print "*************** geno expand *************",  j
			mode = "geno"
			seed_correction(seed_file, chr_name, mode)
			os.system("cp haplotype_expanded.txt haplotype.txt")	
		
			for k in range (0,sub_cycle):
				print "*************** split *************",  k
				mode = "split"
				seed_correction(seed_file, chr_name, mode)
				os.system("cp haplotype_new.txt haplotype.txt")	
		
		for j in range (0,sub_cycle):
			print "*************** split *************",  j
			mode = "split"
			seed_correction(seed_file, chr_name, mode)
			os.system("cp haplotype_new.txt haplotype.txt")	
		
	
	file_name = "haplotype.txt"
	hifi_run(file_name, chr_name)
	hifiAccuCheck("imputed_"+file_name, chr_name)
	
	"""
	seed_dict_bkup = seed_dict.copy()
	print "seed_dict before", len(seed_dict)
	for i in range (0,ref_cycle_number):
		print "########### ref expand cycle #########", i
		file_name = "haplotype_expanded.txt"
		seed_expand_ref()
		seed_recover_extract_ref()
		
		output_file = file_name + "_" + str(i)
		os.system("cp haplotype_expanded.txt " + output_file)
		seed_std_compare(output_file, chr_name)
	print "seed_dict after", len(seed_dict)
	seed_dict = seed_dict_bkup.copy()
	"""
	
	
	#return revised_seed_dict	
	
	
	

def seed_correction(seed_file, chr_name, mode):

	global geno_dict
	global geno_homo_dict
	global geno_hetero_dict
	global hap_std_dict
	global seed_dict
	global seed_homo_dict
	global seed_hetero_dict
	global hap_ref_dict
	
	global seed_title_info
	global ref_title_info
	global geno_title_info
	global seed_file_name
	global number_of_subfile
		
	ref_file_name = "refHaplos.txt"
	ref_tuple = load_raw_data(ref_file_name, raw_data_format)
	ref_title_info = ref_tuple[0]
	hap_ref_dict = ref_tuple[1]

	print "seed_file is :", seed_file
	seed_file_name = seed_file[:seed_file.find('.')].strip()
	seed_tuple = load_seed_data(seed_file)
	seed_title_info = seed_tuple[0]
	seed_dict = seed_tuple[1]

	#print "seed_title_info", seed_title_info
	print "total_seed_number: ", len(seed_dict)
	

	genotype_file = file_path + "genotype_NA10847_"+chr_name+".txt"	# for all
	geno_tuple = load_raw_data(genotype_file, raw_data_format)
	geno_title_info = geno_tuple[0]
	geno_dict = geno_tuple[1]
	print "total_geno_number: ", len(geno_dict)

	hap_std_file = file_path + "ASW_"+chr_name+"_child_hap_refed.txt"	
	hap_std_dict = load_hap_std(hap_std_file)

	group_tuple = group_seed(seed_dict, geno_dict)
	seed_homo_dict = group_tuple[0]
	seed_hetero_dict = group_tuple[1]
	print "seed_homo_dict", len(seed_homo_dict)
	print "seed_hetero_dict", len(seed_hetero_dict)
	
	group_tuple = group_seed(geno_dict, geno_dict)
	geno_homo_dict = group_tuple[0]
	geno_hetero_dict = group_tuple[1]
	print "geno_homo_dict", len(geno_homo_dict)
	print "geno_hetero_dict", len(geno_hetero_dict)
			
	if mode == "split":
		seed_std_compare(seed_file, chr_name)
		seed_error_remove(seed_dict)
		seed_error_remove_extract(seed_dict)
		seed_std_compare("haplotype_new.txt", chr_name)
	elif mode == "recover":
		#revised_seed_dict = seed_error_remove_extract(seed_dict)
		#seed_dict = load_seed_data("haplotype.txt_ori")[1]
		revised_seed_dict = load_seed_data("haplotype_new.txt")[1]
		seed_recover(seed_dict, revised_seed_dict)
		seed_recover_extract(seed_hetero_dict, revised_seed_dict)
		seed_std_compare("haplotype_new.txt", chr_name)
	elif mode == "add":
		#revised_seed_dict = seed_error_remove_extract(seed_dict)
		revised_seed_dict = load_seed_data("haplotype_new.txt")[1]
		seed_recover_extract(seed_hetero_dict, revised_seed_dict)
		seed_std_compare("haplotype_new.txt", chr_name)
	elif mode == "expand":
		seed_std_compare("haplotype.txt", chr_name)
		seed_expand(seed_file)
		same_to_A_dict = seed_std_compare("haplotype_t.txt", chr_name)
		same_to_A_dict = load_seed_data_from_dict(same_to_A_dict)
		same_to_A_dict = dict_add(same_to_A_dict, seed_homo_dict)
		output_revised_seed("haplotype_expanded.txt", same_to_A_dict)

		file_name = "haplotype_expanded.txt"
		hifi_run(file_name, chr_name)
		hifiAccuCheck("imputed_"+file_name, chr_name)
		"""
		seed_tuple = load_seed_data(seed_file)
		seed_title_info = seed_tuple[0]
		seed_dict = seed_tuple[1]

		seed_error_remove(seed_dict)
		seed_error_remove_extract(seed_dict)
		seed_std_compare("haplotype.txt", chr_name)
		
		seed_tuple = load_seed_data(seed_file)
		seed_title_info = seed_tuple[0]
		seed_dict = seed_tuple[1]

		seed_error_remove(seed_dict)
		seed_error_remove_extract(seed_dict)
		seed_std_compare("haplotype.txt", chr_name)
		"""
	elif mode == "ref":
		file_name = "haplotype_expanded.txt"
		seed_expand_ref()
		seed_recover_extract_ref()
		seed_std_compare(file_name, chr_name)
		#hifi_run(file_name, chr_name)
		#hifiAccuCheck("imputed_"+file_name, chr_name)
	elif mode == "refe":
		#seed_expand_ref()
		file_name = "haplotype_expanded.txt"
		seed_recover_extract_ref()
		seed_tuple = seed_std_compare(file_name, chr_name)
		#same_to_A_dict = seed_tuple[0]
		#same_to_A_dict = load_seed_data_from_dict(same_to_A_dict)
		#same_to_A_dict = dict_add(same_to_A_dict, seed_homo_dict)
		#output_revised_seed("haplotype_a.txt", same_to_A_dict)
		"""
		hifi_run(file_name, chr_name)
		hifiAccuCheck("imputed_"+file_name, chr_name)
		
		file_name = "haplotype_a.txt"
		hifi_run(file_name, chr_name)
		hifiAccuCheck("imputed_"+file_name, chr_name)
		"""
	
	elif mode == "mref":
		multple_ref_expand(seed_file, chr_name, mode)
	
	elif mode == "geno":
		seed_expand_geno()
		seed_geno_extract()
		file_name = "haplotype_expanded.txt"
		seed_std_compare(file_name, chr_name)
	elif mode == "genoe":
		seed_geno_extract()
		file_name = "haplotype_expanded.txt"
		seed_std_compare(file_name, chr_name)
	
	elif mode == "combine":
		combine_hifi_seed()
	
		
	else:
		
		"""check error seed"""
		file_name = "haplotype_expanded.txt"
		seed_tuple = seed_std_compare(file_name, chr_name)
		same_to_A_dict = seed_tuple[0]
		same_to_A_dict = load_seed_data_from_dict(same_to_A_dict)
		same_to_B_dict = seed_tuple[1]
		same_to_B_dict = load_seed_data_from_dict(same_to_B_dict)
		#output_revised_seed("haplotype_A.txt", same_to_A_dict)
		#output_revised_seed("haplotype_B.txt", same_to_B_dict)
		position = 409399
		position = 1109818
		position = 128659696
		
		temp_dict = {}
		temp_dict[position] = same_to_B_dict[position]
		print position, same_to_B_dict[position].allele_new
		print position, hap_std_dict[position][2], hap_std_dict[position][3]
		temp_dict[position].allele_new = hap_std_dict[position][2]
		print position, temp_dict[position].allele_new
		#same_to_A_dict = dict_add(same_to_A_dict, temp_dict)
		same_to_A_dict = dict_add(same_to_A_dict, seed_homo_dict)
		
		output_revised_seed("haplotype_one.txt", same_to_A_dict)
		seed_tuple = seed_std_compare("haplotype_one.txt", chr_name)
		
		hifi_run("haplotype_one.txt", chr_name)
		hifi_dict = {}
		hifi_dict = load_hifi_result("imputed_haplotype_one.txt", hifi_dict)
		a = 0
		hifi_B_dict = {}
		hifi_B_dict[position] = same_to_B_dict[position]
		for position, snp in same_to_B_dict.iteritems():
			
			if snp.allele_new != hifi_dict[position].allele_new:
				a += 1
				hifi_B_dict[position] = hifi_dict[position]
				print position, snp.allele_new, hifi_dict[position].allele_new, 
				if position in hap_std_dict:
					print hap_std_dict[position][2], hap_std_dict[position][3]
		print a
		file_name = "refHaplos.txt"
		for position, seed in hifi_B_dict.iteritems():
			calculate_maf(file_name, position)
		#output_revised_seed("haplotype_hifi_B.txt", hifi_B_dict)
		
		#pass
		#compare_std_geno_ref()
		#compare_geno_ref()
		#generate_std_seed(8000)
		"""
		seed_12000 = load_raw_data("haplotype_a.txt", raw_data_format)[1]
		seed_std = load_raw_data("haplotype_std.txt", raw_data_format)[1]
		difference_dict = dict_substract(seed_std, seed_12000)
		output_revised_seed_dict("haplotype_12000_std_dif.txt", difference_dict)
		
		added_dict = dict_add(seed_12000, seed_std)
		output_revised_seed_dict("haplotype_12000_std_added.txt", added_dict)
		
		
		revised_seed_dict = seed_error_remove_extract(seed_dict)
		same_to_B_dict = seed_std_compare("haplotype_new.txt", chr_name)
		print same_to_B_dict
		file_name = "refHaplos.txt"
		
		hifiAccuCheck("imputed_haplotype_1.txt", chr_name)
		for position, seed in same_to_B_dict.iteritems():
			calculate_maf(file_name, position)	
		error_seed_distance(seed_dict, same_to_B_dict)
		error_seed_distance(seed_hetero_dict, same_to_B_dict)
		"""
		
		
		
		"""
		output_revised_seed_without_error(revised_seed_dict, same_to_B_dict)
		seed_std_compare("haplotype_without_error.txt", chr_name)
		output_revised_seed_with_error(revised_seed_dict, same_to_B_dict)
		"""


if __name__=='__main__':
	options = get_args()
	seed_file = options.seedFile
	chr_name = options.chrName
	mode = options.mode	
	#seed_correction(seed_file, chr_name, mode)
	"""
	# for removing error seed
	for i in range (0,1):
		seed_correction(seed_file, chr_name, mode)
		os.system("mkdir -p " + str(i))
		os.system("cp haplotype_new* imputed_haplotype_* haplotype_?.txt " + str(i))
		os.system("cp haplotype_new.txt haplotype.txt")
		
	
	# for adding new seed
	for i in range (0,4):
		seed_correction(seed_file, chr_name, mode)
		os.system("mkdir -p " + str(i))
		os.system("cp haplotype_new* imputed_haplotype_* haplotype_ori_?.txt haplotype_ori_10.txt haplotype_recoved.txt " + str(i))
		#os.system("cp haplotype_recoved.txt haplotype_new.txt")
	"""
	"""
	# for ref expand
	for i in range (0,3):
		print "########### ref #########", i
		mode = "ref"
		seed_correction(seed_file, chr_name, mode)
		os.system("cp haplotype_expanded.txt haplotype.txt")
		for j in range (0,3):
			print "*************** split *************",  j, "ref", i
			mode = "split"
			seed_correction(seed_file, chr_name, mode)
			os.system("cp haplotype_new.txt haplotype.txt")	
		mode = "ref"
	"""
	"""
	for i in range (0,5):
		print "########### ref expand cycle #########", i
		seed_correction(seed_file, chr_name, mode)
		os.system("cp haplotype_expanded.txt haplotype.txt")
	
	for i in range (0,1):
		print "########### geno expand cycle #########", i
		seed_correction(seed_file, chr_name, mode)
		os.system("cp haplotype_new.txt haplotype.txt")
	
	file_name = "haplotype.txt"
	hifi_run(file_name, chr_name)
	hifiAccuCheck("imputed_"+file_name, chr_name)
	"""
	
	for i in range (0,1):
		seed_correction(seed_file, chr_name, mode)
		#multple_ref_expand(seed_file, chr_name, mode)
	

