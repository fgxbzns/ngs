#!/usr/bin/python


# ######################################################################################
# Author: Guoxing Fu
# Impute haplotype from low-depth single chromosome sequencing data
# ######################################################################################

import os, glob, subprocess, random, operator, time, math, copy
from optparse import OptionParser

from tools import *
from seed_std_compare import seed_std_compare
from calculate_maf import calculate_maf
from hifiAccuCheck_v2 import hifiAccuCheck
from cluster import get_cluster
from refMerger_v5 import refMerger
import multiprocessing
from data_dicts import *




class seeds:
	def __init__(self):
		self.rsID = ""
		self.position = 0
		self.allele_ori = ""
		self.allele_new = ""
		self.allele_new_percentage = 0
		self.allele_dict = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
		self.bridge_info = []


class bridge_info:
	def __init__(self):
		self.refID = ""
		self.ref_allele = ""
		self.first_window_size = 0
		self.second_window_size = 0
		self.first_second_size = 0
		self.gap_size = 0


def output_revised_seed(filename, revised_seed_dict):
	seed_new_file = open(currentPath + filename, "w")
	print >> seed_new_file, data_dict.seed_title_info
	revised_seed_sorted_list = sort_dict_by_key(revised_seed_dict)  # need to sort the snps by position
	for snp in revised_seed_sorted_list:
		seed = snp[1]
		line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
		print >> seed_new_file, line
	seed_new_file.close()


def output_revised_seed_dict(filename, revised_seed_dict):
	seed_new_file = open(currentPath + filename, "w")
	print >> seed_new_file, data_dict.seed_title_info
	revised_seed_sorted_list = sort_dict_by_key(revised_seed_dict)  # need to sort the snps by position
	for snp in revised_seed_sorted_list:
		seed = snp[1]
		line = seed[0] + "\t" + str(seed[1]) + "\t" + seed[2]
		print >> seed_new_file, line
	seed_new_file.close()


def output_revised_seed_without_error(revised_seed_dict, same_to_B_dict):
	seed_new_file = open(currentPath + "haplotype_without_error.txt", "w")
	print >> seed_new_file, data_dict.seed_title_info
	for position, seed in revised_seed_dict.iteritems():
		if position not in same_to_B_dict:
			line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
			print >> seed_new_file, line
	seed_new_file.close()


def output_revised_seed_with_error(revised_seed_dict, same_to_B_dict):
	seed_new_file = open(currentPath + "haplotype_with_error.txt", "w")
	print >> seed_new_file, data_dict.seed_title_info
	for position, seed in revised_seed_dict.iteritems():
		if position in same_to_B_dict:
			line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_ori
			print >> seed_new_file, line
	seed_new_file.close()


def hifi_process(file_number, number_of_subfile, hap_subfile_name, geno_subfile_name="genotype.txt",
                 ref_subfile_name="refHaplos.txt"):
	maf_step = float(random.randrange(10, 40)) / (100.0)
	print "maf_step is: ", maf_step
	if file_number < (number_of_subfile - 3):
		#hifi = program_path + "hifi_fu_revise " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " " + str(
		hifi = "./hifi_fu_revise " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " " + str(
			maf_step) + " &"
	else:
		hifi = "./hifi_fu_revise " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " " + str(
			maf_step) + ""
		#hifi = program_path + "hifi_fu_revise " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " " + str(

	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()


def load_hap_qscore(number_of_subfile):
	hifi_dict = {}
	qscore_dict = {}
	for file_number in range(number_of_subfile):
		hifi_subfile_name = "imputed_haplotype_" + str(file_number) + ".txt"
		hifi_dict = load_hifi_result(hifi_subfile_name, hifi_dict)
		qscore_subfile_name = "qscore_haplotype_" + str(file_number) + ".txt"
		qscore_dict = load_qscore_result(qscore_subfile_name, qscore_dict)
	return (hifi_dict, qscore_dict)


def load_window_info(number_of_subfile):
	window_info_dict = {}
	for file_number in range(number_of_subfile):
		window_subfile_name = "window_haplotype_" + str(file_number) + ".txt"
		temp_data_dict = load_raw_data(window_subfile_name)[1]
		for position, elements in temp_data_dict.iteritems():
			try:
				if position not in window_info_dict:
					window_info_dict[position] = {}
				if file_number not in window_info_dict[position]:
					window_info_dict[position][file_number] = elements
				else:
					print "file_number duplicates", file_number
			except:
				# print "error at ", file_name, position, elements
				pass
	print "total snp number in: ", len(window_info_dict)
	return window_info_dict


def seed_error_remove():
	# reduce error seed from ori seed

	number_of_subfile = data_dict.number_of_subfile
	print "seed_hetero_dict new", len(data_dict.seed_hetero_dict)

	seed_hetero_sorted_list = sort_dict_by_key(data_dict.seed_hetero_dict)
	seed_number_ceilling = int(math.ceil(float(len(seed_hetero_sorted_list)) / 100) * 100)
	# print "seed_hetero_number_ceilling: ", seed_number_ceilling
	seed_removed_in_each_subfile = seed_number_ceilling / number_of_subfile
	print "hetero_seed_removed_in_each_subfile: ", seed_removed_in_each_subfile
	#seed_removed_in_last_subfile = int(math.fmod(len(seed_hetero_sorted_list), seed_removed_in_each_subfile))
	#print "seed_removed_in_last_subfile: ", seed_removed_in_last_subfile

	seed_homo_sorted_list = [x for x in data_dict.seed_homo_dict.iteritems()]

	print "seed_hetero_sorted_list", len(seed_hetero_sorted_list)
	for file_number in range(number_of_subfile):
		hap_subfile_name = data_dict.seed_file_name + "_" + str(file_number) + ".txt"
		output_subfile = open(currentPath + hap_subfile_name, "w")
		print >> output_subfile, data_dict.seed_title_info

		seed_hetero_dict_bkup = data_dict.seed_hetero_dict.copy()

		for i in range(int(seed_removed_in_each_subfile)):
			try:
				forward_index = i * number_of_subfile + file_number
				backward_index = (i + 1) * number_of_subfile - file_number
				position = seed_hetero_sorted_list[forward_index][0]
				del data_dict.seed_hetero_dict[position]

				position = seed_hetero_sorted_list[backward_index][0]
				del data_dict.seed_hetero_dict[position]

				#randomly del the second seed in whole range
				random_index = random.randrange((len(seed_hetero_sorted_list) - 1))
				while random_index == (forward_index) or random_index == (backward_index):
					random_index = random.randrange((len(seed_hetero_sorted_list) - 1))

				position = seed_hetero_sorted_list[random_index][0]
				del data_dict.seed_hetero_dict[position]
				""" error, size of the sorted list will change each time after one element is deleted"""
			except:
				pass

		print "seed_hetero_sorted_list new", len(seed_hetero_sorted_list)

		sub_seed_dict = dict_add(data_dict.seed_hetero_dict, data_dict.seed_homo_dict)
		sub_seed_list = sort_dict_by_key(sub_seed_dict)

		print len(sub_seed_list)
		for seed in sub_seed_list:
			line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
			print >> output_subfile, line
		output_subfile.close()
		data_dict.seed_hetero_dict = seed_hetero_dict_bkup.copy()

		hifi_process(file_number, number_of_subfile, hap_subfile_name)
		"""
		process = multiprocessing.Process(target=hifi_mlp, args=(hap_subfile_name))
		process.start()
		process_list.append(process)
	
	for process in process_list:
	    process.join()
		"""


def seed_error_remove_extract():

	number_of_subfile = data_dict.number_of_subfile

	revised_seed_dict = {}
	hifi_dict = data_dict.seed_dict.copy()
	print "hifi_dict initial", len(hifi_dict)

	for file_number in range(number_of_subfile):
		input_subfile_name = "imputed_haplotype_" + str(file_number) + ".txt"
		hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)

	hifi_sorted_list = sort_dict_by_key(hifi_dict)
	for snp in hifi_sorted_list:
		position = snp[0]
		seed = snp[1]
		max_base = keywithmaxval(seed.allele_dict)
		max_value = seed.allele_dict[max_base]
		seed.allele_new_percentage = float(max_value) / float(number_of_subfile)
		if seed.allele_new_percentage * 100 >= 100 and max_base == seed.allele_ori and position in data_dict.seed_dict:
			seed.allele_new = max_base
			revised_seed_dict[position] = seed
		if position in data_dict.seed_hetero_dict and seed.allele_ori != max_base and seed.allele_new_percentage * 100 >= 80:
			hap_std = data_dict.hap_std_dict[position]
			allele_dict = hifi_dict[position].allele_dict
			line = seed.rsID + "\t" + str(seed.position) + "\t" + hap_std[2] + "\t" + \
			       hap_std[3] + "\t" + seed.allele_ori + "\t" + max_base + "\t" + str(
				seed.allele_new_percentage) + "\t" + str(
				allele_dict['A'] + allele_dict['T'] + allele_dict['C'] + allele_dict['G'])

	print "new seed total number", len(revised_seed_dict)
	output_revised_seed("haplotype_error_removed.txt", revised_seed_dict)
	return revised_seed_dict


def seed_recover(seed_dict, revised_seed_dict):

	number_of_subfile = data_dict.number_of_subfile
	"""two usages 1. to recover experiment seed from remove step. 2. to remove errors in expanded seed"""
	removed_seed_dict = dict_substract(data_dict.seed_hetero_dict, revised_seed_dict)

	print "removed_seed_dict", len(removed_seed_dict)
	seed_number_ceilling = int(math.ceil(float(len(removed_seed_dict)) / 100) * 100)
	seed_added_in_each_subfile = seed_number_ceilling / number_of_subfile
	print "hetero_seed_added_in_each_subfile: ", seed_added_in_each_subfile

	removed_seed_list = sort_dict_by_key(removed_seed_dict)

	seed_homo_sorted_list = [x for x in data_dict.seed_homo_dict.iteritems()]

	print "revised_seed_dict", len(revised_seed_dict)
	for file_number in range(number_of_subfile):
		hap_subfile_name = data_dict.seed_file_name + "_" + str(file_number) + ".txt"
		# print "output: ", hap_subfile_name
		output_subfile = open(currentPath + hap_subfile_name, "w")
		print >> output_subfile, data_dict.seed_title_info
		revised_seed_dict_bkup = copy.deepcopy(revised_seed_dict)

		for i in range(int(seed_added_in_each_subfile)):
			try:
				index = i * number_of_subfile + file_number
				position = removed_seed_list[index][0]
				if position not in revised_seed_dict:
					revised_seed_dict[position] = removed_seed_list[index][1]
				else:
					print "seed already in new seed dict"

				random_pos = random.randrange(0, number_of_subfile)
				while random_pos == position:
					random_pos = random.randrange(0, number_of_subfile)
				if random_pos not in revised_seed_dict:
					revised_seed_dict[random_pos] = removed_seed_dict[random_pos]
				else:
					print "seed already in new seed dict"
			except:
				pass

		print "revised_seed_dict new", len(revised_seed_dict)

		revised_seed_list = sort_dict_by_key(revised_seed_dict)

		for seed in revised_seed_list:
			line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
			print >> output_subfile, line
		output_subfile.close()
		revised_seed_dict = copy.deepcopy(revised_seed_dict_bkup)
		hifi_process(file_number, number_of_subfile, hap_subfile_name)


def seed_recover_extract(seed_hetero_dict, revised_seed_dict):

	number_of_subfile = data_dict.number_of_subfile
	recovered_seed_dict = {}
	removed_seed_dict = dict_substract(data_dict.seed_hetero_dict, revised_seed_dict)
	# print "removed_seed_dict", len(removed_seed_dict)
	hifi_dict = {}

	for file_number in range(number_of_subfile):
		input_subfile_name = "imputed_haplotype_" + str(file_number) + ".txt"
		hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)
	print "hifi_dict size: ", len(hifi_dict)
	for position, snp in removed_seed_dict.iteritems():
		if position in hifi_dict:
			seed = hifi_dict[position]
			max_base = keywithmaxval(seed.allele_dict)
			max_value = seed.allele_dict[max_base]
			seed.allele_new_percentage = float(max_value) / float(number_of_subfile)
			#if seed.allele_new_percentage*100 >= 90:	# ref
			if seed.allele_new_percentage >= 0.9 and max_base == data_dict.seed_hetero_dict[position].allele_ori:  # ori
				seed.allele_new = max_base
				recovered_seed_dict[position] = seed

	verify_expanded_seed_by_cluster(recovered_seed_dict)
	revised_seed_dict = dict_add(revised_seed_dict, recovered_seed_dict)
	print "new seed total number", len(recovered_seed_dict)
	output_revised_seed("haplotype_recoved.txt", recovered_seed_dict)
	output_revised_seed("haplotype_expanded.txt", revised_seed_dict)
	return recovered_seed_dict


def seed_expand_qs(seed_file):
	expanded_seed_dict = {}
	hifi_dict = {}
	revised_seed_dict = load_seed_data(seed_file)[1]

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
		if float(qscore_dict[position][3]) == 1.0 and position in data_dict.geno_dict and seed.allele_new in \
				data_dict.geno_dict[position][2]:
			expanded_seed_dict[position] = seed
		elif float(qscore_dict[position][3]) >= 0.9 and position in data_dict.geno_dict and seed.allele_new in \
				data_dict.geno_dict[position][2]:
			a += 1
			expanded_seed_dict[position] = seed

	print "a", a
	revised_seed_dict = dict_add(revised_seed_dict, expanded_seed_dict)

	print "new seed total number", len(expanded_seed_dict)
	output_revised_seed("haplotype_expanded_only.txt", expanded_seed_dict)
	output_revised_seed("haplotype_expanded.txt", revised_seed_dict)
	return revised_seed_dict


def seed_expand_ref_hetero():
	seed_ref_difference_dict = dict_substract(data_dict.hap_ref_dict, data_dict.seed_dict)

	print "seed_ref_difference", len(seed_ref_difference_dict)
	ref_number_ceilling = int(math.ceil(float(len(seed_ref_difference_dict)) / 100) * 100)
	ref_removed_in_each_subfile = ref_number_ceilling / data_dict.number_of_subfile
	print "ref_removed_in_each_subfile: ", ref_removed_in_each_subfile

	seed_ref_difference_list = sort_dict_by_key(seed_ref_difference_dict)
	seed_ref_same_dict = dict_substract(data_dict.hap_ref_dict, seed_ref_difference_dict)

	for file_number in range(number_of_subfile):
		#pos_del_record = open("pos_deledted_"+str(file_number), "w")
		ref_subfile_name = "refHaplos" + "_" + str(file_number) + ".txt"
		ref_subfile = open(currentPath + ref_subfile_name, "w")
		print >> ref_subfile, data_dict.ref_title_info
		seed_ref_difference_dict_bkup = seed_ref_difference_dict.copy()
		geno_dict_bkup = data_dict.geno_dict.copy()

		for i in range(int(ref_removed_in_each_subfile)):
			#for i in range(int(15000)):
			try:
				forward_index = i * number_of_subfile + file_number
				position = seed_ref_difference_list[forward_index][0]
				if position in data_dict.geno_dict:
					del data_dict.geno_dict[position]
				#print >> pos_del_record, seed_ref_difference_list[forward_index][0]		
				del seed_ref_difference_dict[position]

				random_index = random.randrange(0, (len(seed_ref_difference_dict) - 1))
				while random_index == (forward_index):
					random_index = random.randrange(0, (len(seed_ref_difference_dict) - 1))
				position = seed_ref_difference_list[random_index][0]
				if position in data_dict.geno_dict:
					del data_dict.geno_dict[position]
				del seed_ref_difference_dict[position]
			except:
				pass

		print "seed_ref_difference_dict new", len(seed_ref_difference_dict)

		sub_ref_dict = dict_add(seed_ref_difference_dict, seed_ref_same_dict)
		sub_ref_list = sort_dict_by_key(sub_ref_dict)

		for seed in sub_ref_list:
			print >> ref_subfile, list_to_line(seed[1])
		ref_subfile.close()
		seed_ref_difference_dict = seed_ref_difference_dict_bkup.copy()

		# generate new genotype file
		geno_subfile_name = "genotype" + "_" + str(file_number) + ".txt"
		geno_subfile = open(currentPath + geno_subfile_name, "w")
		print >> geno_subfile, data_dict.geno_title_info
		#print "geno_dict reduced", len(geno_dict)
		geno_sorted_list = sort_dict_by_key(data_dict.geno_dict)
		for geno in geno_sorted_list:
			print >> geno_subfile, list_to_line(geno[1])
		geno_subfile.close()
		data_dict.geno_dict = geno_dict_bkup.copy()

		hap_subfile_name = "haplotype" + "_" + str(file_number) + ".txt"

		os.system("cp haplotype.txt " + hap_subfile_name)
		hifi_process(file_number, number_of_subfile, hap_subfile_name, geno_subfile_name, ref_subfile_name)


def load_qscore_result(file_name, qscore_dict):
	data_dict = load_raw_data(file_name, raw_data_format)[1]
	for position, elements in data_dict.iteritems():
		try:
			if position not in qscore_dict:
				qscore_dict[position] = float(elements[3])
			else:
				qscore_dict[position] = qscore_dict[position] + float(elements[3].strip())
		except:
			# print "error at ", file_name, position, elements
			pass
	return qscore_dict


def seed_recover_extract_ref():
	recovered_seed_dict = {}
	seed_ref_difference_dict = dict_substract(data_dict.hap_ref_dict, data_dict.seed_dict)
	hifi_dict, qscore_dict = load_hap_qscore(data_dict.number_of_subfile)

	hifi_sorted_list = sort_dict_by_key(hifi_dict)
	print "hifi_dict size", len(hifi_dict)

	hifi_sorted_pos_list = hifi_dict.keys()

	position_distance = data_dict.ref_position_distance
	expand_range = data_dict.ref_expand_range

	for position, snp in data_dict.seed_hetero_dict.iteritems():
		if position in hifi_sorted_pos_list:
			pos_index = hifi_sorted_pos_list.index(position)
			start_index = (pos_index - expand_range) if pos_index - expand_range >= 0 else 0
			end_index = (pos_index + expand_range) if (pos_index + expand_range) < len(hifi_sorted_list) else len(
				hifi_sorted_list) - 1
			for new_index in range(start_index, end_index):
				if math.fabs(hifi_sorted_list[new_index][0] - position) < position_distance:
					#if True:
					seed = hifi_sorted_list[new_index][1]
					max_base = keywithmaxval(seed.allele_dict)
					max_value = seed.allele_dict[max_base]
					seed.allele_new_percentage = float(max_value) / float(number_of_subfile)
					if seed.allele_new_percentage >= data_dict.allele_new_percentage:  #and position in geno_dict and max_base in geno_dict[position][2]:
						if hifi_sorted_list[new_index][0] in qscore_dict and qscore_dict[position] / float(
								number_of_subfile) >= data_dict.qscore_threshold:
							seed.allele_new = max_base
							recovered_seed_dict[hifi_sorted_list[new_index][0]] = seed

	revised_seed_dict = dict_add(data_dict.seed_dict, recovered_seed_dict)

	print "recovered seed number", len(recovered_seed_dict)
	output_revised_seed("haplotype_recoved.txt", recovered_seed_dict)
	file_name = "haplotype_expanded.txt"
	output_revised_seed(file_name, revised_seed_dict)

	return revised_seed_dict


def seed_recover_extract_ref_1():
	""" no distance restriction """
	recovered_seed_dict = {}
	seed_ref_difference_dict = dict_substract(data_dict.hap_ref_dict, data_dict.seed_dict)
	hifi_dict, qscore_dict = load_hap_qscore(number_of_subfile)
	print "after", len(hifi_dict)
	hifi_dict = dict_substract(hifi_dict, data_dict.seed_dict)
	print "after", len(hifi_dict)
	hifi_sorted_list = sort_dict_by_key(hifi_dict)

	window_info_file = open(currentPath + "window_info_output.txt", "w")
	window_info_dict = load_window_info(number_of_subfile)

	a = 0
	s_to_A = 0
	s_to_B = 0
	for position, snp in hifi_dict.iteritems():
		if position in window_info_dict and len(window_info_dict[position]) >= int(
						number_of_subfile / 2) and position in data_dict.hap_std_dict and \
						data_dict.hap_std_dict[position][2] != 'X':
			window_info = window_info_dict[position]
			print >> window_info_file, position, data_dict.hap_std_dict[position][2], data_dict.hap_std_dict[position][
				3]

			for file_number in range(number_of_subfile):
				if file_number in window_info:
					info = window_info[file_number]
					window_pos = info[9:]
					homo_window = True
					distance_distribution_score = 0
					for w_pos in window_pos:
						distance_distribution_score += (int(w_pos) - position) ** 2
						if int(w_pos) in data_dict.geno_hetero_dict:
							homo_window = False
					distance_distribution_score = math.sqrt(distance_distribution_score) / len(window_pos)

					# temp_hap_list.append(info[4])
					print >> window_info_file, file_number, info[2], info[4], info[6], format(float(info[8]),
					                                                                          "0.4f"), distance_distribution_score, homo_window
					for w_pos in window_pos:
						if int(w_pos) in data_dict.geno_dict:
							print >> window_info_file, w_pos, data_dict.geno_dict[int(w_pos)][2],
					print >> window_info_file, " "
				else:
					print >> window_info_file, "\t", "\t", "\t", "\t"

			temp_hap_list = [info[4] for info in window_info.values()]
			if len(temp_hap_list) >= int(number_of_subfile / 2):
				# remove min and max window size
				temp_size_dict = {file_number: int(info[6]) for file_number, info in window_info.iteritems()}
				del window_info[keywithminval(window_info)]
				del window_info[keywithmaxval(window_info)]

				temp_hap_set = set(temp_hap_list)
				# if qscore_dict[position]/float(len(temp_hap_list)) >= data_dict.qscore_threshold:
				#if len(temp_hap_list) >= int(number_of_subfile/2):
				if True:
					seed = snp
					temp_set_list = list(temp_hap_set)
					if len(temp_hap_set) == 1:
						a += 1
						seed.allele_new = temp_set_list[0][0]
						if seed.allele_new == data_dict.hap_std_dict[position][2]:
							s_to_A += 1
						elif seed.allele_new == data_dict.hap_std_dict[position][3]:
							s_to_B += 1
					elif len(temp_hap_set) == 2:
						seed.allele_new = temp_set_list[0][0] if temp_hap_list.count(
							temp_set_list[0]) >= temp_hap_list.count(temp_set_list[1]) else temp_set_list[1][0]
					recovered_seed_dict[position] = seed

	print "same snp in all hifi", a
	print "s_to_A", s_to_A
	print "s_to_B", s_to_B
	window_info_file.close()

	print "revised_seed_dict seed number", len(recovered_seed_dict)
	verify_expanded_seed_by_cluster(recovered_seed_dict)
	print "revised_seed_dict seed number cluster", len(recovered_seed_dict)

	revised_seed_dict = dict_add(data_dict.seed_dict, recovered_seed_dict)
	"""
	print "****** error seed distance"
	error_seed_distance(revised_seed_dict, recovered_seed_dict)
	"""
	print "recovered seed number", len(recovered_seed_dict)
	# print "new seed total number", len(revised_seed_dict)
	output_revised_seed("haplotype_recoved.txt", recovered_seed_dict)
	file_name = "haplotype_expanded.txt"
	output_revised_seed(file_name, revised_seed_dict)
	"""
	same_to_B_dict = seed_std_compare(file_name, data_dict.chr_name)[1]

	window_info_file = open(currentPath + "window_info_output_B.txt", "w")

	homo_window_total = 0
	for position, snp in same_to_B_dict.iteritems():
		if position in window_info_dict and position in data_dict.hap_std_dict and data_dict.hap_std_dict[position][
			2] != 'X':
			window_info = window_info_dict[position]
			print >> window_info_file, position, data_dict.hap_std_dict[position][2], data_dict.hap_std_dict[position][
				3]

			temp_hap_list = []
			for file_number in range(number_of_subfile):
				if file_number in window_info:
					info = window_info[file_number]
					window_pos = info[9:]
					homo_window = True
					distance_distribution_score = 0
					for w_pos in window_pos:
						distance_distribution_score += (int(w_pos) - position) ** 2
						if int(w_pos) in data_dict.geno_hetero_dict:
							homo_window = False

					distance_distribution_score = math.sqrt(distance_distribution_score)

					temp_hap_list.append(info[4])
					print >> window_info_file, file_number, info[2], info[4], info[6], format(float(info[8]),
					                                                                          "0.4f"), distance_distribution_score, homo_window
					for w_pos in window_pos:
						if int(w_pos) in data_dict.hap_std_dict:
							print >> window_info_file, w_pos, data_dict.hap_std_dict[int(w_pos)][2], \
								data_dict.hap_std_dict[int(w_pos)][3], data_dict.geno_dict[int(w_pos)][2],
						if int(w_pos) in window_info_dict and file_number in window_info_dict[int(w_pos)]:
							print >> window_info_file, window_info_dict[int(w_pos)][file_number][4]
						else:
							print >> window_info_file, " "
					print >> window_info_file, " "

				else:
					print >> window_info_file, "\t", "\t", "\t", "\t"

	window_info_file.close()

	print "homo_window_total", homo_window_total
	"""

	return revised_seed_dict


def remove_single_refID():
	refID_dict = {}

	refID_A_count_dict = {}
	refID_B_count_dict = {}

	refID_list = data_dict.ref_title_info.strip().split()[2:]
	print len(refID_list)

	hifi_dict = {}
	hifi_dict = load_raw_data("imputed_" + data_dict.seed_file)[1]
	print len(hifi_dict)

	window_info_dict = load_raw_data("window_" + data_dict.seed_file)[1]
	print len(window_info_dict)

	for pos, elements in window_info_dict.iteritems():
		window_info = elements[9:]
		match_to_A_index_list = range(2, len(refID_list))
		match_to_B_index_list = range(2, len(refID_list))
		hifi_seq_A = ""
		hifi_seq_B = ""
		for window_pos in window_info:
			window_pos = int(window_pos)
			hifi_seq_A += hifi_dict[window_pos][2]
			hifi_seq_B += hifi_dict[window_pos][3]

			match_to_A_index_list = [index for index in match_to_A_index_list if
			                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][2]]
			match_to_B_index_list = [index for index in match_to_B_index_list if
			                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][3]]
		match_to_A_refID = [refID_list[index] for index in match_to_A_index_list]
		match_to_B_refID = [refID_list[index] for index in match_to_B_index_list]
		for refID in match_to_A_refID:
			refID_A_count_dict[refID] = 0 if refID not in refID_A_count_dict else (refID_A_count_dict[refID] + 1)
		for refID in match_to_B_refID:
			refID_B_count_dict[refID] = 0 if refID not in refID_B_count_dict else (refID_B_count_dict[refID] + 1)
		match_to_A_seq = ""
		match_to_B_seq = ""

		for window_pos in window_info:
			window_pos = int(window_pos)
			if len(match_to_A_index_list) > 0:
				match_to_A_seq += data_dict.hap_ref_dict[window_pos][match_to_A_index_list[0]]
			if len(match_to_B_index_list) > 0:
				match_to_B_seq += data_dict.hap_ref_dict[window_pos][match_to_B_index_list[0]]

		#refID_dict[pos] = (hifi_seq_A, hifi_seq_B, list_to_line(match_to_A_refID), match_to_A_seq, list_to_line(match_to_B_refID), match_to_B_seq)
		refID_dict[pos] = (hifi_seq_A, hifi_seq_B, match_to_A_refID, match_to_A_seq, match_to_B_refID, match_to_B_seq)

	#same_to_B_dict = seed_std_compare("imputed_" + data_dict.seed_file, data_dict.chr_name)[0]

	# output hifi results filtered by number of refID
	hifi_new_output = open("non_one.txt", 'w')

	print >> hifi_new_output, data_dict.seed_title_info
	seed_dict_from_hifi = {}
	for pos in hifi_dict.keys():
		#if (pos in refID_dict and len(refID_dict[pos][2]) <= 50) or (pos in refID_dict and len(refID_dict[pos][4]) <= 50):
		refID_cutoff = 2
		if (pos in refID_dict and len(refID_dict[pos][2]) <= refID_cutoff) or (
						pos in refID_dict and len(refID_dict[pos][4]) <= refID_cutoff):
			pass
		#a = hifi_dict[pos]
		#print >> hifi_new_output, a[0], a[1], a[3], a[2],
		else:
			print >> hifi_new_output, list_to_line(hifi_dict[pos])
			seed_dict_from_hifi[pos] = hifi_dict[pos]
		#print >> hifi_new_output, hifi_dict[pos][0], hifi_dict[pos][1], hifi_dict[pos][2]
	hifi_new_output.close()
	hifiAccuCheck("non_one.txt", chr_name)

	output_revised_seed_dict("seed_from_hifi.txt", seed_dict_from_hifi)


def get_refID():
	log_file = open("refID_log.txt", 'w')

	refID_dict = {}

	refID_A_count_dict = {}
	refID_B_count_dict = {}

	refID_list = data_dict.ref_title_info.strip().split()[2:]
	print len(refID_list)
	# print refID_list

	hifi_dict = {}
	hifi_dict = load_raw_data("imputed_" + data_dict.seed_file)[1]
	print len(hifi_dict)

	window_info_dict = load_raw_data("window_" + data_dict.seed_file)[1]
	print len(window_info_dict)

	for pos, elements in window_info_dict.iteritems():
		window_info = elements[9:]
		match_to_A_index_list = range(2, len(refID_list))
		match_to_B_index_list = range(2, len(refID_list))
		hifi_seq_A = ""
		hifi_seq_B = ""
		for window_pos in window_info:
			window_pos = int(window_pos)
			hifi_seq_A += hifi_dict[window_pos][2]
			hifi_seq_B += hifi_dict[window_pos][3]

			match_to_A_index_list = [index for index in match_to_A_index_list if
			                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][2]]
			match_to_B_index_list = [index for index in match_to_B_index_list if
			                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][3]]
		match_to_A_refID = [refID_list[index] for index in match_to_A_index_list]
		match_to_B_refID = [refID_list[index] for index in match_to_B_index_list]
		for refID in match_to_A_refID:
			refID_A_count_dict[refID] = 0 if refID not in refID_A_count_dict else (refID_A_count_dict[refID] + 1)
		for refID in match_to_B_refID:
			refID_B_count_dict[refID] = 0 if refID not in refID_B_count_dict else (refID_B_count_dict[refID] + 1)
		match_to_A_seq = ""
		match_to_B_seq = ""

		for window_pos in window_info:
			window_pos = int(window_pos)
			if len(match_to_A_index_list) > 0:
				match_to_A_seq += data_dict.hap_ref_dict[window_pos][match_to_A_index_list[0]]
			if len(match_to_B_index_list) > 0:
				match_to_B_seq += data_dict.hap_ref_dict[window_pos][match_to_B_index_list[0]]

		refID_dict[pos] = (hifi_seq_A, hifi_seq_B, match_to_A_refID, match_to_A_seq, match_to_B_refID, match_to_B_seq)

	same_to_B_dict = seed_std_compare("imputed_" + data_dict.seed_file, data_dict.chr_name)[0]

	for pos, items in refID_dict.iteritems():
		print >> log_file, pos, items[0], items[1], len(items[2]), items[3], len(items[4]), items[5]
		#print pos, items[0], items[1], items[2], items[3], items[4], items[5]

	# output hifi results filtered by number of refID
	hifi_new_output = open("non_one.txt", 'w')
	print >> hifi_new_output, data_dict.seed_title_info
	seed_dict_from_hifi = {}
	for pos in hifi_dict.keys():
		if (pos in refID_dict and len(refID_dict[pos][2]) <= 2) or (pos in refID_dict and len(refID_dict[pos][4]) <= 2):
			pass
		#a = hifi_dict[pos]
		#print >> hifi_new_output, a[0], a[1], a[3], a[2],
		else:
			print >> hifi_new_output, list_to_line(hifi_dict[pos])
			seed_dict_from_hifi[pos] = hifi_dict[pos]
		#print >> hifi_new_output, hifi_dict[pos][0], hifi_dict[pos][1], hifi_dict[pos][2]
	hifi_new_output.close()

	output_revised_seed_dict("seed_from_hifi.txt", seed_dict_from_hifi)

	refID_A_count_sorted_list = sort_dict_by_value(refID_A_count_dict)
	refID_B_count_sorted_list = sort_dict_by_value(refID_B_count_dict)

	for i in range(100):
		print >> log_file, refID_A_count_sorted_list[i][0], refID_A_count_sorted_list[i][1], \
			refID_B_count_sorted_list[i][0], refID_B_count_sorted_list[i][1]

	"""" check overlapping of refIDs of adjacent snps"""
	refID_sorted_list = sort_dict_by_key(refID_dict)
	all_zero = 0
	all_zero_dict = {}
	for i in range(len(refID_sorted_list)):
		pos = refID_sorted_list[i][0]
		current_snp_data = refID_sorted_list[i][1]
		if i != 0 and i != len(refID_sorted_list) - 1:
			if pos in same_to_B_dict:
				previous_snp_data = refID_sorted_list[i - 1][1]
				p_refID_A_common = [id for id in current_snp_data[2] if id in previous_snp_data[2]]
				p_refID_B_common = [id for id in current_snp_data[4] if id in previous_snp_data[4]]
				next_snp_data = refID_sorted_list[i + 1][1]
				n_refID_A_common = [id for id in current_snp_data[2] if id in next_snp_data[2]]
				n_refID_B_common = [id for id in current_snp_data[4] if id in next_snp_data[4]]

				print >> log_file, pos, window_info_dict[pos][6], "pre A", len(p_refID_A_common), "pre B", len(
					p_refID_B_common), \
					"next A", len(n_refID_A_common), "next B", len(n_refID_B_common)

				if len(p_refID_A_common) == 0 and len(p_refID_B_common) == 0 and len(n_refID_A_common) == 0 and len(
						n_refID_B_common) == 0:
					all_zero += 1
					all_zero_dict[pos] = ""

	print "all_zero", all_zero
	print "seed_dict_from_hifi before:", len(seed_dict_from_hifi)
	remove_all_zero_dict = dict_substract(seed_dict_from_hifi, all_zero_dict)
	print "remove_all_zero_dict after:", len(remove_all_zero_dict)

	hifi_new_output = open("non_one_remove_all_zero.txt", 'w')
	print >> hifi_new_output, data_dict.seed_title_info
	for pos in remove_all_zero_dict:
		print >> hifi_new_output, list_to_line(hifi_dict[pos])
	hifi_new_output.close()
	seed_std_compare("seed_from_hifi.txt", chr_name)

	print "**************** group here"
	i = 0
	temp_continus_snp = []
	n_refID_A_common = []
	n_refID_B_common = []

	seed_dict_from_linkage = {}
	linkage_size_dict = []

	while i < (len(refID_sorted_list) - 1):
		pos = refID_sorted_list[i][0]
		current_snp_data = refID_sorted_list[i][1]
		next_snp_data = refID_sorted_list[i + 1][1]
		if len(temp_continus_snp) == 0:
			n_refID_A_common = [id for id in current_snp_data[2] if id in next_snp_data[2]]
			print >> log_file, "current_snp_data[2]", current_snp_data[2]
			print >> log_file, "current_snp_data[4]", current_snp_data[4]
			print >> log_file, "next_snp_data[2]", next_snp_data[2]
			print >> log_file, "next_snp_data[4]", next_snp_data[4]

			#print "n_refID_A_common", n_refID_A_common
			n_refID_B_common = [id for id in current_snp_data[4] if id in next_snp_data[4]]
			#print "n_refID_B_common", n_refID_B_common
			#n_refID_A_common = [id for id in current_snp_data[2]] if len(n_refID_A_common) == 0 else n_refID_A_common
			#n_refID_B_common = [id for id in current_snp_data[4]] if len(n_refID_B_common) == 0 else n_refID_B_common
			# to keep a copy of the last common refID of this group. n_refID_A/B_common may become zero in the last check
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
		else:
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
			n_refID_A_common = [id for id in n_refID_A_common if id in next_snp_data[2]]
			n_refID_B_common = [id for id in n_refID_B_common if id in next_snp_data[4]]

		if len(n_refID_A_common) > 0 and len(
				n_refID_B_common) > 0:  # use and to make sure both A and B have common refID
			temp_continus_snp.append(pos)
			i += 1
		else:
			if len(temp_continus_snp) > 0:
				for pos in temp_continus_snp:
					print >> log_file, pos, window_info_dict[pos][6], refID_dict[pos][2], refID_dict[pos][4]

				next_pos = refID_sorted_list[i + 1][0]  # to include the last pos in the group
				print >> log_file, next_pos, window_info_dict[next_pos][6], refID_dict[next_pos][2], \
					refID_dict[next_pos][4]
				print >> log_file, temp_continus_snp, last_common_refID_A, last_common_refID_B

				linkage_size_dict.append(len(temp_continus_snp))

			if len(temp_continus_snp) > 150:
				for pos in temp_continus_snp:
					seed_dict_from_linkage[pos] = list_to_line(hifi_dict[pos])
				#seed_dict_from_linkage.append(list_to_line(hifi_dict[pos]))

			temp_continus_snp = []
			n_refID_A_common = []
			n_refID_B_common = []
			last_common_refID_A = []
			last_common_refID_B = []
			i += 1
			print >> log_file, ""

	linkage_size_dict.sort()
	print "min linkage size:", linkage_size_dict[0]
	print "max linkage size:", linkage_size_dict[-1]

	print "seed_dict_from_hifi before:", len(seed_dict_from_hifi)
	#temp_dict = dict_substract(seed_dict_from_hifi, seed_dict_from_linkage)
	temp_dict = seed_dict_from_linkage
	#temp_dict = dict_add(data_dict.seed_dict, seed_dict_from_linkage)
	print "temp_dict before:", len(temp_dict)

	temp_list = data_dict.seed_dict.keys() + temp_dict.keys()
	temp_list.sort()
	#temp_list = sort_dict_by_key(temp_dict)

	hifi_new_output = open("hifi_from_linkage.txt", 'w')
	print >> hifi_new_output, data_dict.seed_title_info

	for pos in temp_list:
		#print >> hifi_new_output, list_to_line(pos[1])
		#print >> hifi_new_output, data[1]
		if pos in data_dict.seed_dict.keys():
			print >> hifi_new_output, data_dict.seed_dict[pos].rsID, data_dict.seed_dict[pos].position, \
				data_dict.seed_dict[pos].allele_new
		if pos in temp_dict:
			print >> hifi_new_output, temp_dict[pos]

	hifi_new_output.close()
	#hifiAccuCheck("hifi_from_linkage.txt", chr_name)
	seed_std_compare("hifi_from_linkage.txt", chr_name)


def add_seed_by_linkage_longestLD():
	refID_dict = {}

	refID_A_count_dict = {}
	refID_B_count_dict = {}

	refID_list = data_dict.ref_title_info.strip().split()[2:]
	print "refID_list", len(refID_list)
	# print refID_list

	hifi_dict = {}
	#hifi_dict = load_hifi_result("imputed_" + data_dict.seed_file, hifi_dict)
	hifi_dict = load_raw_data("imputed_" + data_dict.seed_file)[1]
	print "hifi_dict", len(hifi_dict)

	window_info_dict = load_raw_data("window_" + data_dict.seed_file)[1]
	print "window_info_dict", len(window_info_dict)

	for pos, elements in window_info_dict.iteritems():
		if elements[2] != 'NN' and elements[4] != 'NN':
			window_info = elements[9:]
			match_to_A_index_list = range(2, len(refID_list))
			match_to_B_index_list = range(2, len(refID_list))
			hifi_seq_A = ""
			hifi_seq_B = ""
			for window_pos in window_info:
				try:
					window_pos = int(window_pos)
					hifi_seq_A += hifi_dict[window_pos][2]
					hifi_seq_B += hifi_dict[window_pos][3]

					match_to_A_index_list = [index for index in match_to_A_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][2]]
					match_to_B_index_list = [index for index in match_to_B_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][3]]
				except:
					pass
					print hifi_seq_A, hifi_seq_B, window_info
					sys.exit(1)
			match_to_A_refID = [refID_list[index] for index in match_to_A_index_list]
			match_to_B_refID = [refID_list[index] for index in match_to_B_index_list]

			for refID in match_to_A_refID:
				refID_A_count_dict[refID] = 0 if refID not in refID_A_count_dict else (refID_A_count_dict[refID] + 1)
			for refID in match_to_B_refID:
				refID_B_count_dict[refID] = 0 if refID not in refID_B_count_dict else (refID_B_count_dict[refID] + 1)
			match_to_A_seq = ""
			match_to_B_seq = ""

			for window_pos in window_info:
				window_pos = int(window_pos)
				if window_pos != 0:
					if len(match_to_A_index_list) > 0 and window_pos in data_dict.hap_ref_dict:
						match_to_A_seq += data_dict.hap_ref_dict[window_pos][match_to_A_index_list[0]]
					if len(match_to_B_index_list) > 0 and window_pos in data_dict.hap_ref_dict:
						match_to_B_seq += data_dict.hap_ref_dict[window_pos][match_to_B_index_list[0]]

			#refID_dict[pos] = (hifi_seq_A, hifi_seq_B, list_to_line(match_to_A_refID), match_to_A_seq, list_to_line(match_to_B_refID), match_to_B_seq)
			refID_dict[pos] = (
				hifi_seq_A, hifi_seq_B, match_to_A_refID, match_to_A_seq, match_to_B_refID, match_to_B_seq)

	print "**************** group here"
	refID_sorted_list = sort_dict_by_key(refID_dict)

	i = 0
	temp_continus_snp = []
	n_refID_A_common = []
	n_refID_B_common = []

	seed_dict_from_linkage = {}
	linkage_size_dict = {}

	snp_ld_length_dict = {}  # to store the LD block size of each snp

	while i < (len(refID_sorted_list) - 1):
		#last_refID_A_common = []
		pos = refID_sorted_list[i][0]
		current_snp_data = refID_sorted_list[i][1]
		next_snp_data = refID_sorted_list[i + 1][1]
		if len(temp_continus_snp) == 0:
			n_refID_A_common = [id for id in current_snp_data[2] if id in next_snp_data[2]]
			n_refID_B_common = [id for id in current_snp_data[4] if id in next_snp_data[4]]
			# to keep a copy of the last common refID of this group. n_refID_A/B_common may become zero in the last check
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
		else:
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
			n_refID_A_common = [id for id in n_refID_A_common if id in next_snp_data[2]]
			n_refID_B_common = [id for id in n_refID_B_common if id in next_snp_data[4]]

		if len(n_refID_A_common) > 0 and len(
				n_refID_B_common) > 0:  # use and to make sure both A and B have common refID
			temp_continus_snp.append(pos)
			i += 1
			last_refID_A_common = n_refID_A_common
		else:
			if len(temp_continus_snp) > 0:
				ld_block_size = len(temp_continus_snp)
				linkage_size_dict[ld_block_size] = (temp_continus_snp, last_refID_A_common)
				for pos in temp_continus_snp:
					if pos not in snp_ld_length_dict:
						snp_ld_length_dict[pos] = ld_block_size
					else:
						print "duplicate pos in snp_ld_length_dict:", pos

			temp_continus_snp = []
			n_refID_A_common = []
			n_refID_B_common = []
			last_common_refID_A = []
			last_common_refID_B = []
			i += 1
		#print >> log_file, ""

	linkage_size_sorted_list = sort_dict_by_key(linkage_size_dict)
	linkage_size_sorted_list.reverse()
	print "max linkage size:", linkage_size_sorted_list[0][0], linkage_size_sorted_list[0][1][1]
	print "min linkage size:", linkage_size_sorted_list[-1][0]

	max_linkage_pos_list = []
	max_linkage_pos_list = linkage_size_sorted_list[0][1][0]

	same_to_A_dict, same_to_B_dict = seed_std_compare("imputed_haplotype.txt", data_dict.chr_name)

	"""
	for list in linkage_size_sorted_list:
		if list[0] > 0:
			#print list[0],#, list[1][0]
			average_maf = 0
			for pos in list[1][0]:
				imputed_seed, seed_frequence = data_dict.hap_ref_allele_frequence_dict[pos][0] \
				if hifi_dict[pos][2] == data_dict.hap_ref_allele_frequence_dict[pos][0][0] else data_dict.hap_ref_allele_frequence_dict[pos][1]
				#print pos, imputed_seed, seed_frequence
				average_maf += seed_frequence
			#max_linkage_pos_list.extend(list[1][0])
			print average_maf/float(list[0]), float(len([pos for pos in list[1][0] if pos in same_to_B_dict]))/float(len(list[1][0]))
	#max_linkage_pos_list.extend(linkage_size_sorted_list[2][1])
	"""
	print "len(max_linkage_pos_list)", len(max_linkage_pos_list)

	ori_seed_pos_list = data_dict.seed_dict.keys()
	print len(ori_seed_pos_list)

	seed_pos_list = data_dict.seed_dict.keys()

	for pos in max_linkage_pos_list:
		if pos not in ori_seed_pos_list:
			seed = seeds()
			seed.rsID = hifi_dict[pos][0]
			seed.position = int(hifi_dict[pos][1])
			seed.allele_ori = hifi_dict[pos][2]
			seed.allele_new = hifi_dict[pos][2]
			data_dict.seed_dict[int(pos)] = seed
		else:
			#print "pos in linked region and in ori seed: ", pos
			pass

	print len(seed_pos_list)

	new_seed_file_name = "haplotype_link.txt"
	output_revised_seed(new_seed_file_name, data_dict.seed_dict)
	same_to_A_dict, same_to_B_dict = seed_std_compare(new_seed_file_name, chr_name)


def add_seed_by_linkage_100():
	"""
	add bridge if length > 100, otherwise add the longest
	:return:
	"""
	refID_dict = {}

	refID_A_count_dict = {}
	refID_B_count_dict = {}

	refID_list = data_dict.ref_title_info.strip().split()[2:]
	print "refID_list", len(refID_list)
	# print refID_list

	hifi_dict = {}
	hifi_dict = load_raw_data("imputed_" + data_dict.seed_file)[1]
	print "hifi_dict", len(hifi_dict)

	window_info_dict = load_raw_data("window_" + data_dict.seed_file)[1]
	print "window_info_dict", len(window_info_dict)

	for pos, elements in window_info_dict.iteritems():
		if elements[2] != 'NN' and elements[4] != 'NN':
			window_info = elements[9:]
			match_to_A_index_list = range(2, len(refID_list))
			match_to_B_index_list = range(2, len(refID_list))
			hifi_seq_A = ""
			hifi_seq_B = ""
			for window_pos in window_info:

				try:
					window_pos = int(window_pos)
					hifi_seq_A += hifi_dict[window_pos][2]
					hifi_seq_B += hifi_dict[window_pos][3]

					match_to_A_index_list = [index for index in match_to_A_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][2]]
					match_to_B_index_list = [index for index in match_to_B_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][3]]
				except:
					pass
					print window_info
					sys.exit(1)
			match_to_A_refID = [refID_list[index] for index in match_to_A_index_list]
			match_to_B_refID = [refID_list[index] for index in match_to_B_index_list]

			for refID in match_to_A_refID:
				refID_A_count_dict[refID] = 0 if refID not in refID_A_count_dict else (refID_A_count_dict[refID] + 1)
			for refID in match_to_B_refID:
				refID_B_count_dict[refID] = 0 if refID not in refID_B_count_dict else (refID_B_count_dict[refID] + 1)
			match_to_A_seq = ""
			match_to_B_seq = ""

			for window_pos in window_info:
				window_pos = int(window_pos)
				if window_pos != 0:
					if len(match_to_A_index_list) > 0 and window_pos in data_dict.hap_ref_dict:
						match_to_A_seq += data_dict.hap_ref_dict[window_pos][match_to_A_index_list[0]]
					if len(match_to_B_index_list) > 0 and window_pos in data_dict.hap_ref_dict:
						match_to_B_seq += data_dict.hap_ref_dict[window_pos][match_to_B_index_list[0]]

			#refID_dict[pos] = (hifi_seq_A, hifi_seq_B, list_to_line(match_to_A_refID), match_to_A_seq, list_to_line(match_to_B_refID), match_to_B_seq)
			refID_dict[pos] = (
				hifi_seq_A, hifi_seq_B, match_to_A_refID, match_to_A_seq, match_to_B_refID, match_to_B_seq)

	print "**************** group here"
	refID_sorted_list = sort_dict_by_key(refID_dict)

	i = 0
	temp_continus_snp = []
	n_refID_A_common = []
	n_refID_B_common = []

	seed_dict_from_linkage = {}
	linkage_size_dict = {}

	snp_ld_length_dict = {}  # to store the LD block size of each snp

	while i < (len(refID_sorted_list) - 1):
		#last_refID_A_common = []
		pos = refID_sorted_list[i][0]
		current_snp_data = refID_sorted_list[i][1]
		next_snp_data = refID_sorted_list[i + 1][1]
		if len(temp_continus_snp) == 0:
			n_refID_A_common = [id for id in current_snp_data[2] if id in next_snp_data[2]]
			n_refID_B_common = [id for id in current_snp_data[4] if id in next_snp_data[4]]
			# to keep a copy of the last common refID of this group. n_refID_A/B_common may become zero in the last check
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
		else:
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
			n_refID_A_common = [id for id in n_refID_A_common if id in next_snp_data[2]]
			n_refID_B_common = [id for id in n_refID_B_common if id in next_snp_data[4]]

		if len(n_refID_A_common) > 0 and len(
				n_refID_B_common) > 0:  # use and to make sure both A and B have common refID
			temp_continus_snp.append(pos)
			i += 1
			last_refID_A_common = n_refID_A_common
		else:
			if len(temp_continus_snp) > 0:
				ld_block_size = len(temp_continus_snp)
				linkage_size_dict[ld_block_size] = (temp_continus_snp, last_refID_A_common)
				for pos in temp_continus_snp:
					if pos not in snp_ld_length_dict:
						snp_ld_length_dict[pos] = ld_block_size
					else:
						print "duplicate pos in snp_ld_length_dict:", pos

			temp_continus_snp = []
			n_refID_A_common = []
			n_refID_B_common = []
			last_common_refID_A = []
			last_common_refID_B = []
			i += 1
		#print >> log_file, ""

	linkage_size_sorted_list = sort_dict_by_key(linkage_size_dict)
	linkage_size_sorted_list.reverse()
	print "max linkage size:", linkage_size_sorted_list[0][0], linkage_size_sorted_list[0][1][1]
	print "min linkage size:", linkage_size_sorted_list[-1][0]

	max_linkage_pos_list = []
	#max_linkage_pos_list = linkage_size_sorted_list[0][1][0]

	for list in linkage_size_sorted_list:
		if list[0] >= 5:
			#print list[0]
			max_linkage_pos_list.extend(list[1][0])
			print "list extended:", len(list[1][0])
		else:
			max_linkage_pos_list.extend(list[1][0])
			print "list extended and stopped:", len(list[1][0])
			break

	# same_to_A_dict, same_to_B_dict = seed_std_compare("imputed_haplotype.txt", data_dict.chr_name)

	print "len(max_linkage_pos_list)", len(max_linkage_pos_list)

	ori_seed_pos_list = data_dict.seed_dict.keys()
	print len(ori_seed_pos_list)

	seed_pos_list = data_dict.seed_dict.keys()
	#seed_pos_list.extend([x for x in max_linkage_pos_list if x not in ori_seed_pos_list])
	#seed_pos_list.sort()

	for pos in max_linkage_pos_list:
		if pos not in ori_seed_pos_list:
			seed = seeds()
			seed.rsID = hifi_dict[pos][0]
			seed.position = int(hifi_dict[pos][1])
			seed.allele_ori = hifi_dict[pos][2]
			seed.allele_new = hifi_dict[pos][2]
			data_dict.seed_dict[int(pos)] = seed
		else:
			#print "pos in linked region and in ori seed: ", pos
			pass

	print len(seed_pos_list)

	new_seed_file_name = "haplotype.txt"
	output_revised_seed(new_seed_file_name, data_dict.seed_dict)
	# same_to_A_dict, same_to_B_dict = seed_std_compare(new_seed_file_name, chr_name)


def add_seed_by_bridge():
	# use bridge to impute missing snps in middle

	# a dict to store ref_ID for each position
	refID_dict = {}

	refID_A_count_dict = {}
	refID_B_count_dict = {}

	# list to store the ref_ID
	refID_list = data_dict.ref_title_info.strip().split()[2:]
	print "refID_list", len(refID_list)
	# print refID_list

	hifi_dict = {}
	#hifi_dict = load_hifi_result("imputed_" + data_dict.seed_file, hifi_dict)
	hifi_dict = load_raw_data("imputed_" + data_dict.seed_file)[1]
	print "hifi_dict", len(hifi_dict)

	window_info_dict = load_raw_data("window_" + data_dict.seed_file)[1]
	print "window_info_dict", len(window_info_dict)

	for pos, elements in window_info_dict.iteritems():
		if elements[2] != 'NN' and elements[4] != 'NN':
			window_info = elements[9:]
			match_to_A_index_list = range(2, len(refID_list))
			match_to_B_index_list = range(2, len(refID_list))
			hifi_seq_A = ""
			hifi_seq_B = ""
			for window_pos in window_info:
				try:
					window_pos = int(window_pos)
					hifi_seq_A += hifi_dict[window_pos][2]
					hifi_seq_B += hifi_dict[window_pos][3]

					match_to_A_index_list = [index for index in match_to_A_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][2]]
					match_to_B_index_list = [index for index in match_to_B_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][3]]
				except:
					pass
					print window_info
					sys.exit(1)
			match_to_A_refID = [refID_list[index] for index in match_to_A_index_list]
			match_to_B_refID = [refID_list[index] for index in match_to_B_index_list]

			for refID in match_to_A_refID:
				refID_A_count_dict[refID] = 0 if refID not in refID_A_count_dict else (refID_A_count_dict[refID] + 1)
			for refID in match_to_B_refID:
				refID_B_count_dict[refID] = 0 if refID not in refID_B_count_dict else (refID_B_count_dict[refID] + 1)
			match_to_A_seq = ""
			match_to_B_seq = ""

			for window_pos in window_info:
				window_pos = int(window_pos)
				if window_pos != 0:
					if len(match_to_A_index_list) > 0 and window_pos in data_dict.hap_ref_dict:
						match_to_A_seq += data_dict.hap_ref_dict[window_pos][match_to_A_index_list[0]]
					if len(match_to_B_index_list) > 0 and window_pos in data_dict.hap_ref_dict:
						match_to_B_seq += data_dict.hap_ref_dict[window_pos][match_to_B_index_list[0]]

			#refID_dict[pos] = (hifi_seq_A, hifi_seq_B, list_to_line(match_to_A_refID), match_to_A_seq, list_to_line(match_to_B_refID), match_to_B_seq)
			refID_dict[pos] = (
				hifi_seq_A, hifi_seq_B, match_to_A_refID, match_to_A_seq, match_to_B_refID, match_to_B_seq)

	print "**************** group here"
	refID_sorted_list = sort_dict_by_key(refID_dict)

	i = 0
	temp_continus_snp = []
	n_refID_A_common = []
	n_refID_B_common = []

	seed_dict_from_linkage = {}
	linkage_size_dict = {}
	linkage_size_list = []

	snp_ld_length_dict = {}  # to store the LD block size of each snp

	while i < (len(refID_sorted_list) - 1):
		#last_refID_A_common = []
		pos = refID_sorted_list[i][0]
		current_snp_data = refID_sorted_list[i][1]
		next_snp_data = refID_sorted_list[i + 1][1]
		if len(temp_continus_snp) == 0:
			n_refID_A_common = [id for id in current_snp_data[2] if id in next_snp_data[2]]
			n_refID_B_common = [id for id in current_snp_data[4] if id in next_snp_data[4]]
			# to keep a copy of the last common refID of this group. n_refID_A/B_common may become zero in the last check
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
		else:
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
			n_refID_A_common = [id for id in n_refID_A_common if id in next_snp_data[2]]
			n_refID_B_common = [id for id in n_refID_B_common if id in next_snp_data[4]]

		if len(n_refID_A_common) > 0 and len(
				n_refID_B_common) > 0:  # use and to make sure both A and B have common refID
			temp_continus_snp.append(pos)
			i += 1
			last_refID_A_common = n_refID_A_common
		else:
			if len(temp_continus_snp) > 0:
				ld_block_size = len(temp_continus_snp)
				linkage_size_dict[ld_block_size] = (temp_continus_snp, last_refID_A_common)
				linkage_size_list.append((ld_block_size, temp_continus_snp, last_refID_A_common))
				for pos in temp_continus_snp:
					if pos not in snp_ld_length_dict:
						snp_ld_length_dict[pos] = ld_block_size
					else:
						print "duplicate pos in snp_ld_length_dict:", pos

			temp_continus_snp = []
			n_refID_A_common = []
			n_refID_B_common = []
			last_common_refID_A = []
			last_common_refID_B = []
			i += 1
		#print >> log_file, ""

	linkage_size_sorted_list = sort_dict_by_key(linkage_size_dict)
	linkage_size_sorted_list.reverse()
	print "max linkage size:", linkage_size_sorted_list[0][0], linkage_size_sorted_list[0][1][1]
	print "min linkage size:", linkage_size_sorted_list[-1][0]

	print "total size in dict", len(linkage_size_dict)
	print "total size in list", len(linkage_size_list)

	print "total size in dict", linkage_size_dict.keys()
	#print "total size in list", linkage_size_list

	max_linkage_pos_list = []
	max_linkage_pos_dict = {}
	max_linkage_pos_list = linkage_size_sorted_list[0][1][0]

	same_to_A_dict, same_to_B_dict = seed_std_compare("imputed_haplotype.txt", data_dict.chr_name)

	linkage_dict = {}  # used to store the size of each bridge

	for data in linkage_size_list:
		if data[0] in linkage_dict:
			linkage_dict[data[0]] += 1
		else:
			linkage_dict[data[0]] = 1
		if data[0] >= 90:
			print data[0]

	# check the accuracy based on bridge size, output data for every bridge
	for data in linkage_size_list:
		bridge_size = data[0]
		correct_seed = 0
		for pos in data[1]:
			if bridge_seed_equalsTo_std(pos, hifi_dict[pos][2]):
				correct_seed += 1
		bridge_accuracy = round(float(correct_seed) / len(data[1]), 3)

	# check the accuracy based on bridge size, group bridge by size
	bridge_accuracy_dict = {}

	for data in linkage_size_list:
		bridge_size = data[0]
		if bridge_size not in bridge_accuracy_dict:
			bridge_accuracy_dict[bridge_size] = (0, 0)
		correct_seed = 0
		for pos in data[1]:
			if bridge_seed_equalsTo_std(pos, hifi_dict[pos][2]):
				correct_seed += 1
		temp_tuple = (
			bridge_accuracy_dict[bridge_size][0] + correct_seed, bridge_accuracy_dict[bridge_size][1] + len(data[1]))
		bridge_accuracy_dict[bridge_size] = temp_tuple

	bridge_accuracy_sorted_list = sort_dict_by_key(bridge_accuracy_dict)

	refID_bridge_dict = {}

	for ref_id in refID_list:
		for data in linkage_size_list:
			if ref_id in data[2]:
				if ref_id not in refID_bridge_dict:
					refID_bridge_dict[ref_id] = []
				refID_bridge_dict[ref_id].append(data)

			# print the number of bridge on each ref
			#for ref_id in refID_bridge_dict.keys():
			#print ref_id, len(refID_bridge_dict[ref_id])
			#pass

	# to get the refID with longest bridge
	refID_bridge_length_dict = {}

	for ref_id, bridge_on_one_ref in refID_bridge_dict.iteritems():
		if ref_id not in refID_bridge_length_dict:
			refID_bridge_length_dict[ref_id] = 0
		for bridge in bridge_on_one_ref:
			#print len(bridge[2])
			refID_bridge_length_dict[ref_id] += len(bridge[1])

	refID_bridge_length_sorted_list = sort_dict_by_key(refID_bridge_length_dict)
	print "length ", refID_bridge_length_sorted_list[0]

	refID_longest_bridge = refID_bridge_length_sorted_list[0][0]

	bridge_dict = {}  # a matrix to store bridge information

	for data in linkage_size_list:
		for pos in data[1]:
			#print data[1]
			if pos not in bridge_dict:
				bridge_dict[pos] = []
			for ref_id in data[2]:
				bridge_dict[pos].append(ref_id)
			#print pos, len(bridge_dict[pos])
	print "bridge_dict ", len(bridge_dict)

	#print bridge pattern
	bridge_dict_sorted_list = sort_dict_by_key(bridge_dict)
	bridge_pos = window_info_dict.keys()
	bridge_pos.sort()
	#print bridge_pos

	print len(bridge_pos)

	correct_pos = 0
	wrong_pos = 0

	selected_gap_pos_list = []

	seed_in_gap_dict = {}  # to store the pos between two bridge anchors

	refID_list_withrs_pos = data_dict.ref_title_info.strip().split()  # used to locate the index of the allele at certain position with refID

	first_window_size_dict = {}
	second_window_size_dict = {}
	gap_size_dict = {}

	for ref_id in refID_bridge_dict.keys():
		imputation_window_in_ref = refID_bridge_dict[ref_id]
		for i in range(len(imputation_window_in_ref) - 1):
			first_window_size = imputation_window_in_ref[i][0]
			second_window_size = imputation_window_in_ref[i + 1][0]
			last_pos_in_first_window = imputation_window_in_ref[i][1][-1]
			first_pos_in_second_window = imputation_window_in_ref[i + 1][1][0]
			size_of_gap_on_bridge = bridge_pos.index(first_pos_in_second_window) - bridge_pos.index(
				last_pos_in_first_window) - 1
			if first_window_size not in first_window_size_dict:
				first_window_size_dict[first_window_size] = 0
			if second_window_size not in second_window_size_dict:
				second_window_size_dict[second_window_size] = 0
			if size_of_gap_on_bridge not in gap_size_dict:
				gap_size_dict[size_of_gap_on_bridge] = 0

	first_window_size_sorted_list = first_window_size_dict.keys()
	#second_window_size_sorted_list = second_window_size_dict.keys()
	gap_size_sorted_list = gap_size_dict.keys()

	window_c = first_window_size_sorted_list[-10]
	bridge_c = gap_size_sorted_list[20]

	#print "first_window_size_sorted_list", first_window_size_sorted_list, window_c
	#print "second_window_size_sorted_list", second_window_size_sorted_list
	#print "gap_size_sorted_list", gap_size_sorted_list, bridge_c

	for ref_id in refID_bridge_dict.keys():
		imputation_window_in_ref = refID_bridge_dict[ref_id]
		if True:
			#if len(imputation_window_in_ref) >= 2:
			for i in range(len(imputation_window_in_ref) - 1):
				first_window_size = imputation_window_in_ref[i][0]
				second_window_size = imputation_window_in_ref[i + 1][0]
				# check the size of the two bridge anchors
				#if True:
				if first_window_size >= window_c or second_window_size >= window_c:
					last_pos_in_first_window = imputation_window_in_ref[i][1][-1]
					first_pos_in_second_window = imputation_window_in_ref[i + 1][1][0]
					size_of_gap_on_bridge = bridge_pos.index(first_pos_in_second_window) - bridge_pos.index(
						last_pos_in_first_window) - 1
					#if True:
					if size_of_gap_on_bridge <= bridge_c:
						#print "gap is ", size_of_gap_on_bridge
						gap_pos_list = bridge_pos[bridge_pos.index(last_pos_in_first_window) + 1: bridge_pos.index(
							first_pos_in_second_window)]
						"""
						# to print the accuracy of gap, and window size of each bridge anchor.
						"""
						correct_snp_in_gap = 0
						total_snp_in_gap = 0  # do not include snp in seed
						for pos in gap_pos_list:
							if True:
								#if pos not in seed_in_gap_dict:
								refID_index = refID_list_withrs_pos.index(ref_id)
								#print "refID_index", refID_index
								ref_allele = data_dict.hap_ref_dict[pos][refID_index]

								if bridge_seed_equalsTo_std(pos, hifi_dict[pos][2]):
									correct_snp_in_gap += 1
								total_snp_in_gap += 1
						#gap_accuracy = round(float(correct_snp_in_gap)/total_snp_in_gap, 3)

						if True:
							#print "gap_info", first_window_size, second_window_size, len(gap_pos_list), gap_accuracy
							#if first_window_size >= 60 and second_window_size >= 60 and gap_accuracy == 0 and len(gap_pos_list) >= 2:

							for gap_pos in gap_pos_list:
								#selected_gap_pos_list.append(gap_pos)
								max_linkage_pos_list.append(gap_pos)
								max_linkage_pos_dict[gap_pos] = ""

							first_w_pos_list = imputation_window_in_ref[i][1]
							second_w_pos_list = imputation_window_in_ref[i + 1][1]

							#print "gap_info", first_window_size, second_window_size, len(gap_pos_list), gap_accuracy
							#print list_to_line(first_w_pos_list) + ";", list_to_line(gap_pos_list) + ";", list_to_line(second_w_pos_list)
							#and (gap_accuracy == 1 or gap_accuracy == 0):
							#print pos, "****", len(gap_pos_list)
							"""
							print "std_data_A",
							for first_pos in first_w_pos_list:
								if first_pos in data_dict.hap_std_dict:
									print data_dict.hap_std_dict[first_pos][2],
							print ";",
							for gap_pos in gap_pos_list:
								if gap_pos in data_dict.hap_std_dict:
									print data_dict.hap_std_dict[gap_pos][2],
							print ";",
							for second_pos in second_w_pos_list:
								if second_pos in data_dict.hap_std_dict:
									print data_dict.hap_std_dict[second_pos][2],
							print ""

							print "std_data_B",
							for first_pos in first_w_pos_list:
								if first_pos in data_dict.hap_std_dict:
									print data_dict.hap_std_dict[first_pos][3],
							print ";",
							for gap_pos in gap_pos_list:
								if gap_pos in data_dict.hap_std_dict:
									print data_dict.hap_std_dict[gap_pos][3],
							print ";",
							for second_pos in second_w_pos_list:
								if second_pos in data_dict.hap_std_dict:
									print data_dict.hap_std_dict[second_pos][3],
							print ""

							print "lowDepth_data",
							for first_pos in first_w_pos_list:
								print hifi_dict[first_pos][2],
							print ";",
							for gap_pos in gap_pos_list:
								print hifi_dict[gap_pos][2],
							print ";",
							for second_pos in second_w_pos_list:
								print hifi_dict[second_pos][2],
							print ""
							"""

						#print len(gap_pos_list), size_of_gap_on_bridge
						#print "***********", len(selected_gap_pos_list)

						#print gap_pos_list
						#max_linkage_pos_list.extend(gap_pos_list)
						print "max_linkage_pos_dict size", len(max_linkage_pos_dict)



						#for pos in max_linkage_pos_dict.keys():
						for pos in max_linkage_pos_list:
							refID_index = refID_list_withrs_pos.index(ref_id)
							#print "refID_index", refID_index
							ref_allele = data_dict.hap_ref_dict[pos][refID_index]
							#print "allele", ref_allele

							if pos not in seed_in_gap_dict:
								seed = seeds()
								seed.rsID = hifi_dict[pos][0]
								seed.position = int(pos)
								seed.allele_ori = hifi_dict[pos][2]
								seed.allele_new = hifi_dict[pos][2]

								bridge_info_t = bridge_info()
								bridge_info_t.refID = ref_id
								bridge_info_t.ref_allele = ref_allele
								bridge_info_t.first_window_size = first_window_size

								bridge_info_t.second_window_size = second_window_size
								#print "first_second_size", bridge_info_t.first_window_size
								#print "second_window_size", bridge_info_t.second_window_size
								bridge_info_t.first_second_size = first_window_size + second_window_size
								#print "first_second_size", bridge_info_t.first_second_size
								bridge_info_t.gap_size = size_of_gap_on_bridge

								seed.bridge_info.append(bridge_info_t)
								seed_in_gap_dict[pos] = seed
								if bridge_seed_equalsTo_std(pos, hifi_dict[pos][2]):
									correct_pos += 1
								else:
									wrong_pos += 1

							for base in seed_in_gap_dict[pos].allele_dict.keys():
								if base == ref_allele:
									seed_in_gap_dict[pos].allele_dict[base] += 1

	print "new seed accuracy", round(float(correct_pos) / (correct_pos + wrong_pos))

	#print "accuracy perc*****", float(correct_pos)/(correct_pos + wrong_pos)
	# check the accuracy based on bridge anchor and river size
	bridge_anchor_accuracy_dict = {}

	for pos, seed in seed_in_gap_dict.iteritems():
		for bridge_info_t in seed.bridge_info:
			#print "bridge_info_t.first_second_size", bridge_info_t.first_second_size
			if bridge_info_t.first_second_size not in bridge_anchor_accuracy_dict:
				bridge_anchor_accuracy_dict[bridge_info_t.first_second_size] = (0, 0)
			if bridge_seed_equalsTo_std(pos, bridge_info_t.ref_allele):
				bridge_anchor_accuracy_dict[bridge_info_t.first_second_size] = (
					bridge_anchor_accuracy_dict[bridge_info_t.first_second_size][0] + 1,
					bridge_anchor_accuracy_dict[bridge_info_t.first_second_size][1] + 1)
			else:
				bridge_anchor_accuracy_dict[bridge_info_t.first_second_size] = (
					bridge_anchor_accuracy_dict[bridge_info_t.first_second_size][0],
					bridge_anchor_accuracy_dict[bridge_info_t.first_second_size][1] + 1)

	bridge_river_accuracy_sorted_list = sort_dict_by_key(bridge_anchor_accuracy_dict)
	"""
	for data in bridge_river_accuracy_sorted_list:
		#print "bridge anchor size, accuracy", data[0], data[1][0], data[1][1], round(float(data[1][0])/data[1][1], 3)
		#print "bridge anchor", data[0], round(float(data[1][0])/data[1][1], 3)
		pass
	"""
	# check the accuracy based on bridge anchor and river size
	river_accuracy_dict = {}

	for pos, seed in seed_in_gap_dict.iteritems():
		for bridge_info_t in seed.bridge_info:
			if bridge_info_t.gap_size not in river_accuracy_dict:
				river_accuracy_dict[bridge_info_t.gap_size] = (0, 0)
			if bridge_seed_equalsTo_std(pos, bridge_info_t.ref_allele):
				river_accuracy_dict[bridge_info_t.gap_size] = (river_accuracy_dict[bridge_info_t.gap_size][0] + 1,
				                                               river_accuracy_dict[bridge_info_t.gap_size][1] + 1)
			else:
				river_accuracy_dict[bridge_info_t.gap_size] = (river_accuracy_dict[bridge_info_t.gap_size][0],
				                                               river_accuracy_dict[bridge_info_t.gap_size][1] + 1)

	river_accuracy_sorted_list = sort_dict_by_key(river_accuracy_dict)
	for data in river_accuracy_sorted_list:
		#print "river size, accuracy", data[0], data[1][0], data[1][1], round(float(data[1][0])/data[1][1], 3)
		pass

	# extend bridge on one ref
	for bridge in refID_bridge_dict[refID_longest_bridge]:
		#print bridge[1]
		if len(bridge[1]) >= 15:
			max_linkage_pos_list.extend(list(bridge[1]))

	"""
	for list in linkage_size_sorted_list:
		if list[0] >= 100:
			#print list[0]
			max_linkage_pos_list.extend(list[1][0])
			print "list extended:", len(list[1][0])
		else:
			max_linkage_pos_list.extend(list[1][0])
			print "list extended and stopped:", len(list[1][0])
			break
	"""

	print "len(max_linkage_pos_list)", len(max_linkage_pos_list)

	ori_seed_pos_list = data_dict.seed_dict.keys()
	print len(ori_seed_pos_list)

	seed_pos_list = data_dict.seed_dict.keys()
	#seed_pos_list.extend([x for x in max_linkage_pos_list if x not in ori_seed_pos_list])
	#seed_pos_list.sort()

	for pos in max_linkage_pos_list:
		if pos not in ori_seed_pos_list:
			seed = seeds()
			seed.rsID = hifi_dict[pos][0]
			seed.position = int(hifi_dict[pos][1])
			seed.allele_ori = hifi_dict[pos][2]
			seed.allele_new = hifi_dict[pos][2]
			data_dict.seed_dict[int(pos)] = seed
		else:
			#print "pos in linked region and in ori seed: ", pos
			pass

	print len(seed_pos_list)

	new_seed_file_name = "haplotype_river.txt"
	output_revised_seed(new_seed_file_name, data_dict.seed_dict)
	same_to_A_dict, same_to_B_dict = seed_std_compare(new_seed_file_name, chr_name)

	hap_bkup = "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict))
	os.system("cp haplotype_river.txt " + hap_bkup)
	os.system("mv " + hap_bkup + " seed_file")


def bridge_seed_equalsTo_std(pos, bridge_seed):
	if pos in data_dict.hap_std_dict:
		std_A = data_dict.hap_std_dict[pos][2]
		std_B = data_dict.hap_std_dict[pos][3]

		if bridge_seed == std_A:
			return True
		elif std_A == "X" or std_B == "X":
			return True
		elif std_A == "N" or std_B == "N":
			return True
		else:
			if (std_A == "A" and std_B == "T") or (std_A == "C" and std_B == "G") or (
							std_A == "T" and std_B == "A") or (std_A == "G" and std_B == "C"):
				return True
			else:
				return False
	else:
		return True


def add_seed_by_linkage_Jan212014():
	refID_dict = {}

	refID_A_count_dict = {}
	refID_B_count_dict = {}

	refID_list = data_dict.ref_title_info.strip().split()[2:]
	print "refID_list", len(refID_list)
	# print refID_list

	hifi_dict = {}
	#hifi_dict = load_hifi_result("imputed_" + data_dict.seed_file, hifi_dict)
	hifi_dict = load_raw_data("imputed_" + data_dict.seed_file)[1]
	print "hifi_dict", len(hifi_dict)

	window_info_dict = load_raw_data("window_" + data_dict.seed_file)[1]
	print "window_info_dict", len(window_info_dict)

	for pos, elements in window_info_dict.iteritems():
		if elements[2] != 'NN' and elements[4] != 'NN':
			window_info = elements[9:]
			match_to_A_index_list = range(2, len(refID_list))
			match_to_B_index_list = range(2, len(refID_list))
			hifi_seq_A = ""
			hifi_seq_B = ""
			for window_pos in window_info:
				try:
					window_pos = int(window_pos)
					hifi_seq_A += hifi_dict[window_pos][2]
					hifi_seq_B += hifi_dict[window_pos][3]

					match_to_A_index_list = [index for index in match_to_A_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][2]]
					match_to_B_index_list = [index for index in match_to_B_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][3]]
				except:
					print window_info
			match_to_A_refID = [refID_list[index] for index in match_to_A_index_list]
			match_to_B_refID = [refID_list[index] for index in match_to_B_index_list]
			#if pos == 21917858:
			#	print match_to_A_refID, match_to_A_refID
			for refID in match_to_A_refID:
				refID_A_count_dict[refID] = 0 if refID not in refID_A_count_dict else (refID_A_count_dict[refID] + 1)
			for refID in match_to_B_refID:
				refID_B_count_dict[refID] = 0 if refID not in refID_B_count_dict else (refID_B_count_dict[refID] + 1)
			match_to_A_seq = ""
			match_to_B_seq = ""

			for window_pos in window_info:
				window_pos = int(window_pos)
				if window_pos != 0:
					if len(match_to_A_index_list) > 0:
						match_to_A_seq += data_dict.hap_ref_dict[window_pos][match_to_A_index_list[0]]
					if len(match_to_B_index_list) > 0:
						match_to_B_seq += data_dict.hap_ref_dict[window_pos][match_to_B_index_list[0]]

			#refID_dict[pos] = (hifi_seq_A, hifi_seq_B, list_to_line(match_to_A_refID), match_to_A_seq, list_to_line(match_to_B_refID), match_to_B_seq)
			refID_dict[pos] = (
				hifi_seq_A, hifi_seq_B, match_to_A_refID, match_to_A_seq, match_to_B_refID, match_to_B_seq)

	print "**************** group here"
	refID_sorted_list = sort_dict_by_key(refID_dict)

	i = 0
	temp_continus_snp = []
	n_refID_A_common = []
	n_refID_B_common = []

	seed_dict_from_linkage = {}
	linkage_size_dict = {}

	snp_ld_length_dict = {}  # to store the LD block size of each snp

	while i < (len(refID_sorted_list) - 1):
		#last_refID_A_common = []
		pos = refID_sorted_list[i][0]
		current_snp_data = refID_sorted_list[i][1]
		next_snp_data = refID_sorted_list[i + 1][1]
		if len(temp_continus_snp) == 0:
			n_refID_A_common = [id for id in current_snp_data[2] if id in next_snp_data[2]]
			n_refID_B_common = [id for id in current_snp_data[4] if id in next_snp_data[4]]
			# to keep a copy of the last common refID of this group. n_refID_A/B_common may become zero in the last check
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
		else:
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
			n_refID_A_common = [id for id in n_refID_A_common if id in next_snp_data[2]]
			n_refID_B_common = [id for id in n_refID_B_common if id in next_snp_data[4]]

		if len(n_refID_A_common) > 0 and len(
				n_refID_B_common) > 0:  # use and to make sure both A and B have common refID
			temp_continus_snp.append(pos)
			i += 1
			last_refID_A_common = n_refID_A_common
		else:
			if len(temp_continus_snp) > 0:
				ld_block_size = len(temp_continus_snp)
				linkage_size_dict[ld_block_size] = (temp_continus_snp, last_refID_A_common)
				for pos in temp_continus_snp:
					if pos not in snp_ld_length_dict:
						snp_ld_length_dict[pos] = ld_block_size
					else:
						print "duplicate pos in snp_ld_length_dict:", pos

			temp_continus_snp = []
			n_refID_A_common = []
			n_refID_B_common = []
			last_common_refID_A = []
			last_common_refID_B = []
			i += 1
		#print >> log_file, ""

	linkage_size_sorted_list = sort_dict_by_key(linkage_size_dict)
	linkage_size_sorted_list.reverse()
	print "max linkage size:", linkage_size_sorted_list[0][0], linkage_size_sorted_list[0][1][1]
	print "min linkage size:", linkage_size_sorted_list[-1][0]

	max_linkage_pos_list = []
	max_linkage_pos_list = linkage_size_sorted_list[0][1][0]

	same_to_A_dict, same_to_B_dict = seed_std_compare("imputed_haplotype.txt", data_dict.chr_name)

	for list in linkage_size_sorted_list:
		if list[0] > 0:
			print list[0],  #, list[1][0]
			average_maf = 0
			for pos in list[1][0]:
				imputed_seed, seed_frequence = data_dict.hap_ref_allele_frequence_dict[pos][0] \
					if hifi_dict[pos][2] == data_dict.hap_ref_allele_frequence_dict[pos][0][0] else \
					data_dict.hap_ref_allele_frequence_dict[pos][1]
				#print pos, imputed_seed, seed_frequence
				average_maf += float(seed_frequence)
			#max_linkage_pos_list.extend(list[1][0])
			#print average_maf/float(list[0]), float(len([pos for pos in list[1][0] if pos in same_to_B_dict]))/float(len(list[1][0]))
	#max_linkage_pos_list.extend(linkage_size_sorted_list[2][1])
	print len(max_linkage_pos_list)

	ori_seed_pos_list = data_dict.seed_dict.keys()
	print len(ori_seed_pos_list)

	seed_pos_list = data_dict.seed_dict.keys()
	#seed_pos_list.extend([x for x in max_linkage_pos_list if x not in ori_seed_pos_list])
	#seed_pos_list.sort()

	for pos in max_linkage_pos_list:
		if pos not in ori_seed_pos_list:
			seed = seeds()
			seed.rsID = hifi_dict[pos][0]
			seed.position = int(hifi_dict[pos][1])
			seed.allele_ori = hifi_dict[pos][2]
			seed.allele_new = hifi_dict[pos][2]
			data_dict.seed_dict[int(pos)] = seed
		else:
			#print "pos in linked region and in ori seed: ", pos
			pass

	print len(seed_pos_list)

	new_seed_file_name = "haplotype.txt"
	output_revised_seed(new_seed_file_name, data_dict.seed_dict)
	same_to_A_dict, same_to_B_dict = seed_std_compare(new_seed_file_name, chr_name)

	os.system("cp haplotype.txt " + "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict)))

	# output information
	#rs#, pos, maf, std_allele,
	error_analysis_dict = {}

	with open("imputed_info.info", 'w') as log_file:
		print >> log_file, "rs#", "pos", "allele_freq_A", "allele_freq_B", "minor_allele_frequence", "hap_std", "imputed_allele", \
			"imputed_allele_frequence", "accuracy", \
			"hap_freq_A", "hap_freq_B", "window_size", \
			"window_distance", "num_x", "x/seed_in_window", "ld_block_size"

		window_info_sorted_list = sort_dict_by_key(window_info_dict)
		print "window_info_dict", len(window_info_dict)

		for data in window_info_sorted_list:
			pos = int(data[0])

			elements = data[1]
			if elements[2] != 'NN' and elements[4] != 'NN':
				rs_id = elements[0]
				#x_symbol = elements[2]
				window_info = elements[9:]
				allele_frequence_data = data_dict.hap_ref_allele_frequence_dict[pos]
				allele_freq_A = allele_frequence_data[0][0] + ":" + str(allele_frequence_data[0][1])
				allele_freq_B = allele_frequence_data[1][0] + ":" + str(allele_frequence_data[1][1])
				minor_allele_frequence = str(allele_frequence_data[0][1]) if allele_frequence_data[0][1] <= \
				                                                             allele_frequence_data[1][1] else str(
					allele_frequence_data[1][1])

				hap_std = data_dict.hap_std_dict[pos][2] if pos in data_dict.hap_std_dict else "NA"
				imputed_allele = hifi_dict[pos][2]
				imputed_allele_frequence = str(allele_frequence_data[0][1]) if imputed_allele == \
				                                                               allele_frequence_data[0][0] else str(
					allele_frequence_data[1][1])

				std_impute_compare = -1
				if hap_std == imputed_allele:
					std_impute_compare = 1
				elif hap_std == "X" or hap_std == "NA":
					std_impute_compare = 2
				else:
					std_impute_compare = 0

				hap_freq_A = len(refID_dict[pos][2])
				hap_freq_B = len(refID_dict[pos][4])
				#refID_dict[pos] = (hifi_seq_A, hifi_seq_B, match_to_A_refID, match_to_A_seq, match_to_B_refID, match_to_B_seq)

				window_size = len(window_info)
				window_distance = int(window_info[-1]) - int(window_info[0])
				ld_block_size = snp_ld_length_dict[pos] if pos in snp_ld_length_dict else 0
				num_x = len([position for position in window_info if
				             int(position) in window_info_dict and window_info_dict[int(position)][2] == "XX"])
				x_over_window_size = format(float(num_x) / window_size, "0.3f")
				#	print data
				print >> log_file, rs_id, pos, allele_freq_A, allele_freq_B, minor_allele_frequence, \
					hap_std, imputed_allele, imputed_allele_frequence, \
					std_impute_compare, hap_freq_A, hap_freq_B, \
					window_size, window_distance, num_x, x_over_window_size, ld_block_size

				target_data = window_distance
				#target_data = x_over_window_size
				if std_impute_compare == 0 or std_impute_compare == 1:
					if target_data not in error_analysis_dict:
						error_analysis_dict[target_data] = [0, 0]
					else:
						if std_impute_compare == 0:
							error_analysis_dict[target_data][0] = error_analysis_dict[target_data][0] + 1
						elif std_impute_compare == 1:
							error_analysis_dict[target_data][1] = error_analysis_dict[target_data][1] + 1
			else:
				#print pos
				pass

	error_analysis_sorted_list = sort_dict_by_key(error_analysis_dict)
	for data in error_analysis_sorted_list:
		data_freq = data[0]
		temp_list = data[1]
		if temp_list[0] != 0 or temp_list[1] != 0:
			print data_freq, temp_list[0], temp_list[1], format(
				float(temp_list[0]) / float(temp_list[0] + temp_list[1]), "0.3f")
		#pass

	# window_distance
	total_num_if_different_freq = len(error_analysis_dict)
	inteval = total_num_if_different_freq / 80
	print "inteval", inteval
	i = 0
	j = 1
	error_data = 0
	correct_data = 0
	#lower_bound = total_num_if_different_freq[0][0]
	for data in error_analysis_sorted_list:
		data_freq = data[0]
		temp_list = data[1]
		if temp_list[0] != 0 or temp_list[1] != 0:
			if i <= j * inteval:
				error_data += temp_list[0]
				correct_data += temp_list[1]
				i += 1
			else:
				#print data_freq, error_data, correct_data, format(float(error_data)/float(error_data + correct_data), "0.3f")
				j += 1
				error_data = 0
				correct_data = 0


def add_seed_by_linkage():
	refID_dict = {}

	refID_A_count_dict = {}
	refID_B_count_dict = {}

	refID_list = data_dict.ref_title_info.strip().split()[2:]
	print "refID_list", len(refID_list)
	# print refID_list

	hifi_dict = {}
	#hifi_dict = load_hifi_result("imputed_" + data_dict.seed_file, hifi_dict)
	hifi_dict = load_raw_data("imputed_" + data_dict.seed_file)[1]
	print "hifi_dict", len(hifi_dict)

	window_info_dict = load_raw_data("window_" + data_dict.seed_file)[1]
	print "window_info_dict", len(window_info_dict)

	for pos, elements in window_info_dict.iteritems():
		if elements[2] != 'NN' and elements[4] != 'NN':
			window_info = elements[9:]
			match_to_A_index_list = range(2, len(refID_list))
			match_to_B_index_list = range(2, len(refID_list))
			hifi_seq_A = ""
			hifi_seq_B = ""
			for window_pos in window_info:
				try:
					window_pos = int(window_pos)
					hifi_seq_A += hifi_dict[window_pos][2]
					hifi_seq_B += hifi_dict[window_pos][3]

					match_to_A_index_list = [index for index in match_to_A_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][2]]
					match_to_B_index_list = [index for index in match_to_B_index_list if
					                         data_dict.hap_ref_dict[window_pos][index] == hifi_dict[window_pos][3]]
				except:
					print window_info
			match_to_A_refID = [refID_list[index] for index in match_to_A_index_list]
			match_to_B_refID = [refID_list[index] for index in match_to_B_index_list]
			#if pos == 21917858:
			#	print match_to_A_refID, match_to_A_refID
			for refID in match_to_A_refID:
				refID_A_count_dict[refID] = 0 if refID not in refID_A_count_dict else (refID_A_count_dict[refID] + 1)
			for refID in match_to_B_refID:
				refID_B_count_dict[refID] = 0 if refID not in refID_B_count_dict else (refID_B_count_dict[refID] + 1)
			match_to_A_seq = ""
			match_to_B_seq = ""

			for window_pos in window_info:
				window_pos = int(window_pos)
				if window_pos != 0:
					if len(match_to_A_index_list) > 0:
						match_to_A_seq += data_dict.hap_ref_dict[window_pos][match_to_A_index_list[0]]
					if len(match_to_B_index_list) > 0:
						match_to_B_seq += data_dict.hap_ref_dict[window_pos][match_to_B_index_list[0]]

			#refID_dict[pos] = (hifi_seq_A, hifi_seq_B, list_to_line(match_to_A_refID), match_to_A_seq, list_to_line(match_to_B_refID), match_to_B_seq)
			refID_dict[pos] = (
				hifi_seq_A, hifi_seq_B, match_to_A_refID, match_to_A_seq, match_to_B_refID, match_to_B_seq)

			"""
			if pos == 21917858:
				print refID_dict[pos]
				for window_pos in window_info:
					window_pos = int(window_pos)
					print hifi_dict[window_pos][2], hifi_dict[window_pos][3],
					print list_to_line(data_dict.hap_ref_dict[window_pos]),
					print " "
			"""
	print "**************** group here"
	refID_sorted_list = sort_dict_by_key(refID_dict)

	i = 0
	temp_continus_snp = []
	n_refID_A_common = []
	n_refID_B_common = []

	seed_dict_from_linkage = {}
	linkage_size_dict = {}

	snp_ld_length_dict = {}  # to store the LD block size of each snp

	while i < (len(refID_sorted_list) - 1):
		#last_refID_A_common = []
		pos = refID_sorted_list[i][0]
		current_snp_data = refID_sorted_list[i][1]
		next_snp_data = refID_sorted_list[i + 1][1]
		if len(temp_continus_snp) == 0:
			n_refID_A_common = [id for id in current_snp_data[2] if id in next_snp_data[2]]
			n_refID_B_common = [id for id in current_snp_data[4] if id in next_snp_data[4]]
			# to keep a copy of the last common refID of this group. n_refID_A/B_common may become zero in the last check
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
		else:
			last_common_refID_A = n_refID_A_common
			last_common_refID_B = n_refID_B_common
			n_refID_A_common = [id for id in n_refID_A_common if id in next_snp_data[2]]
			n_refID_B_common = [id for id in n_refID_B_common if id in next_snp_data[4]]

		if len(n_refID_A_common) > 0 and len(
				n_refID_B_common) > 0:  # use and to make sure both A and B have common refID
			temp_continus_snp.append(pos)
			i += 1
			last_refID_A_common = n_refID_A_common
		else:
			if len(temp_continus_snp) > 0:
				ld_block_size = len(temp_continus_snp)
				linkage_size_dict[ld_block_size] = (temp_continus_snp, last_refID_A_common)
				for pos in temp_continus_snp:
					if pos not in snp_ld_length_dict:
						snp_ld_length_dict[pos] = ld_block_size
					else:
						print "duplicate pos in snp_ld_length_dict:", pos

			temp_continus_snp = []
			n_refID_A_common = []
			n_refID_B_common = []
			last_common_refID_A = []
			last_common_refID_B = []
			i += 1
		#print >> log_file, ""

	linkage_size_sorted_list = sort_dict_by_key(linkage_size_dict)
	linkage_size_sorted_list.reverse()
	print "max linkage size:", linkage_size_sorted_list[0][0], linkage_size_sorted_list[0][1][1]
	print "min linkage size:", linkage_size_sorted_list[-1][0]

	max_linkage_pos_list = []

	for list in linkage_size_sorted_list:
		if list[0] >= 100:
			#print list[0]
			max_linkage_pos_list.extend(list[1][0])
		elif list[0] >= 50:
			print list[0]
			ld_block_score = 8.5
			for pos in list[1][0]:
				if pos in window_info_dict:
					final_score = 0
					elements = window_info_dict[pos]
					if elements[2] != 'NN' and elements[4] != 'NN':
						window_info = elements[9:]
						hap_freq_A = len(refID_dict[pos][2])
						hap_freq_B = len(refID_dict[pos][4])
						window_size = len(window_info)

						hap_freq_A_score = 9 if hap_freq_A >= 120 else 7
						hap_freq_B_score = 9 if hap_freq_B >= 120 else 7
						if window_size >= 200:
							window_size_score = 6
						elif window_size >= 100:
							window_size_score = 9
						else:
							window_size_score = 8
					final_score = ld_block_score * 0.55 + hap_freq_A_score * 0.15 + hap_freq_B_score * 0.15 + window_size_score * 0.15
					print final_score
					#print type(final_score)
					if final_score >= 8.0:
						#print final_score
						max_linkage_pos_list.append(pos)
		else:
			ld_block_score = 7.0
			for pos in list[1][0]:
				if pos in window_info_dict:
					final_score = 0
					elements = window_info_dict[pos]
					if elements[2] != 'NN' and elements[4] != 'NN':
						window_info = elements[9:]
						hap_freq_A = len(refID_dict[pos][2])
						hap_freq_B = len(refID_dict[pos][4])
						window_size = len(window_info)

						hap_freq_A_score = 9 if hap_freq_A >= 120 else 7
						hap_freq_B_score = 9 if hap_freq_B >= 120 else 7
						if window_size >= 200:
							window_size_score = 6
						elif window_size >= 100:
							window_size_score = 9
						else:
							window_size_score = 8
					final_score = ld_block_score * 0.55 + hap_freq_A_score * 0.15 + hap_freq_B_score * 0.15 + window_size_score * 0.15
					print final_score
					#print type(final_score)
					if final_score >= 8.0:
						#print final_score
						max_linkage_pos_list.append(pos)

	ori_seed_pos_list = data_dict.seed_dict.keys()
	print len(ori_seed_pos_list)

	seed_pos_list = data_dict.seed_dict.keys()
	new_seed_dict = {}
	for pos in max_linkage_pos_list:
		if pos not in ori_seed_pos_list:
			seed = seeds()
			seed.rsID = hifi_dict[pos][0]
			seed.position = int(hifi_dict[pos][1])
			seed.allele_ori = hifi_dict[pos][2]
			seed.allele_new = hifi_dict[pos][2]
			data_dict.seed_dict[int(pos)] = seed
		#new_seed_dict[int(pos)] = seed
		else:
			#print "pos in linked region and in ori seed: ", pos
			pass

	print len(seed_pos_list)

	new_seed_file_name = "haplotype_ext.txt"
	output_revised_seed(new_seed_file_name, data_dict.seed_dict)
	#output_revised_seed(new_seed_file_name, new_seed_dict)
	same_to_A_dict, same_to_B_dict = seed_std_compare(new_seed_file_name, chr_name)


# os.system("cp haplotype.txt " + "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict)))

def seed_recover_extract_ref_cluster():
	recovered_seed_dict = {}

	seed_ref_difference_dict = dict_substract(hap_ref_dict, seed_dict)

	hifi_dict, qscore_dict = load_hap_qscore(number_of_subfile)

	hifi_sorted_list = sort_dict_by_key(hifi_dict)
	hifi_sorted_pos_list = []
	j = 0
	while j < len(hifi_sorted_list):
		hifi_sorted_pos_list.append(hifi_sorted_list[j][0])
		j += 1
	print "hifi_sorted_pos_list", len(hifi_sorted_pos_list)

	added_by_cluster = 0
	for num, snp_list in ref_cluster_dict.iteritems():
		for cluster_dict in snp_list:
			for pos_1, ref_1 in cluster_dict.iteritems():
				pos_1 = int(pos_1)
				if pos_1 in hifi_dict and pos_1 not in seed_dict:

					seed = hifi_dict[pos_1]
					max_base = keywithmaxval(seed.allele_dict)
					max_value = seed.allele_dict[max_base]
					seed.allele_new_percentage = float(max_value) / float(number_of_subfile)
					if seed.allele_new_percentage >= 0.70:
						try:
							seed = seeds()
							seed.rsID = ref_1[0].strip()
							seed.position = int(pos_1)
							seed.allele_ori = max_base
							seed.allele_new = max_base
							seed_dict[int(pos_1)] = seed
							added_by_cluster += 1
						except:
							print max_base, "not in ", list_to_line(ref_1)

	output_filename = "haplotype_cluster.txt"
	output_revised_seed(output_filename, seed_dict)
	seed_std_compare(output_filename, chr_name)

	file_name = output_filename
	hifi_run(file_name, chr_name)
	hifiAccuCheck("imputed_" + file_name, chr_name)
	print "added_by_cluster", added_by_cluster


def seed_expand_geno():
	global random_geno_seed_dict

	seed_geno_difference_dict = dict_substract(data_dict.geno_hetero_dict, data_dict.seed_hetero_dict)
	seed_geno_difference_sorted_list = sort_dict_by_key(seed_geno_difference_dict)

	random_geno_seed_dict = {}
	geno_seed_selected_number = 5000
	for i in range(0, geno_seed_selected_number):
		random_index = random.randrange(0, len(seed_geno_difference_dict))
		position = seed_geno_difference_sorted_list[random_index][0]
		while position in random_geno_seed_dict:
			random_index = random.randrange(0, len(seed_geno_difference_dict))
			position = seed_geno_difference_sorted_list[random_index][0]
		random_geno_seed_dict[position] = seed_geno_difference_dict[position]

	geno_seed_added_each_file = geno_seed_selected_number / number_of_subfile
	print "geno_seed_added_each_file: ", geno_seed_added_each_file

	random_geno_seed_sorted_list = sort_dict_by_key(random_geno_seed_dict)

	for geno_allele_index in (0, 1):  # 0 for A, 1 for B
		#process_list = []
		for file_number in range(number_of_subfile):
			hap_subfile_name = data_dict.seed_file_name + "_" + str(geno_allele_index) + "_" + str(file_number) + ".txt"
			output_subfile = open(currentPath + hap_subfile_name, "w")
			print >> output_subfile, data_dict.seed_title_info
			seed_hetero_dict_bkup = data_dict.seed_hetero_dict.copy()
			#print "seed_hetero_sorted_list original", len(data_dict.seed_hetero_dict)

			for i in range(int(geno_seed_added_each_file)):
				try:
					forward_index = i * number_of_subfile + file_number
					position = random_geno_seed_sorted_list[forward_index][0]
					elements = random_geno_seed_sorted_list[forward_index][1]
					if position not in data_dict.seed_hetero_dict:
						seed = seeds()
						seed.rsID = elements[0].strip()
						seed.position = int(elements[1].strip())
						""" take A from geno first """
						seed.allele_ori = elements[2].strip()[geno_allele_index]
						#print position, seed.position, seed.allele_ori, elements[2].strip()
						data_dict.seed_hetero_dict[position] = seed
				except:
					#print "error", forward_index
					pass
			"""last position in seed may not be the same with others this needs to be checked"""
			#print "seed_hetero_sorted_list new", len(data_dict.seed_hetero_dict)

			sub_seed_dict = dict_add(data_dict.seed_hetero_dict, data_dict.seed_homo_dict)
			sub_seed_list = sort_dict_by_key(sub_seed_dict)

			#print "new sub_seed_list size", len(sub_seed_list)
			for seed in sub_seed_list:
				line = seed[1].rsID + "\t" + str(seed[1].position) + "\t" + seed[1].allele_ori
				print >> output_subfile, line
			output_subfile.close()
			data_dict.seed_hetero_dict = seed_hetero_dict_bkup.copy()
			hifi_process(file_number, number_of_subfile, hap_subfile_name)
			"""
			maf_step = float(random.randrange(10, 40))/(100.0)
			print "maf_step is: ", maf_step
			hifi = program_path + "hifi_fu_ref " + hap_subfile_name
			hifi_process = subprocess.Popen(hifi, shell=True)
			hifi_process.wait()
			"""
			"""
			maf_step = float(random.randrange(10, 40))/(100.0)
			geno_file_name="genotype.txt" 
			ref_file_name="refHaplos.txt"
			process = multiprocessing.Process(target=hifi_mlp, args=(hap_subfile_name, geno_file_name, ref_file_name, maf_step))
			process.start()
			process_list.append(process)
	
		for process in process_list:
		    process.join()
			"""


def seed_geno_extract():
	revised_seed_dict = {}

	for geno_allele_index in (0, 1):  # 0 for A, 1 for B

		hifi_dict = {}
		qscore_dict = {}
		# this one has geno_allele_index, cannot use the load_hifi_qscore function
		for file_number in range(number_of_subfile):
			hifi_subfile_name = "imputed_haplotype_" + str(geno_allele_index) + "_" + str(file_number) + ".txt"
			if os.path.exists(hifi_subfile_name):
				hifi_dict = load_hifi_result(hifi_subfile_name, hifi_dict)

			qscore_subfile_name = "qscore_haplotype_" + str(geno_allele_index) + "_" + str(file_number) + ".txt"
			if os.path.exists(hifi_subfile_name):
				qscore_dict = load_qscore_result(qscore_subfile_name, qscore_dict)
		#print "qscore_dict size: ", len(qscore_dict)

		#random_geno_seed_sorted_list = sort_dict_by_key(random_geno_seed_dict)
		hifi_dict_sorted_list = sort_dict_by_key(hifi_dict)

		added_seed = 0
		for snp in hifi_dict_sorted_list:
			#for snp in random_geno_seed_sorted_list:
			position = snp[0]
			if position not in revised_seed_dict and position in hifi_dict and position in data_dict.geno_dict:
				seed = hifi_dict[position]
				max_base = keywithmaxval(seed.allele_dict)
				max_value = seed.allele_dict[max_base]
				seed.allele_new_percentage = float(max_value) / float(number_of_subfile)
				if seed.allele_new_percentage >= 0.90 and max_base == data_dict.geno_dict[position][2][
					geno_allele_index]:
					if position in qscore_dict and qscore_dict[position] / float(number_of_subfile) >= 0.70:
						#if position in hap_std_dict and max_base == hap_std_dict[position][2]:
						#	if True:
						added_seed += 1
						seed.allele_new = max_base
						revised_seed_dict[position] = seed
					#print position, max_base, geno_dict[position][2][geno_allele_index]
					#print "added seed by", geno_allele_index, added_seed

	print "added seed total number", len(revised_seed_dict)

	verify_expanded_seed_by_cluster(revised_seed_dict)
	print "added seed total number after cluster check", len(revised_seed_dict)

	revised_seed_dict = dict_add(data_dict.seed_dict, revised_seed_dict)
	#print "new seed total number", len(revised_seed_dict)
	output_revised_seed("haplotype_expanded.txt", revised_seed_dict)
	return revised_seed_dict


def error_seed_distance(seed_dict, same_to_B_dict):
	seed_sorted_list = sort_dict_by_key(seed_dict)
	seed_pos_list = []
	for j in range(len(seed_sorted_list)):
		seed_pos_list.append(seed_sorted_list[j][0])
	for position, seed in same_to_B_dict.iteritems():
		if position in seed_pos_list:
			index = seed_pos_list.index(position)
			print "added_pos: ", position,
			print seed_sorted_list[index - 1][0], position - seed_sorted_list[index - 1][0],
			print seed_sorted_list[index + 1][0], seed_sorted_list[index + 1][0] - position


def seed_verify():
	experiment_seed_file_name = "haplotype_expanded.txt_1822_100"  #original seed from experiment, hifh accuracy
	middle_expand_seed_file = "haplotype.txt_4247"  #"haplotype.txt_4942_0"	#"haplotype.txt_4574"		# seed expanded from seed_ori
	"""current seed file contains more expanded seed, which will be verified, haplotype.txt_8403_28"""
	# generate geno and ref file from ori seed

	experiment_seed_dict = load_seed_data(experiment_seed_file_name)[1]
	experiment_seed_homo_dict, experiment_seed_hetero_dict = group_seed(experiment_seed_dict, data_dict.geno_dict)
	print "experiment_seed_hetero_dict ori", len(experiment_seed_hetero_dict)
	print "seed_dict ori", len(data_dict.seed_dict)
	seed_remove_experiment_dict = dict_substract(data_dict.seed_dict, experiment_seed_hetero_dict)
	seed_remove_experiment_file_name = "haplotype_expr.txt"
	print "seed_remove_experiment_dict", len(seed_remove_experiment_dict)
	output_revised_seed(seed_remove_experiment_file_name, seed_remove_experiment_dict)

	hifi = program_path + "hifi_fu_ref " + seed_remove_experiment_file_name
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()

	hifi_dict = {}
	hifi_subfile_name = "imputed_haplotype_expr.txt"
	hifi_dict = load_hifi_result(hifi_subfile_name, hifi_dict)

	exp_seed_new_dict = {index: value for index, value in hifi_dict.iteritems() if index in experiment_seed_hetero_dict}
	print "exp_seed_new_dict", len(exp_seed_new_dict)

	changed_seed_dict = {index: value for index, value in exp_seed_new_dict.iteritems() if \
	                     experiment_seed_hetero_dict[index].allele_new != exp_seed_new_dict[index].allele_new \
	                     and data_dict.hap_std_dict[index][2] != "X"}

	print "changed_seed_dict", len(changed_seed_dict)

	for pos, seed in changed_seed_dict.iteritems():
		if pos in data_dict.hap_std_dict:
			print pos, exp_seed_new_dict[pos].allele_new, experiment_seed_hetero_dict[pos].allele_new, \
				data_dict.hap_std_dict[pos][2], data_dict.hap_std_dict[pos][3]

	#print "seed_geno_difference_dict", len(seed_geno_difference_dict)

	#seed_geno_difference_sorted_list = sort_dict_by_key(seed_geno_difference_dict)
	middle_expand_seed_dict = load_seed_data(middle_expand_seed_file)[1]
	middle_seed_homo_dict, middle_seed_hetero_dict = group_seed(middle_expand_seed_dict, data_dict.geno_dict)
	middle_exp_difference_bfchange_dict = dict_substract(data_dict.seed_dict, middle_seed_hetero_dict)
	#i = 0
	print "middle_expand_seed_dict", len(middle_expand_seed_dict)

	changed_seed_list = [x for x in changed_seed_dict.iteritems()]

	for i in range(0, len(changed_seed_list)):
		pos = changed_seed_list[i][0]
		print pos
		middle_expand_seed_dict[pos].allele_new = exp_seed_new_dict[pos].allele_new

	seed_changed_file_name = "haplotype_sch.txt"
	print "middle_expand_seed_dict", len(middle_expand_seed_dict)
	output_revised_seed(seed_changed_file_name, middle_expand_seed_dict)

	hifi = program_path + "hifi_fu_ref " + seed_changed_file_name
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()

	hifi_dict = {}
	hifi_subfile_name = "imputed_haplotype_sch.txt"
	hifi_dict = load_hifi_result(hifi_subfile_name, hifi_dict)

	exp_seed_new_dict_aftchange = {index: value for index, value in hifi_dict.iteritems() if
	                               index in middle_exp_difference_bfchange_dict}
	print "exp_seed_new_dict_aftchange", len(exp_seed_new_dict_aftchange)

	exp_changed_seed_dict = {index: value for index, value in exp_seed_new_dict_aftchange.iteritems() if \
	                         middle_exp_difference_bfchange_dict[index].allele_new != exp_seed_new_dict_aftchange[
		                         index].allele_new \
	                         and data_dict.hap_std_dict[index][2] != "X"}

	print "exp_changed_seed_dict", len(exp_changed_seed_dict)

	for pos, seed in exp_changed_seed_dict.iteritems():
		if pos in data_dict.hap_std_dict:
			print pos, exp_seed_new_dict_aftchange[pos].allele_new, middle_exp_difference_bfchange_dict[pos].allele_new, \
				data_dict.hap_std_dict[pos][2], data_dict.hap_std_dict[pos][3]
			del data_dict.seed_dict[pos]
	output_revised_seed("revised_seed.txt", data_dict.seed_dict)

	seed_std_compare("revised_seed.txt", data_dict.chr_name)

def seed_verify_reverse():

	experiment_seed_file_name = "haplotype_expanded.txt_1822_100"
	middle_expand_seed_file = "haplotype.txt_4247"

	experiment_seed_dict = load_seed_data(experiment_seed_file_name)[1]
	experiment_seed_homo_dict, experiment_seed_hetero_dict = group_seed(experiment_seed_dict, data_dict.geno_dict)
	print "experiment_seed_hetero_dict ori", len(experiment_seed_hetero_dict)

	middle_expand_seed_dict = load_seed_data(middle_expand_seed_file)[1]
	#middle_seed_homo_dict, middle_seed_hetero_dict = group_seed(middle_expand_seed_dict, data_dict.geno_dict)
	print "middle_expand_seed_dict ori", len(middle_expand_seed_dict)

	changed_seed_list = [x for x in experiment_seed_hetero_dict.iteritems()]

	for i in range(0, 200):
		pos = changed_seed_list[i][0]
		#print pos
		del middle_expand_seed_dict[pos]

	seed_remove_experiment_file_name = "haplotype_expr.txt"
	print "seed_remove_experiment_dict", len(middle_expand_seed_dict)
	output_revised_seed(seed_remove_experiment_file_name, middle_expand_seed_dict)

	hifi = program_path + "hifi_fu_ref " + seed_remove_experiment_file_name
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()

	hifi_dict = {}
	hifi_subfile_name = "imputed_haplotype_expr.txt"
	hifi_dict = load_hifi_result(hifi_subfile_name, hifi_dict)

	exp_seed_new_dict = {index: value for index, value in hifi_dict.iteritems() if index in experiment_seed_hetero_dict}
	print "exp_seed_new_dict", len(exp_seed_new_dict)

	changed_seed_dict = {index: value for index, value in exp_seed_new_dict.iteritems() if \
	                     experiment_seed_hetero_dict[index].allele_new != exp_seed_new_dict[index].allele_new \
	                     and data_dict.hap_std_dict[index][2] != "X"}

	print "changed_seed_dict", len(changed_seed_dict)

	for pos, seed in changed_seed_dict.iteritems():
		if pos in data_dict.hap_std_dict:
			print pos, exp_seed_new_dict[pos].allele_new, experiment_seed_hetero_dict[pos].allele_new, \
				data_dict.hap_std_dict[pos][2], data_dict.hap_std_dict[pos][3]

	#print "seed_geno_difference_dict", len(seed_geno_difference_dict)
	"""
	#seed_geno_difference_sorted_list = sort_dict_by_key(seed_geno_difference_dict)
	middle_expand_seed_dict = load_seed_data(middle_expand_seed_file)[1]
	middle_seed_homo_dict, middle_seed_hetero_dict = group_seed(middle_expand_seed_dict, data_dict.geno_dict)
	middle_exp_difference_bfchange_dict = dict_substract(data_dict.seed_dict, middle_seed_hetero_dict)

	print "middle_expand_seed_dict", len(middle_expand_seed_dict)
	
	changed_seed_list = [x for x in changed_seed_dict.iteritems()]
	
	for i in range(0, len(changed_seed_list)):
		pos = changed_seed_list[i][0]
		print pos
		middle_expand_seed_dict[pos].allele_new = exp_seed_new_dict[pos].allele_new

	seed_changed_file_name = "haplotype_sch.txt"
	print "middle_expand_seed_dict", len(middle_expand_seed_dict)
	output_revised_seed(seed_changed_file_name, middle_expand_seed_dict)
		
	hifi = program_path + "hifi_fu_ref " + seed_changed_file_name
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()
	
	hifi_dict = {}
	hifi_subfile_name = "imputed_haplotype_sch.txt"
	hifi_dict = load_hifi_result(hifi_subfile_name, hifi_dict)
	
	exp_seed_new_dict_aftchange = {index:value for index, value in hifi_dict.iteritems() if index in middle_exp_difference_bfchange_dict}
	print "exp_seed_new_dict_aftchange", len(exp_seed_new_dict_aftchange)

	exp_changed_seed_dict = {index:value for index, value in exp_seed_new_dict_aftchange.iteritems() if \
						middle_exp_difference_bfchange_dict[index].allele_new != exp_seed_new_dict_aftchange[index].allele_new \
						and data_dict.hap_std_dict[index][2] != "X"}
	
	print "exp_changed_seed_dict", len(exp_changed_seed_dict)
	
	for pos, seed in exp_changed_seed_dict.iteritems():
		if pos in data_dict.hap_std_dict:
			print pos, exp_seed_new_dict_aftchange[pos].allele_new, middle_exp_difference_bfchange_dict[pos].allele_new, data_dict.hap_std_dict[pos][2], data_dict.hap_std_dict[pos][3]
			del data_dict.seed_dict[pos]
	output_revised_seed("revised_seed.txt", data_dict.seed_dict)	
	
	seed_std_compare("revised_seed.txt", data_dict.chr_name)
	"""

def get_seed_group():
	# find the continuous group with size larger than group_size_threshold
	#seed_sorted_list = sort_dict_by_key(data_dict.seed_dict)
	seed_group_list = []
	group_size_threshold = 20
	temp_group_list = []

	seed_sorted_pos_list = [pos for pos, seed in data_dict.seed_dict.iteritems()]
	seed_sorted_pos_list.sort()
	ref_sorted_pos_list = [pos for pos, seed in data_dict.hap_ref_dict.iteritems()]
	ref_sorted_pos_list.sort()

	for i, pos in enumerate(ref_sorted_pos_list):
		if pos in seed_sorted_pos_list:
			temp_group_list.append(pos)
			if (i + 1) < len(ref_sorted_pos_list) and ref_sorted_pos_list[i + 1] not in seed_sorted_pos_list:
				if len(temp_group_list) >= group_size_threshold:
					seed_group_list.append(temp_group_list)
				temp_group_list = []

	print "seed_group_list", len(seed_group_list)

	group_list = []
	for seed_group in seed_group_list:
		print "*********"
		print seed_group[0], seed_group[-1], len(seed_group)
		for pos in seed_group:
			group_list.append(pos)

	print "total group number in seed ", len(seed_group_list)

	group_seed_dict = {index: value for index, value in data_dict.seed_dict.iteritems() if index in group_list}
	print "total snp number of group in seed ", len(group_seed_dict)

	output_revised_seed("group_seed.txt", group_seed_dict)
	seed_std_compare("group_seed.txt", data_dict.chr_name)

def get_group_list(start_pos, end_pos, dict):
	"""get the homo and hetero snps between start and end points. this step can be improved"""
	temp_dict = {pos: snp for pos, snp in dict.iteritems() if pos >= start_pos and pos <= end_pos}
	temp_sorted_list = sort_dict_by_key(temp_dict)
	return temp_sorted_list

def get_seed_group_index():
	# this is different with get_seed_group(). divide the total snp (std, ref) into small partition
	# and check if the known seed percentage in the partition is higher than the existing_seed_percentage
	"""hetero seed distribution"""
	seed_group_list = []
	temp_group_list = []

	ref_hetero_dict = dict_substract(data_dict.hap_ref_dict, data_dict.geno_homo_dict)
	ref_hetero_total_number = len(ref_hetero_dict)
	print ref_hetero_total_number

	seed_sorted_pos_list = [pos for pos, seed in data_dict.seed_hetero_dict.iteritems()]
	seed_sorted_pos_list.sort()
	ref_sorted_pos_list = [pos for pos, seed in ref_hetero_dict.iteritems()]
	ref_sorted_pos_list.sort()

	ref_hetero_in_each_partition = data_dict.seed_group_window_size
	partition_number = int(math.ceil(float(ref_hetero_total_number) / ref_hetero_in_each_partition))
	print "partition_number: ", partition_number
	ref_hetero_in_last_subfile = int(math.fmod(ref_hetero_total_number, ref_hetero_in_each_partition))
	ref_hetero_in_last_subfile = ref_hetero_in_each_partition if ref_hetero_in_last_subfile == 0 else ref_hetero_in_last_subfile
	print "ref_hetero_in_last_subfile: ", ref_hetero_in_last_subfile

	for p_number in range(partition_number):
		for i in range(int(ref_hetero_in_each_partition)):
			index = p_number * ref_hetero_in_each_partition + i
			if index < ref_hetero_total_number:
				pos = ref_sorted_pos_list[index]
				if pos in seed_sorted_pos_list:
					temp_group_list.append(pos)
		seed_group_list.append(temp_group_list)
		temp_group_list = []

	print "seed_group_list", len(seed_group_list)

	geno_subfile_name = "genotype_group.txt"
	geno_subfile = open(currentPath + geno_subfile_name, "w")
	print >> geno_subfile, data_dict.geno_title_info

	ref_subfile_name = "ref_group.txt"
	ref_subfile = open(currentPath + ref_subfile_name, "w")
	print >> ref_subfile, data_dict.ref_title_info

	output_file = open(currentPath + "seed_group_index_percentage.txt", "w")
	print >> output_file, "pos", "number", "percentge"
	#print len(seed_group_list)

	possible_seed = 0

	ref_pos_percentage_dict = {}
	for i, seed_group in enumerate(seed_group_list):
		group_size = len(seed_group)
		seed_percentage_in_ref = 0
		if i == len(seed_group_list) - 1:
			seed_percentage_in_ref = float(group_size) / float(ref_hetero_in_last_subfile)
		#print group_size, seed_percentage_in_ref
		else:
			seed_percentage_in_ref = float(group_size) / float(ref_hetero_in_each_partition)
		#print group_size, seed_percentage_in_ref
		for j in range(int(ref_hetero_in_each_partition)):
			index = i * ref_hetero_in_each_partition + j
			if index < ref_hetero_total_number:
				pos = ref_sorted_pos_list[index]
				#ref_pos_percentage_dict[pos] = seed_percentage_in_ref
				if seed_percentage_in_ref >= data_dict.existing_seed_percentage:
					print >> output_file, pos, group_size, seed_percentage_in_ref
				""" #hetero seed only """
				if seed_percentage_in_ref >= data_dict.existing_seed_percentage:
					if pos in data_dict.seed_dict:
						ref_pos_percentage_dict[pos] = data_dict.seed_dict[pos]
					if pos in data_dict.geno_dict:
						print >> geno_subfile, list_to_line(data_dict.geno_dict[pos])
					if pos in data_dict.hap_ref_dict:
						print >> ref_subfile, list_to_line(data_dict.hap_ref_dict[pos])
		"""
		# add homo seed
		if seed_percentage_in_ref >= data_dict.existing_seed_percentage:
			possible_seed += 1	
			start_pos = seed_group[0]
			end_pos = seed_group[-1]
			#print start_pos, end_pos
			seed_list = get_group_list(start_pos, end_pos, data_dict.seed_dict)
			for snp in seed_list:
				pos = snp[0]
				ref_pos_percentage_dict[pos] = data_dict.seed_dict[pos]
			geno_list = get_group_list(start_pos, end_pos, data_dict.geno_dict)
			for snp in geno_list:
				pos = snp[0]
				print >> geno_subfile, list_to_line(data_dict.geno_dict[pos])
			ref_list = get_group_list(start_pos, end_pos, data_dict.hap_ref_dict)
			for snp in ref_list:
				pos = snp[0]
				print >> ref_subfile, list_to_line(data_dict.hap_ref_dict[pos])
		"""

	last_pos = max(data_dict.seed_dict.keys())
	if last_pos not in ref_pos_percentage_dict:
		ref_pos_percentage_dict[last_pos] = data_dict.seed_dict[last_pos]
		print >> geno_subfile, list_to_line(data_dict.geno_dict[last_pos])
		print >> ref_subfile, list_to_line(data_dict.hap_ref_dict[last_pos])

	output_file.close()
	geno_subfile.close()
	ref_subfile.close()

	print "possible_seed", possible_seed

	hap_subfile = "haplotype_group.txt"
	output_revised_seed(hap_subfile, ref_pos_percentage_dict)
	seed_std_compare(hap_subfile, data_dict.chr_name)

	hifi = program_path + "hifi_fu_revise " + hap_subfile + " " + geno_subfile_name + " " + ref_subfile_name
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()
	hifiAccuCheck("imputed_" + hap_subfile, chr_name)
	#seed_std_compare("imputed_"+hap_subfile, chr_name)

	extended_seed_dict = load_seed_data("imputed_" + hap_subfile)[1]
	data_dict.seed_dict = dict_add(data_dict.seed_dict, extended_seed_dict)
	file_name = "haplotype_expanded.txt"
	output_revised_seed(file_name, data_dict.seed_dict)
	seed_std_compare("haplotype.txt", chr_name)
	same_to_A_dict, same_to_B_dict = seed_std_compare(file_name, chr_name)
	os.system("cp haplotype_expanded.txt haplotype.txt")
	os.system("cp haplotype.txt " + "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict)))

def calculate_seed_group_accuracy():
	# divide the pos in ref into partitions. Each partition has its snp pos list,
	# percentage of seed, accuracy of hifi imputed snps.

	ref_group_list = []
	temp_group_list = []

	seed_sorted_pos_list = data_dict.seed_dict.keys()
	seed_sorted_pos_list.sort()
	ref_sorted_pos_list = data_dict.hap_ref_dict.keys()
	ref_sorted_pos_list.sort()
	ref_total_number = len(ref_sorted_pos_list)
	print "ref_total_number", ref_total_number

	ref_in_each_partition = data_dict.seed_group_window_size
	partition_number = int(math.ceil(float(ref_total_number) / ref_in_each_partition))
	print "partition_number: ", partition_number
	ref_hetero_in_last_subfile = int(math.fmod(ref_total_number, ref_in_each_partition))
	ref_hetero_in_last_subfile = ref_in_each_partition if ref_hetero_in_last_subfile == 0 else ref_hetero_in_last_subfile
	print "ref_hetero_in_last_subfile: ", ref_hetero_in_last_subfile

	for p_number in range(partition_number):
		for i in range(int(ref_in_each_partition)):
			index = p_number * ref_in_each_partition + i
			if index < ref_total_number:
				pos = ref_sorted_pos_list[index]
				#if pos in data_dict.seed_dict.keys():
				temp_group_list.append(pos)
		ref_group_list.append(temp_group_list)
		temp_group_list = []

	print "ref_group_list", len(ref_group_list)
	print ref_group_list[0]
	print ref_group_list[1]
	print ref_group_list[-2]
	print ref_group_list[-1]

	hifi_dict = {}
	hifi_dict = load_hifi_result("imputed_haplotype.txt", hifi_dict)

	seed_file_name = "haplotype.txt"
	same_to_A_dict, same_to_B_dict = seed_std_compare(seed_file_name, data_dict.chr_name)
	hifi_file_name = "imputed_haplotype.txt"
	same_to_AB_dict, AT_GC_dict = hifiAccuCheck(hifi_file_name, data_dict.chr_name)

	cutoff = 0.9
	number = 0
	a_perc = 0

	seed_percentage_dict = {}
	for ref_group in ref_group_list:
		#print ref_group
		if len(ref_group) > 0:
			snp_from_seed = [pos for pos in ref_group if pos in data_dict.seed_dict]
			seed_percentage = round(float(len(snp_from_seed)) / len(ref_group), 2)

			snp_from_hifi = [pos for pos in ref_group if
			                 pos in hifi_dict and pos not in AT_GC_dict and pos not in same_to_B_dict]
			if len(snp_from_hifi) > 0:
				snp_in_A = [pos for pos in snp_from_hifi if pos in same_to_AB_dict]
				snp_in_A_percentage = round(float(len(snp_in_A)) / len(snp_from_hifi), 4)
				#snp_in_B = [pos for pos in snp_from_hifi if pos in same_to_B_dict]
				if seed_percentage > 0.9 and seed_percentage <= 1:
					#print seed_percentage, len(snp_in_A), len(snp_from_hifi), snp_in_A_percentage
					number += 1
					a_perc += snp_in_A_percentage
				if seed_percentage in seed_percentage_dict:
					seed_percentage_dict[seed_percentage].append(snp_in_A_percentage)
				else:
					seed_percentage_dict[seed_percentage] = []
					seed_percentage_dict[seed_percentage].append(snp_in_A_percentage)
	print number, a_perc, a_perc / number
	seed_percentage_ordered_list = sort_dict_by_key(seed_percentage_dict)
	for data in seed_percentage_ordered_list:
		perc, sna_A_list = data[0], data[1]
		print perc, sum(sna_A_list) / len(sna_A_list)
	print ">=0.9", (sum(seed_percentage_dict[0.9]) + sum(seed_percentage_dict[1.0])) / (
		len(seed_percentage_dict[0.9]) + len(seed_percentage_dict[1.0]))


def combine_hifi_seed(input_prefix, ori_seed_file):
	"""
	combine hifi seed from mutiple run based on consistence
	"""

	revised_seed_dict = {}
	hifi_dict = {}
	#cycle_number = data_dict.cycle_number
	cycle_number = 2
	ori_seed_dict = load_seed_data(ori_seed_file)[1]

	for i in range(1, cycle_number + 1):
		input_subfile_name = input_prefix + str(i) + ".txt"
		hifi_dict = load_hifi_result(input_subfile_name, hifi_dict)
	print "hifi_dict seed total number", len(hifi_dict)

	for position, snp in hifi_dict.iteritems():
		if True:
			#if position not in ori_seed_dict:
			seed = hifi_dict[position]
			max_base = keywithmaxval(seed.allele_dict)
			max_value = seed.allele_dict[max_base]
			seed.allele_new_percentage = float(max_value) / cycle_number
			if seed.allele_new_percentage >= 0.9:
				seed.allele_new = max_base
				revised_seed_dict[position] = seed
			else:
				#print seed.allele_new_percentage 
				pass

	print "consistent seed total number", len(revised_seed_dict)
	#revised_seed_dict = dict_add(ori_seed_dict, revised_seed_dict)
	print "new seed total number", len(revised_seed_dict)
	output_filename = "haplotype_combined.txt"
	output_revised_seed(output_filename, revised_seed_dict)

	same_to_A_dict, same_to_B_dict = seed_std_compare(output_filename, data_dict.chr_name)
	os.system("cp " + output_filename + "_haplotype.txt")
	os.system("cp " + output_filename + "_haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict)))


def combine_seed(input_file_1, input_file_2):
	"""
	combine seed from various tests based on consistence
	"""

	revised_seed_dict = {}
	hifi_dict = {}

	hifi_dict = load_hifi_result(input_file_1, hifi_dict)
	print "hifi_dict seed total number 1", len(hifi_dict)
	hifi_dict = load_hifi_result(input_file_2, hifi_dict)
	print "hifi_dict seed total number 2", len(hifi_dict)

	for position, snp in hifi_dict.iteritems():
		if True:
			seed = hifi_dict[position]
			max_base = keywithmaxval(seed.allele_dict)
			max_value = seed.allele_dict[max_base]
			seed.allele_new_percentage = float(max_value) / 2
			if seed.allele_new_percentage >= 0.9:
				seed.allele_new = max_base
				revised_seed_dict[position] = seed
			else:
				#print seed.allele_new_percentage
				pass

	print "consistent seed total number", len(revised_seed_dict)
	#revised_seed_dict = dict_add(ori_seed_dict, revised_seed_dict)
	print "new seed total number", len(revised_seed_dict)
	output_filename = "haplotype_combined.txt"
	output_revised_seed(output_filename, revised_seed_dict)

	same_to_A_dict, same_to_B_dict = seed_std_compare(output_filename, data_dict.chr_name)
	os.system("cp " + output_filename + "_haplotype.txt")
	os.system("cp " + output_filename + "_haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict)))


def verify_expanded_seed_by_cluster(revised_seed_dict):

	for pos, revised_seed in revised_seed_dict.iteritems():
		keep_this_seed = True
		if pos in data_dict.cluster_pos_dict:
			revised_seed_allele = revised_seed.allele_new
			revised_seed_allele_maf_num = get_maf(pos, revised_seed_allele)
			cluster_list = data_dict.ref_cluster_dict[revised_seed_allele_maf_num]
			for cluster_dict in cluster_list:
				if keep_this_seed and pos in cluster_dict:
					for pos_2, ref in cluster_dict.iteritems():
						if pos_2 in data_dict.seed_dict:
							seed_allele = data_dict.seed_dict[pos_2].allele_new
							seed_allele_maf_num = get_maf(pos_2, seed_allele)
							if revised_seed_allele_maf_num != seed_allele_maf_num:
								keep_this_seed = False
								print "revised: ", pos, revised_seed_allele, revised_seed_allele_maf_num
								print "revised: ", pos_2, seed_allele, seed_allele_maf_num
		if not keep_this_seed:
			del data_dict.revised_seed_dict[pos]

def get_maf(position, allele):
	alleles = hap_ref_dict[position]
	alleles = alleles[2:]
	#print "total allele: ", len(alleles)
	unique_alleles = set(alleles)
	n_alleles = len(unique_alleles)
	if n_alleles == 0 or n_alleles > 2:
		print "error in: ", position
		sys.exit(1)
	else:
		#print position
		for ref_allele in unique_alleles:
			#print allele, alleles.count(allele), format((float(alleles.count(allele))/float(len(alleles))), "0.3f")
			if allele == ref_allele:
				maf_num = alleles.count(allele)
	if maf_num != 0:
		return maf_num
	else:
		print "revised allele not in ref_allele please check", position
		sys.exit(1)

def update_cluster():
	print "****update_cluster running*****"
	revised_seed_dict = {}

	ref_file_name = "refHaplos.txt"
	global ref_cluster_dict
	maf_upper_bound = 0.5
	maf_lower_bound = 0.01
	ref_cluster_dict = get_cluster(ref_file_name, maf_upper_bound, maf_lower_bound)

	added_by_cluster = 0
	for num, snp_list in ref_cluster_dict.iteritems():
		for cluster_dict in snp_list:
			#cluster_sorted_list = sort_dict_by_key(cluster_dict)

			for pos_1, ref_1 in cluster_dict.iteritems():
				if int(pos_1) in data_dict.seed_dict:
					seed_1 = data_dict.seed_dict[int(pos_1)].allele_ori
					try:
						index_1 = ref_1.index(seed_1)
						for pos_2, ref_2 in cluster_dict.iteritems():
							if pos_2 != pos_1 and int(pos_2) not in data_dict.seed_dict:
								cluster_seed = ref_2[index_1]
								#print pos_2, cluster_seed, " not in seed"
								seed = seeds()
								seed.rsID = ref_2[0].strip()
								seed.position = int(pos_2)
								seed.allele_ori = cluster_seed
								seed.allele_new = cluster_seed
								revised_seed_dict[int(pos_2)] = seed
								added_by_cluster += 1
					except:
						#print seed_1, "not in ", list_to_line(ref_1)
						pass

	print "new seed total number", len(revised_seed_dict)
	output_filename = "haplotype_added_by_cluster.txt"
	output_revised_seed(output_filename, revised_seed_dict)
	#if len(revised_seed_dict) > 0:
	#	seed_std_compare(output_filename, chr_name)

	revised_seed_dict = dict_add(data_dict.seed_dict, revised_seed_dict)
	output_filename = "haplotype_expanded.txt"
	output_revised_seed(output_filename, data_dict.seed_dict)

	same_to_A_dict, same_to_B_dict = seed_std_compare(output_filename, data_dict.chr_name)
	os.system(
		"cp haplotype_expanded.txt " + "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict)))

	print "added_by_cluster", added_by_cluster


def clean_up():
	os.system(
		"rm refHaplos_?.txt qscore_haplotype_* pos_deledted_* imputed_haplotype_* genotype_* haplotype_?_* haplotype_?.txt")

def load_data_dicts(seed_file, chr_name):
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

	global ref_cycle_number
	ref_cycle_number = 5

	#global maf_step
	maf_step = 0.1

	ref_file_name = "refHaplos.txt"
	ref_title_info, hap_ref_dict = load_raw_data(ref_file_name, raw_data_format)

	#global ref_cluster_dict
	#ref_cluster_dict = get_cluster(ref_file_name)

	print "seed_file is :", seed_file
	seed_file_name = seed_file[:seed_file.find('.')].strip()
	seed_title_info, seed_dict = load_seed_data(seed_file)

	#print "seed_title_info", seed_title_info
	print "total_seed_number: ", len(seed_dict)

	genotype_file = file_path + "genotype_NA10847_" + chr_name + ".txt"  # for all
	geno_title_info, geno_dict = load_raw_data(genotype_file, raw_data_format)
	print "total_geno_number: ", len(geno_dict)

	hap_std_file = file_path + "ASW_" + chr_name + "_child_hap_refed.txt"
	hap_std_dict = load_hap_std(hap_std_file)

	seed_homo_dict, seed_hetero_dict = group_seed(seed_dict, geno_dict)

	print "seed_homo_dict", len(seed_homo_dict)
	print "seed_hetero_dict", len(seed_hetero_dict)

	geno_homo_dict, geno_hetero_dict = group_seed(geno_dict, geno_dict)

	print "geno_homo_dict", len(geno_homo_dict)
	print "geno_hetero_dict", len(geno_hetero_dict)


def output_genohomo(filename):
	revised_seed_dict = data_dict.seed_homo_dict
	seed_new_file = open(currentPath + filename, "w")

	a = 0
	for pos, seed in data_dict.seed_hetero_dict.iteritems():
		if a < 500:
			#revised_seed_dict[pos] = [seed.rsID, seed.position, seed.allele_new]
			revised_seed_dict[pos] = seed
			a += 1

	print >> seed_new_file, data_dict.seed_title_info
	revised_seed_sorted_list = sort_dict_by_key(revised_seed_dict)  # need to sort the snps by position
	for snp in revised_seed_sorted_list:
		seed = snp[1]
		#line = seed[0] + "\t" + str(seed[1]) + "\t" + seed[2][0]
		line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
		print >> seed_new_file, line
	seed_new_file.close()

def ref_bridge_corssover(seed_file, chr_name, mode):
	# used to generate the 100-cycle figure, use LD, does not include ref remove
	sub_cycle = data_dict.cycle_number
	haplotype_file = "haplotype.txt"

	os.system("cp haplotype.txt haplotype_ori.txt")

	record_file = open(data_dict.record_file_name, "w")
	print >> record_file, "id", "total hetero", "A", "B", "B%"
	i = 1
	while i <= 5:

		os.system("cp haplotype.txt haplotype_ori.txt")

		# ref
		remPercent = float(random.randrange(20, 50)) / 100
		print "remPercent", remPercent
		haplotype_file = "haplotype.txt"
		refMerger(haplotype_file, chr_name, remPercent)
		hifi_run(haplotype_file, data_dict.chr_name)
		hifi_result = "imputed_haplotype.txt"
		refMerger(hifi_result, chr_name, 0)
		seed_std_compare(haplotype_file, data_dict.chr_name)
		ref_temp_file = "haplotype_ref_" + str(i) + ".txt"
		os.system("cp haplotype.txt " + ref_temp_file)

		os.system("cp haplotype_ori.txt haplotype.txt")
		for j in range(3):
			hifi_run(haplotype_file, data_dict.chr_name)

			mode = "100"
			print "########### linkage expand #########", j
			seed_correction(seed_file, chr_name, mode)

			same_to_A_dict, same_to_B_dict = seed_std_compare(haplotype_file, data_dict.chr_name)
			seed_same_to_A = len(same_to_A_dict)
			seed_same_to_B = len(same_to_B_dict)
			B_in_hetero = round((float(seed_same_to_B) / float(seed_same_to_A + seed_same_to_B)) * 100, 2)
			print >> record_file, "LD", i, seed_same_to_A + seed_same_to_B, seed_same_to_A, seed_same_to_B, B_in_hetero

			hap_bkup = "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict))
			os.system("cp haplotype.txt " + hap_bkup)
			os.system("mv " + hap_bkup + " seed_file")

		bridge_temp_file = "haplotype_bridge_" + str(i) + ".txt"
		os.system("cp haplotype.txt " + bridge_temp_file)

		combine_seed(ref_temp_file, bridge_temp_file)

		os.system("cp haplotype_combined.txt haplotype.txt")

		i += 1

def final_1(seed_file, chr_name, mode):
	sub_cycle = data_dict.cycle_number
	haplotype_file = "haplotype.txt"

	os.system("cp haplotype.txt haplotype_ori.txt")

	record_file = open(data_dict.record_file_name, "w")
	print >> record_file, "id", "total hetero", "A", "B", "B%"

	remPercent = 0
	print "remPercent", remPercent
	haplotype_file = "haplotype.txt"
	refMerger(haplotype_file, chr_name, remPercent)
	hifi_run(haplotype_file, data_dict.chr_name)

	mode = "100"
	print "########### 100 expand #########",
	seed_correction(seed_file, chr_name, mode)

	same_to_A_dict, same_to_B_dict = seed_std_compare(haplotype_file, data_dict.chr_name)
	seed_same_to_A = len(same_to_A_dict)
	seed_same_to_B = len(same_to_B_dict)
	B_in_hetero = round((float(seed_same_to_B)/float(seed_same_to_A + seed_same_to_B))*100, 2)
	print >> record_file, "100", seed_same_to_A+seed_same_to_B, seed_same_to_A, seed_same_to_B, B_in_hetero

	hap_bkup = "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict))
	os.system("cp haplotype.txt " + hap_bkup)
	os.system("mv " + hap_bkup + " seed_file")

	remPercent = 0
	refMerger(haplotype_file, chr_name, remPercent)

	# mode = "remove"
	# print "########### error remove #########"
	#
	# seed_correction(seed_file, chr_name, mode)
	# os.system("cp haplotype_error_removed.txt " + haplotype_file)
	#
	# same_to_A_dict, same_to_B_dict = seed_std_compare(haplotype_file, data_dict.chr_name)
	# seed_same_to_A = len(same_to_A_dict)
	# seed_same_to_B = len(same_to_B_dict)
	#
	# hap_bkup = "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict))
	# os.system("cp haplotype.txt " + hap_bkup)
	#
	# os.system("mv " + hap_bkup + " seed_file")
	#
	# B_in_hetero = round((float(seed_same_to_B)/float(seed_same_to_A + seed_same_to_B))*100, 2)
	# print >> record_file, "REV", seed_same_to_A + seed_same_to_B, seed_same_to_A, seed_same_to_B, B_in_hetero
	#
	# print "########### recover expand cycle #########"
	# mode = "recover"
	# seed_correction(seed_file, chr_name, mode)
	# os.system("cp haplotype_expanded.txt haplotype.txt")
	# same_to_A_dict, same_to_B_dict = seed_std_compare("haplotype.txt", data_dict.chr_name)
	# os.system("cp haplotype.txt " + "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict)))


	i = 1
	while i <= 1:

		for j in range(3):
			remPercent = 0 if j == 0 else float(random.randrange(20, 40)) / 100
			"""
			print "remPercent", remPercent
			haplotype_file = "haplotype.txt"
			refMerger(haplotype_file, chr_name, remPercent)
			hifi_run(haplotype_file, data_dict.chr_name)

			mode = "linkage"

			print "########### linkage expand #########", i
			seed_correction(seed_file, chr_name, mode)

			same_to_A_dict, same_to_B_dict = seed_std_compare(haplotype_file, data_dict.chr_name)
			seed_same_to_A = len(same_to_A_dict)
			seed_same_to_B = len(same_to_B_dict)
			B_in_hetero = round((float(seed_same_to_B)/float(seed_same_to_A + seed_same_to_B))*100, 2)
			print >> record_file, "LD", i, seed_same_to_A+seed_same_to_B, seed_same_to_A, seed_same_to_B, B_in_hetero

			hap_bkup = "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict))
			os.system("cp haplotype.txt " + hap_bkup)
			os.system("mv " + hap_bkup + " seed_file")

			i += 1
			"""
		"""
		remPercent = 0
		refMerger(haplotype_file, chr_name, remPercent)

		mode = "remove"
		print "########### error remove #########", i

		seed_correction(seed_file, chr_name, mode)
		os.system("cp haplotype_error_removed.txt " + haplotype_file)

		same_to_A_dict, same_to_B_dict = seed_std_compare(haplotype_file, data_dict.chr_name)
		seed_same_to_A = len(same_to_A_dict)
		seed_same_to_B = len(same_to_B_dict)

		hap_bkup = "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict))
		os.system("cp haplotype.txt " + hap_bkup)

		os.system("mv " + hap_bkup + " seed_file")

		B_in_hetero = round((float(seed_same_to_B)/float(seed_same_to_A + seed_same_to_B))*100, 2)
		print >> record_file, "REV", i, seed_same_to_A + seed_same_to_B, seed_same_to_A, seed_same_to_B, B_in_hetero
		"""
		i += 1

		remove_single_refID()

	record_file.close()


#def seed_correction(seed_file, chr_name, mode):
def seed_correction(seed_file, geno_file, ref_file, chr_name, mode):

	#global number_of_subfile

	global data_dict
	data_dict = data_dicts()
	data_dict.seed_file = seed_file
	data_dict.geno_file_name = geno_file
	data_dict.ref_file_name = ref_file
	data_dict.chr_name = chr_name

	data_dict.seed_file = seed_file
	data_dict.load_data_dicts()

	if mode == "overall":
		pass
	elif mode == "test":
		print data_dict.seed_file
		data_dict.geno_file_name

	elif mode == "crossover":
		ref_bridge_corssover(seed_file, chr_name, mode)

	elif mode == "bridge":
		add_seed_by_bridge()

	elif mode == "final":
		final_1(seed_file, chr_name, mode)

	elif mode == "refid":
		get_refID()
	elif mode == "fid":
		remove_single_refID()
	elif mode == "linkage":
		#add_seed_by_linkage()
		add_seed_by_linkage_longestLD()
	elif mode == "100":
		add_seed_by_linkage_100()
	elif mode == "ref_f":
		data_dict.load_ref_allele_frequence()
	elif mode == "std_perc":
		# calculate the percentage of hetero snp in std haplotype
		std_homo, std_hetero = group_seed(data_dict.hap_std_dict, data_dict.geno_dict)
		print len(std_homo), len(std_hetero)
		print "hetero_perc", round(float(len(std_hetero)) / len(data_dict.hap_std_dict), 2)
	elif mode == "remove":
		seed_std_compare(seed_file, chr_name)
		seed_error_remove()
		seed_error_remove_extract()
		same_to_A_dict, same_to_B_dict = seed_std_compare("haplotype_error_removed.txt", chr_name)
		os.system("cp haplotype_error_removed.txt " + "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(
			len(same_to_B_dict)))
	elif mode == "reme":
		seed_error_remove_extract()
		seed_std_compare("haplotype_error_removed.txt", chr_name)
	elif mode == "recover":
		revised_seed_dict = load_seed_data("haplotype_error_removed.txt")[1]
		seed_recover(data_dict.seed_dict, revised_seed_dict)
		seed_recover_extract(data_dict.seed_hetero_dict, revised_seed_dict)
		same_to_A_dict, same_to_B_dict = seed_std_compare("haplotype_expanded.txt", chr_name)
		os.system(
			"cp haplotype_expanded.txt " + "haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict)))

	elif mode == "recovere":
		revised_seed_dict = load_seed_data("haplotype_ori.txt")[1]
		seed_recover_extract(data_dict.seed_hetero_dict, revised_seed_dict)
		seed_std_compare("haplotype_expanded.txt", chr_name)
	elif mode == "expand":
		file_name = "haplotype.txt"
		#hifi_run(file_name, chr_name)
		hifiAccuCheck("imputed_" + file_name, chr_name)
		seed_std_compare("haplotype.txt", chr_name)
		seed_expand_qs(seed_file)
		seed_std_compare("haplotype_expanded.txt", chr_name)

	elif mode == "ref":
		file_name = "haplotype_expanded.txt"
		seed_expand_ref_hetero()
		seed_recover_extract_ref()
		seed_std_compare(file_name, chr_name)
	elif mode == "refe":
		#seed_expand_ref()
		file_name = "haplotype_expanded.txt"
		seed_recover_extract_ref()
		same_to_A_dict, same_to_B_dict = seed_std_compare(file_name, chr_name)
		os.system("cp " + file_name + " haplotype.txt_" + str(len(same_to_A_dict)) + "_" + str(len(same_to_B_dict)))

	elif mode == "rcluster":
		#seed_expand_ref()	
		seed_recover_extract_ref_cluster()

	elif mode == "genohomo":
		seed_file = "geno_homo_seed.txt"
		output_genohomo(seed_file)
		seed_std_compare(seed_file, data_dict.chr_name)
		hifi = program_path + "hifi_fu_revise " + seed_file
		hifi_process = subprocess.Popen(hifi, shell=True)
		hifi_process.wait()
		hifiAccuCheck("imputed_" + seed_file, chr_name)

	elif mode == "combine_hifi":
		#input_prefix = "haplotype_error_removed.txt_"
		#input_prefix = "seed_from_hifi.txt_"
		input_prefix = "imputed_haplotype_"

		combine_hifi_seed(input_prefix, "haplotype.txt")
	elif mode == "combine":
		input_file_1 = "haplotype_ref.txt"
		input_file_2 = "haplotype.txt_4220_866"
		combine_seed(input_file_1, input_file_2)

	elif mode == "seed_v":
		seed_verify_reverse()
	elif mode == "hifi":
		hap_subfile_name = seed_file
		seed_std_compare(seed_file, chr_name)
		geno_subfile_name = "genotype.txt"
		ref_subfile_name = "refHaplos.txt"
		hifi = program_path + "hifi_fu_revise " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " " + str(
			0.1)
		hifi_process = subprocess.Popen(hifi, shell=True)
		hifi_process.wait()
		hifiAccuCheck("imputed_" + seed_file, chr_name)

	elif mode == "sgroup":
		#get_seed_group()
		get_seed_group_index()
	elif mode == "city_accu":
		calculate_seed_group_accuracy()
	else:
		clean_up()


def get_args():
	desc = " Impute haplotype from low-depth single chromosome sequencing data"
	usage = "lowdepth -i seedFile -g genotypeFile -r referenceFile -c chr"
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--seed", type="string", dest="seedFile", help="Input Seed File name", default="null")
	parser.add_option("-g", "--geno", type="string", dest="genotype", help="Input geno file name", default="null")
	parser.add_option("-r", "--ref", type="string", dest="reference", help="Input ref file name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode", help="split or combine", default="null")
	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.seedFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options


if __name__ == '__main__':
	options = get_args()
	seed_file = options.seedFile
	geno_file = options.genotype
	ref_file = options.reference
	chr_name = options.chrName
	mode = options.mode

	start_time = time.time()
	seed_correction(seed_file, geno_file, ref_file, chr_name, mode)
	print "run time is: ", round((time.time() - start_time), 3), "s"
