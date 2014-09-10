#!/usr/bin/python

# calculate the seed percentage and hifi error in a centain range.


import os, glob, subprocess, random, operator, time
from optparse import OptionParser
from tools import *

class data_class():
	def __init__(self):
		self.chr_name = ""

		self.seed_file_name = ""
		self.seed_file_prefix = ""
		self.seed_title_info = ""
		self.seed_dict = {}

		self.hifi_result_file = ""
		self.hifi_result_dict = {}

		self.hap_std_dict = {}

		self.hap_std_file_name = ""

		self.partition_size = 1000

		#self.hap_std_file_name = file_path + "ASW_" + self.chr_name + "_child_hap_refed.txt"  # 454,solid NA10847

		#self.hap_std_file_name = file_path + "NA12878_hap_new_refed.txt"	# simulation data hg18 chr6



def partition():

	data.partition_number = int(math.ceil(float(len(data.hifi_result_dict))/data.partition_size))
	#print "data.partition_number", data.partition_number
	#ref_removed_in_each_subfile = ref_number_ceilling/data_dict.number_of_subfile
	seed_removed_in_last_subfile = int(math.fmod(len(data.hifi_result_dict), data.partition_size))

	current_list_pos = 0

	pos_list = data.hifi_result_dict.keys()
	pos_list.sort()
	hifi_result_size = len(pos_list)

	#while current_list_pos < hifi_result_size:
	for i in range(data.partition_number):
		temp_result_dict = {}
		for j in range(data.partition_size):

			if (i * data.partition_size + j) < hifi_result_size:
				current_list_pos = pos_list[i * data.partition_size + j]
				#print "i * data.partition_size + j", i, j, i * data.partition_size + j, current_list_pos
				temp_result_dict[current_list_pos] = data.hifi_result_dict[current_list_pos]

		#seed_percentage = round(seed_number_in_partition/len(temp_result_dict), 2)
		#print "partition: ", i, len(temp_result_dict)
		seed_percentage, hifi_accuracy_rate = seed_calculation(temp_result_dict)
		print "partition: ", i, len(temp_result_dict), seed_percentage, hifi_accuracy_rate

def seed_calculation(temp_result_dict):
	# to calculate seed percentage and hifi errors

	seed_same_to_A = 0
	seed_same_to_B = 0
	seed_same_to_AB = 0
	seed_X = 0
	seed_N = 0
	seed_not_AB = 0

	seed_number = 0
	hifi_error = 0

	temp_result_dict_size = len(temp_result_dict)
	#temp_result_pos_list = temp_result_dict.keys()
	#temp_result_pos_list.sort()
	for position in temp_result_dict.keys():

		# check seed
		if position in data.seed_dict and position in data.hap_std_dict:
			#seed_number += 1

			line_ref = data.hap_std_dict[position]
			elements_ref = line_ref
			ref_A = elements_ref[2].strip()
			ref_B = elements_ref[3].strip()
			line_seed = data.seed_dict[position]
			elements_seed = line_seed
			seed_A = elements_seed[2].strip()
			#seed_A = elements_seed[3].strip() # for solid data 4 and 6, chr from mother
			if seed_A == ref_A:
				if seed_A == ref_B:
					seed_same_to_AB += 1
					#seed_number += 1
				else:
					seed_same_to_A += 1
					seed_number += 1
			elif seed_A == ref_B:
				seed_same_to_B += 1
			elif ref_A == "X" or ref_B == "X":
				seed_X += 1
				seed_number += 1
			elif ref_A == "N" or ref_B == "N":
				seed_N += 1
				seed_number += 1
			else:
				seed_not_AB += 1

		if position in data.hap_std_dict:
			line_ref = data.hap_std_dict[position]
			elements_ref = line_ref
			ref_A = elements_ref[2].strip()
			ref_B = elements_ref[3].strip()
			line_hifi = data.hifi_result_dict[position]
			elements_hifi = line_hifi
			hifi_A = elements_hifi[2].strip()
			hifi_B = elements_hifi[3].strip()


			if hifi_A == ref_A:
				if hifi_B == ref_B:
					hifi_AB_pos = 0
				else:
					hifi_A_pos = 0
			elif hifi_B == ref_B:
				hifi_B_pos = 0
				#hifi_error += 1
			elif (ref_A == "X" or ref_B == "X"):
				hifi_X_pos = 0
				#hifi_error += 1
			elif (ref_A == "N" or ref_B == "N"):
				hifi_N_pos = 0
				#hifi_error += 1
			else:
				if (ref_A == "A" and ref_B == "T") or (ref_A == "C" and ref_B == "G") or (
							ref_A == "T" and ref_B == "A") or (ref_A == "G" and ref_B == "C"):
					#hifi_error += 1
					pass
				else:
					hifi_error += 1

	#print "seed_number, hifi_error", seed_number, hifi_error
	return round(float(seed_number)/temp_result_dict_size*100, 4), round(float((len(temp_result_dict)-hifi_error))/temp_result_dict_size*100, 2)

def seed_distribution():
	#data.chr_name = "chr5"
	#data.seed_file_name = "song_1_prem_chr5_sorted_rmsk_indel_1_called_seed.txt"

	#data.seed_file_name = "NA12878_hg18ch6_A_0.5x_0.04er_indel_0_called_seed.txt"
	#data.hifi_result_file = "imputed_haplotype.txt"

	data.hap_std_file_name = file_path + "ASW_" + data.chr_name + "_child_hap_refed.txt"
	#data.hap_std_file_name = file_path + "NA12878_hap_new_refed.txt"

	data.seed_dict = load_raw_data(data.seed_file_name)[1]
	data.hifi_result_dict = load_raw_data(data.hifi_result_file)[1]
	data.hap_std_dict = load_raw_data(data.hap_std_file_name)[1]
	data.hap_std_dict = removeN(data.hap_std_dict)

	hap_std_total_number = len(data.hap_std_dict)
	hifi_result_total_number = len(data.hifi_result_dict)

	partition()



def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
	parser.add_option("-i", "--seed", type="string", dest="hifiSeed", help="Input seed file Name", default="null")
	parser.add_option("-r", "--result", type="string", dest="hifiResult", help="Input result file Name", default="null")
	(options, args) = parser.parse_args()
	"""
	if options.chrName == "null" or options.hifiSeed == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	"""
	return options


if __name__ == '__main__':
	options = get_args()
	#path = "/home/guoxing/node1/disk2/solid/song_1/prem_rmsk_indel/seed_distribution/"
	#os.chdir(path)

	start_time = time.time()
	global data
	data = data_class()
	data.chr_name = options.chrName
	data.seed_file_name = options.hifiSeed

	data.hifi_result_file = options.hifiResult

	seed_distribution()

	#print "run time is: ", round((time.time() - start_time), 3), "s"
