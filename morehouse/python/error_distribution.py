#!/usr/bin/python

import os, glob, subprocess, random, operator, time
from optparse import OptionParser
from tools import *
import hifiAccuCheck_v2 as hifiAccuCheck
import seed_std_compare as seed_std_compare


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

		self.hap_std_file_name = file_path + "ASW_" + self.chr_name + "_child_hap_refed.txt"  # 454,solid NA10847

	#self.hap_std_file_name = file_path + "NA12878_hap_new_refed.txt"	# simulation data hg18 chr6


def seed_distribution():
	data.chr_name = "chr5"
	data.seed_file_name = "song_1_prem_chr5_sorted_rmsk_indel_1_called_seed.txt"
	data.hifi_result_file = "imputed_haplotype.txt"

	data.hap_std_file_name = file_path + "ASW_" + data.chr_name + "_child_hap_refed.txt"

	data.seed_dict = load_raw_data(data.seed_file_name)[1]
	data.hifi_result_dict = load_raw_data(data.hifi_result_file)[1]
	data.hap_std_dict = load_raw_data(data.hap_std_file_name)[1]
	data.hap_std_dict = removeN(data.hap_std_dict)

	hap_std_total_number = len(data.hap_std_dict)
	hifi_result_total_number = len(data.hifi_result_dict)
	"""
	print hap_std_total_number, hifi_result_total_number
	same_to_A_dict_hs, same_to_B_dict_hs, same_to_AB_dict_hs, not_same_to_AB_dict_hs, \
	same_position_dict_hs, different_position_dict_hs, AT_GC_dict_hs = hifiAccuCheck.compare_std_result(
		data.hifi_result_dict, data.hap_std_dict)

	print len(same_to_AB_dict_hs)

	same_to_A_dict_ss, same_to_B_dict_ss = seed_std_compare.seed_std_compare(data.seed_file_name, data.chr_name)

	print len(same_to_A_dict_ss)
	print len(same_to_B_dict_ss)
	"""
	output_distribution()

def output_distribution():
	same_position_total_number = 0
	different_position_total_number = 0
	same_AB_total_number = 0
	same_A_total_number = 0
	same_B_total_number = 0
	not_same_AB_total_number = 0

	seed_same_to_A = 0
	seed_same_to_B = 0
	seed_same_to_AB = 0
	seed_X = 0
	seed_N = 0
	seed_not_AB = 0

	error_distribution_output_file_name = "hifi_error_distribution.txt"
	error_distribution_output_file = open(currentPath + error_distribution_output_file_name, "w")
	print >> error_distribution_output_file, "position \t seed_AB \t seed_A \t seed_B \t seed_X \t seed_N \t seed_other \t hifi_AB \t hifi_A \t hifi_B \t hifi_X \t hifi_N \t hifi_error"

	other_axis_value = "20 \t"
	seed_correct_axis_value = "40 \t"
	seed_error_axis_value = "60 \t"
	hifi_error_axis_value = "80 \t"

	snp_hap_hifi_dict_sorted_list = sort_dict_by_key(data.hap_std_dict)

	for snp in snp_hap_hifi_dict_sorted_list:
		position = snp[0]
		seed_AB_pos = "\t"
		seed_A_pos = "\t"
		seed_B_pos = "\t"
		seed_X_pos = "\t"
		seed_N_pos = "\t"
		seed_other_pos = "\t"

		hifi_AB_pos = "\t"
		hifi_A_pos = "\t"
		hifi_B_pos = "\t"
		hifi_X_pos = "\t"
		hifi_N_pos = "\t"
		hifi_other_pos = "\t"

		# check hifi results
		if position in data.hifi_result_dict:
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
					hifi_AB_pos = other_axis_value
				#same_AB_total_number += 1
				else:
					hifi_A_pos = other_axis_value
				#same_A_total_number += 1
			elif hifi_B == ref_B:  # error = B + other	assume A is the selected chr. need to update this for solid data song_4 song_6
				hifi_B_pos = other_axis_value
			#hifi_other_pos = other_axis_value
			#same_B_total_number += 1
			elif ref_A == "X" or ref_B == "X":
				hifi_X_pos = other_axis_value
			elif ref_A == "N" or ref_B == "N":
				hifi_N_pos = other_axis_value
			else:
				hifi_other_pos = hifi_error_axis_value
				not_same_AB_total_number += 1

			# check hifi seeds, these position need to be in ref too
			if position in data.seed_dict:
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
						seed_AB_pos = seed_correct_axis_value
					else:
						seed_same_to_A += 1
						seed_A_pos = seed_correct_axis_value
						seed_AB_pos = seed_correct_axis_value  # include homo seed in correct seed
				elif seed_A == ref_B:
					seed_same_to_B += 1
					seed_B_pos = seed_error_axis_value
				elif ref_A == "X" or ref_B == "X":
					seed_X += 1
					seed_X_pos = other_axis_value
				elif ref_A == "N" or ref_B == "N":
					seed_N += 1
					seed_N_pos = other_axis_value
				else:
					seed_not_AB += 1
					seed_other_pos = other_axis_value
			#if hifi_A != hifi_B: # keep hete points only
			#keep all points
			print >> error_distribution_output_file, str(
				position) + "\t", seed_AB_pos, seed_A_pos, seed_B_pos, seed_X_pos, seed_N_pos, seed_other_pos, hifi_AB_pos, hifi_A_pos, hifi_B_pos, hifi_X_pos, hifi_N_pos, hifi_other_pos

	print "seed_same_to_A", seed_same_to_A
	print "seed_same_to_B", seed_same_to_B
	print "homo seed", seed_same_to_AB
	print "seed_not_AB", seed_not_AB


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
	'''
	data.chr_name = options.chrName
	data.seed_file_name = options.hifiSeed
	data.hifi_result_file = options.hifiResult
	'''
	seed_distribution()

	print "run time is: ", round((time.time() - start_time), 3), "s"
'''


snp_hap_ref_dict = {}
snp_hap_hifi_dict = {}
snp_hap_seed_dict = {}
snp_hap_ref_total_number = 0
snp_hap_ref_total_number_withoutXN = 0
snp_hap_hifi_total_number = 0
snp_hap_seed_total_number = 0

# Reading options
usage = "usage: %prog [options] arg1"
parser = OptionParser(usage=usage)
parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
parser.add_option("-i", "--seed", type="string", dest="hifiSeed", help="Input seed file Name", default="null")

(options, args) = parser.parse_args()
chr_name = options.chrName
seed_file = options.hifiSeed

#hap_ref_file_name = "NA12878_hap_new_refed.txt"	# simulation data
hap_ref_file_name = "ASW_" + chr_name + "_child_hap_refed.txt"

hap_hifi_file_name = "imputed_haplotype.txt_4663"

inputFile_hap_ref = open(file_path + hap_ref_file_name, "r")
inputFile_hap_hifi = open(currentPath + hap_hifi_file_name, "r")
inputFile_hap_seed = open(currentPath + seed_file, "r")

x_number = 0
n_number = 0
a_number = 0
t_number = 0
c_number = 0
g_number = 0

for line in inputFile_hap_ref:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		#if elements[2].strip() != "N" and elements[2].strip() != "X" and elements[3].strip() != "N" and elements[3].strip() != "X":
		if elements[2].strip() != "N" and elements[3].strip() != "N":
			position = elements[1].strip()
			try:
				position = int(position)
				snp_hap_ref_dict[position] = line.strip()
			except ValueError:
				print position
				print line
			snp_hap_ref_total_number += 1

for line in inputFile_hap_hifi:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		if True:
			position = elements[1].strip()
			try:
				position = int(position)
				snp_hap_hifi_dict[position] = line.strip()
			except ValueError:
				print position
				print line
			snp_hap_hifi_total_number += 1

for line in inputFile_hap_seed:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		if True:
			position = elements[1].strip()
			try:
				position = int(position)
				snp_hap_seed_dict[position] = line.strip()
			except ValueError:
				print position
				print line
			snp_hap_seed_total_number += 1

print "snp_hap_ref_total_number", snp_hap_ref_total_number
print "snp_hap_hifi_total_number", snp_hap_hifi_total_number
print "snp_hap_seed_total_number", snp_hap_seed_total_number

same_position_total_number = 0
different_position_total_number = 0
same_AB_total_number = 0
same_A_total_number = 0
same_B_total_number = 0
not_same_AB_total_number = 0

seed_same_to_A = 0
seed_same_to_B = 0
seed_same_to_AB = 0
seed_X = 0
seed_N = 0
seed_not_AB = 0

error_distribution_output_file_name = "hifi_error_distribution.txt"
error_distribution_output_file = open(currentPath + error_distribution_output_file_name, "w")
#print >> error_distribution_output_file, "position seed_AB seed_A seed_B seed_X seed_N seed_other hifi_AB hifi_A hifi_B hifi_X hifi_N hifi_other"
print >> error_distribution_output_file, "position \t seed_AB \t seed_A \t seed_B \t seed_X \t seed_N \t seed_other \t hifi_AB \t hifi_A \t hifi_B \t hifi_X \t hifi_N \t hifi_error"
#print >> error_distribution_output_file, position, seed_AB_pos, seed_A_pos, seed_B_pos, seed_X_pos, seed_N_pos, seed_other_pos, hifi_AB_pos, hifi_A_pos, hifi_B_pos, hifi_X_pos, hifi_N_pos, hifi_other_pos

other_axis_value = "20 \t"
seed_correct_axis_value = "40 \t"
seed_error_axis_value = "60 \t"
hifi_error_axis_value = "80 \t"

snp_hap_hifi_dict_sorted_list = [x for x in snp_hap_hifi_dict.iteritems()]
snp_hap_hifi_dict_sorted_list.sort(key=lambda x: x[0])  # sort by key

#for position, line_hifi in snp_hap_hifi_dict.iteritems():
for snp in snp_hap_hifi_dict_sorted_list:

	position = snp[0]
	seed_AB_pos = "\t"
	seed_A_pos = "\t"
	seed_B_pos = "\t"
	seed_X_pos = "\t"
	seed_N_pos = "\t"
	seed_other_pos = "\t"

	hifi_AB_pos = "\t"
	hifi_A_pos = "\t"
	hifi_B_pos = "\t"
	hifi_X_pos = "\t"
	hifi_N_pos = "\t"
	hifi_other_pos = "\t"

	# check hifi results
	if position in snp_hap_ref_dict:
		line_ref = snp_hap_ref_dict[position]
		elements_ref = line_ref.strip().split()
		ref_A = elements_ref[2].strip()
		ref_B = elements_ref[3].strip()
		line_hifi = snp_hap_hifi_dict[position]
		elements_hifi = line_hifi.strip().split()
		hifi_A = elements_hifi[2].strip()
		hifi_B = elements_hifi[3].strip()

		if hifi_A == ref_A:
			if hifi_B == ref_B:
				hifi_AB_pos = other_axis_value
			#same_AB_total_number += 1
			else:
				hifi_A_pos = other_axis_value
			#same_A_total_number += 1
		elif hifi_B == ref_B:  # error = B + other	assume A is the selected chr. need to update this for solid data song_4 song_6
			hifi_B_pos = other_axis_value
		#hifi_other_pos = other_axis_value
		#same_B_total_number += 1
		elif ref_A == "X" or ref_B == "X":
			hifi_X_pos = other_axis_value
		elif ref_A == "N" or ref_B == "N":
			hifi_N_pos = other_axis_value
		else:
			hifi_other_pos = hifi_error_axis_value
			not_same_AB_total_number += 1

		# check hifi seeds, these position need to be in ref too
		if position in snp_hap_seed_dict:
			line_ref = snp_hap_ref_dict[position]
			elements_ref = line_ref.strip().split()
			ref_A = elements_ref[2].strip()
			ref_B = elements_ref[3].strip()
			line_seed = snp_hap_seed_dict[position]
			elements_seed = line_seed.strip().split()
			seed_A = elements_seed[2].strip()
			#seed_A = elements_seed[3].strip() # for solid data 4 and 6, chr from mother
			if seed_A == ref_A:
				if seed_A == ref_B:
					seed_same_to_AB += 1
					seed_AB_pos = seed_correct_axis_value
				else:
					seed_same_to_A += 1
					seed_A_pos = seed_correct_axis_value
					seed_AB_pos = seed_correct_axis_value  # include homo seed in correct seed
			elif seed_A == ref_B:
				seed_same_to_B += 1
				seed_B_pos = seed_error_axis_value
			elif ref_A == "X" or ref_B == "X":
				seed_X += 1
				seed_X_pos = other_axis_value
			elif ref_A == "N" or ref_B == "N":
				seed_N += 1
				seed_N_pos = other_axis_value
			else:
				seed_not_AB += 1
				seed_other_pos = other_axis_value
		#if hifi_A != hifi_B: # keep hete points only
		#keep all points
		print >> error_distribution_output_file, str(
			position) + "\t", seed_AB_pos, seed_A_pos, seed_B_pos, seed_X_pos, seed_N_pos, seed_other_pos, hifi_AB_pos, hifi_A_pos, hifi_B_pos, hifi_X_pos, hifi_N_pos, hifi_other_pos

print "seed_same_to_A", seed_same_to_A
print "seed_same_to_B", seed_same_to_B
print "homo seed", seed_same_to_AB
print "seed_not_AB", seed_not_AB

inputFile_hap_ref.close()
inputFile_hap_hifi.close()
inputFile_hap_seed.close()
error_distribution_output_file.close()

snp_hap_seed_dict_sorted_list = [x for x in snp_hap_seed_dict.iteritems()]
snp_hap_seed_dict_sorted_list.sort(key=lambda x: x[0])  # sort by key

'''