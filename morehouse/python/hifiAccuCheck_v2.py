#!/usr/bin/python

# location /home/guoxing/tool/morehouse


import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

currentPath = os.getcwd() + '/'
raw_data_format = "list"
file_path = "/home/guoxing/disk2/solid/common_files/"


snp_hap_ori_dict = {}
hifi_result_dict = {}
hap_std_total_number = 0
hap_std_total_number_withoutXN = 0
hifi_result_total_number = 0

def list_to_line(list):
	line = ""
	for a in list:
		line += str(a).strip() + "\t"
	return line.strip()


def load_raw_data(file_name, raw_data_format="list"):
	title_info = ""
	data = {}
	with open(file_name, "r") as fp:
		for line in fp:
			if line.startswith("rsID"):
				title_info = list_to_line(line.strip().split())
			else:
				elements = line.strip().split()
				try:
					# convert the position to int for sorting
					if raw_data_format == "list":
						data[int(elements[1])] = elements
					elif raw_data_format == "string":
						data[int(elements[1])] = line.strip()
				except:
					#print "error in ", line, file_name
					pass
	return (title_info, data)


def removeN(hifi_std_dict):
	temp_dict = {}
	for position, elements in hifi_std_dict.iteritems():
		if elements[2].strip() != "N" and elements[3].strip() != "N":
			temp_dict[position] = hifi_std_dict[position]
	return temp_dict

def allele_similarity(hifi_result_dict, hifi_std_dict):
	same_to_A = 0
	same_to_B = 0
	for position, elements_hifi in hifi_result_dict.iteritems():
		if position in hifi_std_dict:
			hifi_A = elements_hifi[2].strip()
			hifi_B = elements_hifi[3].strip()
			elements_std = hifi_std_dict[position]
			std_A = elements_std[2].strip()
			std_B = elements_std[3].strip()
			if hifi_A != hifi_B:
				if hifi_A == std_A:
					same_to_A += 1
				if hifi_A == std_B:
					same_to_B += 1
	#print "similarity, same_to_A", same_to_A
	#print "similarity, same_to_B", same_to_B
	if same_to_A >= same_to_B:
		print "same_to_A", same_to_A, "total", same_to_A + same_to_B, "A percentage", float(same_to_A)/(same_to_A + same_to_B+1)
		return "similar_to_A"
	else:
		print "same_to_B", same_to_B, "total", same_to_A + same_to_B, "B percentage", float(same_to_B)/(same_to_A + same_to_B+1)
		return "similar_to_B"

def compare_std_result(hifi_result_dict, hifi_std_dict):
	same_to_A_dict = {}
	same_to_B_dict = {}
	same_to_X_dict = {}
	same_to_N_dict = {}
	same_to_AB_dict = {}
	AT_GC_dict = {}
	not_same_to_AB_dict = {}
	same_position_dict = {}
	different_position_dict = {}
	hifi_result_x_dict = {}
	std_x_dict = {}

	same_position_total_number = 0
	different_position_total_number = 0
	similarity = allele_similarity(hifi_result_dict, hifi_std_dict)

	for position, elements_hifi in hifi_result_dict.iteritems():
		if position in hifi_std_dict:
			hifi_A = elements_hifi[2].strip()
			hifi_B = elements_hifi[3].strip()
			#if hifi_A != 'X' and hifi_B != 'X' and hifi_A != 'N' and hifi_B != 'N':
			if True:
				elements_std = hifi_std_dict[position]
				std_A = elements_std[2].strip()
				std_B = elements_std[3].strip()

				if similarity == "similar_to_B":  # for solid data 4 and 6, the hifi seed is from mother (B)
					hifi_A, hifi_B = hifi_B, hifi_A
				# the hifi seed is from father, A
				if hifi_A == std_A:  #A is A
					if hifi_B == std_B:
						same_to_AB_dict[position] = elements_hifi
					else:
						same_to_A_dict[position] = elements_hifi
				elif hifi_B == std_B:
					same_to_B_dict[position] = elements_hifi
				elif std_A == "X" or std_B == "X" or std_A == "N" or std_B == "N":
					same_to_AB_dict[position] = elements_hifi
				else:
					if (std_A == "A" and std_B == "T") or (std_A == "C" and std_B == "G") or (
							std_A == "T" and std_B == "A") or (std_A == "G" and std_B == "C"):
						AT_GC_dict[position] = elements_hifi
						hifi_result_x_dict[position] = elements_hifi
					elif hifi_A == "X" or hifi_B == "X":
						#same_to_AB_dict[position] = elements_hifi
						pass
					else:
						not_same_to_AB_dict[position] = elements_hifi
				same_position_dict[position] = elements_hifi

			if hifi_A == "X" or hifi_B == "X":
				hifi_result_x_dict[position] = elements_hifi
			if std_A == "X" or std_B == "X" or std_A == "N" or std_B == "N":
				std_x_dict[position] = elements_hifi

			#same_position_total_number += 1
		else:
			#different_position_total_number += 1
			different_position_dict[position] = elements_hifi

	return (
	same_to_A_dict, same_to_B_dict, same_to_AB_dict, not_same_to_AB_dict, same_position_dict, different_position_dict,
	AT_GC_dict, hifi_result_x_dict, std_x_dict)


def output_dict(file_name, dict, hifi_std_dict):
	output_file = open(currentPath + file_name, "w")
	for position, elements in dict.iteritems():
		print >> output_file, list_to_line(elements), "std:", list_to_line(hifi_std_dict[position])
	output_file.close()


def seperate_homo_hetero(same_to_AB_dict):
	homo_dict = {}
	hetero_dict = {}
	for position, snp in same_to_AB_dict.iteritems():
		#if snp[2] != 'X' and snp[3] != 'X' and snp[2] != 'N' and snp[3] != 'N':
		if snp[2] == snp[3]:
			homo_dict[position] = snp
		else:
			hetero_dict[position] = snp
	return (homo_dict, hetero_dict)


def hifiAccuCheck(hifi_result_file, chr_name):
	hap_std_file_name = file_path + "ASW_" + chr_name + "_child_hap_refed.txt"  # 454,solid NA10847
	#hap_std_file_name = file_path + "NA12878_hap_new_refed.txt"	# simulation data hg18 chr6
	#hap_std_file_name = "standardhaplotype.txt"         # chr1 NA11919_A
	#hap_std_file_name = std_name

	hifi_std_dict = load_raw_data(hap_std_file_name, raw_data_format)[1]
	hap_std_total_number = len(hifi_std_dict)

	hifi_std_dict = removeN(hifi_std_dict)

	hifi_result_dict = load_raw_data(hifi_result_file, raw_data_format)[1]
	hifi_result_total_number = len(hifi_result_dict)

	print "hap_std_total_number", hap_std_total_number
	print "hifi_result_total_number", hifi_result_total_number

	same_to_A_dict, same_to_B_dict, same_to_AB_dict, not_same_to_AB_dict, same_position_dict, different_position_dict, \
	AT_GC_dict, hifi_result_x_dict, std_x_dict = compare_std_result(hifi_result_dict, hifi_std_dict)

	"""
	output_dict("same_to_A_dict.txt", same_to_A_dict, hifi_std_dict)
	output_dict("same_to_B_dict.txt", same_to_B_dict, hifi_std_dict)
	output_dict("not_same_to_AB_dict.txt", not_same_to_AB_dict, hifi_std_dict)
	output_dict("different_position_dict.txt", different_position_dict, hifi_std_dict)
	"""
	same_A_total_number = len(same_to_A_dict)
	same_B_total_number = len(same_to_B_dict)
	same_AB_total_number = len(same_to_AB_dict)
	not_same_AB_total_number = len(not_same_to_AB_dict)
	same_position_total_number = len(same_position_dict)
	different_position_total_number = len(different_position_dict)
	AT_GC_dict_number = len(AT_GC_dict)
	hifi_result_x_number = len(hifi_result_x_dict)
	std_x_number = len(std_x_dict)

	pencentage_in_common = round(float(same_position_total_number) / hifi_result_total_number * 100, 3)
	error_rate = round(float(not_same_AB_total_number)/(hifi_result_total_number)*100, 3)

	accuracy = round((same_A_total_number + same_B_total_number + same_AB_total_number)/float(same_position_total_number - hifi_result_x_number - AT_GC_dict_number), 4)
	#accuracy = round((same_A_total_number + same_B_total_number + same_AB_total_number)/float(hifi_result_total_number - AT_GC_dict_number - hifi_result_x_number), 4)

	same_AB_homo, same_AB_hetero = seperate_homo_hetero(same_to_AB_dict)
	not_same_AB_homo, not_same_AB_hetero = seperate_homo_hetero(not_same_to_AB_dict)
	different_AB_homo, different_AB_hetero = seperate_homo_hetero(different_position_dict)
	AT_GC_homo, AT_GC_hetero = seperate_homo_hetero(AT_GC_dict)
	imputed_homo, imputed_hetero = seperate_homo_hetero(hifi_result_dict)
	same_position_homo, same_position_hetero = seperate_homo_hetero(same_position_dict)
	#hetero_accuracy = round(float(len(same_AB_hetero))/float(len(same_position_hetero)-len(AT_GC_hetero))*100, 3)
	#homo_accuracy = round(float(len(same_AB_homo))/float(len(same_position_homo)-len(AT_GC_homo))*100, 3)

	#hetero_accuracy = round(float(same_A_total_number + same_B_total_number + len(same_AB_hetero))/float(len(same_AB_hetero)+len(not_same_AB_hetero) \
	#														-AT_GC_dict_number )*100, 3)
	hetero_accuracy = round(
		float(same_A_total_number + same_B_total_number + len(same_AB_hetero)) / float(len(same_position_hetero) \
		                                                                               - AT_GC_dict_number) * 100, 2)
	#homo_accuracy = round(float(len(same_AB_homo)) / float(len(same_position_homo)) * 100, 3)

	coverage = round(float(same_position_total_number)/(hifi_result_total_number)*100, 2)
	coverage = round(float(same_position_total_number)/(hap_std_total_number)*100, 2)

	print "same_position_total_number", same_position_total_number
	print "different_position_total_number", different_position_total_number
	print "same_AB_total_number", same_AB_total_number
	print "same_A_total_number", same_A_total_number
	print "same_B_total_number", same_B_total_number
	print "AT_GC_dict_number", AT_GC_dict_number
	print "not_same_AB_total_number", not_same_AB_total_number
	print "pencentage in common", pencentage_in_common
	#print "coverage", coverage

	#print "len(same_AB_hetero)", len(same_AB_hetero)
	#print "len(same_position_hetero)", len(same_position_hetero)
	#print "len(AT_GC_hetero)", len(AT_GC_hetero)
	#print "hetero_accuracy", hetero_accuracy
	#print "homo_accuracy", homo_accuracy
	print "hifi_result_x_number", hifi_result_x_number
	print "std_x_number", std_x_number

	print "accuracy", accuracy

	accuracy_output_file_name = "hifi_accuracy.txt"
	accuracy_output_file = open(currentPath + accuracy_output_file_name, "w")
	print >> accuracy_output_file, "hifi result file: ", currentPath + hifi_result_file
	print >> accuracy_output_file, "hifi std file: ", hap_std_file_name
	print >> accuracy_output_file, "hap_std_total_number: ", hap_std_total_number
	print >> accuracy_output_file, "hifi_result_total_number: ", hifi_result_total_number
	print >> accuracy_output_file, "same_position_total_number: ", same_position_total_number
	print >> accuracy_output_file, "different_position_total_number: ", different_position_total_number
	print >> accuracy_output_file, "same_A_total_number: ", same_A_total_number
	print >> accuracy_output_file, "same_B_total_number: ", same_B_total_number
	print >> accuracy_output_file, "same_AB_total_number: ", same_AB_total_number
	print >> accuracy_output_file, "not_same_AB_total_number: ", not_same_AB_total_number
	print >> accuracy_output_file, "pencentage in common: ", pencentage_in_common
	print >> accuracy_output_file, "accuracy: ", accuracy

	accuracy_output_file.close()
	"""
	# record data
	data_record_file_name = "solid_process_4.txt"
	data_record_file = open(data_record_path + data_record_file_name, "a")
	print >> data_record_file, currentPath, same_position_total_number, (same_A_total_number+same_B_total_number+not_same_AB_total_number), same_AB_total_number, hifi_result_total_number, accuracy
	data_record_file.close()
	cmd = "grep hifi_data hifi_accuracy.txt >> " + data_record_path + data_record_file_name
	#print cmd
	#os.system(cmd)
	"""
	return same_to_AB_dict, AT_GC_dict


def hifiAccuCheck_file(hifi_result_file, hap_std_file_name):

	hifi_std_dict = load_raw_data(hap_std_file_name, raw_data_format)[1]
	hap_std_total_number = len(hifi_std_dict)

	hifi_std_dict = removeN(hifi_std_dict)

	hifi_result_dict = load_raw_data(hifi_result_file, raw_data_format)[1]
	hifi_result_total_number = len(hifi_result_dict)

	print "hap_std_total_number", hap_std_total_number
	print "hifi_result_total_number", hifi_result_total_number

	same_to_A_dict, same_to_B_dict, same_to_AB_dict, not_same_to_AB_dict, same_position_dict, different_position_dict, \
	AT_GC_dict, hifi_result_x_dict, std_x_dict = compare_std_result(hifi_result_dict, hifi_std_dict)

	for pos in not_same_to_AB_dict.keys():
		#print pos, hifi_result_dict[pos], hifi_std_dict[pos]
		pass


	"""
	output_dict("same_to_A_dict.txt", same_to_A_dict, hifi_std_dict)
	output_dict("same_to_B_dict.txt", same_to_B_dict, hifi_std_dict)
	output_dict("not_same_to_AB_dict.txt", not_same_to_AB_dict, hifi_std_dict)
	output_dict("different_position_dict.txt", different_position_dict, hifi_std_dict)
	"""
	same_A_total_number = len(same_to_A_dict)
	same_B_total_number = len(same_to_B_dict)
	same_AB_total_number = len(same_to_AB_dict)
	not_same_AB_total_number = len(not_same_to_AB_dict)
	same_position_total_number = len(same_position_dict)
	different_position_total_number = len(different_position_dict)
	AT_GC_dict_number = len(AT_GC_dict)
	hifi_result_x_number = len(hifi_result_x_dict)
	std_x_number = len(std_x_dict)

	pencentage_in_common = round(float(same_position_total_number) / hifi_result_total_number * 100, 3)
	error_rate = round(float(not_same_AB_total_number)/(hifi_result_total_number)*100, 3)

	accuracy = round((same_A_total_number + same_B_total_number + same_AB_total_number)/float(same_position_total_number - AT_GC_dict_number), 4)
	#accuracy = round((same_A_total_number + same_B_total_number + same_AB_total_number)/float(hifi_result_total_number - AT_GC_dict_number - hifi_result_x_number), 4)

	same_AB_homo, same_AB_hetero = seperate_homo_hetero(same_to_AB_dict)
	not_same_AB_homo, not_same_AB_hetero = seperate_homo_hetero(not_same_to_AB_dict)
	different_AB_homo, different_AB_hetero = seperate_homo_hetero(different_position_dict)
	AT_GC_homo, AT_GC_hetero = seperate_homo_hetero(AT_GC_dict)
	imputed_homo, imputed_hetero = seperate_homo_hetero(hifi_result_dict)
	same_position_homo, same_position_hetero = seperate_homo_hetero(same_position_dict)
	#hetero_accuracy = round(float(len(same_AB_hetero))/float(len(same_position_hetero)-len(AT_GC_hetero))*100, 3)
	#homo_accuracy = round(float(len(same_AB_homo))/float(len(same_position_homo)-len(AT_GC_homo))*100, 3)

	#hetero_accuracy = round(float(same_A_total_number + same_B_total_number + len(same_AB_hetero))/float(len(same_AB_hetero)+len(not_same_AB_hetero) \
	#													-AT_GC_dict_number )*100, 3)
	#hetero_accuracy = round(
	#	float(same_A_total_number + same_B_total_number + len(same_AB_hetero)) / float(len(same_position_hetero) \
	#	                                                                               - AT_GC_dict_number) * 100, 2)
	#homo_accuracy = round(float(len(same_AB_homo)) / float(len(same_position_homo)) * 100, 3)

	coverage = round(float(same_position_total_number)/(hifi_result_total_number)*100, 2)
	coverage = round(float(same_position_total_number)/(hap_std_total_number)*100, 2)

	print "same_position_total_number", same_position_total_number
	print "different_position_total_number", different_position_total_number
	print "same_AB_total_number", same_AB_total_number
	print "same_A_total_number", same_A_total_number
	print "same_B_total_number", same_B_total_number
	print "AT_GC_dict_number", AT_GC_dict_number
	print "not_same_AB_total_number", not_same_AB_total_number
	print "pencentage in common", pencentage_in_common
	#print "coverage", coverage

	#print "len(same_AB_hetero)", len(same_AB_hetero)
	#print "len(same_position_hetero)", len(same_position_hetero)
	#print "len(AT_GC_hetero)", len(AT_GC_hetero)
	#print "hetero_accuracy", hetero_accuracy
	#print "homo_accuracy", homo_accuracy
	print "hifi_result_x_number", hifi_result_x_number
	print "std_x_number", std_x_number

	print "accuracy", accuracy

	accuracy_output_file_name = "hifi_accuracy.txt"
	accuracy_output_file = open(currentPath + accuracy_output_file_name, "w")
	print >> accuracy_output_file, "hifi result file: ", currentPath + hifi_result_file
	print >> accuracy_output_file, "hifi std file: ", hap_std_file_name
	print >> accuracy_output_file, "hap_std_total_number: ", hap_std_total_number
	print >> accuracy_output_file, "hifi_result_total_number: ", hifi_result_total_number
	print >> accuracy_output_file, "same_position_total_number: ", same_position_total_number
	print >> accuracy_output_file, "different_position_total_number: ", different_position_total_number
	print >> accuracy_output_file, "same_A_total_number: ", same_A_total_number
	print >> accuracy_output_file, "same_B_total_number: ", same_B_total_number
	print >> accuracy_output_file, "same_AB_total_number: ", same_AB_total_number
	print >> accuracy_output_file, "not_same_AB_total_number: ", not_same_AB_total_number
	print >> accuracy_output_file, "pencentage in common: ", pencentage_in_common
	print >> accuracy_output_file, "accuracy: ", accuracy

	accuracy_output_file.close()
	"""
	# record data
	data_record_file_name = "solid_process_4.txt"
	data_record_file = open(data_record_path + data_record_file_name, "a")
	print >> data_record_file, currentPath, same_position_total_number, (same_A_total_number+same_B_total_number+not_same_AB_total_number), same_AB_total_number, hifi_result_total_number, accuracy
	data_record_file.close()
	cmd = "grep hifi_data hifi_accuracy.txt >> " + data_record_path + data_record_file_name
	#print cmd
	#os.system(cmd)
	"""
	return same_to_AB_dict, AT_GC_dict


def get_args():
	desc = "Compare seed and std hap, to check purity of seed"
	usage = "seed_std_compare -i seed_file -c chr#"
	parser = OptionParser(usage=usage)
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="chr11")
	parser.add_option("-i", "--imputed", type="string", dest="hifiResult", help="Input hifiResult file Name",
	                  default="null")
	parser.add_option("-s", "--s", type="string", dest="std_name", help="Input std Name", default="")


	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.hifiResult == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options


if __name__ == '__main__':
	options = get_args()
	chr_name = options.chrName
	global std_name
	std_name = options.std_name
	hifi_result_file = options.hifiResult
	hifiAccuCheck(hifi_result_file, chr_name)
