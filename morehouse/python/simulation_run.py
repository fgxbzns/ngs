#!/usr/bin/python

# location /home/guoxing/tool/morehouse

"""July 24 2013, do not count At, CG snps  """

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser
from tools import *
import data_dicts as data_dicts

from data_dicts import *
from refMerger_v5 import refMerger
from hifiAccuCheck_v2 import hifiAccuCheck, compare_std_result, seperate_homo_hetero
from snpPick_solid import snpPick
from seed_std_compare import seed_std_compare
from seed_correction_v4 import seed_correction


simulation_path = "/home/guoxing/disk2/simulation_data/"
solid_path = "/home/guoxing/disk2/solid/"

dept_list = ['4.0']
#dept_list = ['0.5', '1.0', '2.0', '4.0', '8.0']
error_rate_list = ['0.005', '0.01', '0.02', '0.04']
#error_rate_list = ['0.005']
chr_name = "chr9"


def hifi_test(hap_subfile_name="haplotype.txt", geno_subfile_name="genotype.txt", ref_subfile_name="refHaplos.txt", maf_step=0.1):
	hifi = program_path + "hifi_fu_revise " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " " + str(maf_step) + ""    # make sure the other hifi processes are finished
	#hifi = program_path + "hifi_fu_revise " + seed_input_file
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()

def hifi_test(seed_input_file):
	#hifi = program_path + "hifi_fu_revise " + hap_subfile_name + " " + geno_subfile_name + " " + ref_subfile_name + " " + str(maf_step) + ""    # make sure the other hifi processes are finished
	hifi = program_path + "hifi_fu_revise " + seed_input_file
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()

def simulation_run():
	for depth in dept_list:
		depth_path = simulation_path + depth + 'x/'
		print "depth_path", depth_path
		#os.system("cd " + depth_path)
		for error_rate in error_rate_list:
			#combined_seed_file = depth_path + "NA12878_hg18ch6_A_" + depth + "x_" + error_rate + "er_indel_0_combined_seed.txt"
			#print "combined_seed_file", combined_seed_file
			#refMerger(combined_seed_file, chr_name, 0)
			#seed_file = depth_path + "NA12878_hg18ch6_A_" + depth + "x_" + error_rate + "er_indel_0_combined_seed.txt"
			#print "seed_file", seed_file
			#hap_seed_name = "haplotype_" + depth + "x_" + error_rate + "er.txt"
			#os.system("cp " + "haplotype.txt" + " " + hap_seed_name)
			#hap_seed_name = "haplotype.txt"
			#hifi_test(hap_seed_name)
			error_path = depth_path + "NA12878_hg18ch6_A_" + depth + "x_" + error_rate + "er/"
			#os.system("cd " + error_path)
			#os.system("pwd")
			hifiAccuCheck(error_path + "imputedhaplotype.txt", chr_name)


def depth_cutoff():
	for i in range(1, 2):
		i = 11
		data_path = solid_path + "song_" + str(i) + "/prem_rmsk_indel/"
		os.chdir(data_path)
		os.system("pwd")
		sam_file = "song_" + str(i) + "_prem_" + solid_chr[i] + "_sorted_rmsk_indel.sam"
		print "**************** song_" + str(i)
		print "sam_file", sam_file
		for depth_threshold in range(0, 3):
			snpPick(data_path + sam_file, depth_threshold, solid_chr[i])
			sam_file_name = "song_" + str(i) + "_prem_" + solid_chr[i] + "_sorted_rmsk_indel_"
			seed_std_compare(data_path + sam_file_name + str(depth_threshold) + "_called_seed.txt", solid_chr[i])
			combined_seed_file = data_path + sam_file_name + str(depth_threshold) + "_combined_seed.txt"
			refMerger(combined_seed_file, solid_chr[i], 0)
			hifi_test("haplotype.txt")
			hifiAccuCheck(data_path + "imputed_haplotype.txt", solid_chr[i])


"""
def output_revised_seed(filename, selected_seed_dict):
    seed_new_file = open(currentPath + filename, "w")
    print >> seed_new_file, "rsID    position    NA10847-F       NA10847-M"
    selected_seed_sorted_list = sort_dict_by_key(selected_seed_dict)     # need to sort the snps by position
    for snp in selected_seed_sorted_list:
        seed = snp[1]
        line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
        print >> seed_new_file, line
    seed_new_file.close()
"""


def generate_std_seed_run(seed_number, add_range):
	chr_name = "chr9"
	hap_std_dict = load_seed_data(file_path + "ASW_" + chr_name + "_child_hap_refed.txt")[1]
	hap_std_list = sort_dict_by_key(hap_std_dict)

	genotype_file = file_path + "genotype_NA10847_" + chr_name + ".txt"
	geno_dict = load_raw_data(genotype_file)[1]
	seed_homo_dict, seed_hetero_dict = group_seed(hap_std_dict, geno_dict)

	seed_hetero_list = sort_dict_by_key(seed_hetero_dict)

	selected_seed_dict = {}

	if add_range == "begining":
		# add seed from begining
		for i in range(seed_number - 1):
			selected_seed_dict[seed_hetero_list[i][0]] = seed_hetero_list[i][1]
	elif add_range == "middle":
		# add seed in the middle
		hetero_total = len(seed_hetero_list)
		middle_point = hetero_total / 2
		for i in range(seed_number - 1):
			selected_seed_dict[seed_hetero_list[middle_point + i][0]] = seed_hetero_list[middle_point + i][1]
	elif add_range == "end":
		# add seed from end
		seed_hetero_list.reverse()
		for i in range(1, seed_number):
			selected_seed_dict[seed_hetero_list[i][0]] = seed_hetero_list[i][1]
	elif add_range == "random":
		# randomly adding seed
		i = 0
		while i < seed_number - 1:
			random_index = random.randrange(0, (len(seed_hetero_list) - 1))
			while seed_hetero_list[random_index][0] in selected_seed_dict or seed_hetero_list[random_index][
				1].allele_new == "N" \
					or seed_hetero_list[random_index][1].allele_new == "X":
				#while seed_hetero_list[random_index][0] in selected_seed_dict:
				random_index = random.randrange(0, (len(seed_hetero_list) - 1))
			selected_seed_dict[seed_hetero_list[random_index][0]] = seed_hetero_list[random_index][1]
			i += 1

	# always add the last snp into seed, hifi requirement
	selected_seed_dict[hap_std_list[-1][0]] = hap_std_list[-1][1]

	file_name = "haplotype_std.txt"
	output_revised_seed(file_name, selected_seed_dict)
	seed_std_compare(file_name, chr_name)
	refMerger(file_name, chr_name, 0)
	file_name = "haplotype.txt"
	hifi_run(file_name, chr_name)
	hifiAccuCheck("imputed_" + file_name, chr_name)


def rareSNP_error_rate():

	#data.update_ref_dict()
	#data.update_hap_ref_allele_frequence_dict()
	#data.output_hap_ref_allele_frequence_dict()
	data.read_hap_ref_allele_frequence_dict()
	print len(data.hap_ref_allele_frequence_dict)

	maf_result_pos_dict = {}
	maf_range_high = 1.00
	maf_range_low = 0.99
	for pos, list in data.hap_ref_allele_frequence_dict.iteritems():
		if list[0][1] >= maf_range_low and list[0][1] < maf_range_high:
			maf_result_pos_dict[pos] = ""
	print len(maf_result_pos_dict)

	#hifi_result_file = "imputed_haplotype.txt"
	hifi_result_file = "non_one.txt"
	hifi_result_dict = load_raw_data(hifi_result_file)[1]
	hifi_result_total_number = len(hifi_result_dict)
	maf_result_dict = {pos: value for pos, value in hifi_result_dict.iteritems() if pos in maf_result_pos_dict}
	maf_result_total_number = len(maf_result_dict)

	chr_name = "chr9"
	data.chr_name = chr_name
	data.geno_dict = data.update_geno_dict()

	genotype_file = file_path + "genotype_NA10847_" + chr_name + ".txt"
	geno_dict = load_raw_data(genotype_file)[1]

	print len(geno_dict)
	maf_result_homo_dict, maf_result_hetero_dict = group_seed(maf_result_dict, geno_dict)

	print "maf hetero %", round(float(len(maf_result_hetero_dict))/maf_result_total_number*100,3)

	hap_std_file_name = file_path + "ASW_" + chr_name + "_child_hap_refed.txt"
	#data.hap_std_file_name = file_path + "NA12878_hap_new_refed.txt"
	hap_std_dict = load_raw_data(hap_std_file_name)[1]

	same_to_A_dict, same_to_B_dict, same_to_AB_dict, not_same_to_AB_dict, same_position_dict, different_position_dict, \
	AT_GC_dict = compare_std_result(maf_result_dict, hap_std_dict)

	same_A_total_number = len(same_to_A_dict)
	same_B_total_number = len(same_to_B_dict)
	same_AB_total_number = len(same_to_AB_dict)
	not_same_AB_total_number = len(not_same_to_AB_dict)
	same_position_total_number = len(same_position_dict)
	different_position_total_number = len(different_position_dict)
	AT_GC_dict_number = len(AT_GC_dict)

	pencentage_in_common = format(float(same_position_total_number) / float(maf_result_total_number) * 100, "0.3f")
	accuracy = round(float(same_A_total_number + same_B_total_number + same_AB_total_number) / float(
		same_position_total_number - AT_GC_dict_number) * 100, 3)

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
		                                                                               - AT_GC_dict_number) * 100, 3)
	homo_accuracy = round(float(len(same_AB_homo)) / float(len(same_position_homo)) * 100, 3)

	print "same_position_total_number", same_position_total_number
	print "different_position_total_number", different_position_total_number
	print "same_AB_total_number", same_AB_total_number
	print "same_A_total_number", same_A_total_number
	print "same_B_total_number", same_B_total_number
	print "AT_GC_dict_number", AT_GC_dict_number
	print "not_same_AB_total_number", not_same_AB_total_number
	print "pencentage in common", pencentage_in_common
	print "accuracy", accuracy
	#print "len(same_AB_hetero)", len(same_AB_hetero)
	#print "len(same_position_hetero)", len(same_position_hetero)
	#print "len(AT_GC_hetero)", len(AT_GC_hetero)
	print "hetero_accuracy", hetero_accuracy
	print "homo_accuracy", homo_accuracy

def output_files(file_name, title_info, dict):
	outpuf_file = open(currentPath + file_name, "w")	# for hifi
	print >> outpuf_file, title_info
	sorted_list = sort_dict_by_key(dict)
	for element in sorted_list:
		print >> outpuf_file, element[1][0], element[1][1], element[1][2]
	outpuf_file.close()

def genotype_missing(remPercent):

	ori_genotype_file_name = "genotype_ori.txt"
	geno_title_info, geno_dict = load_raw_data(ori_genotype_file_name)
	ori_geno_total = len(geno_dict)
	print "ori_geno_total", ori_geno_total

	geno_sorted_list = sort_dict_by_key(geno_dict)
	if remPercent > 0:
		print "remPercent is", remPercent
		for i in range(int(remPercent*ori_geno_total)):
			if len(geno_sorted_list) > 1:
				random_index = random.randrange(0, (len(geno_sorted_list)-1))
				position = geno_sorted_list[random_index][0]
				if position in geno_dict:
					del geno_dict[int(position)]
					del geno_sorted_list[random_index]

			geno_sorted_list = sort_dict_by_key(geno_dict)

	print "new geno total", len(geno_dict)
	output_files("genotype.txt", geno_title_info, geno_dict)
	os.system("cp genotype.txt genotype_" + str(remPercent) + ".txt")
	chr_name = "chr9"

	file_name = "haplotype.txt"
	hifi_run(file_name, chr_name)
	mode = "filterbyrefid"
	seed_correction(file_name, chr_name, mode)
	hifiAccuCheck("imputed_" + file_name, chr_name)
	print "filterbyrefid accuracy"
	hifiAccuCheck("non_one.txt", chr_name)

if __name__ == '__main__':
	start_time = time.time()
	#options = get_args()
	global data
	data = data_dicts()
	#chr_name = options.chrName
	#hifi_result_file = options.hifiResult
	#hifiAccuCheck(hifi_result_file, chr_name)
	#simulation_run()
	#depth_cutoff()
	add_range = "middle"
	#add_range = "random"
	#generate_std_seed_run(2500, add_range)
	#rareSNP_error_rate()
	for i in (0.8, 0.9):
		genotype_missing(i)

	print "run time is: ", round((time.time() - start_time), 3), "s"
