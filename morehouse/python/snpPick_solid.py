#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for solid data May 09 2013
# June 22, 2013, added purity. >= 0.90 max_allele_number will be kept as seed. A. B, X, homo will be kept as seeds.
# quality score check. ord('!')-33 > 13

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *
from hifiAccuCheck_v2 import hifiAccuCheck
from seed_std_compare import seed_std_compare

quality_score_threshold = 13
max_allele_percentage_threshold = 0.6

# A from Father, B from Mother
class snps:
	def __init__(self):
		self.rsID = ""
		self.position = 0
		self.A = ""
		self.B = ""
		self.covered_reads_list = []
		self.allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
		self.consistence = True
		self.max_allele = ""
		self.max_allele_number = 0
		self.max_allele_percentage = 0

# class to store reads from sam file
class read:
	def __init__(self, qName, flag, rName, start_position, read_sequence, quality_score_sequence, read_length, covered_snp):
		self.qName = qName
		self.flag = flag
		self.rName = rName
		self.start_position = start_position
		self.read_sequence = read_sequence
		self.quality_score_sequence = quality_score_sequence
		self.read_length = read_length
		self.covered_snp = covered_snp
		
def get_args():
	desc="variation call"
	usage = "snpPick_solid -i haplotypeFile -c chr#" 
	parser = OptionParser(usage = usage)
	parser.add_option("-i", "--haplotype", type="string", dest="haplotypeFile",help = "Input File Name", default="null")
	parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
	parser.add_option("-d", "--threshold", type="string", dest="threshold",help = "Input the depth threshold", default="3")
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="chr11")
	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.samFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

def load_hap_std(file_name):
	hap_std_dict = {}
	title_haplotype = load_raw_data(file_name, raw_data_format)[0]
	data_dict = load_raw_data(file_name, raw_data_format)[1]
	for position, elements in data_dict.iteritems():
		if True:
			try:
				snp = snps()
				snp.rsID = elements[0]
				snp.position = int(position)
				snp.A = elements[2].strip()
				snp.B = elements[3].strip()
				hap_std_dict[position] = snp
			except ValueError:
				print "error in ", file_name, position
	return (title_haplotype, hap_std_dict)	


# start time
start = time.time()		

options = get_args()
haplotype_file = options.haplotypeFile
sam_file = options.samFile
depth_threshold = int(options.threshold)
chr_name = options.chrName

#haplotype_file = "NA12878_hap_new_refed.txt"	# simulation data chr6
haplotype_file = "ASW_"+chr_name+"_child_hap_refed.txt"	# for all
genotype_file = "genotype_NA10847_"+chr_name+".txt"	# for all

	
sam_file_name = sam_file[:(len(sam_file)-4)] + "_" + str(depth_threshold)

print "haplotype file: ", haplotype_file
print "genotype_file : ", genotype_file
print "sam file: ", sam_file_name








data_record_file_name = sam_file_name + "_data_record.txt"
data_record_file = open(currentPath + data_record_file_name, "w")



hap_std_tuple = load_hap_std(file_path + haplotype_file)
title_haplotype = hap_std_tuple[0]
hap_std_dict = hap_std_tuple[1]

genotype_file = "genotype_NA10847_"+chr_name+".txt"	# for all
geno_dict = load_raw_data(file_path + genotype_file, raw_data_format)[1]


reads_list=[]
insert_size = 0


def variant_call(sam_file, hap_std_dict):
	inputfile_sum = open(currentPath + sam_file, "r")
	sam_line_first = inputfile_sum.readline() # the first read line in a pair
	total_reads_num = 0
	covered_snp_total_number = 0

	# solid data, single end. no insert size
	while sam_line_first!='':
		if not sam_line_first.startswith("@"):
			total_reads_num += 1	
			elements_first = sam_line_first.strip().split()
			try:
				read_ID_first = elements_first[0].strip()
				rName_first = elements_first[2].strip()
			except:
				print "error in first read:", sam_line_first						
			# first read
			qName_first = elements_first[0].strip()
			flag_first = elements_first[1].strip()
			start_position_first = int(elements_first[3].strip())
			read_sequence_first = elements_first[9].strip()
			read_length_first = len(read_sequence_first)
			quality_score_sequence_first = elements_first[10].strip()
			i = 0
			while i < read_length_first:
				if (start_position_first+i) in hap_std_dict:							
					covered_snp = read_sequence_first[i]			# ith position is the covered snp
					quality_score_symbol = quality_score_sequence_first[i]
					if (rName_first == chr_name) and (not covered_snp == 'N') and ((ord(quality_score_symbol)-33) > quality_score_threshold):	# check quality_score
						covered_snp_total_number += 1
						hap_std_dict[start_position_first+i].covered_reads_list.append(read(qName_first, flag_first, rName_first, \
						start_position_first, read_sequence_first, quality_score_sequence_first, read_length_first, covered_snp))					
						for base, value in hap_std_dict[start_position_first+i].allele_dict.iteritems():
							if base == covered_snp:
								hap_std_dict[start_position_first+i].allele_dict[base] += 1	
				i += 1			
		sam_line_first = inputfile_sum.readline()	
	inputfile_sum.close()
	hap_std_sorted_list = sort_dict_by_key(hap_std_dict)
	return (total_reads_num, covered_snp_total_number, hap_std_sorted_list)

variant_call_tulpe = variant_call(sam_file, hap_std_dict)
total_reads_num = variant_call_tulpe[0]
covered_snp_total_number = variant_call_tulpe[1]
hap_std_sorted_list = variant_call_tulpe[2]

print "haplotye list size is: ", len(hap_std_dict)	
print "total_reads_num", total_reads_num

def output_coverage_info(hap_std_sorted_list):
	outputFile_reads = open(currentPath + sam_file_name + "_reads.txt", "w")
	outputFile_reads.write("SNP position \t Depth \n")
	outputfile_allele = open(currentPath + sam_file_name+"_allele.txt", "w")
	outputfile_allele.write("Chromosome \t position \t Total Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")
	for snp_data in hap_std_sorted_list:
		snp = snp_data[1]
		if len(snp.covered_reads_list) > depth_threshold:
			outputfile_allele.write(chr_name+"\t"+str(snp.position)+"\t"+str(len(snp.covered_reads_list))+"\t"+str(snp.allele_dict['A'])	\
									+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+"\n")
			outputFile_reads.write("@_" + snp.rsID + "\t" + str(snp.position) + "\t" + snp.A  + "\t" + snp.B  + "\t" + str(len(snp.covered_reads_list)) + "\n")
			for reads in snp.covered_reads_list:
				outputFile_reads.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" \
										+ reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
	outputFile_reads.close()
	outputfile_allele.close()

def update_max_allele_number(max_allele_number, dict):
	if max_allele_number not in dict:
		dict[max_allele_number] = 1
	else:
		dict[max_allele_number] += 1

def output_ABX_files(file_name, title_info, data_dict):
	file = open(currentPath + file_name, "w")
	print >> file, title_info
	for position, snp in data_dict.iteritems():
		print >> file, snp.rsID, snp.position, snp.A, snp.B, snp.max_allele, snp.max_allele_number, snp.max_allele_percentage
		for reads in snp.covered_reads_list:
			print >> file, reads.qName, reads.flag, reads.rName, reads.start_position, reads.covered_snp, reads.read_sequence, reads.quality_score_sequence
	file.close()
		
def output_seed_file(file_name, title_info, data_dict):
	file = open(currentPath + file_name, "w")
	print >> file, title_info
	for position, snp in data_dict.iteritems():
		print >> file, snp.rsID, snp.position, snp.max_allele
	file.close()

# correct errors in seed, for testing purpose
hifi_pure_corrected_file_name = sam_file_name + "_hifi_corrected.txt"
hifi_pure_corrected_file = open(currentPath + hifi_pure_corrected_file_name, "w")		
hifi_pure_corrected_file.write(title_haplotype + "\n")

distribution_file = open(currentPath + sam_file_name + "_distribution.txt", "w")
distribution_file.write("rsID \t phys_position \t snp.A	\t snp.B \t distribution \n")


pure_total = 0
not_in_genotype = 0
A_dict = {}
B_dict = {}
X_dict = {}
not_ABX_dict = {}
not_in_genotype_dict_num = {}

called_seed_dict = {}
same_to_A_dict = {}
same_to_B_dict = {}
same_to_AB_dict = {}
same_to_X_dict = {}
same_to_N_dict = {}
not_in_genotype_dict = {}

for snp_data in hap_std_sorted_list:
	snp = snp_data[1]
	position = snp.position
	if len(snp.covered_reads_list) > depth_threshold:	
		max_allele = keywithmaxval(snp.allele_dict)
		max_allele_number = snp.allele_dict[max_allele]
		snp.max_allele = max_allele 
		snp.max_allele_number = max_allele_number 
		snp.max_allele_percentage = float(max_allele_number)/float(len(snp.covered_reads_list)) 

		if snp.max_allele_percentage >= max_allele_percentage_threshold and max_allele_number > depth_threshold:
			pure_total += 1
			# check genotype to remove called base that does not in genotype
			if max_allele in geno_dict[position][2]:	
				called_seed_dict[position] = snp			
				if max_allele == snp.A or max_allele == snp.B:	#check genotype
					if max_allele == snp.A and max_allele != snp.B:
						same_to_A_dict[position] = snp
						distribution_file.write(snp.rsID + "\t" + str(position) + "\t" + snp.A + "\t" + snp.B + "\t" + "A" + "\n")
						update_max_allele_number(max_allele_number, A_dict)
						hifi_pure_corrected_file.write(snp.rsID + "\t" + str(position) + "\t" + max_allele + "\n")			
					if max_allele == snp.B and max_allele != snp.A:
						same_to_B_dict[position] = snp
						distribution_file.write(snp.rsID + "\t" + str(position) + "\t" + snp.A + "\t" + snp.B + "\t" + "B" + "\n")					
						update_max_allele_number(max_allele_number, B_dict)				
						hifi_pure_corrected_file.write(snp.rsID + "\t" + str(position) + "\t" + snp.A + "\n")				
					if max_allele == snp.B and max_allele == snp.A:
						same_to_AB_dict[position] = snp
						distribution_file.write(snp.rsID + "\t" + str(position) + "\t" + snp.A + "\t" + snp.B + "\t" + "" + "\n")					
						hifi_pure_corrected_file.write(snp.rsID + "\t" + str(position) + "\t" + max_allele + "\n")
				elif snp.A == "X" or snp.B == "X":	# keep these reads too
					same_to_X_dict[position] = snp
					update_max_allele_number(max_allele_number, X_dict)
				else:
					same_to_N_dict[position] = snp
					update_max_allele_number(max_allele_number, not_ABX_dict)									
			else:	# not in genotype, sequencing error???
				not_in_genotype_dict[position] = snp
				update_max_allele_number(max_allele_number, not_in_genotype_dict_num)

hifi_pure_corrected_file.close()
"""
title_info = "rsID \t phys_position \t std_A \t std_B \t max_allele \t max_allele_number \t max_allele_percentage"
file_name = sam_file_name + "_A.txt"
output_ABX_files(file_name, title_info, same_to_A_dict)
file_name = sam_file_name + "_B.txt"
output_ABX_files(file_name, title_info, same_to_B_dict)
file_name = sam_file_name + "_AB.txt"
output_ABX_files(file_name, title_info, same_to_AB_dict)
file_name = sam_file_name + "_X.txt"
output_ABX_files(file_name, title_info, same_to_X_dict)
file_name = sam_file_name + "_notABX.txt"
output_ABX_files(file_name, title_info, same_to_N_dict)
file_name = sam_file_name + "_notInGenotype.txt"
output_ABX_files(file_name, title_info, not_in_genotype_dict)
"""

called_seed_file_name = sam_file_name + "_called_seed.txt"
title_info = title_haplotype
output_seed_file(called_seed_file_name, title_info, called_seed_dict)

called_seed_total = len(called_seed_dict)
same_to_A = len(same_to_A_dict)
same_to_B = len(same_to_B_dict)
same_to_AB = len(same_to_AB_dict)
X_AB = len(same_to_X_dict)
not_ABX = len(same_to_N_dict)
not_in_genotype = len(not_in_genotype_dict)

A_in_hetero = format((float(same_to_A)/float(same_to_A + same_to_B))*100, "0.2f")
B_in_hetero = format((float(same_to_B)/float(same_to_A + same_to_B))*100, "0.2f")
correct_rate_in_hetero = A_in_hetero if A_in_hetero >= B_in_hetero else B_in_hetero

A_in_all = format((float(same_to_A)/float(called_seed_total))*100, "0.2f")
B_in_all = format((float(same_to_B)/float(called_seed_total))*100, "0.2f")
correct_rate_in_all_seed = A_in_all if A_in_all >= B_in_all else B_in_all

hete_seed_rate = format(float(same_to_A + same_to_B)/float(called_seed_total)*100, "0.2f")				

print "called_seed_dict", called_seed_total
print "same_to_A", same_to_A
print "same_to_B", same_to_B
print "same_to_AB: ", same_to_AB
print "X_AB", X_AB
print "not_ABX", not_ABX

print "A_in_hetero", A_in_hetero, "%"
print "B_in_hetero", B_in_hetero, "%"
print "correct_rate_in_hetero: ", correct_rate_in_hetero
print "correct_rate_in_all_seed: ", correct_rate_in_all_seed
print "(A+B)/seed: ", hete_seed_rate
print "not_in_genotype", not_in_genotype
print "pure_total (include alleles pure but may not in genotype)", pure_total

print >> data_record_file, "haplotype file: ", haplotype_file
print >> data_record_file, "haplotye list size is:", len(hap_std_dict)
print >> data_record_file, "genotype_file file: ", genotype_file
print >> data_record_file, "sam file: ", sam_file
print >> data_record_file, "reads list size is: ", total_reads_num

print >> data_record_file, "called_seed_dict: ", called_seed_total
print >> data_record_file, "same_to_A: ", same_to_A
print >> data_record_file, "same_to_B: ", same_to_B
print >> data_record_file, "same_to_AB: ", same_to_AB
print >> data_record_file, "X_AB: ", X_AB
print >> data_record_file, "not_ABX: ", not_ABX
print >> data_record_file, "A_in_hetero: ", A_in_hetero, "%"
print >> data_record_file, "B_in_hetero: ", B_in_hetero, "%"
print >> data_record_file, "correct_rate_in_hetero: ", correct_rate_in_hetero
print >> data_record_file, "correct_rate_in_all_seed: ", correct_rate_in_all_seed
print >> data_record_file, "hete_seed_rate: ", hete_seed_rate
print >> data_record_file, "not_in_genotype: ", not_in_genotype
print >> data_record_file, "pure_total (include alleles pure but may not in genotype): ", pure_total

print >> data_record_file, "all", "sam_file", "chr", "depth_threshold", "same_to_A", "same_to_B", "correct_rate_in_hetero", "same_to_AB", "X_AB", "not_ABX", "called_seed_total", "correct_rate_in_all_seed", "hete_seed_rate", "pure_total", "not_in_genotype"
print >> data_record_file, "data", sam_file, chr_name, depth_threshold, same_to_A, same_to_B, correct_rate_in_hetero, same_to_AB, X_AB, not_ABX, called_seed_total, correct_rate_in_all_seed, hete_seed_rate, pure_total, not_in_genotype

# correct errors in seed, add homo from genotype to seed. for testing purpose
geno_sorted_list = sort_dict_by_key(geno_dict)

hifi_pure_corrected_with_homo_file_name = sam_file_name + "_hifi_homo.txt"
hifi_pure_corrected_with_homo_file = open(currentPath + hifi_pure_corrected_with_homo_file_name, "w")		
hifi_pure_corrected_with_homo_file.write(title_haplotype + "\n")


group_tuple = group_seed(geno_dict, geno_dict)
geno_homo_dict = group_tuple[0]
geno_hetero_dict = group_tuple[1]
print len(geno_homo_dict)
print len(geno_hetero_dict)

group_tuple = group_seed(called_seed_dict, geno_dict)
geno_homo_dict = group_tuple[0]
geno_hetero_dict = group_tuple[1]
print len(geno_homo_dict)
print len(geno_hetero_dict)
"""?????"""

for geno_data in geno_sorted_list:
	position = geno_data[0]
	geno = geno_data[1][2]
	if geno[0] == geno[1]:
		hifi_pure_corrected_with_homo_file.write(geno_data[1][0] + "\t" + str(position) + "\t" + geno[0] + "\n")	
	elif position in hap_std_dict:
		snp = hap_std_dict[position]
		if len(snp.covered_reads_list) > depth_threshold:	
			max_allele = keywithmaxval(snp.allele_dict)
			max_allele_number = snp.allele_dict[max_allele]	
			pure = False
			if float(max_allele_number)/float(len(snp.covered_reads_list)) >= max_allele_percentage_threshold:
				pure = True

			if max_allele_number > depth_threshold:
				if max_allele == snp.A or max_allele == snp.B:	#check genotype
					if max_allele == snp.A and max_allele != snp.B:							
						hifi_pure_corrected_with_homo_file.write(snp.rsID + "\t" + str(position) + "\t" + max_allele + "\n")				
					if max_allele == snp.B and max_allele != snp.A:								
						#hifi_pure_corrected_with_homo_file.write(snp.rsID + "\t" + str(position) + "\t" + snp.A + "\n")		# error_corrected
						hifi_pure_corrected_with_homo_file.write(snp.rsID + "\t" + str(position) + "\t" + max_allele + "\n")			
				elif snp.A == "X" or snp.B == "X":	# keep these reads too
					hifi_pure_corrected_with_homo_file.write(snp.rsID + "\t" + str(position) + "\t" + max_allele + "\n")			
				else:
					pass	
					
						
	
hifi_pure_corrected_with_homo_file.close()		

# caulculate the coverage distribution for each snp
coverage_distribution_file = open(currentPath + sam_file_name + "_coverage_distribution.txt", "w")
coverage_distribution_file.write("coverage depth \t A \t B \t X \t not_ABX \t percentage \n")

for key, value in A_dict.iteritems():
	B_value = 0
	X_value = 0
	notABX_value = 0
	coverage_distribution_file.write(str(key) + "\t" + str(value) + "\t")
	if key in B_dict: 
		coverage_distribution_file.write(str(B_dict[key]) + "\t" )
		B_value = B_dict[key]
	else: coverage_distribution_file.write("\t" )
	if key in X_dict: 
		coverage_distribution_file.write(str(X_dict[key]) + "\t" )
		X_value = X_dict[key]
	else: coverage_distribution_file.write("\t" )
	if key in not_ABX_dict: 
		coverage_distribution_file.write(str(not_ABX_dict[key]) + "\t" )
		notABX_value = not_ABX_dict[key]
	else: coverage_distribution_file.write("\t" )
	coverage_distribution_file.write(format(float(value)/float(value+B_value+X_value+notABX_value), "0.3f") + "\n" )
coverage_distribution_file.close()



end = time.time()
run_time = str(format((end - start), "0.3f"))
print "run time is: " + run_time + "s"


data_record_file.write("run time is: " + run_time + "s \n")
data_record_file.close()

seed_std_compare(called_seed_file_name, chr_name)
seed_std_compare(hifi_pure_corrected_with_homo_file_name, chr_name)





# compare with ori hap file			
"""
unchanged_snp_file_name = sam_file_name + "_unchanged_snp.txt"
unchanged_snp_file = open(currentPath + unchanged_snp_file_name, "w")		
unchanged_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

changed_snp_file_name = sam_file_name + "_changed_snp.txt"
changed_snp_file = open(currentPath + changed_snp_file_name, "w")		
changed_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

# need to update for A or B 	
for snp_data in hap_std_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) > depth_threshold:	
		consisitent = True
		for reads in snp.covered_reads_list:
			if snp.A != reads.covered_snp:
			#if snp.B != reads.covered_snp:
				consisitent = False
		if consisitent:
			unchanged_snp_file.write(str(snp.position) + "\t" + snp.A + "\t" + str(len(snp.covered_reads_list)) + "\t" +"\n")
			#unchanged_snp_file.write(str(snp.position) + "\t" + snp.B + "\t" + str(len(snp.covered_reads_list)) + "\t" +"\n")
		if not consisitent:
			changed_snp_file.write(str(snp.position) + "\t" + snp.A + "\t" + str(len(snp.covered_reads_list)) +"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+ "\n")
			#changed_snp_file.write(str(snp.position) + "\t" + snp.B + "\t" + str(len(snp.covered_reads_list)) +"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+ "\n")

unchanged_snp_file.close()
changed_snp_file.close()
"""


# check symmetric and asymmetric
"""
same_A_file = open(currentPath + sam_file_name + "_same_A.txt", "w")
same_A_file.write("rsID \t phys_position \t NA12878_A \t selected_base \t reads	\n")

pair_A_file = open(currentPath + sam_file_name + "_pair_A.txt", "w")
pair_A_file.write("rsID \t phys_position \t NA12878_A \t selected_base \t reads	\n")

different_A_file = open(currentPath + sam_file_name + "_different_A.txt", "w")
different_A_file.write("rsID \t phys_position \t NA12878_A \t selected_base \t reads \n")

same_A = 0
pair_A = 0
different_A = 0


for snp_data in hap_std_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) > depth_threshold:	
		max_allele = keywithmaxval(snp.allele_dict)
		max_allele_number = snp.allele_dict[max_allele]

		pure = True
		for base in base_list:
			if base != max_allele:
				if snp.allele_dict[base] != 0:
					pure = False
		if pure:
			pure_total += 1
			hap_base = snp.A
			if hap_base == "A":
				if max_allele == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_allele == 'T':
					pair_A += 1
					pair_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B +"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						pair_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				else:
					different_A += 1
					different_A_file.positionID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						different_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")		
			elif hap_base == "T":
				if max_allele == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_allele == 'A':
					pair_A += 1
					pair_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list)) + "\n")
					for reads in snp.covered_reads_list:
						pair_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				else:
					different_A += 1
					different_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						different_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
			elif hap_base == "C":
				if max_allele == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_allele == 'G':
					pair_A += 1
					pair_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list)) + "\n")
					for reads in snp.covered_reads_list:
						pair_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				else:
					different_A += 1
					different_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B +"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						different_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
			elif hap_base == "G":
				if max_allele == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A  +"\t"+ snp.B+"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_allele == 'C':
					pair_A += 1
					pair_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B +"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						pair_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				else:
					different_A += 1
					different_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B + "\t" + str(len(snp.covered_reads_list)) + "\n")
					for reads in snp.covered_reads_list:
						different_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")	
			#else:
				#print hap_base
			


print "same A", same_A
print "symmetric to A is", pair_A
print "asymmetric to A is", different_A

data_record_file.write("same A is: " + str(same_A) + "\n")
data_record_file.write("symmetric to A is: " + str(pair_A) + "\n")
data_record_file.write("asymmetric to A is: " + str(different_A) + "\n")
"""
