#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for solid data May 09 2013

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

program_path = "/home/guoxing/tool/morehouse"
file_path = "/home/guoxing/disk2/solid/common_files/"

currentPath = os.getcwd() + '/'

quality_score_dict = { '!':0, '\"':1, '#':2, '$':3, '%':4, '&':5, '\'':6, '(':7, 
					')':8, '*':9, '+':10, ',':11, '-':12, '.':13 }

#print quality_score_dict['\'']

def getTotalBaseNum(fileName):
	totalBase = 0
	f = open(currentPath+fileName, "r")
	for line in f:
		if not line.startswith('>'):
			totalBase += len(line.strip())
	return totalBase
	f.close()

# A from Father, B from Mother
class snp:
	def __init__(me, rsID, position, A, B, covered_reads_list, allele_dict):
		me.rsID = rsID
		me.position = position
		me.A = A
		me.B = B
		me.covered_reads_list = covered_reads_list
		me.allele_dict = allele_dict

# class to store reads from sam file
class read:
	def __init__(me, qName, flag, rName, start_position, read_sequence, quality_score_sequence, read_length, covered_snp):
		me.qName = qName
		me.flag = flag
		me.rName = rName
		me.start_position = start_position
		me.read_sequence = read_sequence
		me.quality_score_sequence = quality_score_sequence
		me.read_length = read_length
		me.covered_snp = covered_snp
	
def keywithmaxval(dict):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value """  
     v=list(dict.values())
     k=list(dict.keys())
     return k[v.index(max(v))]

# start time
start = time.time()		
		
# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--haplotype", type="string", dest="haplotypeFile",help = "Input File Name", default="null")
parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
parser.add_option("-d", "--threshold", type="string", dest="threshold",help = "Input the depth threshold", default="3")
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="chr11")



(options, args) = parser.parse_args()

haplotype_file = options.haplotypeFile
sam_file = options.samFile
depth_threshold = int(options.threshold)
chr_name = options.chrName

#if haplotype_file == "null":
	#print "Please input the haplotype file name"

haplotype_file = "NA12878_hap_new_refed.txt"	# simulation data chr6
#haplotype_file = "ASW_chr11_child_hap_refed.txt"	# solid song_5_2 chr11
haplotype_file = "ASW_"+chr_name+"_child_hap_refed.txt"	# for all

#genotype_file = "genotype_NA10847_chr11.txt"	# solid song_5_2 chr11
genotype_file = "genotype_NA10847_"+chr_name+".txt"	# for all

if sam_file == "null":
	print "Please input the sam file name"
	

sam_file_name = sam_file[:(len(sam_file)-4)] + "_" + str(depth_threshold)

print "haplotype file: ", haplotype_file
print "genotype_file file: ", genotype_file
print "sam file: ", sam_file_name


sam_path = currentPath

inputFile_hap = open(file_path + haplotype_file, "r")
inputFile_geno = open(file_path + genotype_file, "r")
inputFile_sam = open(currentPath + sam_file, "r")

outputFile_reads = open(currentPath + sam_file_name + "_reads.txt", "w")
outputFile_reads.write("SNP position \t Depth \n")

outputFile_allele = open(currentPath + sam_file_name+"_allele.txt", "w")
outputFile_allele.write("Chromosome \t position \t Total Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

data_record_file_name = sam_file_name + "_data_record.txt"
data_record_file = open(currentPath + data_record_file_name, "w")

# get haplotype info
snp_dict={}
title_haplotype = ""

for line in inputFile_hap:
	if line.startswith("rsID"):
		title_haplotype = line.strip()
	if not line.startswith("rsID"):
		elements = line.strip().split()
		rsID = elements[0].strip()
		position = int(elements[1].strip())
		A = elements[2].strip()
		B = elements[3].strip()
		covered_reads_list = []
		allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
		snp_dict[position] = snp(rsID, position, A, B, covered_reads_list, allele_dict)

inputFile_hap.close()

# get genotype info
geno_dict = {}

for line in inputFile_geno:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		#rsID = elements[0].strip()
		position = int(elements[1].strip())
		alleles = elements[2].strip()
		geno_dict[position] = alleles

inputFile_geno.close()

reads_list=[]
insert_size = 0
"""
hap_homo_file = open(currentPath + "hap_homo.txt", "w")
hap_hete_file = open(currentPath + "hap_hete.txt", "w")

for snp in snp_list:
	if snp.A == snp.B:
		hap_homo_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\n")
	else:
		hap_hete_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\n")
	
hap_homo_file.close()
hap_hete_file.close()
"""
sam_line_first = inputFile_sam.readline() # the first read line in a pair


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
		if True:
			if True:						
				if True: 
					# first read
					qName_first = elements_first[0].strip()
					flag_first = elements_first[1].strip()
					start_position_first = int(elements_first[3].strip())
					read_sequence_first = elements_first[9].strip()
					read_length_first = len(read_sequence_first)
					quality_score_sequence_first = elements_first[10].strip()
					i = 0
					while i < read_length_first:
						if (start_position_first+i) in snp_dict:							
							covered_snp = read_sequence_first[i]			# ith position is the covered snp
							quality_score_symbol = quality_score_sequence_first[i]
							if (rName_first == chr_name) and (not covered_snp == 'N') and (not quality_score_symbol in quality_score_dict):	# check quality_score
								covered_snp_total_number += 1
								snp_dict[start_position_first+i].covered_reads_list.append(read(qName_first, flag_first, rName_first, \
								start_position_first, read_sequence_first, quality_score_sequence_first, read_length_first, covered_snp))					
								if covered_snp == 'A':
									snp_dict[start_position_first+i].allele_dict['A'] += 1
								elif covered_snp == 'T':
									snp_dict[start_position_first+i].allele_dict['T'] += 1
								elif covered_snp == 'C':
									snp_dict[start_position_first+i].allele_dict['C'] += 1
								elif covered_snp == 'G':
									snp_dict[start_position_first+i].allele_dict['G'] += 1			
						i += 1
			
	sam_line_first = inputFile_sam.readline()

inputFile_sam.close()

print "haplotye list size is: ", len(snp_dict)	
print "total_reads_num", total_reads_num

snp_sorted_list = [x for x in snp_dict.iteritems()] 
snp_sorted_list.sort(key=lambda x: x[0]) # sort by key

# output origial haplotype homo and hete alleles
"""
hap_homo_file = open(currentPath + sam_file_name + "_orihap_homoSNP.txt", "w")
hap_hete_file = open(currentPath + sam_file_name + "_orihap_heteSNP.txt", "w")

#rs_different_number = 0
oriHap_homoSnp_number = 0
"""
for snp_data in snp_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) > depth_threshold:
		outputFile_allele.write(chr_name+"\t"+str(snp.position)+"\t"+str(len(snp.covered_reads_list))+"\t"+str(snp.allele_dict['A'])	\
								+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+"\n")
		outputFile_reads.write("@_" + snp.rsID + "\t" + str(snp.position) + "\t" + snp.A  + "\t" + snp.B  + "\t" + str(len(snp.covered_reads_list)) + "\n")
		for reads in snp.covered_reads_list:
			outputFile_reads.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" \
									+ reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
		"""
		if snp.A == snp.B:
			oriHap_homoSnp_number += 1
			rs_different = False
			for reads in snp.covered_reads_list:
				if snp.A != reads.covered_snp:
					rs_different = True
			if rs_different:
				#rs_different_number += 1
				hap_homo_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\n")
				for reads in snp.covered_reads_list:
					if snp.A != reads.covered_snp:
						hap_homo_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) \
						+ "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
		else:
			hap_hete_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\n")
			for reads in snp.covered_reads_list:
				hap_hete_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" \
				+ reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
		"""
"""
hap_homo_file.close()
hap_hete_file.close()
"""

# compare with ori hap file			
"""
unchanged_snp_file_name = sam_file_name + "_unchanged_snp.txt"
unchanged_snp_file = open(currentPath + unchanged_snp_file_name, "w")		
unchanged_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

changed_snp_file_name = sam_file_name + "_changed_snp.txt"
changed_snp_file = open(currentPath + changed_snp_file_name, "w")		
changed_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

# need to update for A or B 	
for snp_data in snp_sorted_list:
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


# prepare the haplotype file for hifi

hifi_pure_file_name = sam_file_name + "_hifi.txt"
hifi_pure_file = open(currentPath + hifi_pure_file_name, "w")		
hifi_pure_file.write(title_haplotype)

"""
hifi_max_file_name = sam_file_name + "_hifi_max.txt"
hifi_max_file = open(currentPath + hifi_max_file_name, "w")		
hifi_max_file.write("rsID \t phys_position \t NA12878_A \t	NA12878_B \n")
"""

base_list = ["A", "T", "C", "G"]

"""
hete_A_max_file = open(currentPath + sam_file_name + "_hete_max_A.txt", "w")
hete_A_max_file.write("rsID \t phys_position \t NA12878_A \t	selected_base \n")

hete_B_max_file = open(currentPath + sam_file_name + "_hete_max_B.txt", "w")
hete_B_max_file.write("rsID \t phys_position \t NA12878_B \t	selected_base \n")
"""
# hete allele, pure means all the variation calls get the same base
hete_A_pure_file = open(currentPath + sam_file_name + "_A.txt", "w")
hete_A_pure_file.write("rsID \t phys_position \t A \t	selected_base \n")

hete_B_pure_file = open(currentPath + sam_file_name + "_B.txt", "w")
hete_B_pure_file.write("rsID \t phys_position \t A \t	selected_base \n")

hete_X_pure_file = open(currentPath + sam_file_name + "_X.txt", "w")
hete_X_pure_file.write("rsID \t phys_position \t A \t NA12878_B \t selected_base \n")

hete_notABX_pure_file = open(currentPath + sam_file_name + "_notABX.txt", "w")
hete_notABX_pure_file.write("rsID \t phys_position \t A \t	selected_base \n")

not_in_genotype_file = open(currentPath + sam_file_name + "_notInGenotype.txt", "w")
not_in_genotype_file.write("rsID \t phys_position \t A \t B	\t selected_base \t depth \ genotype \n")

distribution_file = open(currentPath + sam_file_name + "_distribution.txt", "w")
distribution_file.write("rsID \t phys_position \t snp.A	\t snp.B \t distribution \n")

max_hete = 0
pure_hete = 0

common_total = 0
same_to_A = 0
same_to_B = 0
same_to_AB = 0
X_AB = 0
not_ABX = 0
pure_total = 0
hifi_seed = 0
not_in_genotype = 0
A_dict = {}
B_dict = {}
X_dict = {}
not_ABX_dict = {}
not_in_genotype_dict = {}

for snp_data in snp_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) > depth_threshold:	
		max_base = keywithmaxval(snp.allele_dict)
		max_value = snp.allele_dict[max_base]
		#hifi_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\t" + str(max_value) + "\n")
		#hifi_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")
		"""
		if snp.A != snp.B:		#hete
			max_hete += 1
			if max_base == snp.A:
				hete_A_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + max_base + "\t" + str(max_value) + "\n")
			if max_base == snp.B:
				hete_B_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.B + "\t" + max_base + "\t" + str(max_value) + "\n")
		"""
		pure = True
		for base in base_list:
			if base != max_base:
				if snp.allele_dict[base] != 0:
					pure = False
		if pure and max_value > depth_threshold:
			pure_total += 1
			# check genotype to remove called base that does not in genotype
			if max_base in geno_dict[snp.position]:	
				hifi_seed += 1	
				hifi_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")	
				
				if max_base == snp.A or max_base == snp.B:	#check genotype
					#hifi_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")
					common_total += 1
					if max_base == snp.A and max_base != snp.B:
						same_to_A += 1
						#hifi_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")
						hete_A_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\t" +  max_base + "\t" + str(max_value) + "\n")
						distribution_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\t" + "A" + "\n")
						if max_value not in A_dict:
							A_dict[max_value] = 1
						else:
							A_dict[max_value] += 1					
					if max_base == snp.B and max_base != snp.A:
						same_to_B += 1
						hete_B_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\t" + max_base + "\t" + str(max_value) + "\n")
						distribution_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\t" + "B" + "\n")					
						if max_value not in B_dict:
							B_dict[max_value] = 1
						else:
							B_dict[max_value] += 1				
						for reads in snp.covered_reads_list:
							hete_B_pure_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) \
							+ "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")					
					if max_base == snp.B and max_base == snp.A:
						same_to_AB += 1
						#hifi_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")
						distribution_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\t" + "" + "\n")
				elif snp.A == "X" or snp.B == "X":	# keep these reads too
					X_AB += 1
					#hifi_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")
					hete_X_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t"+ snp.B + "\t" + max_base + "\t" + str(max_value) + "\n")
					if max_value not in X_dict:
						X_dict[max_value] = 1
					else:
						X_dict[max_value] += 1
				else:
					not_ABX += 1
					hete_notABX_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t"+ snp.B + "\t" + max_base + "\t" + str(max_value) +"\n")
					if max_value not in not_ABX_dict:
						not_ABX_dict[max_value] = 1
					else:
						not_ABX_dict[max_value] += 1	
				
					
			else:	# not in genotype, sequencing error???
				not_in_genotype += 1
				
				not_in_genotype_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t"+ snp.B + "\t" + max_base + "\t" + str(max_value) + "\t" + geno_dict[snp.position] + "\n")
				if max_value not in not_in_genotype_dict:
					not_in_genotype_dict[max_value] = 1
				else:
					not_in_genotype_dict[max_value] += 1
				
				for reads in snp.covered_reads_list:
						not_in_genotype_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
		

correct_rate = format((float(same_to_A)/float(same_to_A + same_to_B))*100, "0.2f")
A_seed_rate = format(float(same_to_A)/float(hifi_seed)*100, "0.2f")
hete_seed_rate = format(float(same_to_A + same_to_B)/float(hifi_seed)*100, "0.2f")				

print "common_total", common_total
print "same_to_A", same_to_A
print "same_to_B", same_to_B
print "correct rate: ", correct_rate
print "homo: ", same_to_AB
print "X_AB", X_AB
print "not_ABX", not_ABX
print "hifi_seed", hifi_seed
print "A/seed: ", A_seed_rate
print "(A+B)/seed: ", hete_seed_rate
print "not_in_genotype", not_in_genotype
print "pure_total", pure_total


#data_record_file.write("pure and in genotype" + "\n")
#data_record_file.write("total in genotype is: " + str(common_total) + "\n")
print >> data_record_file, "haplotype file: ", haplotype_file
print >> data_record_file, "haplotye list size is:", str(len(snp_dict))
print >> data_record_file, "genotype_file file: ", genotype_file
print >> data_record_file, "sam file: ", sam_file
print >> data_record_file, "reads list size is: ", str(total_reads_num)

data_record_file.write("same to snp.A: " + str(same_to_A) + "\n")
data_record_file.write("same to snp.B: " + str(same_to_B) + "\n")
data_record_file.write("correct rate: " + str(correct_rate) + "\n")
data_record_file.write("homo: " + str(same_to_AB) + "\n")
data_record_file.write("X_AB: " + str(X_AB) + "\n")
data_record_file.write("not_ABX: " + str(not_ABX) + "\n")
data_record_file.write("total seed for hifi: " + str(hifi_seed) + "\n")


print >> data_record_file, "A/seed: ", A_seed_rate
print >> data_record_file, "(A+B)/seed: ", hete_seed_rate

data_record_file.write("not_in_genotype: " + str(not_in_genotype) + "\n")
data_record_file.write("pure_total (include alleles pure but may not in genotype): " + str(pure_total) + "\n")

#data_record_file.write("All \t A \t B \t A/(A+B) \t \homo \t X \t N \t hifi_seed \t A/seed \t (A+B)/seed \t pure_total \t not_ingeno \n")
print >> data_record_file, "all", "sam_file", "chr", "depth_threshold", "same_to_A", "same_to_B", "correct_rate", "same_to_AB", "X_AB", "not_ABX", "hifi_seed", "A_seed_rate", "hete_seed_rate", "pure_total", "not_in_genotype"
print >> data_record_file, "data", sam_file, chr_name, depth_threshold, same_to_A, same_to_B, correct_rate, same_to_AB, X_AB, not_ABX, hifi_seed, A_seed_rate, hete_seed_rate, pure_total, not_in_genotype

hifi_pure_file.close()
"""
#hifi_max_file.close()
#print "max_hete", max_hete
#print "pure_hete", pure_hete

hete_A_max_file.close()
hete_B_max_file.close()
"""
hete_A_pure_file.close()
hete_B_pure_file.close()
hete_X_pure_file.close()
hete_notABX_pure_file.close()

# caulculate the coverage distribution for each snp
coverage_distribution_file = open(currentPath + sam_file_name + "_coverage_distribution.txt", "w")
coverage_distribution_file.write("coverage depth \t A \t B \t X \t not_ABX \t percentage \n")

print "A: ", A_dict			
print "B: ", B_dict
print "X: ", X_dict
print "not_ABX: ", not_ABX_dict
print "not_in_genotype_dict: ", not_in_genotype_dict


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


for snp_data in snp_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) > depth_threshold:	
		max_base = keywithmaxval(snp.allele_dict)
		max_value = snp.allele_dict[max_base]

		pure = True
		for base in base_list:
			if base != max_base:
				if snp.allele_dict[base] != 0:
					pure = False
		if pure:
			pure_total += 1
			hap_base = snp.A
			if hap_base == "A":
				if max_base == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_base == 'T':
					pair_A += 1
					pair_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B +"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						pair_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				else:
					different_A += 1
					different_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						different_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")		
			elif hap_base == "T":
				if max_base == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_base == 'A':
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
				if max_base == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A +"\t"+ snp.B+ "\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_base == 'G':
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
				if max_base == hap_base:
					same_A += 1
					same_A_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A  +"\t"+ snp.B+"\t" + str(len(snp.covered_reads_list))  + "\n")
					for reads in snp.covered_reads_list:
						same_A_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
				elif max_base == 'C':
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

end = time.time()
run_time = str(format((end - start), "0.3f"))
print "run time is: " + run_time + "s"


data_record_file.write("run time is: " + run_time + "s \n")

data_record_file.close()

outputFile_reads.close()
outputFile_allele.close()

