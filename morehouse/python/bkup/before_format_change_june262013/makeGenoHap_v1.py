#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# generate genotype and haplotye for hifi

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

file_path = "/home/guoxing/morehouse_new_hap/fils/"
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--haplotype", type="string", dest="haplotypeFile",help = "Input haplotype File Name", default="null")
parser.add_option("-g", "--genotype", type="string", dest="genotypeFile",help = "Input genotype File Name", default="null")
parser.add_option("-r", "--ref", type="string", dest="refFile",help = "Input ref File Name", default="null")

(options, args) = parser.parse_args()

haplotype_file = options.haplotypeFile
genotype_file = options.genotypeFile
ref_file = options.refFile

if haplotype_file == "null":
	print "Please input the haplotype file name"

print haplotype_file

if genotype_file == "null":
	print "Please input the genotype file name"

genotype_input_file_name = "genotype_NA12878.txt"
print genotype_input_file_name
	
if ref_file == "null":
	print "Please input the ref file name"

ref_file_name = "refHaplos.txt"
print ref_file_name

snp_dic = {}		# store info from ref
tag_info = ""
total_snp_number = 0
common_snp_number = 0

ref_input_file = open(currentPath + ref_file_name, "r")

for line in ref_input_file:
		if line.startswith("rsID"):
			if tag_info == "":
				tag_info = line.strip()		
		if not line.startswith("rsID"):
			elements = line.strip().split()
			position = elements[1].strip()				
			try:
				position = int(position)
				snp_dic[position] = line.strip()
			except ValueError:
				print line
			
ref_input_file.close()

# dict can not be sorted. need to convert it to a list first
snp_sorted_list = [x for x in snp_dic.iteritems()] 

snp_sorted_list.sort(key=lambda x: x[0]) # sort by key
#snp_sorted_list.sort(key=lambda x: x[1]) # sort by value

# keep only the snps that are available in ref data
genotype_input_file = open(currentPath + genotype_input_file_name, "r")

genotype_output_file_name = "genotype.txt"
genotype_output_file = open(currentPath + genotype_output_file_name, "w")
genotype_output_file.write("rsID \t phys_position \t NA12878 \n")

genotype_dict = {}

for line in genotype_input_file:
		if not line.startswith("rsID"):
			elements = line.strip().split()
			position = elements[1].strip()							
			try:
				position = int(position)
				genotype_dict[position] = line.strip()
			except ValueError:
				print "error in ", line

# keep only the snps that are available in ref data
haplotype_input_file = open(currentPath + haplotype_file, "r")

haplotype_output_file_name = "haplotype.txt"
haplotype_output_file = open(currentPath + haplotype_output_file_name, "w")
haplotype_output_file.write("rsID \t phys_position \t NA12878_A \t NA12878_B \n")

haplotype_dict = {}

for line in haplotype_input_file:
		if not line.startswith("rsID"):
			elements = line.strip().split()
			position = elements[1].strip()							
			try:
				position = int(position)
				haplotype_dict[position] = line.strip()
			except ValueError:
				print "error in ", line

haplotype_sorted_list = [x for x in haplotype_dict.iteritems()] 
haplotype_sorted_list.sort(key=lambda x: x[0]) # sort by key

				
last_haplotype_position = haplotype_sorted_list[-1][0]
print "last_haplotype_position", last_haplotype_position

# if the base_info contains the same amount of elements as the tag_info does, it is common	
for snp in snp_sorted_list:
	# assume last_haplotype_position is the smallest in all three. Make sure the last position in three files is the same. Limited by current version of hifi.
	if int(snp[0]) <= int(last_haplotype_position):  
		# outputFile.write(snp[1] + "\n")
		if snp[0] in genotype_dict:
			genotype_output_file.write(genotype_dict[snp[0]] + "\n")
		if snp[0] in haplotype_dict:
			haplotype_output_file.write(haplotype_dict[snp[0]] + "\n")
		common_snp_number += 1

print "total snp number in ref", len(snp_dic)
print "common_snp_number", common_snp_number
#print tag_info

#outputFile.close()
genotype_input_file.close()
genotype_output_file.close()

haplotype_input_file.close()
haplotype_output_file.close()

