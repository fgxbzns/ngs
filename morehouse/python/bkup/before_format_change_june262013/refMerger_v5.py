#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# prepare geno, hap, ref files for hifi. 
# This is to remove hap of father and mother from the ref file. 

# not finished. 

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

file_path = "/home/guoxing/disk2/solid/common_files/"

currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--haplotype", type="string", dest="haplotypeFile",help = "Input File Name", default="null")
parser.add_option("-s", "--genotype", type="string", dest="samFile",help = "Input File Name", default="null")
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="chr11")

(options, args) = parser.parse_args()

haplotype_file = options.haplotypeFile
#genotype_file = options.samFile
chr_name = options.chrName

#chr_name = "chr11"

snp_dic = {}
tag_info = ""
total_snp_number = 0
common_snp_number = 0


for infile in glob.glob(os.path.join(file_path,"*"+chr_name+"_???.phased")):	# add * for chrX
	file_name = infile[(infile.find("hapmap3")):].strip()
	inputFile = open(file_path + file_name, "r")
	print file_name
	for line in inputFile:
		if line.startswith("rsID"):
			if tag_info == "":
				tag_info = line.strip()
			else:
				elements = line.strip().split()
				first_tag = elements[2].strip()
				tag_info += "\t" + line[line.find(first_tag):].strip()
				
				
		elements = line.strip().split()
		try:
			ID_index = elements.index(person_ID)
			print ID_index
		except:
			print person_ID, "not found in file"
			break
				
				
				
				
						
		if not line.startswith("rsID"):
			elements = line.strip().split()
			position = elements[1].strip()				
			try:
				position = int(position)
				if position not in snp_dic:
					snp_dic[position] = line.strip()
				else:
					elements = line.strip().split()
					first_base = elements[2].strip()
					snp_dic[position] += "\t" + line[line.find(first_base):].strip()
			except ValueError:
				print file_name, position	
				print line
			
	inputFile.close()

outputFile = open(currentPath + "refHaplos.txt", "w")	# for hifi
outputFile.write(tag_info + "\n")

total_element_number = len(tag_info.strip().split())
print "total_element_number", total_element_number

# dict can not be sorted. need to convert it to a list first
snp_sorted_list = [x for x in snp_dic.iteritems()] 

snp_sorted_list.sort(key=lambda x: x[0]) # sort by key
#snp_sorted_list.sort(key=lambda x: x[1]) # sort by value

# keep only the snps that are available in ref data
#genotype_input_file_name = "genotype_NA12878_chr6.txt"	# for simulation data
#genotype_input_file_name = "genotype_NA10847_chr11.txt"	# for solid song_5_2 chr11
genotype_input_file_name = "genotype_NA10847_"+chr_name+".txt"	
#genotype_input_file_name = "genotype_NA12878_"+chr_name+".txt"	

print genotype_input_file_name
genotype_input_file = open(file_path + genotype_input_file_name, "r")

genotype_output_file_name = "genotype.txt"
genotype_output_file = open(currentPath + genotype_output_file_name, "w")

title_genotype = ""

genotype_dict = {}

for line in genotype_input_file:
	if line.startswith("rsID"):
		title_genotype = line.strip()
	if not line.startswith("rsID"):
		elements = line.strip().split()
		position = elements[1].strip()							
		try:
			position = int(position)
			genotype_dict[position] = line.strip()
		except ValueError:
			print "error in ", line

genotype_output_file.write(title_genotype + "\n")

# keep only the snps that are available in ref data
#haplotype_input_file_name = "haplotype_NA12878.txt"
#haplotype_input_file_name = "NA12878_hap_new.txt"
#haplotype_input_file_name = "NA12878_hg18ch6_A_4.0x_0.01er_hifi_pure.txt"


if haplotype_file == "null":
	print "Please input the haplotype file name"


print haplotype_file
haplotype_input_file = open(currentPath + haplotype_file, "r")

#haplotype_output_file_name = haplotype_file[:haplotype_file.find(".")]+"_refed.txt"	# for making hap_refed file
haplotype_output_file_name = "haplotype.txt" # for hifi
haplotype_output_file = open(currentPath + haplotype_output_file_name, "w")
#haplotype_output_file.write("rsID \t phys_position \t NA12878_A \t NA12878_B \n")

title_haplotype = ""

haplotype_dict = {}

for line in haplotype_input_file:
	if line.startswith("rsID"):
		title_haplotype = line.strip()
	if not line.startswith("rsID"):
		elements = line.strip().split()
		position = elements[1].strip()							
		try:
			position = int(position)
			haplotype_dict[position] = line.strip()
		except ValueError:
			print "error in ", line

haplotype_output_file.write(title_haplotype + "\n")

haplotype_sorted_list = [x for x in haplotype_dict.iteritems()] 
haplotype_sorted_list.sort(key=lambda x: x[0]) # sort by key

				
last_haplotype_position = haplotype_sorted_list[-1][0]
print "last_haplotype_position", last_haplotype_position

# if the base_info contains the same amount of elements as the tag_info does, it is common	
for snp in snp_sorted_list:
	# assume last_haplotype_position is the smallest in all three. Make sure the last position in three files is the same. Limited by current version of hifi.
	if len(snp[1].strip().split()) == total_element_number and int(snp[0]) <= int(last_haplotype_position):  
		outputFile.write(snp[1] + "\n")
		if snp[0] in genotype_dict:
			genotype_output_file.write(genotype_dict[snp[0]] + "\n")
		if snp[0] in haplotype_dict:
			haplotype_output_file.write(haplotype_dict[snp[0]] + "\n")
		common_snp_number += 1

print "total snp number in ref", len(snp_dic)
print "common_snp_number", common_snp_number
#print tag_info

outputFile.close()
genotype_input_file.close()
genotype_output_file.close()

haplotype_input_file.close()
haplotype_output_file.close()

