#!/usr/bin/python

# location /home/guoxing/tool/morehouse

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

file_path = "/home/guoxing/disk2/solid/common_files/"
NA10847_F12146_M12239_path = file_path + "NA10847_F12146_M12239/"
data_record_path = "/home/guoxing/disk2/solid/common_files/data_record/"

currentPath = os.getcwd() + '/'

snp_hap_ori_dict = {}
snp_hap_hifi_dict = {}
snp_hap_ori_total_number = 0
snp_hap_ori_total_number_withoutXN = 0
snp_hap_hifi_total_number = 0

# for NA10847 only
genotype_father_ID = "NA12146"
genotype_mather_ID = "NA12239"

genotype_father_dict = {}
genotype_mather_dict = {}

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="chr11")
(options, args) = parser.parse_args()
chr_name = options.chrName

genotype_father_file_name = "genotype_" + genotype_father_ID + "_" + chr_name + ".txt" 
inputFile_genotype_father = open(NA10847_F12146_M12239_path + genotype_father_file_name, "r")

genotype_mather_file_name = "genotype_" + genotype_mather_ID + "_" + chr_name + ".txt" 
inputFile_genotype_mather = open(NA10847_F12146_M12239_path + genotype_mather_file_name, "r")

hete_in_father = 0
hete_in_mather = 0
hete_in_child = 0


for line in inputFile_genotype_father:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		position = elements[1].strip()	
		try:
			position = int(position)
			genotype_father_dict[position] = elements[2].strip()
			if elements[2].strip()[0] != elements[2].strip()[1]:
				hete_in_father += 1
		except ValueError:
			print position	
			print line

for line in inputFile_genotype_mather:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		position = elements[1].strip()	
		try:
			position = int(position)
			genotype_mather_dict[position] = elements[2].strip()
			if elements[2].strip()[0] != elements[2].strip()[1]:
				hete_in_mather += 1
		except ValueError:
			print position	
			print line

print "genotype_father_file_name", genotype_father_file_name
print "genotype_mather_file_name", genotype_mather_file_name
print "hete_in_father", hete_in_father
print "hete_in_mather", hete_in_mather


#hap_ori_file_name = "NA12878_hap_new_refed.txt"	# simulation data
hap_ori_file_name = "ASW_"+chr_name+"_child_hap_refed.txt"	
#hap_ori_file_name = "ASW_chr5_child_hap_refed.txt"	# solid song_1 chr5

hap_hifi_file_name = "imputedhaplotype.txt"

inputFile_hap_ori = open(file_path + hap_ori_file_name, "r")
inputFile_hap_hifi = open(currentPath + hap_hifi_file_name, "r")

accuracy_output_file_name = "hifi_accuracy.txt"
accuracy_output_file = open(currentPath + accuracy_output_file_name, "w")

x_number = 0
n_number = 0
a_number = 0
t_number = 0
c_number = 0
g_number = 0

for line in inputFile_hap_ori:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		#if elements[2].strip() != "N" and elements[2].strip() != "X" and elements[3].strip() != "N" and elements[3].strip() != "X":
		if elements[2].strip() != "N" and elements[3].strip() != "N":
			if elements[2].strip() != elements[3].strip():
				hete_in_child += 1
			position = elements[1].strip()	
			try:
				position = int(position)
				snp_hap_ori_dict[position] = line.strip()
			except ValueError:
				print file_name, position	
				print line
			snp_hap_ori_total_number += 1
		"""
		if elements[2].strip() == "N":
			n_number += 1
		if elements[2].strip() == "X":
			x_number += 1
		if elements[2].strip() == "A":
			a_number += 1
		if elements[2].strip() == "T":
			t_number += 1
		if elements[2].strip() == "C":
			c_number += 1
		if elements[2].strip() == "G":
			g_number += 1
	
print "n_number", n_number
print "x_number", x_number
print "a_number", a_number
print "t_number", t_number
print "c_number", c_number
print "g_number", g_number
"""			

print "hete_in_child", hete_in_child

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

print "snp_hap_ori_total_number", snp_hap_ori_total_number
print "snp_hap_hifi_total_number", snp_hap_hifi_total_number

same_position_total_number = 0
different_position_total_number = 0
same_AB_total_number = 0
same_A_total_number = 0
same_B_total_number = 0
not_same_AB_total_number = 0
triple_heterozygous_total_number = 0

difference_output_file_name = "hifi_difference.txt"
difference_output_file = open(currentPath + difference_output_file_name, "w")

for position, line_hifi in snp_hap_hifi_dict.iteritems():
	# check triple heterozygous, the position in father, mother and child are all heterozygous
	
	
	if position in snp_hap_ori_dict:
		triple_heterozygous = False
		if position in genotype_father_dict and position in genotype_mather_dict:
			line_ori = snp_hap_ori_dict[position]
			elements_ori = line_ori.strip().split()
			ori_A = elements_ori[2].strip()
			ori_B = elements_ori[3].strip()
			father_A = genotype_father_dict[position][0]
			father_B = genotype_father_dict[position][1]
			mather_A = genotype_mather_dict[position][0]
			mather_B = genotype_mather_dict[position][1]
			line_ori = snp_hap_ori_dict[position]
			elements_hifi = line_hifi.strip().split()
			hifi_A = elements_hifi[2].strip()
			hifi_B = elements_hifi[3].strip()
			#if hifi_A != hifi_B:
			#	print >> accuracy_output_file, hifi_A, hifi_B, father_A, father_B, mather_A, mather_B
			if hifi_A != hifi_B and father_A != father_B and mather_A != mather_B:
				triple_heterozygous = True
				triple_heterozygous_total_number += 1
				print >> accuracy_output_file, hifi_A, hifi_B, father_A, father_B, mather_A, mather_B
		father_A = genotype_father_dict[position][0]
		father_B = genotype_father_dict[position][1]
		mather_A = genotype_mather_dict[position][0]
		mather_B = genotype_mather_dict[position][1]
		if not triple_heterozygous:
			line_ori = snp_hap_ori_dict[position]
			elements_hifi = line_hifi.strip().split()
			hifi_A = elements_hifi[2].strip()
			hifi_B = elements_hifi[3].strip()
			elements_ori = line_ori.strip().split()
			ori_A = elements_ori[2].strip()
			ori_B = elements_ori[3].strip()
			
			# the hifi seed is from father, A
			if hifi_A == ori_A:	#A is A
				if hifi_B == ori_B:
					same_AB_total_number += 1	# same AB
				else:
					same_A_total_number += 1	# same A
			elif hifi_B == ori_B:
				same_B_total_number += 1	# same B
			elif ori_A == "X" or ori_B == "X" or ori_A == "N" or ori_B == "N":
				same_AB_total_number += 1	# same AB
			else:
				not_same_AB_total_number += 1
				print >> difference_output_file, hifi_A, hifi_B, ori_A, ori_B, father_A, father_B, mather_A, mather_B
			same_position_total_number += 1
			"""
			# for solid data 4 and 6, the hifi seed is from mother (B)
			if hifi_A == ori_B:	#A is B
				if hifi_B == ori_A:
					same_AB_total_number += 1	# same AB
				else:
					same_A_total_number += 1	# same A
			elif hifi_B == ori_A:
				same_B_total_number += 1	# same B
			elif ori_A == "X" or ori_B == "X" or ori_A == "N" or ori_B == "N":
				same_AB_total_number += 1	# same AB
			else:
				not_same_AB_total_number += 1
			same_position_total_number += 1
			"""
	else:
		different_position_total_number += 1
		#difference_output_file.write(line_hifi + "\n")

		
print "same_position_total_number", same_position_total_number
print "different_position_total_number", different_position_total_number

print "same_AB_total_number", same_AB_total_number
print "same_A_total_number", same_A_total_number
print "same_B_total_number", same_B_total_number
print "not_same_AB_total_number", not_same_AB_total_number
print "triple_heterozygous_total_number", triple_heterozygous_total_number

pencentage_in_common = format(float(same_position_total_number)/float(snp_hap_hifi_total_number)*100, "0.2f")
accuracy = format(float(same_A_total_number + same_B_total_number + same_AB_total_number)/float(same_position_total_number)*100, "0.2f")	

print "pencentage in common", float(same_position_total_number)/float(snp_hap_hifi_total_number)
print "accuracy", float(same_AB_total_number)/float(same_position_total_number)


accuracy_output_file.write("snp_hap_ori_total_number: " + str(snp_hap_ori_total_number) + "\n")
accuracy_output_file.write("snp_hap_hifi_total_number: " + str(snp_hap_hifi_total_number) + "\n")
accuracy_output_file.write("same_position_total_number: " + str(same_position_total_number) + "\n")
accuracy_output_file.write("different_position_total_number: " + str(different_position_total_number) + "\n")
accuracy_output_file.write("same_AB_total_number: " + str(same_AB_total_number) + "\n")
accuracy_output_file.write("same_A_total_number: " + str(same_A_total_number) + "\n")
accuracy_output_file.write("same_B_total_number: " + str(same_B_total_number) + "\n")
accuracy_output_file.write("not_same_AB_total_number: " + str(not_same_AB_total_number) + "\n")
accuracy_output_file.write("triple_heterozygous_total_number: " + str(triple_heterozygous_total_number) + "\n")
accuracy_output_file.write("pencentage in common: " + str(pencentage_in_common) + "\n")	
accuracy_output_file.write("accuracy: " + str(accuracy) + "\n")

inputFile_hap_ori.close()
inputFile_hap_hifi.close()
accuracy_output_file.close()

difference_output_file.close()

"""
# record data
data_record_file_name = "solid_process_4.txt"
data_record_file = open(data_record_path + data_record_file_name, "a")
print >> data_record_file, "hifi_data", same_position_total_number, (same_A_total_number+same_B_total_number+not_same_AB_total_number), snp_hap_hifi_total_number, accuracy
data_record_file.close()
cmd = "grep hifi_data hifi_accuracy.txt >> " + data_record_path + data_record_file_name
print cmd
os.system(cmd)
"""
"""
data_record_file = open(data_record_path + data_record_file_name, "a")
print >> data_record_file, ""
data_record_file.close()
"""
