#!/usr/bin/python

# location /home/guoxing/tool/morehouse

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

file_path = "/home/guoxing/morehouse_new_hap/fils/"
currentPath = os.getcwd() + '/'

class snp_hap_ori:
	def __init__(me, rsID, position, A, B):
		me.rsID = rsID
		me.position = position
		me.A = A
		me.B = B

class snp_hap_hifi:
	def __init__(me, rsID, position, A, B):
		me.rsID = rsID
		me.position = position
		me.A = A
		me.B = B

snp_hap_ori_dict = {}
snp_hap_hifi_dict = {}
snp_hap_ori_total_number = 0
snp_hap_ori_total_number_withoutXN = 0
snp_hap_hifi_total_number = 0

hap_ori_file_name = "NA12878_hap_new_refed.txt"
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
		if elements[2].strip() != "N" and elements[2].strip() != "X" and elements[3].strip() != "N" and elements[3].strip() != "X":
		#if True:
			position = elements[1].strip()	
			try:
				position = int(position)
				snp_hap_ori_dict[position] = line.strip()
			except ValueError:
				print file_name, position	
				print line
			snp_hap_ori_total_number += 1
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
"""		
print "n_number", n_number
print "x_number", x_number

print "a_number", a_number
print "t_number", t_number
print "c_number", c_number
print "g_number", g_number
"""			
for line in inputFile_hap_hifi:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		if True:
		#if elements[2].strip() != "N" and elements[2].strip() != "X" and elements[3].strip() != "N" and elements[3].strip() != "X":
			position = elements[1].strip()	
			try:
				position = int(position)
				snp_hap_hifi_dict[position] = line.strip()
			except ValueError:
				print file_name, position	
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

difference_output_file_name = "hifi_difference.txt"
difference_output_file = open(currentPath + difference_output_file_name, "w")

for position, line_hifi in snp_hap_hifi_dict.iteritems():
	if position in snp_hap_ori_dict:
		line_ori = snp_hap_ori_dict[position]
		elements_hifi = line_hifi.strip().split()
		elements_ori = line_ori.strip().split()
		# same AB
		if elements_hifi[2].strip() == elements_ori[2].strip():
			if elements_hifi[3].strip() == elements_ori[3].strip():
				same_AB_total_number += 1	# same AB
			else:
				same_A_total_number += 1	# same A
		elif elements_hifi[3].strip() == elements_ori[3].strip():
			same_B_total_number += 1	# same B
		else:
			not_same_AB_total_number += 1
		same_position_total_number += 1
		#difference_output_file.write(line_hifi + "\t" + line_ori + "\n")
	else:
		different_position_total_number += 1
		difference_output_file.write(line_hifi + "\n")

		
print "same_position_total_number", same_position_total_number
print "different_position_total_number", different_position_total_number

print "same_AB_total_number", same_AB_total_number
print "same_A_total_number", same_A_total_number
print "same_B_total_number", same_B_total_number
print "not_same_AB_total_number", not_same_AB_total_number

print "pencentage in common", float(same_position_total_number)/float(snp_hap_hifi_total_number)
print "accuracy", float(same_AB_total_number)/float(same_position_total_number)


accuracy_output_file.write("snp_hap_ori_total_number: " + str(snp_hap_ori_total_number) + "\n")
accuracy_output_file.write("snp_hap_hifi_total_number: " + str(snp_hap_hifi_total_number) + "\n")
accuracy_output_file.write("same_AB_total_number: " + str(same_AB_total_number) + "\n")
accuracy_output_file.write("same_A_total_number: " + str(same_A_total_number) + "\n")
accuracy_output_file.write("same_B_total_number: " + str(same_B_total_number) + "\n")
accuracy_output_file.write("not_same_AB_total_number: " + str(not_same_AB_total_number) + "\n")
accuracy_output_file.write("pencentage in common: " + str(float(same_position_total_number)/float(snp_hap_hifi_total_number)) + "\n")			
accuracy_output_file.write("accuracy: " + str(float(same_A_total_number + same_B_total_number + same_AB_total_number)/float(same_position_total_number)) + "\n")

inputFile_hap_ori.close()
inputFile_hap_hifi.close()
accuracy_output_file.close()

difference_output_file.close()
