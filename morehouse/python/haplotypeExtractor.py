#!/usr/bin/python
import os, glob, subprocess, random, operator
from optparse import OptionParser

# location /home/guoxing/tool/morehouse

# replace X with N. Keep N

projectPath = '/home/guoxing/morehouse/'
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--file", type="string", dest="fileName",help = "Input File Name", default="null")
(options, args) = parser.parse_args()

fileName = options.fileName

NA_number = "NA12878"
NA_A = NA_number+"_A"
NA_B = NA_number+"_B"
chromsome = "6"

#file specific
chromsome_position = 1
rs_number_position = 0
pos_position = 2
A_position = 3
B_position = 4

inputFile = open(currentPath+fileName, "r")
outputFile = open(currentPath+NA_number+"_hg18ch"+chromsome+"_hap.txt", "w")
XFile = open(currentPath+NA_number+"_hg18ch"+chromsome+"_X.txt", "w")
recordFile = open(currentPath+NA_number+"_hg18ch"+chromsome+"_record.txt", "w")
outputFile.write("rsID \t"+ "phys_position \t" + NA_A + "\t" + NA_B + "\t" + "\n");	

# First line
line = inputFile.readline()
#elements = line.split()
#print elements[rs_number_position].strip()
#print elements[pos_position].strip()
#print elements[A_position].strip()
total = 0
no_rs = 0
x_snp = 0

for line in inputFile:
	#if line.startswith(chromsome):
	elements = line.split()
	if 	elements[chromsome_position].strip() == chromsome:
		rs_number = elements[rs_number_position].strip()
		phys_position = elements[pos_position].strip()
		A_char = elements[A_position].strip()
		B_char = elements[B_position].strip()
		if (A_char == "X" or B_char == "X") and elements[rs_number_position].strip() != ".":
			XFile.write(rs_number+"\t"+phys_position+"\t" +A_char +"\t"+ B_char+"\t"+"\n")
			x_snp += 1
		if A_char == "X":
			A_char = "N"
		if B_char == "X":
			B_char = "N"
		if elements[rs_number_position].strip() != ".":  # remove snp without a rs number
			outputFile.write(rs_number+"\t"+phys_position+"\t" +A_char +"\t"+ B_char+"\t"+"\n")	
		total += 1
		if elements[rs_number_position].strip() == ".":
			no_rs += 1

inputFile.close()
outputFile.close()
XFile.close()

no_rs_percentage = no_rs*100/total
x_snp_percentage = x_snp*100/total

recordFile.write("Total snp number:"+"\t"+str(total)+"\n")
recordFile.write("snp without rs number:"+"\t"+str(no_rs)+"\n")
recordFile.write("snp without rs percentage:"+"\t"+str(no_rs_percentage)+"% \n")
recordFile.write("snp with X number:"+"\t"+str(x_snp)+"\n")
recordFile.write("snp with X percentage:"+"\t"+str(x_snp_percentage)+"% \n")


print "Total snp number: ", total
print "snp without rs number: ", no_rs
print "snp without rs percentage: ", no_rs_percentage	
print "snp with X number: ", no_rs
print "snp with X percentage: ", x_snp_percentage	



