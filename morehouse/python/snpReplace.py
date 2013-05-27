#!/usr/bin/python
import os, glob, subprocess, random, operator
from optparse import OptionParser

currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--chrfile", type="string", dest="ChrFileName",help = "Input Chromosome File", default="null")
parser.add_option("-p", "--hapfile", type="string", dest="haplotypeFileName",help = "Input Haplotype File", default="null")
(options, args) = parser.parse_args()

ChrFileName = options.ChrFileName
haplotypeFileName = options.haplotypeFileName

#A from Father, B from Mother
class snp:
	def __init__(me, position, A, B):
		me.position = position
		me.A = A
		me.B = B

ChrFileName = "hg18chr6_mar12.fa"
haplotypeFileName = "NA12878_hap_new_refed.txt"

#dataSetName = haplotypeFileName[:haplotypeFileName.find('.')-4]

dataSetName = "NA12878_hg18ch6"

output_A = dataSetName+"_A.fa"
output_B = dataSetName+"_B.fa"

inputFile_chr = open(currentPath+ChrFileName, "r")
inputFile_hap = open(currentPath+haplotypeFileName, "r")

outputFile_A = open(currentPath+output_A, "w")
outputFile_A.write(">"+dataSetName+"_A \n")
outputFile_B = open(currentPath+output_B, "w")
outputFile_B.write(">"+dataSetName+"_B \n")

outputFile_record = open(currentPath+dataSetName+"_record.txt", "w")
outputFile_record.write("chr_position \t"+"hap_position \t"+"Original \t"+"NA12878_A \t"+"NA12878_B"+"\n")

# get haplotype info
snpList=[]

for line in inputFile_hap:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		position = int(elements[1].strip())
		A = elements[2].strip()
		B = elements[3].strip()
		snpList.append(snp(position, A, B))
		
# this file is already sorted	
snpList.sort(key=operator.attrgetter('position'))	
	
current_position = 0
snp_index = 0
snp_position = int(snpList[snp_index].position)
		
for line in inputFile_chr:
	if not line.startswith(">"):
		line = line.strip()	
		for base in line:
			current_position  += 1
			if current_position != snp_position:
				outputFile_A.write(base)
				outputFile_B.write(base)
			if current_position == snp_position:
				outputFile_record.write(str(current_position)+"\t"+str(snp_position)+"\t"+base+"\t"+snpList[snp_index].A+"\t"+snpList[snp_index].B+"\n")
				outputFile_A.write(snpList[snp_index].A)
				outputFile_B.write(snpList[snp_index].B)
				if snp_index < (len(snpList)-1):
					snp_index += 1
					snp_position = int(snpList[snp_index].position)
		outputFile_A.write("\n")
		outputFile_B.write("\n")

inputFile_chr.close()
inputFile_hap.close()
outputFile_A.close()
outputFile_B.close()

print "size of snpList: ",len(snpList)


