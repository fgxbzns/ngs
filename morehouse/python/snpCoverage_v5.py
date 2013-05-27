#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for simulation NGS data

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

file_path = "/home/guoxing/morehouse_new_hap/fils/"
currentPath = os.getcwd() + '/'

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
	def __init__(me, qName, flag, rName, start_position, read_sequence, base):
		me.qName = qName
		me.flag = flag
		me.rName = rName
		me.start_position = start_position
		me.read_sequence = read_sequence
		me.base = base
	
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

(options, args) = parser.parse_args()

haplotype_file = options.haplotypeFile
sam_file = options.samFile

reads_length = 100

#if haplotype_file == "null":
	#print "Please input the haplotype file name"

haplotype_file = "NA12878_hap_new_refed.txt"

if sam_file == "null":
	print "Please input the sam file name"
	#sam_file = ".sam"

sam_file_name = sam_file[:(len(sam_file)-4)]

print "haplotype file: ", haplotype_file
print "sam file: ", sam_file_name


inputFile_hap = open(file_path + haplotype_file, "r")
inputFile_sam = open(currentPath + sam_file, "r")
"""
outputFile_sorted_sam = open(currentPath + sam_file_name + "_sorted_sam.txt", "w")
outputFile_sorted_sam.write("qname \t flag \t rname \t position \t sequence Depth \n")
"""
outputFile_reads = open(currentPath + sam_file_name + "_reads.txt", "w")
outputFile_reads.write("SNP position \t Depth \n")

outputFile_allele = open(currentPath + sam_file_name+"_allele.txt", "w")
outputFile_allele.write("Chromosome \t position \t Total Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

data_record_file_name = sam_file_name + "_data_record.txt"
data_record_file = open(currentPath + data_record_file_name, "w")

# get haplotype info
snp_list=[]

for line in inputFile_hap:
	if not line.startswith("rsID"):
		elements = line.strip().split()
		rsID = elements[0].strip()
		position = int(elements[1].strip())
		A = elements[2].strip()
		B = elements[3].strip()
		covered_reads_list = []
		allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
		snp_list.append(snp(rsID, position, A, B, covered_reads_list, allele_dict))
		
# this file is already sorted	
snp_list.sort(key=operator.attrgetter('position'))	

reads_list=[]

for line in inputFile_sam:
	if not line.startswith("@"):
		elements = line.strip().split()
		try:
			rName = elements[2].strip()
			gap_info = elements[8].strip()
		except:
			print line
		if (rName == "chr6") and (not gap_info == "0"):  # remove reads without chr name
			qName = elements[0].strip()
			flag = elements[1].strip()
			start_position = int(elements[3].strip())
			read_sequence = elements[9].strip()
			base = ""
			reads_list.append(read(qName, flag, rName, start_position, read_sequence, base))

reads_list.sort(key=operator.attrgetter('start_position'))	

# keep a record of sorted sam reads
#for reads in reads_list:
#	outputFile_sorted_sam.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.read_sequence + "\n")

print "haplotye list size is: ", len(snp_list)
print "reads list size is: ", len(reads_list)

current_read = 0
reads_index = 0

for snp in snp_list:
	reads_index = current_read
	if reads_list[reads_index].start_position <= snp.position:
		if reads_list[reads_index].start_position <= (snp.position - reads_length):
			while reads_list[reads_index].start_position <= (snp.position - reads_length):
				reads_index += 1
				current_read += 1
		while reads_list[reads_index].start_position <= snp.position and reads_list[reads_index].start_position > (snp.position - reads_length):
			current_reads = reads_list[reads_index]
			base = current_reads.read_sequence[int(snp.position)-int(current_reads.start_position)]
			if not base == 'N':
				snp.covered_reads_list.append(read(current_reads.qName, current_reads.flag, current_reads.rName, current_reads.start_position, current_reads.read_sequence, base))					
				if base == 'A':
					snp.allele_dict['A'] += 1
				elif base == 'T':
					snp.allele_dict['T'] += 1
				elif base == 'C':
					snp.allele_dict['C'] += 1
				elif base == 'G':
					snp.allele_dict['G'] += 1									
			reads_index += 1
						
for snp in snp_list:
	if len(snp.covered_reads_list) > 0:
		outputFile_allele.write("chr6"+"\t"+str(snp.position)+"\t"+str(len(snp.covered_reads_list))+"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+"\n")
		outputFile_reads.write("@_" + str(snp.position) + "\t" + snp.A  + "\t" + snp.B  + "\t" + str(len(snp.covered_reads_list)) + "\n")
		for reads in snp.covered_reads_list:
			outputFile_reads.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.base + "\t" + reads.read_sequence + "\n")

"""
# compare with ori hap file

unchanged_snp_file_name = sam_file_name + "_unchanged_snp.txt"
unchanged_snp_file = open(currentPath + unchanged_snp_file_name, "w")		
unchanged_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

changed_snp_file_name = sam_file_name + "_changed_snp.txt"
changed_snp_file = open(currentPath + changed_snp_file_name, "w")		
changed_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

# need to update for A or B 	
for snp in snp_list:
	if len(snp.covered_reads_list) > 0:	
		consisitent = True
		for reads in snp.covered_reads_list:
			if snp.A != reads.base:
			#if snp.B != reads.base:
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

hifi_pure_file_name = sam_file_name + "_hifi_pure.txt"
hifi_pure_file = open(currentPath + hifi_pure_file_name, "w")		
hifi_pure_file.write("rsID \t phys_position \t NA12878_A \t	NA12878_B \n")
"""
hifi_max_file_name = sam_file_name + "_hifi_max.txt"
hifi_max_file = open(currentPath + hifi_max_file_name, "w")		
hifi_max_file.write("rsID \t phys_position \t NA12878_A \t	NA12878_B \n")
"""
base_list = ["A", "T", "C", "G"]

for snp in snp_list:
	if len(snp.covered_reads_list) > 0:	
		max_base = keywithmaxval(snp.allele_dict)
		#max_value = snp.allele_dict[max_base]
		#hifi_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\t" + str(max_value) + "\n")
		#hifi_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")
		pure = True
		for base in base_list:
			if base != max_base:
				if snp.allele_dict[base] != 0:
					pure = False
		if pure:
			#hifi_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\t" + str(max_value) + "\n")
			hifi_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")


hifi_pure_file.close()
#hifi_max_file.close()
	
end = time.time()
run_time = str(end - start)
run_time = run_time[:(run_time.find('.')+3)]
print "run time is: " + run_time + "s"

data_record_file.write("sam file is: " + sam_file + "\n")
data_record_file.write("haplotype file is: " + haplotype_file + "\n")
data_record_file.write("haplotye list size is: " + str(len(snp_list)) + "\n")
data_record_file.write("reads list size is: " + str(len(reads_list)) + "\n")
data_record_file.write("run time is: " + run_time + "s \n")
data_record_file.write("\n")

data_record_file.close()

inputFile_hap.close()
inputFile_sam.close()
outputFile_reads.close()
outputFile_allele.close()
#outputFile_sorted_sam.close()
