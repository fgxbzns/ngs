#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# this is for real NGS data

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

program_path = "/home/guoxing/tool/morehouse"
file_path = "/home/guoxing/morehouse_new_hap/fils/"
currentPath = os.getcwd() + '/'

# to import other python files
#from home.guoxing.tool.morehouse.qualityScoreDict import quality_score_dict
"""
import sys
sys.path.append(program_path)
import quality_score_dict
"""


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

(options, args) = parser.parse_args()

haplotype_file = options.haplotypeFile
sam_file = options.samFile

# define variables

first_read_length = 35
second_read_length = 37
insert_size_lower_bond = 100
insert_size_upper_bond = 1000
depth_cutoff = 1



#reads_length = 100

#if haplotype_file == "null":
	#print "Please input the haplotype file name"

haplotype_file = "NA12878_hap_new_refed.txt"
#haplotype_file = "nbt_1739_S4_chr6.txt"
#haplotype_file = "NA12878_hg18ch6_hap.txt"


if sam_file == "null":
	print "Please input the sam file name"
	sam_file = "SRR077303_18_10000.txt"
	sam_file = "SRR077303_18.sam"
	

sam_file_name = sam_file[:(len(sam_file)-4)]

print "haplotype file: ", haplotype_file
print "sam file: ", sam_file_name


sam_path = "/home/guoxing/disk2/depth_real_data/"

inputFile_hap = open(file_path + haplotype_file, "r")
#inputFile_sam = open(currentPath + sam_file, "r")
inputFile_sam = open(sam_path + sam_file, "r")
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
snp_dict={}

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
		snp_dict[position] = snp(rsID, position, A, B, covered_reads_list, allele_dict)

reads_list=[]
insert_size = 0

hap_homo_file = open(currentPath + "hap_homo.txt", "w")
hap_hete_file = open(currentPath + "hap_hete.txt", "w")

for snp in snp_list:
	if snp.A == snp.B:
		hap_homo_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\n")
	else:
		hap_hete_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\n")
	
hap_homo_file.close()
hap_hete_file.close()

sam_line_first = inputFile_sam.readline() # the first read line in a pair


total_reads_num = 0
XA_total = 0
reads_within_size_limit = 0
covered_snp_total_number = 0

while sam_line_first!='':
	if not sam_line_first.startswith("@"):
		total_reads_num += 1	
		elements_first = sam_line_first.strip().split()
		try:
			read_ID_first = elements_first[0].strip()
			rName_first = elements_first[2].strip()
			insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for second read is negative
		except:
			print "error in first read:", sam_line_first
		if (insert_size_first >= insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):
			reads_within_size_limit += 1
			multiple_maping_first = False
			# check multiple maping for first read
			try:
				XA = elements_first[21].strip()
				#print XA
				XA_total += 1
				multiple_maping_first = True
			except:
				multiple_maping_first = False

			# if the first read is within insert size limit, check the second read
			sam_line_second = inputFile_sam.readline()
			total_reads_num += 1
			elements_second = sam_line_second.strip().split()
			multiple_maping_second = False
			try:
				read_ID_second = elements_second[0].strip()
				rName_second = elements_second[2].strip()
				insert_size_second = abs(int(elements_second[8].strip()))			#  insert_size for second read is negative
			except:
				print "error in second read:", sam_line_second
			if read_ID_first == read_ID_second:		# check if the two reads belong to the same pair
				try:
					XA = elements_second[21].strip()
					#print XA
					XA_total += 1
					multiple_maping_second = True
				except:
					multiple_maping_second = False
				if (not multiple_maping_first)	or (not multiple_maping_second): # keep the pair as long as one read is not multiple mapping
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
							
							#print "start_position", start_position
							#print "snp_position", snp_dict[start_position+i].position
							covered_snp = read_sequence_first[i]			# ith position is the covered snp
							quality_score_symbol = quality_score_sequence_first[i]
							if (rName_first == 'chr6')and (not covered_snp == 'N') and (not quality_score_symbol in quality_score_dict):	# check quality_score
								covered_snp_total_number += 1
								snp_dict[start_position_first+i].covered_reads_list.append(read(qName_first, flag_first, rName_first, start_position_first, read_sequence_first, quality_score_sequence_first, read_length_first, covered_snp))					
								if covered_snp == 'A':
									snp_dict[start_position_first+i].allele_dict['A'] += 1
								elif covered_snp == 'T':
									snp_dict[start_position_first+i].allele_dict['T'] += 1
								elif covered_snp == 'C':
									snp_dict[start_position_first+i].allele_dict['C'] += 1
								elif covered_snp == 'G':
									snp_dict[start_position_first+i].allele_dict['G'] += 1			
						i += 1
					
					# second read
					qName_second = elements_second[0].strip()
					flag_second = elements_second[1].strip()
					start_position_second = int(elements_second[3].strip())
					read_sequence_second = elements_second[9].strip()
					read_length_second = len(read_sequence_second)
					quality_score_sequence_second = elements_second[10].strip()
					i = 0
					while i < read_length_second:
						if (start_position_second+i) in snp_dict:
							
							#print "start_position", start_position
							#print "snp_position", snp_dict[start_position+i].position
							covered_snp = read_sequence_second[i]			# ith position is the covered snp
							quality_score_symbol = quality_score_sequence_second[i]
							if (rName_second == 'chr6') and (not covered_snp == 'N') and (not quality_score_symbol in quality_score_dict):
								covered_snp_total_number += 1
								snp_dict[start_position_second+i].covered_reads_list.append(read(qName_second, flag_second, rName_second, start_position_second, read_sequence_second, quality_score_sequence_second, read_length_second, covered_snp))					
								if covered_snp == 'A':
									snp_dict[start_position_second+i].allele_dict['A'] += 1
								elif covered_snp == 'T':
									snp_dict[start_position_second+i].allele_dict['T'] += 1
								elif covered_snp == 'C':
									snp_dict[start_position_second+i].allele_dict['C'] += 1
								elif covered_snp == 'G':
									snp_dict[start_position_second+i].allele_dict['G'] += 1
							#else:
							#	print quality_score_dict[quality_score_symbol]	
						i += 1
			
			else:
				print "first and second read ID do not match", read_ID_first, read_ID_second			
			
			
	sam_line_first = inputFile_sam.readline()


print "haplotye list size is: ", len(snp_dict)	
print "total_reads_num", total_reads_num
print "reads_within_size_limit", reads_within_size_limit
print "XA_total", XA_total
print "covered_snp_total_number", covered_snp_total_number

"""

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
			covered_snp = ""
			reads_list.append(read(qName, flag, rName, start_position, read_sequence, covered_snp))

reads_list.sort(key=operator.attrgetter('start_position'))	

# keep a record of sorted sam reads
#for reads in reads_list:
#	outputFile_sorted_sam.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.read_sequence + "\n")
"""


snp_sorted_list = [x for x in snp_dict.iteritems()] 
snp_sorted_list.sort(key=lambda x: x[0]) # sort by key

hap_homo_file = open(currentPath + "hap_homo.txt", "w")
hap_hete_file = open(currentPath + "hap_hete.txt", "w")

rs_different_number = 0
homo_number = 0

for snp_data in snp_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) >= depth_cutoff:
		outputFile_allele.write("chr6???"+"\t"+str(snp.position)+"\t"+str(len(snp.covered_reads_list))+"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+"\n")
		outputFile_reads.write("@_" + str(snp.position) + "\t" + snp.A  + "\t" + snp.B  + "\t" + str(len(snp.covered_reads_list)) + "\n")
		for reads in snp.covered_reads_list:
			outputFile_reads.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
		if snp.A == snp.B:
			homo_number += 1
			#hap_homo_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\n")
			rs_different = False
			for reads in snp.covered_reads_list:
				if snp.A != reads.covered_snp:
					rs_different = True
					#hap_homo_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
			if rs_different:
				rs_different_number += 1
				hap_homo_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\n")
				for reads in snp.covered_reads_list:
					if snp.A != reads.covered_snp:
						hap_homo_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")
		else:
			hap_hete_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + snp.B + "\n")
			for reads in snp.covered_reads_list:
				hap_hete_file.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\t" + reads.quality_score_sequence + "\n")


print "rs_different_number", rs_different_number
print "homo_number", homo_number

hap_homo_file.close()
hap_hete_file.close()
			


"""
for position, snp in snp_dict.iteritems():
	snp = snp_dict[position]
	if len(snp.covered_reads_list) > 0:
		outputFile_allele.write("chr6???"+"\t"+str(snp.position)+"\t"+str(len(snp.covered_reads_list))+"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+"\n")
		outputFile_reads.write("@_" + str(snp.position) + "\t" + snp.A  + "\t" + snp.B  + "\t" + str(len(snp.covered_reads_list)) + "\n")
		for reads in snp.covered_reads_list:
			outputFile_reads.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\n")
"""



"""						
for snp in snp_list:
	if len(snp.covered_reads_list) > 0:
		outputFile_allele.write("chr6"+"\t"+str(snp.position)+"\t"+str(len(snp.covered_reads_list))+"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+"\n")
		outputFile_reads.write("@_" + str(snp.position) + "\t" + snp.A  + "\t" + snp.B  + "\t" + str(len(snp.covered_reads_list)) + "\n")
		for reads in snp.covered_reads_list:
			outputFile_reads.write(reads.qName + "\t" + reads.flag + "\t" + reads.rName + "\t" + str(reads.start_position) + "\t" + reads.covered_snp + "\t" + reads.read_sequence + "\n")


"""
# compare with ori hap file

unchanged_snp_file_name = sam_file_name + "_unchanged_snp.txt"
unchanged_snp_file = open(currentPath + unchanged_snp_file_name, "w")		
unchanged_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

changed_snp_file_name = sam_file_name + "_changed_snp.txt"
changed_snp_file = open(currentPath + changed_snp_file_name, "w")		
changed_snp_file.write("SNP position \t SNP_A \t Depth \t A_depth \t T_depth \t C_depth \t G_depth \n")

# need to update for A or B 	
for snp_data in snp_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) >= depth_cutoff:	
		consisitent = True
		for reads in snp.covered_reads_list:
			#if snp.A != reads.covered_snp:
			if snp.B != reads.covered_snp:
				consisitent = False
		if consisitent:
			#unchanged_snp_file.write(str(snp.position) + "\t" + snp.A + "\t" + str(len(snp.covered_reads_list)) + "\t" +"\n")
			unchanged_snp_file.write(str(snp.position) + "\t" + snp.B + "\t" + str(len(snp.covered_reads_list)) + "\t" +"\n")
		if not consisitent:
			#changed_snp_file.write(str(snp.position) + "\t" + snp.A + "\t" + str(len(snp.covered_reads_list)) +"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+ "\n")
			changed_snp_file.write(str(snp.position) + "\t" + snp.B + "\t" + str(len(snp.covered_reads_list)) +"\t"+str(snp.allele_dict['A'])+"\t"+str(snp.allele_dict['T'])+"\t"+str(snp.allele_dict['C'])+"\t"+str(snp.allele_dict['G'])+ "\n")

unchanged_snp_file.close()
changed_snp_file.close()



# prepare the haplotype file for hifi

hifi_pure_file_name = sam_file_name + "_hifi_pure.txt"
hifi_pure_file = open(currentPath + hifi_pure_file_name, "w")		
hifi_pure_file.write("rsID \t phys_position \t NA12878_A \t	NA12878_B \n")




hifi_max_file_name = sam_file_name + "_hifi_max.txt"
hifi_max_file = open(currentPath + hifi_max_file_name, "w")		
hifi_max_file.write("rsID \t phys_position \t NA12878_A \t	NA12878_B \n")


base_list = ["A", "T", "C", "G"]


hete_A_max_file = open(currentPath + sam_file_name + "_hete_max_A.txt", "w")
hete_A_max_file.write("rsID \t phys_position \t NA12878_A \t	selected_base \n")

hete_B_max_file = open(currentPath + sam_file_name + "_hete_max_B.txt", "w")
hete_B_max_file.write("rsID \t phys_position \t NA12878_B \t	selected_base \n")

hete_A_pure_file = open(currentPath + sam_file_name + "_hete_pure_A.txt", "w")
hete_A_pure_file.write("rsID \t phys_position \t NA12878_A \t	selected_base \n")

hete_B_pure_file = open(currentPath + sam_file_name + "_hete_pure_B.txt", "w")
hete_B_pure_file.write("rsID \t phys_position \t NA12878_B \t	selected_base \n")

max_hete = 0
pure_hete = 0


com_total = 0
com_A = 0
com_B = 0
com_AB = 0

for snp_data in snp_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) >= depth_cutoff:	
		max_base = keywithmaxval(snp.allele_dict)
		max_value = snp.allele_dict[max_base]
		hifi_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\t" + str(max_value) + "\n")
		#hifi_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")
		
		if snp.A != snp.B:		#hete
			max_hete += 1
			if max_base == snp.A:
				hete_A_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + max_base + "\t" + str(max_value) + "\n")
			if max_base == snp.B:
				hete_B_max_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.B + "\t" + max_base + "\t" + str(max_value) + "\n")
		
		pure = True
		for base in base_list:
			if base != max_base:
				if snp.allele_dict[base] != 0:
					pure = False
		if pure:
			#hifi_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\t" + str(max_value) + "\n")
			if max_base == snp.A or max_base == snp.B:
				com_total += 1
				if max_base == snp.A and max_base != snp.B:
					com_A += 1
				if max_base == snp.B and max_base != snp.A:
					com_B += 1
				if max_base == snp.B and max_base == snp.A:
					com_AB += 1
				hifi_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + max_base + "\n")
			if snp.A != snp.B:		#hete
				pure_hete += 1
				if max_base == snp.A:
					hete_A_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.A + "\t" + max_base + "\t" + str(max_value) + "\n")
				if max_base == snp.B:
					hete_B_pure_file.write(snp.rsID + "\t" + str(snp.position) + "\t" + snp.B + "\t" + max_base + "\t" + str(max_value) + "\n")


print "com_total", com_total
print "com_A", com_A
print "com_B", com_B
print "com_AB", com_AB

hifi_pure_file.close()
hifi_max_file.close()
	


print "max_hete", max_hete
print "pure_hete", pure_hete

hete_A_max_file.close()
hete_B_max_file.close()
hete_A_pure_file.close()
hete_B_pure_file.close()


same_A_file = open(currentPath + sam_file_name + "_same_A.txt", "w")
same_A_file.write("rsID \t phys_position \t NA12878_A \t selected_base \t reads	\n")

pair_A_file = open(currentPath + sam_file_name + "_pair_A.txt", "w")
pair_A_file.write("rsID \t phys_position \t NA12878_A \t selected_base \t reads	\n")

different_A_file = open(currentPath + sam_file_name + "_different_A.txt", "w")
different_A_file.write("rsID \t phys_position \t NA12878_A \t selected_base \t reads \n")

same_A = 0
pair_A = 0
different_A = 0
pure_total = 0

for snp_data in snp_sorted_list:
	snp = snp_data[1]
	if len(snp.covered_reads_list) >= depth_cutoff:	
		max_base = keywithmaxval(snp.allele_dict)
		max_value = snp.allele_dict[max_base]

		pure = True
		for base in base_list:
			if base != max_base:
				if snp.allele_dict[base] != 0:
					pure = False
		if pure:
			pure_total += 1
			hap_base = snp.B
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
			

print "pure_total", pure_total
print "same_A", same_A
print "pair_A", pair_A
print "different_A", different_A

data_record_file.write("pure_total is: " + str(pure_total) + "\n")
data_record_file.write("same_A is: " + str(same_A) + "\n")
data_record_file.write("pair_A is: " + str(pair_A) + "\n")
data_record_file.write("different_A is: " + str(different_A) + "\n")




end = time.time()
run_time = str(end - start)
run_time = run_time[:(run_time.find('.')+3)]
print "run time is: " + run_time + "s"

data_record_file.write("sam file is: " + sam_file + "\n")
data_record_file.write("haplotype file is: " + haplotype_file + "\n")
data_record_file.write("haplotye list size is: " + str(len(snp_dict)) + "\n")
data_record_file.write("reads list size is: " + str(total_reads_num) + "\n")
data_record_file.write("run time is: " + run_time + "s \n")
data_record_file.write("\n")

data_record_file.write("rs_different_number is: " + str(rs_different_number) + "s \n")
data_record_file.write("homo_number is: " + str(homo_number) + "s \n")



data_record_file.close()

inputFile_hap.close()
inputFile_sam.close()

outputFile_reads.close()

outputFile_allele.close()
#outputFile_sorted_sam.close()

