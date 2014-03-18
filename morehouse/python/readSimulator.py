#!/usr/bin/python

# NGS reads simulator

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

projectPath = '/home/guoxing/project/'
currentPath = os.getcwd() + '/'

def getTotalBaseNum(fileName):
	totalBase = 0
	f = open(currentPath+fileName, "r")
	for line in f:
		if not line.startswith('>'):
			totalBase += len(line.strip())
	return totalBase
	f.close()

def convertBasePair(base):
	if base == "A":
		base = "T"
	elif base == "T":
		base = "A"
	elif base == "C":
		base = "G"
	elif base == "G":
		base = "C"
	return base
		
# 3'-5' second read
def convertSecondRead_1(second_read_coverage_sequence):
	temp_sequence = ""
	for base in second_read_coverage_sequence:
		temp_sequence += convertBasePair(base)
	return temp_sequence

# 5'-3' second read	
def convertSecondRead_2(second_read_coverage_sequence):
	temp_sequence = ""
	for base in second_read_coverage_sequence:
		temp_sequence += convertBasePair(base)
	return temp_sequence[::-1]

# introduce error into the reads
def replaceAllele(reads):
	allele_selector = random.randrange(0, len(reads))
	base_list = ["A", "T", "C", "G"]
	base_selector = random.randrange(0, len(base_list))
	while reads[allele_selector] == base_list[base_selector]:
		base_selector = random.randrange(0, len(base_list))
	reads = reads[:allele_selector] + base_list[base_selector] + reads[(allele_selector+1):]
	#print "allele_selector", allele_selector
	#print "base_selector", base_selector
	#print reads
	return reads

# this could be used for single end simulation
def introduceError(error_number, paired_end_first_read_list, paired_end_second_read_list):
	for error in range (0, error_number):
		#print "error", error
		list_selector = random.randrange(0, 2)
		#print "list_selector", list_selector
		if list_selector == 0:  # choose first read list
			#print "list_selector first", list_selector
			reads_selector = random.randrange(0, len(paired_end_first_read_list))
			print "reads_selector first", reads_selector
			paired_end_first_read_list[reads_selector] = replaceAllele(paired_end_first_read_list[reads_selector])
		elif list_selector == 1: # choose second read list
			#print "list_selector sec", list_selector
			reads_selector = random.randrange(0, len(paired_end_second_read_list))
			paired_end_second_read_list[reads_selector] = replaceAllele(paired_end_second_read_list[reads_selector]) 	


# start time
start = time.time()		
		
# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--file", type="string", dest="chromosomeFile",help = "Input File Name", default="null")
parser.add_option("-p", "--percentage", type="float", dest="percentage",help = "Chromosome coverage percentage", default="0.3")
parser.add_option("-l", "--readLength", type="int", dest="readLength",help = "read Length", default="100")
parser.add_option("-g", "--gap", type="int", dest="gap", help = "gap size", default="200")
parser.add_option("-d", "--depth", type="float", dest="depth",help = "depth", default="0.1")
parser.add_option("-e", "--errorRate", type="float", dest="errorRate",help = "error rate", default="0.005")
(options, args) = parser.parse_args()

chromosome_file = options.chromosomeFile
coveragePercentage = options.percentage
read_length = options.readLength
gap = options.gap
depth = options.depth
error_rate = options.errorRate

#input file
if chromosome_file == "null":
	chromosome_file = "NA12878_hg18ch6_A.fa"
chromosome_file_name = chromosome_file[:chromosome_file.find('.')] + "_" + str(depth) + "x" + "_" + str(error_rate) + "er" 
chromosome_file_input = open(currentPath + chromosome_file, "r")

#output files
fragment_file_name = chromosome_file_name + "_fragment_" + str(coveragePercentage)+".fa"
fragment_file_output = open(currentPath + fragment_file_name, "w")

single_end_file_name = chromosome_file_name + "_single_end.fa"

paired_end_file_name_first_read = chromosome_file_name + "_paired_first.fa"
paired_end_file_name_second_read = chromosome_file_name + "_paired_second.fa"

paired_end_file_name_first_read_no_error = chromosome_file_name + "_paired_end_first_no_error.fa"
paired_end_file_name_second_read_no_error = chromosome_file_name + "_paired_end_second_no_error.fa"

errors_file_name = chromosome_file_name + "_errors.fa"

data_record_file_name = chromosome_file_name + "_data_record.txt"

lower_bound = 500
upper_bound = 800
sequence = ""
total_base_number = getTotalBaseNum(chromosome_file)
covered_base_number = int(total_base_number * coveragePercentage)
window_number = covered_base_number / ((lower_bound + upper_bound) / 2)
window_size = total_base_number / window_number

current_position = 0

window_ini = current_position
window_end = current_position + window_size
window_sequenc = ""

fragment_size = random.randrange(lower_bound,upper_bound)
fragment_ini = random.randrange(window_ini,(window_end-fragment_size))
fragment_end = fragment_ini + fragment_size
fragment_sequenc = ""


fragment_ini = random.randrange(0,(window_size-fragment_size))
fragment_end = fragment_ini + fragment_size
fragment_total_number = 0
fragment_list = []

fragment_leftover_sequence = ""

lineNumber = 0

# skip title line
line = chromosome_file_input.readline()

while line != "":
	line = chromosome_file_input.readline().strip()
	#print line
	current_position += len(line)
	if current_position <= window_end:
		window_sequenc += line
		#line = inputFile.readline().strip()
		#current_position += len(line)
	
	if current_position > window_end:
		fragment_sequence = window_sequenc[fragment_ini:fragment_end].upper()   # convert all letters to uppercase
		fragment_total_number += 1
		fragment_list.append(fragment_sequence)
		fragment_file_output.write("fragment ID is: "+str(fragment_total_number)+"  size is : \t"+str(len(fragment_sequence))+"\n")
		fragment_file_output.write(fragment_sequence+"\n")		
		fragment_leftover_sequence = window_sequenc[fragment_end:] 
		fragment_sequenc = ""
		window_sequenc = fragment_leftover_sequence
		fragment_leftover_sequence_size = len(fragment_leftover_sequence)
		fragment_leftover_sequence = ""
		window_end += window_size
		fragment_ini = random.randrange(0,(fragment_leftover_sequence_size+window_size-fragment_size))
		fragment_size = random.randrange(lower_bound,upper_bound)
		fragment_end = fragment_ini + fragment_size	

chromosome_file_input.close()
fragment_file_output.close()

total_reads_number = int(depth * total_base_number * coveragePercentage / (2 * read_length))
single_end_reads_number = total_reads_number
paired_end_reads_number_each_read = total_reads_number / 2

total_error_number = int(total_reads_number * read_length * error_rate)

print "total_reads_number", total_reads_number
print "paired_end_reads_number_each_read", paired_end_reads_number_each_read

paired_end_first_read_list = []
paired_end_second_read_list = []

# files containg reads without errors
#paired_end_first_read_no_rerror_output = open(currentPath + paired_end_file_name_first_read_no_error, "w")
#paired_end_second_read_no_rerror_output = open(currentPath + paired_end_file_name_second_read_no_error, "w")

quality_score_sequence = ""
for quality_score in range(0, read_length):
	 quality_score_sequence += "R"
#print quality_score_sequence


#for reads_ID <= paired_end_reads_number_each_read:
for reads_ID in range(1, (paired_end_reads_number_each_read+1)):
	first_read_start_point = 0
	first_read_end_point = 0
	first_read_sequence = ""

	second_read_start_point = 0
	second_read_end_point = 0
	second_read_coverage_sequence = ""  # store the covered DNA sequence of second read
	second_read_sequence = ""

	random_fragment_number = random.randrange(0,fragment_total_number)
	random_fragment_sequence = fragment_list[random_fragment_number]
	first_read_start_point_select_range = len(random_fragment_sequence) - 2 * read_length  - gap
	first_read_start_point = random.randrange(0,first_read_start_point_select_range)
	first_read_end_point = first_read_start_point + read_length
	first_read_sequence = random_fragment_sequence[first_read_start_point:first_read_end_point]
	
	# keep a record of no error reads
	#paired_end_first_read_no_rerror_output.write("@"+chromosome_file_name +"_seq1_" + str(reads_ID) + "\t length = " + str(read_length) +"\n")
	#paired_end_first_read_no_rerror_output.write(first_read_sequence + "\n")
	#paired_end_first_read_no_rerror_output.write("+"+chromosome_file_name +"_seq1_" + str(reads_ID) + "\t length = " + str(read_length) +"\n")
	#paired_end_first_read_no_rerror_output.write(quality_score_sequence + "\n")
	paired_end_first_read_list.append(first_read_sequence)

	second_read_start_point = first_read_end_point + gap + read_length
	second_read_end_point = second_read_start_point - read_length
	second_read_coverage_sequence = random_fragment_sequence[second_read_end_point:second_read_start_point]  # opposite direction, end first, then second
	second_read_sequence = convertSecondRead_2(second_read_coverage_sequence)
	
	# keep a record of no error reads
	#paired_end_second_read_no_rerror_output.write("@"+chromosome_file_name +"_seq1_" + str(reads_ID) + "\t length = " + str(read_length) + "\n")
	#paired_end_second_read_no_rerror_output.write(second_read_sequence + "\n")
	#paired_end_second_read_no_rerror_output.write("+"+chromosome_file_name +"_seq1_" + str(reads_ID) + "\t length = " + str(read_length) + "\n")
	#paired_end_second_read_no_rerror_output.write(quality_score_sequence + "\n")
	paired_end_second_read_list.append(second_read_sequence)

#paired_end_first_read_no_rerror_output.close()
#paired_end_second_read_no_rerror_output.close()


paired_end_file_name_first_read_error_record = chromosome_file_name + "_paired_end_first_error_record.txt"
paired_end_file_name_second_read_error_record = chromosome_file_name + "_paired_end_second_error_record.txt"

paired_end_file_name_first_read_error_record_output = open(currentPath + paired_end_file_name_first_read_error_record, "w")
paired_end_file_name_first_read_error_record_output.write("reads ID \t allele_position \t From \t To \n")

paired_end_file_name_second_read_error_record_output = open(currentPath + paired_end_file_name_second_read_error_record, "w")
paired_end_file_name_second_read_error_record_output.write("reads ID \t allele_position \t From \t To \n")

paired_end_first_read_error_dic = {'': []}
paired_end_second_read_error_dic = {'': []}

base_list = ["A", "T", "C", "G"]

for error in range (0, total_error_number):
	list_selector = random.randrange(0, 2)
	if list_selector == 0:  # choose first read list
		reads_selector = random.randrange(0, len(paired_end_first_read_list))
		if reads_selector not in paired_end_first_read_error_dic:
			paired_end_first_read_error_dic[reads_selector]=[]
			temp_reads = paired_end_first_read_list[reads_selector]
			allele_selector = random.randrange(0, len(temp_reads))
			paired_end_first_read_error_dic[reads_selector].append(allele_selector)
			base_selector = random.randrange(0, len(base_list))
			while temp_reads[allele_selector] == base_list[base_selector]:  # make sure it won't be the same base
				base_selector = random.randrange(0, len(base_list))
			#print "before: ", paired_end_first_read_list[reads_selector]
			#print "allele_selector", allele_selector
			#print "base_selector", base_selector			
			paired_end_first_read_list[reads_selector] = temp_reads[:allele_selector] + base_list[base_selector] + temp_reads[(allele_selector+1):]
			paired_end_file_name_first_read_error_record_output.write(str((reads_selector+1))+"\t"+str(allele_selector)+"\t"+temp_reads[allele_selector]+"\t"+base_list[base_selector]+"\n")
			#print "after: ", paired_end_first_read_list[reads_selector]
		else:
			temp_reads = paired_end_first_read_list[reads_selector]
			allele_selector = random.randrange(0, len(temp_reads))
			while allele_selector in paired_end_first_read_error_dic[reads_selector]: # choose a different allele if it is already changed
				allele_selector = random.randrange(0, len(temp_reads))
			paired_end_first_read_error_dic[reads_selector].append(allele_selector)
			base_selector = random.randrange(0, len(base_list))
			while temp_reads[allele_selector] == base_list[base_selector]:  # make sure it won't be the same base
				base_selector = random.randrange(0, len(base_list))
				
			#print "before: ", paired_end_first_read_list[reads_selector]
			#print "allele_selector", allele_selector
			#print "base_selector", base_selector			
			paired_end_first_read_list[reads_selector] = temp_reads[:allele_selector] + base_list[base_selector] + temp_reads[(allele_selector+1):]
			paired_end_file_name_first_read_error_record_output.write(str((reads_selector+1))+"\t"+str(allele_selector)+"\t"+temp_reads[allele_selector]+"\t"+base_list[base_selector]+"\n")
			#print "after: ", paired_end_first_read_list[reads_selector]	
			
		#paired_end_first_read_list[reads_selector] = replaceAllele(paired_end_first_read_list[reads_selector])
		#paired_end_file_name_first_read_error_record_output.write(str((reads_selector+1)) + "\n")
	
	
	elif list_selector == 1: # choose second read list
		reads_selector = random.randrange(0, len(paired_end_second_read_list))
		if reads_selector not in paired_end_second_read_error_dic:
			paired_end_second_read_error_dic[reads_selector]=[]
			temp_reads = paired_end_second_read_list[reads_selector]
			allele_selector = random.randrange(0, len(temp_reads))
			paired_end_second_read_error_dic[reads_selector].append(allele_selector)
			base_selector = random.randrange(0, len(base_list))
			while temp_reads[allele_selector] == base_list[base_selector]:  # make sure it won't be the same base
				base_selector = random.randrange(0, len(base_list))
			#print "before: ", paired_end_second_read_list[reads_selector]
			#print "allele_selector", allele_selector
			#print "base_selector", base_selector			
			paired_end_second_read_list[reads_selector] = temp_reads[:allele_selector] + base_list[base_selector] + temp_reads[(allele_selector+1):]
			paired_end_file_name_second_read_error_record_output.write(str((reads_selector+1))+"\t"+str(allele_selector)+"\t"+temp_reads[allele_selector]+"\t"+base_list[base_selector]+"\n")
			#print "after: ", paired_end_second_read_list[reads_selector]
		else:
			temp_reads = paired_end_first_read_list[reads_selector]
			allele_selector = random.randrange(0, len(temp_reads))
			while allele_selector in paired_end_second_read_error_dic[reads_selector]: # choose a different allele if it is already changed
				allele_selector = random.randrange(0, len(temp_reads))
			paired_end_second_read_error_dic[reads_selector].append(allele_selector)
			base_selector = random.randrange(0, len(base_list))
			while temp_reads[allele_selector] == base_list[base_selector]:  # make sure it won't be the same base
				base_selector = random.randrange(0, len(base_list))
			paired_end_second_read_list[reads_selector] = temp_reads[:allele_selector] + base_list[base_selector] + temp_reads[(allele_selector+1):]
			paired_end_file_name_second_read_error_record_output.write(str((reads_selector+1))+"\t"+str(allele_selector)+"\t"+temp_reads[allele_selector]+"\t"+base_list[base_selector]+"\n")
		
		
		#paired_end_second_read_list[reads_selector] = replaceAllele(paired_end_second_read_list[reads_selector])
		#paired_end_file_name_second_read_error_record_output.write(str((reads_selector+1)) + "\n")

paired_end_file_name_first_read_error_record_output.close()
paired_end_file_name_second_read_error_record_output.close()

# files containg all reads with errors
paired_end_first_read_output = open(currentPath + paired_end_file_name_first_read, "w")
paired_end_second_read_output = open(currentPath + paired_end_file_name_second_read, "w")

for reads_ID in range(1, (paired_end_reads_number_each_read+1)):
	paired_end_first_read_output.write("@"+chromosome_file_name +"_seq1_" + str(reads_ID) + "\t length = " + str(read_length) +"\n")
	paired_end_first_read_output.write(paired_end_first_read_list[reads_ID-1] + "\n")
	paired_end_first_read_output.write("+"+chromosome_file_name +"_seq1_" + str(reads_ID) + "\t length = " + str(read_length) +"\n")
	paired_end_first_read_output.write(quality_score_sequence + "\n")

	paired_end_second_read_output.write("@"+chromosome_file_name +"_seq1_" + str(reads_ID) + "\t length = " + str(read_length) + "\n")
	paired_end_second_read_output.write(paired_end_second_read_list[reads_ID-1] + "\n")
	paired_end_second_read_output.write("+"+chromosome_file_name +"_seq1_" + str(reads_ID) + "\t length = " + str(read_length) + "\n")
	paired_end_second_read_output.write(quality_score_sequence + "\n")

paired_end_first_read_output.close()
paired_end_second_read_output.close()

end = time.time()
run_time = str(end - start)
run_time = run_time[:(run_time.find('.')+3)]
print chromosome_file_name + " run time: " + run_time + "s"

# output data record file
data_record_file = open(currentPath + data_record_file_name, "w")
data_record_file.write("Chromosome coverage percentage: " + str(coveragePercentage*100) + "% \n")
data_record_file.write("read_length: " + str(read_length) + "\n")
data_record_file.write("gap: " + str(gap) + "\n")
data_record_file.write("depth: " + str(depth) + "\n")
data_record_file.write("error_rate: " + str(error_rate*100) + "% \n")
data_record_file.write("total run time: " + run_time + "s \n")
data_record_file.write("\n")

data_record_file.write("total base number in " + chromosome_file + ": " + str(total_base_number) + "\n")
data_record_file.write("total fragment number is: " + str(fragment_total_number) + "\n")
data_record_file.write("total_reads_number: " + str(total_reads_number) + "\n")
data_record_file.write("total reads number in each pair-end file: " + str(paired_end_reads_number_each_read) + "\n")
data_record_file.write("total_error_number: " + str(total_error_number) + "\n")
data_record_file.write("\n")

data_record_file.write("Fragment example. \n")
data_record_file.write("first_read_start_point_select_range: " + str(first_read_start_point_select_range) + "\n")
data_record_file.write("random_fragment_sequence length: " + str(len(random_fragment_sequence)) + "\n")
data_record_file.write("random_fragment_sequence: " + random_fragment_sequence + "\n")
data_record_file.write("\n")

data_record_file.write("first read start position: " + str(first_read_start_point) + "\n")
data_record_file.write("first read end position: " + str(first_read_end_point) + "\n")
data_record_file.write("first read length: " + str(len(first_read_sequence)) + "\n")
data_record_file.write("first read sequence: " + first_read_sequence + "\n")
data_record_file.write("\n")

data_record_file.write("second read start position: " + str(second_read_start_point) + "\n")
data_record_file.write("second read end position: " + str(second_read_end_point) + "\n")
data_record_file.write("second read length: " + str(len(second_read_coverage_sequence)) + "\n")
data_record_file.write("second read covered chromosome sequence: \n" + second_read_coverage_sequence + "\n")
data_record_file.write("second read sequence: \n" + convertSecondRead_2(second_read_coverage_sequence) + "\n")
data_record_file.write("\n")

data_record_file.close()

os.system("rm -rf " + chromosome_file_name)
os.system("mkdir " + chromosome_file_name)
os.system("mv "+chromosome_file_name + "_* " + chromosome_file_name + "/")

depth_folder_name = chromosome_file[:chromosome_file.find('.')] + "_" + str(depth) + "x"
os.system("mkdir " + depth_folder_name)
os.system("mv "+chromosome_file_name + "* " + depth_folder_name + "/")

chromosome_folder_name = chromosome_file[:chromosome_file.find('.')]
os.system("mkdir " + chromosome_folder_name)
os.system("mv "+depth_folder_name + "* /" + chromosome_folder_name + "/")
