#!/usr/bin/python

# run bwa

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

ref_path = "/home/guoxing/disk2/UCSC_hg18_index_lm/"
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--inputfile", type="string", dest="inputfile",help = "Input File Name", default="null")
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="all")
(options, args) = parser.parse_args()

# start time
start = time.time()	

fastq_file = options.inputfile
fastq_file_name = fastq_file[:(fastq_file.find('fastq')-1)].strip()
chr_name = options.chrName

print "fastq_file_name: ", fastq_file_name

#ref_file_name = "hg18chr.fa"
ref_file_name = ref_path + "hg18chr.fa"
if not chr_name == "all":
	ref_file_name = ref_path + "hg18chr_" + chr_name + ".fa"

sai_file = fastq_file_name + "_" + chr_name + ".sai"
sam_file = fastq_file_name + "_" + chr_name + ".sam"

bwa_sai = "bwa" + " aln " + ref_file_name + " " + fastq_file + " > " + sai_file
print bwa_sai

bwa_sai_Process = subprocess.Popen(bwa_sai, shell=True)
bwa_sai_Process.wait()

bwa_sam = "bwa" + " samse " + ref_file_name + " " + sai_file + " " + fastq_file + " > " + sam_file
print bwa_sam

bwa_sam_Process = subprocess.Popen(bwa_sam, shell=True)
bwa_sam_Process.wait()
"""
bwa_long = "bwa" + " bwasw " + ref_file_name + " " + fastq_file + " > " + sam_file
print bwa_long

bwa_long_Process = subprocess.Popen(bwa_long, shell=True)
bwa_long_Process.wait()
"""

end = time.time()
run_time = str(format((end - start), "0.3f"))
print "run time is: " + run_time + "s"

data_record_file_name = fastq_file_name + "_BWAtime_" + ".txt"
data_record_file = open(currentPath + data_record_file_name, "a")
print >> data_record_file, "run time is: " + run_time + "s"
data_record_file.close()

"""

def get_args():
	desc = "remove primer"
	usage = "primerRemove -i input_fastq" 
	parser = OptionParser(usage=usage, description=desc) 
	parser.add_option("-i", "--inputFile", type="string", dest="inputFile",help = "Input File Name", default="null")
	
	
	(options, args) = parser.parse_args()
	if options.inputFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	global input_file_name
	global primer_list

	options = get_args()
	input_file_name = options.inputFile
	
	start_time = time.time()

	elapse_time = time.time() - start_time
	print "elapse_time is: " + str(format(elapse_time, "0.3f")) + "s"
"""

