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

(options, args) = parser.parse_args()

fastq_file = options.inputfile
fastq_file_name = fastq_file[:(fastq_file.find('fastq')-1)].strip()

print "fastq_file_name: ", fastq_file_name

ref_file_name = ref_path + "hg18chr.fa"
#ref_file_name = "hg18chr.fa"

sai_file = fastq_file_name + ".sai"
sam_file = fastq_file_name + ".sam"
"""
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


