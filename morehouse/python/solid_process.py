#!/usr/bin/python

# run bwa

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

program_path = "/home/guoxing/tool/morehouse/python/"
other_path = "/home/guoxing/tool/morehouse/other/"
samtools_path = other_path + "samtools-0.1.18/"
ref_path = "/home/guoxing/disk2/UCSC_hg18_index_lm/"
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 

parser.add_option("-c", "--cFile", type="string", dest="csfastaFile",help = "Input csfasta Name", default="null")
parser.add_option("-q", "--qFile", type="string", dest="qualFile",help = "Input qual Name", default="null")
parser.add_option("-o", "--output", type="string", dest="fastqFile",help = "output fastq file ", default="null")

(options, args) = parser.parse_args()

csfasta_file_name = options.csfastaFile
qual_file_name = options.qualFile
fastq_file_name = options.fastqFile

# solid2fastq
solid2fastq = program_path + "solid2fastq.py " + " -c " + csfasta_file_name + " -q " + qual_file_name + " -o " + fastq_file_name
print "solid2fastq runing"
print solid2fastq
solid2fastq_Process = subprocess.Popen(solid2fastq, shell=True)
solid2fastq_Process.wait()

# primer remove
print "primerRemove runing"
priRem_input_file_name = fastq_file_name + ".fastq"
primerRemove = program_path + "primerRemove.py " + " -i " + priRem_input_file_name
print primerRemove
primerRemove_Process = subprocess.Popen(primerRemove, shell=True)
primerRemove_Process.wait()

# BWA 
print "BWA runing"
bwa_input_file = fastq_file_name + "_prem.fastq"
bwa_input_file_name = bwa_input_file[:bwa_input_file.find('.')].strip()
print "bwa_input_file: ", bwa_input_file

ref_file_name = ref_path + "hg18chr.fa"

sai_file = bwa_input_file_name + ".sai"
sam_file = bwa_input_file_name + ".sam"

bwa_sai = "bwa" + " aln " + ref_file_name + " " + bwa_input_file + " > " + sai_file
print bwa_sai
bwa_sai_Process = subprocess.Popen(bwa_sai, shell=True)
bwa_sai_Process.wait()

bwa_sam = "bwa" + " samse " + ref_file_name + " " + sai_file + " " + bwa_input_file + " > " + sam_file
print bwa_sam
bwa_sam_Process = subprocess.Popen(bwa_sam, shell=True)
bwa_sam_Process.wait()

"""
# for long reads data, like 454
bwa_long = "bwa" + " bwasw " + ref_file_name + " " + fastq_file + " > " + sam_file
print bwa_long

bwa_long_Process = subprocess.Popen(bwa_long, shell=True)
bwa_long_Process.wait()
"""
#sam_file = "song_1_prem.sam"
#bwa_input_file_name = "song_1_prem"

# grep
print "grep chr from sam file"
grep = "grep chr " + sam_file + " > " + bwa_input_file_name + "_chr.sam"
print grep
grep_Process = subprocess.Popen(grep, shell=True)
grep_Process.wait()

# chr percentage
print "chrPercentage runing"
chrPercentage = program_path + "chrPercentage.py " + " -i " + bwa_input_file_name + "_chr.sam"
print chrPercentage
chrPercentage_Process = subprocess.Popen(chrPercentage, shell=True)
chrPercentage_Process.wait()

# sort the sam file
print "sort runing"
sort_input = bwa_input_file_name + "_chr"
sam2bam = samtools_path + "samtools view -bS " + sort_input + ".sam > " + sort_input + ".bam"
print sam2bam
sam2bam_Process = subprocess.Popen(sam2bam, shell=True)
sam2bam_Process.wait()

sortbam = samtools_path + "samtools sort " + sort_input + ".bam " + sort_input + "_sorted"
print sortbam
sortbam_Process = subprocess.Popen(sortbam, shell=True)
sortbam_Process.wait()

bam2sam = samtools_path + "samtools view " + sort_input + "_sorted.bam > " + sort_input + "_sorted.sam"
print bam2sam
bam2sam_Process = subprocess.Popen(bam2sam, shell=True)
bam2sam_Process.wait()

chr_name = ""
chr_percentage_file = open(currentPath + bwa_input_file_name + "_chr.sam_percentage",'r')
line_1 = chr_percentage_file.readline()
line_2 = chr_percentage_file.readline()
print "chr with highest percentage: ", line_2.strip()
elements = line_2.strip().split()
chr_name = elements[0].strip()
chr_percentage_file.close()

# grep chr#
print "grep " + chr_name +" from sam file"
grep = "grep " + chr_name + " " + sort_input + "_sorted.sam > " + bwa_input_file_name + "_" + chr_name + "_sorted.sam"
print grep
grep_Process = subprocess.Popen(grep, shell=True)
grep_Process.wait()

# repeat masker remove
repeatRemove_input = bwa_input_file_name + "_" + chr_name + "_sorted.sam"
repeatRemove = program_path + "repeatRemove_sorted.py -s " + repeatRemove_input
print repeatRemove
repeatRemove_Process = subprocess.Popen(repeatRemove, shell=True)
repeatRemove_Process.wait()

# chr percentage
print "chrPercentage runing"
chrPercentage = program_path + "chrPercentage.py " + " -i " + bwa_input_file_name + "_" + chr_name + "_sorted_rmsk.sam"
print chrPercentage
chrPercentage_Process = subprocess.Popen(chrPercentage, shell=True)
chrPercentage_Process.wait()




