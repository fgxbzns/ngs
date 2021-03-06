#!/usr/bin/python

# run bwa

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

program_path = "/home/guoxing/disk2/ngs/morehouse/python/"
other_path = "/home/guoxing/disk2/ngs/morehouse/other/"
samtools_path = other_path + "samtools-0.1.18/"
ref_path = "/home/guoxing/disk2/UCSC_hg18_index_lm/"
data_record_path = "/home/guoxing/disk2/solid/common_files/data_record/"

currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
"""
parser.add_option("-c", "--cFile", type="string", dest="csfastaFile",help = "Input csfasta Name", default="null")
parser.add_option("-q", "--qFile", type="string", dest="qualFile",help = "Input qual Name", default="null")
parser.add_option("-o", "--output", type="string", dest="fastqFile",help = "output fastq file ", default="null")
"""
parser.add_option("-s", "--rmskFile", type="string", dest="rmskFile",help = "Input rmsk Name", default="null")
parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="chr11")
parser.add_option("-d", "--depth", type="string", dest="depth",help = "Input depth threshold", default="3")


(options, args) = parser.parse_args()
"""
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
"""
"""
# priRem_input_file_name = "song_1_prep.fastq"
priRem_input_file_name = options.bwaFile
bwa_input_file = priRem_input_file_name
bwa_input_file_name = priRem_input_file_name[:priRem_input_file_name.find('.')].strip()

"""
"""
# BWA 
print "BWA runing"
bwa_input_file = priRem_input_file_name
bwa_input_file_name = priRem_input_file_name[:priRem_input_file_name.find('.')].strip()
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
"""
# for long reads data, like 454
bwa_long = "bwa" + " bwasw " + ref_file_name + " " + fastq_file + " > " + sam_file
print bwa_long

bwa_long_Process = subprocess.Popen(bwa_long, shell=True)
bwa_long_Process.wait()
"""
#sam_file = "song_1_prem.sam"
#bwa_input_file_name = "song_1_prem"

#bwa_input_file_name = options.bwaFile
#sam_file = bwa_input_file_name + ".sam"
"""
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
"""
"""
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

"""

rmsk_input_file = options.rmskFile
chr_name = options.chrName
depth = options.depth

# priRem_input_file_name = "song_1_prep.fastq"
#priRem_input_file_name = options.bwaFile
#bwa_input_file = priRem_input_file_name
rmsk_input_file_name = rmsk_input_file[:rmsk_input_file.find('.')].strip()
"""
# sam process
print "sam_process runing"
sam_process = program_path + "sam_process.py " + " -s " + rmsk_input_file
print sam_process
sam_process_Process = subprocess.Popen(sam_process, shell=True)
sam_process_Process.wait()
"""
# snpPick_solid
print "snpPick_solid runing"
#snpPick_solid = program_path + "snpPick_solid.py " + " -s " + rmsk_input_file_name + ".sam -c " + chr_name

#snpPick_solid = program_path + "snpPick_solid.py " + " -s " + rmsk_input_file_name + "_indel.sam -c " + chr_name + " -d " + depth		#sam_processed
snpPick_solid = program_path + "snpPick_solid.py " + " -s " + rmsk_input_file_name + ".sam -c " + chr_name + " -d " + depth		#sam_processed
print snpPick_solid
snpPick_solid_Process = subprocess.Popen(snpPick_solid, shell=True)
snpPick_solid_Process.wait()


data_record_file_name = "solid_process_4.txt"
try:
	data_record_file = open(data_record_path + data_record_file_name, "r")
except:
	data_record_file = open(data_record_path + data_record_file_name, "w")
	print >> data_record_file, "all", "sam_file", "chr", "depth_threshold", "same_to_A", "same_to_B", "correct_rate", \
	"same_to_AB", "X_AB", "not_ABX", "hifi_seed", "A_seed_rate", "hete_seed_rate", "pure_total", "not_in_genotype", "hifi", "error", "total", "accuracy"
data_record_file = open(data_record_path + data_record_file_name, "a")
"""
print >> data_record_file, "==================================================="
print >> data_record_file, "time is: " + str(time.time())
print >> data_record_file, "current path is: " + currentPath
print >> data_record_file, "chr_name is: " + chr_name
print >> data_record_file, "depth is: " + depth
"""

data_record_file.close()
cmd = "grep data " + rmsk_input_file_name + "_"+depth+"_data_record.txt >> " + data_record_path + data_record_file_name
print cmd
os.system(cmd)

# refMerger
print "refMerger runing"
refMerger = program_path + "refMerger_v4.py " + " -i " + rmsk_input_file_name + "_"+depth+"_hifi.txt -c " + chr_name
print refMerger
refMerger_Process = subprocess.Popen(refMerger, shell=True)
refMerger_Process.wait()

# hifi
print "hifi runing"
hifi = program_path + "hifi &"
print hifi
hifi_Process = subprocess.Popen(hifi, shell=True)
hifi_Process.wait()


# hifiAccuCheck
print "hifiAccuCheck runing"
hifiAccuCheck = program_path + "hifiAccuCheck_v2.py -c " + chr_name
print hifiAccuCheck
hifiAccuCheck_Process = subprocess.Popen(hifiAccuCheck, shell=True)
hifiAccuCheck_Process.wait()



