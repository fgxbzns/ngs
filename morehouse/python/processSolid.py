#!/usr/bin/python

# location /home/guoxing/tool/morehouse

# pre-process solid file for solid2fastaq. To remove "1" in csfasta file

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

program_path = "/home/guoxing/tool/morehouse"
file_path = "/home/guoxing/morehouse_new_hap/fils/"
currentPath = os.getcwd() + '/'


# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--input", type="string", dest="csfastaFile",help = "Input File Name", default="null")

(options, args) = parser.parse_args()

csfasta_File = options.csfastaFile

inputFile_csfasta = open(currentPath + csfasta_File, "r")
outputFile_csfasta = open(currentPath + csfasta_File + "_processed.csfasta", "w")

for line in inputFile_csfasta:
	if line.startswith(">"):
		elements = line.strip().split()
		read_ID = elements[0].strip()	
		outputFile_csfasta.write(read_ID + "\n")
	else:
		outputFile_csfasta.write(line + "\n")
		
inputFile_csfasta.close()
outputFile_csfasta.close()
		
