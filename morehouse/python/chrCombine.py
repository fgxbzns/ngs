#!/usr/bin/python

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

currentPath = os.getcwd() + '/'

outputFile = open(currentPath + "hg18chr_gx.fa", "w")

i = 1
name_list = []
while i <= 22:
	name_list.append(i)
	i += 1

name_list.append('X')
name_list.append('Y')

for name in name_list:
	file_name = 'chr' + str(name) + ".fa"
	print "current file: ", file_name
	inputFile = open(file_name, "r")
	for line in inputFile:
		outputFile.write(line.strip() + "\n")
	inputFile.close()

"""
for infile in glob.glob(os.path.join(currentPath,'*.fa')):
	file_name = infile[:(infile.find('.'))].strip() + ".fa"
	print "current file: ", file_name
	inputFile = open(file_name, "r")
	for line in inputFile:
		outputFile.write(line.strip() + "\n")
"""	

outputFile.close()
