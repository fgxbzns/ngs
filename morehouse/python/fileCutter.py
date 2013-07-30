#!/usr/bin/python

# take certain lines of a file

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

currentPath = os.getcwd() + '/'

def getTotalBaseNum(fileName):
	totalBase = 0
	f = open(currentPath+fileName, "r")
	for line in f:
		if not line.startswith('>'):
			totalBase += len(line.strip())
	return totalBase
	f.close()

def wccount(file_name):
    out = subprocess.Popen(['wc', '-l', file_name],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-i", "--inputfile", type="string", dest="inputfile",help = "Input File Name", default="null")
parser.add_option("-l", "--lineNumber", type="string", dest="lineNumber",help = "lines to take", default="null")
parser.add_option("-s", "--startLine", type="string", dest="startLine",help = "start line", default="null")
parser.add_option("-e", "--endLine", type="string", dest="endLine",help = "end line", default="null")
parser.add_option("-t", "--totalLine", dest="totalLine",help = "total number of lines", default=False)
parser.add_option("-p", "--percentage", type="string", dest="percentage",help = "percentage of total lines", default="null")

(options, args) = parser.parse_args()

inputfile = options.inputfile
lineNumber = options.lineNumber
startLine = options.startLine
endLine = options.endLine
totalLine = options.totalLine
percentage = options.percentage

input_file_name = inputfile[:inputfile.find('.')]
#output_file_name = input_file_name + "_" + str(lineNumber) + ".txt"

inputFile_ori = open(currentPath + inputfile, "r")
#outputFile = open(currentPath + output_file_name, "w")

if lineNumber != "null":
	lineNumber = int(lineNumber)
	output_file_name = input_file_name + "_" + str(lineNumber) + ".txt"
	outputFile = open(currentPath + output_file_name, "w")
	i = 0
	while i < lineNumber:
		line = inputFile_ori.readline().strip()
		outputFile.write(line + "\n")
		i += 1
	outputFile.close()

if startLine != "null" and endLine != "null":
	startLine = int(startLine)
	endLine = int(endLine)
	output_file_name = input_file_name + "_" + str(startLine) + "_" + str(endLine) + ".txt"
	outputFile = open(currentPath + output_file_name, "w")
	i = 0
	while i < startLine:
		inputFile_ori.readline()
		i += 1
	while i >= startLine and i <= endLine:
		line = inputFile_ori.readline().strip()
		outputFile.write(line + "\n")
		i += 1
	outputFile.close()

if totalLine:
	"""
	i = 0
	for line in inputFile_ori:
		i += 1
	print "total line is: ", i
	"""
	print "total line is: ", wccount(inputfile)
	
if percentage != "null":
	percentage = float(percentage)
	output_file_name = input_file_name + "_" + str(percentage) +  ".txt"
	outputFile = open(currentPath + output_file_name, "w")
	total_line = 0
	for line in inputFile_ori:
		total_line += 1
	lineNumber = int(total_line * percentage)
	i = 0
	while i < lineNumber:
		line = inputFile_ori.readline().strip()
		outputFile.write(line + "\n")
		i += 1
	outputFile.close()
	

inputFile_ori.close()
#outputFile.close()
