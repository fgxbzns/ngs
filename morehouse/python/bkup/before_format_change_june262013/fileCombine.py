#!/usr/bin/python

import os, glob, subprocess, random, operator, time
from optparse import OptionParser

class info:
	def __init__(me, file_name, data):
		me.file_name = file_name
		me.data = data

class tag:
	def __init__(me, file_name, tag_info, tag_blank):
		me.file_name = file_name
		me.tag_info = tag_info
		me.tag_blank = tag_blank

currentPath = os.getcwd() + '/'
file_name_list = []
snp_dic = {'': []}

for infile in glob.glob(os.path.join(currentPath,'tag*.txt')):
	file_name = infile[(infile.find('tag')):].strip()
	inputFile_tag = open(currentPath + file_name, "r")
	for line in inputFile_tag:
		if not line.startswith("chr19"):
			elements = line.strip().split()
			rsID = elements[2].strip()
			tag_info = line[line.find(rsID):].strip()
			tag_blank = ""
			temp_info = tag_info.strip().split()
			for a in temp_info:
				tag_blank += "\t"
			file_name_list.append(tag(file_name, tag_info, tag_blank))
		if line.startswith("chr19"):
			elements = line.strip().split()
			position = int(elements[1].strip())
			if position not in snp_dic:
				snp_dic[position]=[]		
			rsID = elements[2].strip()
			data = line[line.find(rsID):].strip()
			snp_dic[position].append(info(file_name, data))
	inputFile_tag.close()

outputFile = open(currentPath + "output.txt", "w")
outputFile.write("chr19 \t position")

print file_name_list[0].tag_info
print len(file_name_list[0].tag_info)
print file_name_list[0].tag_blank
print len(file_name_list[0].tag_blank)

print len(snp_dic)
print len(file_name_list)

for tag in file_name_list:
	outputFile.write("\t" + tag.tag_info)
	#print tag.file_name
outputFile.write("\n")
"""
for position, info_list in snp_dic.iteritems():
	outputFile.write("chr19 \t" + str(position))
	for info in info_list:
		for tag in file_name_list:
			if tag.file_name == info.file_name:
				outputFile.write("&" + info.data)
			else:
				outputFile.write("&" + tag_blank)
	outputFile.write("\n")
"""

for position, info_list in snp_dic.iteritems():
	outputFile.write("chr19 \t" + str(position) + "\t")
	for info in info_list:
		exist = False
		for tag in file_name_list:
			if tag.file_name == info.file_name:
				outputFile.write(info.data)
				exist = True
			if not exist:
			#print "no exist"
				outputFile.write(tag_blank)
	outputFile.write("\n")

outputFile.close()
