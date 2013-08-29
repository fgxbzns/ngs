#!/usr/bin/python
#######################################################################################
# find enzyme cut seq in dna seq
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser

import re
from tools import *

def enzyme_search(enzyme_seq, DNA_seq):
	matchObj = re.search( enzyme_seq, DNA_seq, re.M|re.I)
	if matchObj:
	   #print "find in", DNA_seq, matchObj.group()
	   return True
	else:
	   #print "No match!!"
	   return False

def seq_convert(enzyme_seq):
	temp_seq = ""
	symbol_dict = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'U':'U', 'W':'[AT]', 'S':'[CG]', \
				'M':'[AC]', 'K':'[GT]', 'R':'[AG]', 'Y':'[CT]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'}	
	for letter in enzyme_seq:
		temp_seq += symbol_dict[letter]
	return temp_seq

def load_enzyme_seq(file_name):
	en_list = []
	fp2 = open(file_name,'r')
	for line in fp2:
		elements = line.strip().split()
		try:
			enzyme_name = elements[0]
			enzyme_cut = seq_convert(elements[1])
			en_list.append((enzyme_name, enzyme_cut))
		except:
			print "error in: ", line
	fp2.close()
	print en_list
	return en_list

def load_dna_seq(file_name):
	fp1=open(file_name, 'r')
	seq_list = []
	for line in fp1:
	    elements = line.strip().split()
	    pos=elements[0]
	    seq=elements[1]
	    recseq=elements[2]
	    seq_list.append((pos,seq,recseq))
	fp1.close()
	return seq_list

def check_dna_seq(en_list, seq_list):
	output3=open('output.txt','w') 
	for pos, seq, recseq in seq_list:
	    for enzyme_name, enzyme_cut in en_list:
	    	cut_in_seq = enzyme_search(enzyme_cut, seq)
	    	cut_in_recseq = enzyme_search(enzyme_cut, recseq)
	    	if (cut_in_seq and not cut_in_recseq) or (not cut_in_seq and cut_in_recseq):
	            print >> output3, pos, seq, recseq, enzyme_name, enzyme_cut
	output3.close()

def get_args():
	desc="calculate_maf"
	usage = "" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-e", "--en", type="string", dest="en_name",help = "Input enzyme file name", default="null")
	parser.add_option("-s", "--seq", type="string", dest="seq_name",help = "Input seq file name", default="null")
	(options, args) = parser.parse_args()
	if options.en_name == "null" or options.seq_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__=='__main__':
	options = get_args()
	en_file_name = options.en_name
	seq_file_name = options.seq_name
	"""
	enzyme_seq = "[CG][CA]"
	DNA_seq = "CCCGGTCCGACAAAAATGA"
	enzyme_search(enzyme_seq, DNA_seq)
	enzyme_seq = "GACNNNNRTGA"
	enzyme_seq = seq_convert(enzyme_seq)
	print enzyme_seq
	enzyme_search(enzyme_seq, DNA_seq)
	"""
	check_dna_seq(load_enzyme_seq(en_file_name), load_dna_seq(seq_file_name))

