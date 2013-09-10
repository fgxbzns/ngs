#!/usr/bin/python
#######################################################################################
# Common tools
#######################################################################################
from __future__ import division
import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser
from ctypes import *
#from hifiAccuCheck_v2 import hifiAccuCheck


"""cannot import other files"""

file_path = "/home/guoxing/disk2/solid/common_files/"
program_path = "/home/guoxing/disk2/ngs/morehouse/python/"
bash_path = "/home/guoxing/disk2/ngs/bash/"
data_record_path = "/home/guoxing/disk2/solid/common_files/data_record/"
currentPath = os.getcwd() + '/'

raw_data_format = "list"

quality_score_dict = { '!':0, '\"':1, '#':2, '$':3, '%':4, '&':5, '\'':6, '(':7, 
					')':8, '*':9, '+':10, ',':11, '-':12, '.':13 }

def usage():
	print "%s [seed_file] [chr]" % sys.argv[0]

def load_raw_data(file_name, raw_data_format = "list"):		
	title_info = ""
	data = {}
	fp = open(file_name, "r")
	for line in fp:
		if line.startswith("rsID"):
			title_info = line.strip()
		else:
			elements = line.strip().split()
			try:
				# convert the position to int for sorting
				if raw_data_format == "list":
					data[int(elements[1])] = elements
				elif raw_data_format == "string":
					data[int(elements[1])] = line.strip()
			except ValueError:
				#print "error in ", line
				pass
	fp.close()
	return (title_info, data)

def hifi_run(file_name, chr_name):
	hifi = program_path + "hifi_fu_ref " + file_name + " genotype.txt refHaplos.txt 0.10"
	hifi_process = subprocess.Popen(hifi, shell=True)
	hifi_process.wait()
	#hifiAccuCheck("imputed_" + file_name, chr_name)

def hifi_mlp(hap_file_name="haplotype.txt", geno_file_name="genotype.txt", ref_file_name="refHaplos.txt", MAFSTEP = 0.1):
	"""
	add extern "C" in front of the function in cpp, otherwise it will be claimed as undefined.
	gcc hifi_fu_python.cpp -fPIC -shared -o libhifi_fu.so
	
	"""
	process = cdll.LoadLibrary("/home/guoxing/disk2/ngs/morehouse/other/libhifi_fu.so")
	process.hifi_m.argtypes = [c_char_p, c_char_p, c_char_p, c_float]
	process.hifi_m(hap_file_name, geno_file_name, ref_file_name, MAFSTEP)
	#hifiAccuCheck("imputed_" + hap_file_name, "chr9")

# remove the snps with "N" in std hap
def removeN(hifi_std_dict):
	temp_dict = {}
	for position, elements in hifi_std_dict.iteritems():
		if elements[2].strip() != "N" and elements[3].strip() != "N":
			temp_dict[position] = hifi_std_dict[position]
	return temp_dict

def list_to_line(list):
	line = ""
	for a in list:
		line += a + "\t"
	return line
	"""
	a = ['a', 'b', 'c']
	print "".join(a)
	"""

def dict_substract(large_dict, small_dict):
	return {index:value for index, value in large_dict.iteritems() if index not in small_dict}

def dict_add(revised_seed_dict, recovered_seed_dict):
	for position, seed in recovered_seed_dict.iteritems():
		if position not in revised_seed_dict:
			revised_seed_dict[position] = seed
		else: 
			#print position
			pass
	return revised_seed_dict

def sort_dict_by_key(input_dict):
	sorted_list = []
	sorted_list = [x for x in input_dict.iteritems()] 
	sorted_list.sort(key=lambda x: x[0]) # sort by key
	#sorted_list.sort(key=lambda x: x[1]) # sort by value
	return sorted_list

def wccount(file_name):
    out = subprocess.Popen(['wc', '-l', file_name],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])
   
def getTotalBaseNum(fileName):
	totalBase = 0
	f = open(currentPath+fileName, "r")
	for line in f:
		if not line.startswith('>'):
			totalBase += len(line.strip())
	return totalBase
	f.close()   
 
def keywithmaxval(dict):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value """  
     v=list(dict.values())
     k=list(dict.keys())
     return k[v.index(max(v))]

# group the seed into homo and hetero groups
def group_seed(seed_dict, geno_dict):
	seed_homo_dict = {}
	seed_hetero_dict = {}
	for position, snp in seed_dict.iteritems():
		if position in geno_dict:	# pos in seed may not be in geno
			geno_allele = geno_dict[position][2]
			if geno_allele[0] == geno_allele[1]:
				seed_homo_dict[position] = snp
			else:
				seed_hetero_dict[position] = snp
		else:
			seed_hetero_dict[position] = snp
	return (seed_homo_dict, seed_hetero_dict)
"""
def time():
	start = time.time()	
	end = time.time()
	run_time = str(format((end - start), "0.3f"))
	print "run time is: " + run_time + "s"
"""
# group the seed into homo and hetero groups
def group_seed(seed_dict, geno_dict):
	seed_homo_dict = {}
	seed_hetero_dict = {}
	for position, snp in seed_dict.iteritems():
		if position in geno_dict:	# pos in seed may not be in geno
			geno_allele = geno_dict[position][2]
			if geno_allele[0] == geno_allele[1]:
				seed_homo_dict[position] = snp
			else:
				seed_hetero_dict[position] = snp
		else:
			seed_hetero_dict[position] = snp
	return (seed_homo_dict, seed_hetero_dict)
