#!/usr/bin/python
#######################################################################################
# Common tools
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

file_path = "/home/guoxing/disk2/solid/common_files/"
program_path = "/home/guoxing/disk2/ngs/morehouse/python/"
data_record_path = "/home/guoxing/disk2/solid/common_files/data_record/"
currentPath = os.getcwd() + '/'


def usage():
	print "%s [seed_file] [chr]" % sys.argv[0]

def load_raw_data(file_name):		
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
				data[int(elements[1])] = elements
			except ValueError:
				#print "error in ", line
				pass
	fp.close()
	return (title_info, data)

def sort_dict_by_key(input_dict):
	sorted_list = []
	sorted_list = [x for x in input_dict.iteritems()] 
	sorted_list.sort(key=lambda x: x[0]) # sort by key
	return sorted_list

def wccount(filename):
    out = subprocess.Popen(['wc', '-l', filename],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])
    
def keywithmaxval(dict):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value """  
     v=list(dict.values())
     k=list(dict.keys())
     return k[v.index(max(v))]
