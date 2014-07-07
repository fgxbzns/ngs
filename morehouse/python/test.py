#!/usr/bin/python
#######################################################################################
# check rough similarity of two chr
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy, ctypes
from optparse import OptionParser

from ctypes import *


import os
import threading
import multiprocessing

import re
from tools import *
from hifiAccuCheck_v2 import hifiAccuCheck


def load_chr(file_name):
	chr_seq = ""
	fp = open(file_name, "r")
	for line in fp:
		chr_seq += line.strip()
	fp.close()
	return chr_seq

def select_seq(chr_seq, point_num):
	seq_dict = {}
	total_num = len(chr_seq)
	
	number_ceilling = int(math.ceil(float(total_num)/100)*100)
	num_in_each_part = number_ceilling/point_num
	
	for i in range(point_num):
		seq_dict[i*num_in_each_part] = chr_seq[i*num_in_each_part:(i*num_in_each_part+200)]
	return seq_dict

def search_seq(chr_seq, seq_dict):
	pass

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

"""
# worker function
def worker(sign, lock):
    lock.acquire()
    print(sign, os.getpid())
    lock.release()

# Main
print('Main:',os.getpid())

# Multi-thread
record = []
lock  = threading.Lock()
for i in range(5):
    thread = threading.Thread(target=hifi_process)
    thread.start()
    record.append(thread)

for thread in record:
    thread.join()
"""



def offer(queue):
	queue.put



if __name__=='__main__':
	#options = get_args()
	#file_name = options.ref_name
	"""
	enzyme_seq = "[CG][CA]"
	DNA_seq = "CCCGGTCCGACAAAAATGA"
	enzyme_search(enzyme_seq, DNA_seq)
	enzyme_seq = "GACNNNNRTGA"
	enzyme_seq = seq_convert(enzyme_seq)
	print enzyme_seq
	enzyme_search(enzyme_seq, DNA_seq)
	"""
	#hifi = cdll.LoadLibrary("/home/guoxing/disk2/ngs/morehouse/other/libhifi_fu.so")
	#so = 
	#int argc, char** argv
	#hifi.main.argtypes = [c_int,c_char_p]
	
	#hifi.main(1, "haplotype.txt")
	hap_file_name="haplotype.txt"
	chr_name = "chr9"
	"""
	hifi_process(MAFSTEP = 0.1)
	hifiAccuCheck("imputed_"+hap_file_name, chr_name)
	hifi_process(MAFSTEP = 0.5)
	hifiAccuCheck("imputed_"+hap_file_name, chr_name)
	"""
	# Multi-process
	record = []
	lock = multiprocessing.Lock()
	for i in range(5):
	    process = multiprocessing.Process(target=hifi_mlp, )
	    process.start()
	    record.append(process)
	
	for process in record:
	    process.join()
	
	print "abdcd"
	if os.path.exists("haplotype.txt"):
		print "a"





#!/usr/bin/env python

import sys
import os

if len(sys.argv) != 3:
    print "usage: %s fastq bcfile" % sys.argv[0]
    sys.exit(1)

outpath=os.path.dirname(sys.argv[2])    
seq={}
for x in open(sys.argv[2],'r'):
    v=x.rstrip('\r\n').split(' ')
    seq[v[1]]=[0,v[0],open(os.path.join(outpath,v[0]+'.txt'),'w')]

with open(sys.argv[1],'r') as fi:
    k=0
    while True:
        k+=1
        if k%1000==0:
            print 'Processing sequence '+str(k)
        row1=fi.readline()
        if not row1:
            break
        row2=fi.readline()
        row3=fi.readline()
        row4=fi.readline()
        for i in range(5,16):
            s=row2[:i]
            if s in seq:
                entry=seq[s]
                entry[0]+=1
                entry[2].write(''.join([row1,row2,row3,row4]))
                break

with open(os.path.join(outpath,'count.txt'),'w') as fo:
    for key,val in seq.items():
        fo.write('%s\t%s\t%d\n' % (val[1],key,val[0]))
        val[2].close()


"""
	song_3
	~/disk2/solid/song_3/prem_rmsk_indel/refM_test/2280/mref$
    song_5
    /home/guoxing/disk2/solid/song_5/prem_rmsk_indel/hifi_1/7766_45/7819_48
"""


sudo apt-get install r-cran-rcpp -y
sudo apt-get install r-cran-dbi -y





"""
ERROR: dependency ‘Rcpp’ is not available for package ‘plyr’
* removing ‘/home/guoxing/R/x86_64-pc-linux-gnu-library/2.14/plyr’
ERROR: dependency ‘DBI’ is not available for package ‘RSQLite’
* removing ‘/home/guoxing/R/x86_64-pc-linux-gnu-library/2.14/RSQLite’
ERROR: dependencies ‘plyr’, ‘Rcpp’ are not available for package ‘reshape2’
* removing ‘/home/guoxing/R/x86_64-pc-linux-gnu-library/2.14/reshape2’
ERROR: dependency ‘plyr’ is not available for package ‘scales’
* removing ‘/home/guoxing/R/x86_64-pc-linux-gnu-library/2.14/scales’
ERROR: dependency ‘plyr’ is not available for package ‘reshape’
* removing ‘/home/guoxing/R/x86_64-pc-linux-gnu-library/2.14/reshape’
ERROR: dependencies ‘plyr’, ‘reshape2’, ‘scales’ are not available for package ‘ggplot2’
* removing ‘/home/guoxing/R/x86_64-pc-linux-gnu-library/2.14/ggplot2’
ERROR: dependencies ‘RSQLite’, ‘reshape’, ‘ggplot2’ are not available for package ‘cummeRbund’
* removing ‘/home/guoxing/R/x86_64-pc-linux-gnu-library/2.14/cummeRbund’
"""