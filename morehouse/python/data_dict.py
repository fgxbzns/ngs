#!/usr/bin/python
#######################################################################################
# Common tools
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

seed_dict = {}
geno_dict = {}
hap_ref_dict = {}   
hap_std_dict = {} 
hifi_result_dict = {}  

class seeds:
    def __init__(self):
        self.rsID = ""
        self.position = 0
        self.allele_ori = ""
        self.allele_new = ""
        self.allele_new_percentage = 0
        self.allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}

def load_seed_data(file_name):
    seed_dict = {}
    data_tuple = load_raw_data(file_name)
    seed_title_info = data_tuple[0]
    data_dict = data_tuple[1]
    for position, elements in data_dict.iteritems():
        seed = seeds()
        seed.rsID = elements[0].strip()
        seed.position = int(elements[1].strip())
        seed.allele_ori = elements[2].strip()
        seed_dict[position] = seed
    return (seed_title_info, seed_dict)

def load_hap_std(file_name):
    hap_std_dict = {}
    data_dict = load_raw_data(file_name)[1]
    for position, elements in data_dict.iteritems():
        if elements[2].strip() != "N" and elements[3].strip() != "N":
            try:
                position = elements[1].strip()
                position = int(position)
                hap_std_dict[position] = elements
            except ValueError:
                print "error in ", file_name, position
    return hap_std_dict    