#!/usr/bin/python
#######################################################################################
# Common tools
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser
from tools import *


class data_dicts:
    
    def __init__(self):
        self.chr_name = ""
        
        self.seed_file = "haplotype.txt"
        self.seed_file_name = ""
        self.seed_title_info = ""
        self.seed_dict = {}
        self.seed_homo_dict = {}
        self.seed_hetero_dict = {}
        
        self.geno_file_name = "genotype.txt"
        self.geno_title_info = ""
        self.geno_dict = {}
        self.geno_homo_dict = {}
        self.geno_hetero_dict = {}
        
        self.ref_file_name = "refHaplos.txt"
        self.ref_title_info = ""
        self.hap_ref_dict = {}
        
        self.hap_std_dict = {}
        self.ref_cluster_dict = {}
        self.cluster_pos_dict = {}
          
        self.number_of_subfile = 10
        self.ref_cycle_number = 3
        
        self.maf_upper_bound = 0.5
        self.maf_lower_bound = 0.3
        self.ref_cluster_dict = {}
        self.cluster_pos_dict = {}

    def update_seed_dict(self):
        #print "seed_file is :", self.seed_file
        self.seed_file_name = self.seed_file[:seed_file.find('.')].strip()
        self.seed_title_info, self.seed_dict = load_seed_data(self.seed_file)
        if len(self.geno_dict) == 0:
            update_geno_dict(self)
        self.seed_homo_dict, self.seed_hetero_dict = group_seed(self.seed_dict, self.geno_dict)
        print "total_seed_number: ", len(self.seed_dict)
        print "seed_homo_dict", len(self.seed_homo_dict)
        print "seed_hetero_dict", len(self.seed_hetero_dict)
        
    def update_geno_dict(self):
        genotype_file = file_path + "genotype_NA10847_" + self.chr_name + ".txt"
        self.geno_title_info, self.geno_dict = load_raw_data(genotype_file)
        self.geno_homo_dict, self.geno_hetero_dict = group_seed(self.geno_dict, self.geno_dict)
        print "total_seed_number: ", len(self.geno_dict)
        print "geno_homo_dict", len(self.geno_homo_dict)
        print "geno_hetero_dict", len(self.geno_hetero_dict)
        
    def update_ref_dict(self):
        self.ref_title_info, self.hap_ref_dict = load_raw_data(self.ref_file_name)

    def update_hap_std_dict(self):  
        hap_std_file = file_path + "ASW_"+chr_name+"_child_hap_refed.txt"    
        self.hap_std_dict = load_hap_std(hap_std_file)
        print "total_hap_std_dict_number: ", len(self.hap_std_dict)
    
    def update_ref_cluster_dict(self):
        self.ref_cluster_dict = get_cluster(self.ref_file_name, self.maf_upper_bound, self.maf_lower_bound)
        for maf_num, cluster_list in self.ref_cluster_dict.iteritems():
            for cluster_dict in cluster_list:
                for pos, ref in cluster_dict.iteritems():
                    self.cluster_pos_dict[pos] = ref
        print "cluster_pos_dict: ", len(self.cluster_pos_dict)
    
    def load_data_dicts(self):
        update_seed_dict(self)
        update_geno_dict(self)
        update_ref_dict(self)
        update_hap_std_dict(self)
        update_ref_cluster_dict(self)
        
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
    data_tuple = load_raw_data(file_name, raw_data_format)
    seed_title_info = data_tuple[0]
    data_dict = data_tuple[1]
    for position, elements in data_dict.iteritems():
        seed = seeds()
        seed.rsID = elements[0].strip()
        seed.position = int(elements[1].strip())
        seed.allele_ori = elements[2].strip()
        seed.allele_new = elements[2].strip()
        seed_dict[position] = seed
    return (seed_title_info, seed_dict)

def load_seed_data_from_dict(data_dict):
    seed_dict = {}
    for position, elements in data_dict.iteritems():
        seed = seeds()
        seed.rsID = elements[0].strip()
        seed.position = int(elements[1].strip())
        seed.allele_ori = elements[2].strip()
        seed.allele_new = elements[2].strip()
        seed_dict[position] = seed
    return seed_dict

def load_hap_std(file_name):
    hap_std_dict = {}
    data_dict = load_raw_data(file_name, raw_data_format)[1]
    for position, elements in data_dict.iteritems():
        """ ??? N X """
        if elements[2].strip() != "N" and elements[3].strip() != "N": 
            try:
                position = elements[1].strip()
                position = int(position)
                hap_std_dict[position] = elements
            except ValueError:
                print "error in ", file_name, position
    return hap_std_dict    

def load_hifi_result(file_name, hifi_dict):
    data_dict = load_raw_data(file_name, raw_data_format)[1]
    for position, elements in data_dict.iteritems():        
        try:
            if position not in hifi_dict:
                seed = seeds()
                seed.rsID = elements[0]
                seed.position = int(elements[1])
                hifi_dict[position] = seed            
            # need to modify this if seed is from second snp column, mother side
            hifi_dict[position].allele_new = elements[2].strip()
            # update the number of the allele from hifi result
            for base, value in hifi_dict[position].allele_dict.iteritems():
                if base == hifi_dict[position].allele_new:
                    hifi_dict[position].allele_dict[base] += 1
        except:
                #print "error at ", file_name, position, elements
                pass
    print "total snp number in : ", file_name, len(data_dict) 
    return hifi_dict

# group the seed into homo and hetero groups
def group_seed(seed_dict, geno_dict):
    seed_homo_dict = {}
    seed_hetero_dict = {}
    for position, snp in seed_dict.iteritems():
        if position in geno_dict:    # pos in seed may not be in geno
            geno_allele = geno_dict[position][2]
            if geno_allele[0] == geno_allele[1]:
                seed_homo_dict[position] = snp
            else:
                seed_hetero_dict[position] = snp
        else:
            seed_hetero_dict[position] = snp
    return (seed_homo_dict, seed_hetero_dict)

def output_revised_seed(filename, revised_seed_dict):
    seed_new_file = open(currentPath + filename, "w")
    print >> seed_new_file, seed_title_info
    revised_seed_sorted_list = sort_dict_by_key(revised_seed_dict)     # need to sort the snps by position
    for snp in revised_seed_sorted_list:
        seed = snp[1]
        line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
        print >> seed_new_file, line
    seed_new_file.close()
    
def output_revised_seed_dict(filename, revised_seed_dict):
    seed_new_file = open(currentPath + filename, "w")
    print >> seed_new_file, seed_title_info
    revised_seed_sorted_list = sort_dict_by_key(revised_seed_dict)     # need to sort the snps by position
    for snp in revised_seed_sorted_list:
        seed = snp[1]
        line = seed[0] + "\t" + str(seed[1]) + "\t" + seed[2]
        print >> seed_new_file, line
    seed_new_file.close()
    
def output_revised_seed_without_error(revised_seed_dict, same_to_B_dict):
    seed_new_file = open(currentPath + "haplotype_without_error.txt", "w")
    print >> seed_new_file, seed_title_info
    for position, seed in revised_seed_dict.iteritems():
        if position not in same_to_B_dict:
            line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
            print >> seed_new_file, line
    seed_new_file.close()

def output_revised_seed_with_error(revised_seed_dict, same_to_B_dict):
    seed_new_file = open(currentPath + "haplotype_with_error.txt", "w")
    print >> seed_new_file, seed_title_info
    for position, seed in revised_seed_dict.iteritems():
        if position in same_to_B_dict:
            line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_ori
            print >> seed_new_file, line
    seed_new_file.close()

if __name__=='__main__':

    data_dict = data_dicts()
    data_dict.chr_name = "chr9"
    data_dict.load_data_dicts()

