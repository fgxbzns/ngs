#!/usr/bin/python
#######################################################################################
# Common tools
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser
from tools import *
from cluster import get_cluster



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
        self.hap_ref_major_allele_dict = {}
        self.hap_ref_allele_frequence_dict = {}

        self.record_file_name = "seed_correction_record.txt"
        
        self.hap_std_dict = {}
        self.ref_cluster_dict = {}
        self.cluster_pos_dict = {}
          
        self.number_of_subfile = 10
        self.cycle_number = 3
        
        self.maf_upper_bound = 0.5
        self.maf_lower_bound = 0.3
        self.ref_cluster_dict = {}
        self.cluster_pos_dict = {}
        
        # parameters for ref extract
        #self.ref_position_distance = 500000
       # self.ref_expand_range = 5
        self.ref_position_distance = 500000
        self.ref_expand_range = 10
        self.remPercent = 0.6
        self.allele_new_percentage = 0.80
        self.qscore_threshold = 0.70
        
        # parameters for seed group
        self.seed_group_window_size = 20
        self.existing_seed_percentage = 0.95

    def update_seed_dict(self):
        #print "seed_file is :", self.seed_file
        self.seed_file_name = self.seed_file[:self.seed_file.find('.')].strip()
        self.seed_title_info, self.seed_dict = load_seed_data(self.seed_file)
        if len(self.geno_dict) == 0:
            self.update_geno_dict()
        self.seed_homo_dict, self.seed_hetero_dict = group_seed(self.seed_dict, self.geno_dict)
        print "total_seed_number: ", len(self.seed_dict)
        print "seed_homo_dict", len(self.seed_homo_dict)
        print "seed_hetero_dict", len(self.seed_hetero_dict)
        
    def update_geno_dict(self):
        genotype_file = file_path + "genotype_NA10847_" + self.chr_name + ".txt"
        self.geno_title_info, self.geno_dict = load_raw_data(genotype_file)
        self.geno_homo_dict, self.geno_hetero_dict = group_seed(self.geno_dict, self.geno_dict)
        print "total_geno_number: ", len(self.geno_dict)
        print "geno_homo_dict", len(self.geno_homo_dict)
        print "geno_hetero_dict", len(self.geno_hetero_dict)
        
    def update_ref_dict(self):
        self.ref_title_info, self.hap_ref_dict = load_raw_data(self.ref_file_name)
        print "total_ref_number: ", len(self.hap_ref_dict)

    def update_hap_std_dict(self):  
        hap_std_file = file_path + "ASW_"+self.chr_name+"_child_hap_refed.txt"    
        self.hap_std_dict = load_hap_std(hap_std_file)
        print "total_hap_std_dict_number: ", len(self.hap_std_dict)
    
    def update_ref_cluster_dict(self):
        self.ref_cluster_dict = get_cluster(self.ref_file_name, self.maf_upper_bound, self.maf_lower_bound)
        for maf_num, cluster_list in self.ref_cluster_dict.iteritems():
            for cluster_dict in cluster_list:
                for pos, ref in cluster_dict.iteritems():
                    self.cluster_pos_dict[pos] = ref
        print "cluster_pos_dict: ", len(self.cluster_pos_dict)
    
    def update_ref_major_allele(self):
        # prepare the list with major alleles
        maf_dict = {}
        for pos, snp in self.hap_ref_dict.iteritems():
            alleles = snp[2:]
            unique_alleles = set(alleles)
            if '-' in unique_alleles:
	            unique_alleles.remove('-')
            n_alleles = len(unique_alleles)
            major_allele = ""
            if n_alleles == 0 or n_alleles > 2:
                print "maf error in ref: ", pos, unique_alleles
                #sys.exit(1)
            else:
                maf_temp_list = []
                for ref_allele in unique_alleles:
                    maf_temp_list.append(ref_allele, alleles.count(ref_allele))
                major_allele = maf_temp_list[0][0] if maf_temp_list[0][1] >= maf_temp_list[1][1] else maf_temp_list[1][0]
            if pos not in maf_dict:
                maf_dict[pos] = major_allele
            else:
                print "duplicated pos", pos
                pass
            #maf_dict[maf_num].append(pos)

        print len(maf_dict)
        self.hap_ref_major_allele_dict = maf_dict
    
    def update_hap_ref_allele_frequence_dict(self):
        maf_dict = {}
        for pos, snp in self.hap_ref_dict.iteritems():
            alleles = snp[2:]
            number_of_alleles = len(alleles)
            unique_alleles = set(alleles)
            if '-' in unique_alleles:
	            unique_alleles.remove('-')
            n_alleles = len(unique_alleles)
            major_allele = ""
            if n_alleles == 0 or n_alleles > 2:
                print "maf error in ref: ", pos, unique_alleles
                sys.exit(1)
            else:      
                maf_temp_list = []
                for ref_allele in unique_alleles:
                    maf_temp_list.append((ref_allele, format((float(alleles.count(ref_allele))/number_of_alleles), "0.3f")))
                if pos not in maf_dict:
                    maf_dict[pos] = maf_temp_list
                else:
                    print "duplicated pos", pos
                    pass
        print len(maf_dict)
        self.hap_ref_allele_frequence_dict = maf_dict
    
    def output_hap_ref_allele_frequence_dict(self):
       file_name = "hap_ref_allele_frequence_dict.txt"
       ref_allele_frequence_sorted_list = sort_dict_by_key(self.hap_ref_allele_frequence_dict)
       with open(file_name, "w") as output_file:
           for data in ref_allele_frequence_sorted_list:
                pos = data[0]
                allele_frequence_list = data[1]
                try:
                    if allele_frequence_list[0][1] >= allele_frequence_list[1][1]:
                        print >> output_file, pos, allele_frequence_list[0][0], allele_frequence_list[0][1], allele_frequence_list[1][0], allele_frequence_list[1][1]
                    else:
                        print >> output_file, pos, allele_frequence_list[1][0], allele_frequence_list[1][1], allele_frequence_list[0][0], allele_frequence_list[0][1] 
                except:
                    print "output_hap_ref_allele_frequence_dict", pos, allele_frequence_list
    
    def read_hap_ref_allele_frequence_dict(self):
        maf_dict = {}
        file_name = "hap_ref_allele_frequence_dict.txt"
        with open(file_name, "r") as input_file:
            for line in input_file:
                elements = line.strip().split()
                pos = int(elements[0])
                temp_list = []
                temp_list.append((elements[1], float(elements[2])))
                temp_list.append((elements[3], float(elements[4])))
                maf_dict[pos] = temp_list
        self.hap_ref_allele_frequence_dict = maf_dict
                    
    def load_data_dicts(self):
        self.update_geno_dict()
        self.update_seed_dict()
        self.update_ref_dict()
        self.update_hap_std_dict()
        #self.update_ref_cluster_dict()
    
    def load_seed_geno_ref(self):
        self.update_geno_dict()
        self.update_seed_dict()
        self.update_ref_dict()
        
    def load_ref_allele_frequence(self):
        self.update_hap_ref_allele_frequence_dict()
        self.output_hap_ref_allele_frequence_dict()
        #self.read_hap_ref_allele_frequence_dict()
        #pass
      
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
    seed_title_info, temp_data_dict = load_raw_data(file_name, raw_data_format)
    for position, elements in temp_data_dict.iteritems():
        seed = seeds()
        seed.rsID = elements[0].strip()
        seed.position = int(elements[1].strip())
        seed.allele_ori = elements[2].strip()
        seed.allele_new = elements[2].strip()
        seed_dict[position] = seed
    return (seed_title_info, seed_dict)

def load_seed_data_from_dict(temp_data_dict):
    seed_dict = {}
    for position, elements in temp_data_dict.iteritems():
        seed = seeds()
        seed.rsID = elements[0].strip()
        seed.position = int(elements[1].strip())
        seed.allele_ori = elements[2].strip()
        seed.allele_new = elements[2].strip()
        seed_dict[position] = seed
    return seed_dict

def load_hap_std(file_name):
    hap_std_dict = {}
    temp_data_dict = load_raw_data(file_name, raw_data_format)[1]
    for position, elements in temp_data_dict.iteritems():
        """ ??? N X """
        if elements[2].strip() != "N" and elements[3].strip() != "N": 
            try:
                position = int(elements[1].strip())
                hap_std_dict[position] = elements
            except ValueError:
                print "error in ", file_name, position
    return hap_std_dict    

def load_hifi_result(file_name, hifi_dict):
    temp_data_dict = load_raw_data(file_name)[1]
    for position, elements in temp_data_dict.iteritems():        
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
    print "total snp number in : ", file_name, len(temp_data_dict) 
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

if __name__=='__main__':

    data_dict = data_dicts()
    data_dict.chr_name = "chr9"
    data_dict.load_data_dicts()
    data_dict.load_seed_geno_ref()














def compare_std_hap_ref():
    ref_file_name = "refHaplos.txt"
    hap_ref_dict = load_raw_data(ref_file_name, raw_data_format)[1]

    std_x_dict = {}
    std_not_in_ref_dict = {}
    ref_homo = 0

    for position, snp in hap_std_dict.iteritems():
        if position in hap_ref_dict:
            exist_in_ref = False
            std_A = snp[2]
            std_B = snp[3]
            alleles = hap_ref_dict[position]
            alleles = alleles[2:]
            unique_alleles = list(set(alleles))
            n_alleles = len(unique_alleles)
            if n_alleles == 0 or n_alleles > 2:
                print "error in: ", position
                sys.exit(1)
            else:
                try:
                    if n_alleles == 1:
                        ref_homo += 1
                    if std_A == std_B:
                        if n_alleles == 2:
                            if std_A == unique_alleles[0] or std_A == unique_alleles[1]:
                                exist_in_ref = True
                        if n_alleles == 1:
                            if std_A == unique_alleles[0]:
                                exist_in_ref = True
                    else:
                        if n_alleles == 2:
                            if (std_A == unique_alleles[0] and std_B == unique_alleles[1]) or (std_A == unique_alleles[1] and std_B == unique_alleles[0]):
                                exist_in_ref = True
                            if n_alleles == 1:    # hetero_std, homo_ref
                                pass
                    if not exist_in_ref:
                        if std_A == 'X':
                            std_x_dict[position] = list_to_line(snp)
                        else:
                            std_not_in_ref_dict[position] = list_to_line(snp)
                            print position, unique_alleles, n_alleles, std_A, std_B, geno_dict[position][2]
                except:
                    #print position, unique_alleles, n_alleles
                    pass
    print "std_x_dict: ", len(std_x_dict)
    print "std_not_in_ref_dict: ", len(std_not_in_ref_dict)
    print "ref_homo", ref_homo

def compare_geno_ref():
    ref_file_name = "refHaplos.txt"
    hap_ref_dict = load_raw_data(ref_file_name, raw_data_format)[1]

    std_x_dict = {}
    std_n_dict = {}
    std_not_in_ref_dict = {}
    ref_homo = 0

    for position, snp in geno_dict.iteritems():
        if position in hap_ref_dict:
            exist_in_ref = False
            std_A = snp[2][0]
            std_B = snp[2][1]
            alleles = hap_ref_dict[position]
            alleles = alleles[2:]
            unique_alleles = list(set(alleles))
            n_alleles = len(unique_alleles)
            if n_alleles == 0 or n_alleles > 2:
                print "error in: ", position
                sys.exit(1)
            else:
                try:
                    if n_alleles == 1:
                        ref_homo += 1
                    if std_A == std_B:
                        if n_alleles == 2:
                            if std_A == unique_alleles[0] or std_A == unique_alleles[1]:
                                exist_in_ref = True
                        if n_alleles == 1:
                            if std_A == unique_alleles[0]:
                                exist_in_ref = True
                    else:
                        if n_alleles == 2:
                            if (std_A == unique_alleles[0] and std_B == unique_alleles[1]) or (std_A == unique_alleles[1] and std_B == unique_alleles[0]):
                                exist_in_ref = True
                            if n_alleles == 1:    # hetero_std, homo_ref
                                pass
                    if not exist_in_ref:
                        if std_A == 'X':
                            std_x_dict[position] = list_to_line(snp)
                        if std_A == 'N':
                            std_n_dict[position] = list_to_line(snp)
                        else:
                            std_not_in_ref_dict[position] = list_to_line(snp)
                            print position, unique_alleles, n_alleles, std_A, std_B, geno_dict[position][2]
                except:
                    #print position, unique_alleles, n_alleles
                    pass
    print "std_x_dict: ", len(std_x_dict)
    print "std_n_dict: ", len(std_n_dict)
    print "std_not_in_ref_dict: ", len(std_not_in_ref_dict)
    print "ref_homo", ref_homo

def output_revised_seed(filename, selected_seed_dict):
    seed_new_file = open(currentPath + filename, "w")
    print >> seed_new_file, "rsID    position    NA10847-F       NA10847-M"
    selected_seed_sorted_list = sort_dict_by_key(selected_seed_dict)     # need to sort the snps by position
    for snp in selected_seed_sorted_list:
        seed = snp[1]
        line = seed.rsID + "\t" + str(seed.position) + "\t" + seed.allele_new
        print >> seed_new_file, line
    seed_new_file.close()

def generate_std_seed(seed_number):
        chr_name = "chr1"
        hap_std_dict = load_seed_data(file_path+"ASW_"+chr_name+"_child_hap_refed.txt")[1]
        hap_std_list = sort_dict_by_key(hap_std_dict)
        
        genotype_file = file_path + "genotype_NA10847_" + chr_name + ".txt"
        geno_dict = load_raw_data(genotype_file)[1]
        seed_homo_dict, seed_hetero_dict = group_seed(hap_std_dict, geno_dict)

        seed_hetero_list = sort_dict_by_key(seed_hetero_dict)

        selected_seed_dict = {}
        i = 0
        while i < seed_number-1:
            random_index = random.randrange(0,(len(seed_hetero_list)-1))
            while seed_hetero_list[random_index][0] in selected_seed_dict:
                random_index = random.randrange(0,(len(seed_hetero_list)-1))
            selected_seed_dict[seed_hetero_list[random_index][0]] = seed_hetero_list[random_index][1]
            i += 1
        # always add the last snp into seed, hifi requirement
        selected_seed_dict[hap_std_list[-1][0]] = hap_std_list[-1][1]

        file_name = "haplotype.txt"
        output_revised_seed(file_name, selected_seed_dict)
        """
        seed_std_compare(file_name, chr_name)
        refMerger(file_name, chr_name, 0)
        hifi_run(file_name, chr_name)
        hifiAccuCheck("imputed_"+file_name, chr_name)
        """

