#!/usr/bin/python
#######################################################################################
# hifi v2
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *
from calculate_maf import calculate_maf
from cluster import get_cluster


class hifi:
	
	def __init__(self):
	
		self.seed_file_name = "haplotype.txt"
		self.seed_file_prefix = ""
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
		self.ref_dict = {}
		
		self.hap_std_dict = {}
		self.ref_cluster_dict = {}
		self.cluster_pos_dict = {}
		
		self.hgr_pos = []
		self.h_xn_dict = {}
		self.maf_dict = {}
		
		self.pos_to_impute = []
		
		self.window_initial_size = 16
		self.window_size_upper_bound = 500
		self.window_size_lower_bound = 10
		  
	def update_seed_dict(self):
		# print "seed_file is :", self.seed_file
		self.seed_file_prefix = self.seed_file_name[:self.seed_file_name.find('.')].strip()
		self.seed_title_info, self.seed_dict = load_raw_data(self.seed_file_name)
		if len(self.geno_dict) == 0:
			self.update_geno_dict()
		self.seed_homo_dict, self.seed_hetero_dict = group_seed(self.seed_dict, self.geno_dict)
		print "total_seed_number: ", len(self.seed_dict)
		print "seed_homo_dict", len(self.seed_homo_dict)
		print "seed_hetero_dict", len(self.seed_hetero_dict)
		
	def update_geno_dict(self):
		# genotype_file = file_path + "genotype_NA10847_" + self.chr_name + ".txt"
		self.geno_title_info, self.geno_dict = load_raw_data(self.geno_file_name)
		self.geno_homo_dict, self.geno_hetero_dict = group_seed(self.geno_dict, self.geno_dict)
		print "total_geno_number: ", len(self.geno_dict)
		print "geno_homo_dict", len(self.geno_homo_dict)
		print "geno_hetero_dict", len(self.geno_hetero_dict)
		
	def update_ref_dict(self):
		self.ref_title_info, self.ref_dict = load_raw_data(self.ref_file_name)
		print "total_ref_number: ", len(self.ref_dict)
	
	def update_hap_std_dict(self):  
		hap_std_file = file_path + "ASW_" + self.chr_name + "_child_hap_refed.txt"	
		self.hap_std_dict = load_hap_std(hap_std_file)
		print "total_hap_std_dict_number: ", len(self.hap_std_dict)
	"""
	def update_ref_cluster_dict(self):
		self.ref_cluster_dict = get_cluster(self.ref_file_name, self.maf_upper_bound, self.maf_lower_bound)
		for maf_num, cluster_list in self.ref_cluster_dict.iteritems():
			for cluster_dict in cluster_list:
				for pos, ref in cluster_dict.iteritems():
					self.cluster_pos_dict[pos] = ref
		print "cluster_pos_dict: ", len(self.cluster_pos_dict)
	"""
	def load_data_dicts(self):
		self.update_geno_dict()
		self.update_seed_dict()
		self.update_ref_dict()
		self.update_hap_std_dict()
		# self.update_ref_cluster_dict()
	
	def load_seed_geno_ref(self):
		self.update_geno_dict()
		self.update_seed_dict()
		self.update_ref_dict()
	
	def merge_pos(self):
		"""
		This step is already done in refMerger, all positions in hap
		and geno are in ref.		
		
		put all in one list, choose unique, sort
		
		"""
		self.hgr_pos.extend(self.seed_dict.keys())
		self.hgr_pos.extend(self.geno_dict.keys())
		self.hgr_pos.extend(self.ref_dict.keys())
		#print len(self.hgr_pos)
		#print len(self.ref_dict)
		#print len(list(set(self.hgr_pos)))
		self.hgr_pos = list(set(self.hgr_pos))
		self.hgr_pos.sort()
		#print len(self.hgr_pos)
	

	def hg_merge(self):
		
		"""xx must be hetero"""
		
		h_xn_dict = {}
		hap_geno_discrepancy_dict = {}
		
		for pos in self.hgr_pos:
			if pos not in self.seed_dict and pos not in self.geno_dict:
				h_xn_dict[pos] = ('N', 'N')
			elif pos not in self.seed_dict and pos in self.geno_dict:
				if self.geno_dict[pos][2][0] == self.geno_dict[pos][2][1]:
					h_xn_dict[pos] = (self.geno_dict[pos][2][0], self.geno_dict[pos][2][0])
				else:
					h_xn_dict[pos] = ('X', 'X')
			elif pos in self.seed_dict and pos not in self.geno_dict:
				h_xn_dict[pos] = (self.seed_dict[pos][2], 'N')
			elif pos in self.seed_dict and pos in self.geno_dict:
				if self.geno_dict[pos][2][0] == self.geno_dict[pos][2][1]:
					if self.seed_dict[pos][2] == self.geno_dict[pos][2][0]:
						h_xn_dict[pos] = (self.geno_dict[pos][2][0], self.geno_dict[pos][2][0])
					else:
						h_xn_dict[pos] = (self.geno_dict[pos][2][0], self.geno_dict[pos][2][0])
						hap_geno_discrepancy_dict[pos] = (self.seed_dict[pos][2], self.geno_dict[pos][2])
				else:
					if self.seed_dict[pos][2] == self.geno_dict[pos][2][0]:
						h_xn_dict[pos] = (self.geno_dict[pos][2][0], self.geno_dict[pos][2][1])
					elif self.seed_dict[pos][2] == self.geno_dict[pos][2][1]:
						h_xn_dict[pos] = (self.geno_dict[pos][2][1], self.geno_dict[pos][2][0])
					else:
						h_xn_dict[pos] = ('X', 'X')
						hap_geno_discrepancy_dict[pos] = (self.seed_dict[pos][2], self.geno_dict[pos][2])
			else:
				print "geno AN or NA", pos
		
		#print len(hap_geno_discrepancy_dict)				
		self.h_xn_dict = h_xn_dict
		self.output_dict("h_xn.txt", self.h_xn_dict)
	
	def output_dict(self, filename, data_dict):
		seed_new_file = open(currentPath + filename, "w")
		print >> seed_new_file, self.seed_title_info
		revised_seed_sorted_list = sort_dict_by_key(data_dict)
		for snp in revised_seed_sorted_list:
			pos = snp[0]
			seed = snp[1]
			print >> seed_new_file, pos, seed[0], seed[1]
		seed_new_file.close()			
	
	def get_maf(self):
		# calculate maf for each snp in ref
		# maf is represented by number here
		maf_dict = {}
		for pos, snp in self.ref_dict.iteritems():
			alleles = snp[2:]
			unique_alleles = set(alleles)
			n_alleles = len(unique_alleles)
			if n_alleles == 0 or n_alleles > 2:
				print "maf error in ref: ", position
				sys.exit(1)
			else:
				maf_temp_list = []
				for ref_allele in unique_alleles:
					maf_temp_list.append(alleles.count(ref_allele))
				maf_num = min(maf_temp_list)	
			if maf_num not in maf_dict:
				maf_dict[maf_num] = []
			maf_dict[maf_num].append(pos)
		"""
		for pos, list in maf_dict.iteritems():
			print pos
			for pos_1 in list:
				print pos_1,
			print ""
		"""
		print len(maf_dict)
		self.maf_dict = maf_dict
	
	def combine_msg(self, maf_step):
		#combine_maf_layer > maf_step, know_seed and geno_homo
		
		#maf_num_list = [pos for pos in self.maf_dict.keys() if pos >= maf_step]
		#print len(maf_num_list)

		pos_to_impute = []
		for pos, pos_list in self.maf_dict.iteritems():
			if pos >= maf_step:
				pos_to_impute.extend(pos_list)
		"""
		print len(pos_to_impute)
		pos_to_impute.extend(self.seed_dict.keys())
		print len(pos_to_impute)
		pos_to_impute.extend(self.geno_homo_dict.keys())
		print len(pos_to_impute)
		"""
		pos_to_impute.extend(self.seed_dict.keys())
		pos_to_impute.extend(self.geno_homo_dict.keys())
		pos_to_impute = list(set(pos_to_impute))
		pos_to_impute.sort()
		print len(pos_to_impute)
		self.pos_to_impute = pos_to_impute
	
	def to_impute_window(self, impute_type):
		window_half = self.window_initial_size/2
		pos_total = len(self.pos_to_impute)
		window_center = 0
		window_start = 0
		window_end = 0
		total_impute_win_num = 0
		for i, pos in enumerate(self.pos_to_impute):
			if self.h_xn_dict[pos][0] == impute_type or self.h_xn_dict[pos][1] == impute_type:
				window_center = self.pos_to_impute[i]
				window_start = i - window_half if i - window_half >= 0 else 0
				window_end = i + window_half if i + window_half < pos_total else (pos_total - 1)
				if window_start == window_end:
					print "window size is zero"
					sys.exit(1)
				
				impute_window_pos_list = []
				for i in range (window_start, window_end+1):
					position = self.pos_to_impute[i]
					h_xn = self.h_xn_dict[position]
					if h_xn[0] != 'N' and h_xn[1] != 'N' and h_xn[0] != 'X' and h_xn[1] != 'X':
						impute_window_pos_list.append(position)
				impute_window_pos_list.append(window_center)
				"""
				prepare the impute window, if the size is smaller than lower bound,
				add more pos from either top or bottom
				"""
				print "impute_window_pos_list size: ", len(impute_window_pos_list)
				expand_direction = "up"
				while len(impute_window_pos_list) < self.window_size_lower_bound:
					if expand_direction == "up":
						if window_start -1 != 0:
							impute_window_pos_list.append(self.pos_to_impute[window_start-1])
							#print "up", self.pos_to_impute[window_start-1]
						expand_direction = "down"
					elif expand_direction == "down":
						if window_end +1 != pos_total:
							impute_window_pos_list.append(self.pos_to_impute[window_end+1])
						expand_direction = "up"
						#print "down", self.pos_to_impute[window_end+1]
				total_impute_win_num += 1
				
				match_to_A = 0
				match_to_B = 0
				shrink_point = "top"
				impute_window_pos_list_size = len(impute_window_pos_list)	
				while len(impute_window_pos_list) >= self.window_size_lower_bound and len(impute_window_pos_list) <= self.window_size_upper_bound \
				 and match_to_A != 1 and match_to_B != 1:
					#while match_to_A != 1 and match_to_B != 1:
						match_to_A, match_to_B = self.impute_X(impute_window_pos_list, window_center)
						
						if shrink_point == "top":
							del impute_window_pos_list[0]
							shrink_point = "bottom"
						elif shrink_point == "bottom":
							del impute_window_pos_list[-1]
							shrink_point = "top"
				
				
				#sys.exit(1)
				
		print "total_impute_win_num", total_impute_win_num
		 
			
	def impute_X(self, impute_window_pos_list, window_center):
		impute_window_pos_list.sort()
		total_haplotype_number = len(self.ref_title_info.split()) - 2
		print "total_haplotype_number", total_haplotype_number
		ref_list = []
		ref_x_list = []
		
		window_center_index = impute_window_pos_list.index(window_center)
		print "window_center_index", window_center_index
		
		for i in range(2, total_haplotype_number):
			temp_line = ""
			#temp_line_x = ""
			for pos in impute_window_pos_list:
				temp_line += self.ref_dict[pos][i]
				"""
				if pos == window_center:
					temp_line_x += 'X'
				else:
					temp_line_x += self.ref_dict[pos][i]
				"""
			if temp_line not in ref_list:
				ref_list.append(temp_line)
				#ref_x_list.append(temp_line_x)
		
		for ref in ref_list:
			ref_x_list.append(ref[:window_center_index] + 'X' + ref[window_center_index+1:])
		
		#print ref_list
		#print ref_x_list
		
		hap_xn_A = ""
		hap_xn_B = ""		
		for pos in impute_window_pos_list:
			hap_xn_A += self.h_xn_dict[pos][0]
			hap_xn_B += self.h_xn_dict[pos][1]
		print hap_xn_A
		print hap_xn_B
		
		
		match_to_A = 0
		match_to_B = 0
		print hap_xn_A[window_center_index], hap_xn_B[window_center_index]
		"""
		if hap_xn_A[window_center_index] != 'X' and hap_xn_B[window_center_index] == 'X':
			match_to_A = ref_list.count(hap_xn_A)
			match_to_B = ref_x_list.count(hap_xn_B)
		elif hap_xn_A[window_center_index] == 'X' and hap_xn_B[window_center_index] == 'X':
			match_to_A = ref_x_list.count(hap_xn_A)
			match_to_B = ref_x_list.count(hap_xn_B)
		"""
		match_to_A = ref_x_list.count(hap_xn_A)
		match_to_B = ref_x_list.count(hap_xn_B)
		
		"""
		for ref in ref_list:
			print ref,
			if ref == hap_xn_A:
				print ref, "match"
				#break
		"""
			
		print match_to_A, match_to_B	
		return match_to_A, match_to_B
	
	def run(self):
		self.load_seed_geno_ref()	
		self.merge_pos()
		self.hg_merge()
		self.get_maf()
		self.combine_msg(186)
		self.to_impute_window('X')

def get_args():
	desc = "Compare seed and std hap, to check purity of seed"

	usage = "seed_std_compare -i seed_file -c chr#" 
	parser = OptionParser(usage=usage, description=desc) 
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
	parser.add_option("-i", "--seed", type="string", dest="hifiSeed", help="Input seed file Name", default="null")
	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.hifiSeed == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	"""
	options = get_args()
	chr_name = options.chrName
	seed_input_file = options.hifiSeed	
	seed_std_compare(seed_input_file, chr_name)
	"""
	hifi = hifi()
	hifi.run()











