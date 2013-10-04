#!/usr/bin/python
#######################################################################################
# hifi v2
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *
from calculate_maf import calculate_maf
from cluster import get_cluster
from hifiAccuCheck_v2 import hifiAccuCheck

from seed_std_compare import seed_std_compare



class hifi:
	
	def __init__(self):
		self.chr_name = "chr9"
	
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
		self.h_xn_imputed_dict = {}
		self.maf_dict = {}
		
		self.pos_to_impute = []
		
		self.window_initial_size = 16
		self.window_size_upper_bound = 100
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
		print "ref_title_info: ", len(self.ref_title_info.split())
	
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
		#self.output_dict("h_xn.txt", self.h_xn_dict)
		#seed_std_compare("h_xn.txt", self.chr_name)
		
	
	def output_dict(self, filename, data_dict):
		seed_new_file = open(currentPath + filename, "w")
		print >> seed_new_file, self.seed_title_info
		revised_seed_sorted_list = sort_dict_by_key(data_dict)
		for snp in revised_seed_sorted_list:
			pos = snp[0]
			seed = snp[1]
			if seed[0] != 'X' and seed[1] != 'X' and seed[0] != 'N' and seed[0] != 'N': 
				print >> seed_new_file, self.ref_dict[pos][0], pos, seed[0], seed[1]
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
		#pos_to_impute.extend(self.geno_homo_dict.keys())
		pos_to_impute = list(set(pos_to_impute))
		pos_to_impute.sort()
		print len(pos_to_impute)
		self.pos_to_impute = pos_to_impute
	
	def to_impute_window(self, impute_type):

		window_half = self.window_initial_size/2
		pos_total = len(self.pos_to_impute)

		imputed_snp_this_round = 1
		#if True:
		while imputed_snp_this_round != 0:
			imputed_snp_this_round = 0
		
			window_center = 0
			window_start = 0
			window_end = 0
			total_impute_win_num = 0
			xx_number = 0
			
			
			for i, pos in enumerate(self.pos_to_impute):
				if self.h_xn_dict[pos][0] == impute_type or self.h_xn_dict[pos][1] == impute_type:
					xx_number += 1
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
					#print "impute_window_pos_list size: ", len(impute_window_pos_list)
					expand_direction = "up"
					while len(impute_window_pos_list) < self.window_size_lower_bound:
						if expand_direction == "up":
							if window_start -1 != 0:
								window_start = window_start - 1
								impute_window_pos_list.append(self.pos_to_impute[window_start])
								#print "up", self.pos_to_impute[window_start-1]
							expand_direction = "down"
						elif expand_direction == "down":
							if window_end +1 != pos_total:
								window_end += 1
								impute_window_pos_list.append(self.pos_to_impute[window_end])
							expand_direction = "up"
							#print "down", self.pos_to_impute[window_end+1]
					total_impute_win_num += 1
					
					match_to_A = 0
					match_to_B = 0
					solution_size = 0
					solution_list = []
					shrink_point = "top"
					expand_direction = "up"
					
					previous_round_solution_size = 0
					impute_cycle = 0
					impute_cycle_limit = 100
					
					impute_window_pos_list_size = len(impute_window_pos_list)	
					while impute_window_pos_list_size >= self.window_size_lower_bound and impute_window_pos_list_size <= self.window_size_upper_bound \
					 and solution_size != 1 and impute_cycle < impute_cycle_limit:
							impute_cycle += 1
						#while match_to_A != 1 and match_to_B != 1:
							solution_list = self.impute_X(impute_window_pos_list, window_center, impute_type)
							 
							solution_size = len(solution_list)
							
							# compare previous solution current solution
							if previous_round_solution_size == 2 and solution_size == 0:
								# remove the snp added in last round, skip it then add another one to check for solution
								# impute_window_pos_list is sorted, so need to check which one is added last time
								index_to_remove = 0 if expand_direction == "down" else -1
								#print "remove 2 to 0", index_to_remove, impute_window_pos_list[index_to_remove]
								del impute_window_pos_list[index_to_remove]
								solution_size = 2
									
							elif previous_round_solution_size == 0 and solution_size == 2:
								if expand_direction == "down":
									if window_start -1 != 0:
										window_start = window_start - 1
									else:
										expand_direction == "up"
										window_end += 1
								elif expand_direction == "up":
									if window_end +1 != pos_total:
										window_end += 1
									else:
										expand_direction == "down"
										window_start = window_start - 1
		
							#print "solution_size", solution_size
							if solution_size == 1:
								#print window_center, solution_list
								self.h_xn_dict[window_center] = (solution_list[0][0], solution_list[0][1])
								self.h_xn_imputed_dict[window_center] = (solution_list[0][0], solution_list[0][1], impute_window_pos_list)
								imputed_snp_this_round += 1
								#if impute_type == 'N':
								#	print window_center, solution_list
								#print "impute_window_pos_list size: ", impute_window_pos_list_size
							elif solution_size == 0:
								if shrink_point == "top":
									window_start = window_start + 1
									del impute_window_pos_list[0]
									shrink_point = "bottom"
								elif shrink_point == "bottom":
									window_end = window_end - 1
									del impute_window_pos_list[-1]
									shrink_point = "top"
							#elif solution_size > 2:
							#	solution_size == 2
							elif solution_size >= 2:
								if expand_direction == "up":
									if window_start -1 != 0:
										window_start = window_start - 1
										impute_window_pos_list.append(self.pos_to_impute[window_start])
										#print "up", self.pos_to_impute[window_start-1]
									expand_direction = "down"
								elif expand_direction == "down":
									if window_end +1 != pos_total:
										window_end += 1
										impute_window_pos_list.append(self.pos_to_impute[window_end])
									expand_direction = "up"
							else:
								print "more than 2 solutions at", window_center
		
							previous_round_solution_size = solution_size
							if window_center not in impute_window_pos_list:
								impute_window_pos_list.append(window_center)
							impute_window_pos_list_size = len(impute_window_pos_list)
								
				#sys.exit(1)
			print "imputed_snp_this_round", imputed_snp_this_round		
		#print "total_impute_win_num", total_impute_win_num
		#print "xx_number", xx_number

	def get_matched_window(self, impute_window_pos_list, window_center, hap_type):
		# find the window that matches A or B
		hap_index = 0 if hap_type == 'A' else 1
		total_haplotype_number = len(self.ref_title_info.split()) - 2
		
		match_to_A_index_list_new = range(2, total_haplotype_number)
		for pos in impute_window_pos_list:
			if pos != window_center:
				match_to_A_index_list_new = [index for index in match_to_A_index_list_new if self.ref_dict[pos][index] == self.h_xn_dict[pos][hap_index]]		
		
		match_to_A_list = []
		ref_x_list = []
		
		for index in match_to_A_index_list_new:
			temp_ref = ""
			for pos in impute_window_pos_list:
				temp_ref += self.ref_dict[pos][index]
			if temp_ref not in match_to_A_list:
				match_to_A_list.append(temp_ref)
					
		#print "match_to_A_index_list_new", match_to_A_index_list_new
		#print "match_to_A_index_list_new", len(match_to_A_list)
		return match_to_A_list
		
	def impute_X(self, impute_window_pos_list, window_center, impute_type):
		
		impute_window_pos_list.sort()
		window_center_index = impute_window_pos_list.index(window_center)
		
		match_to_A_list = self.get_matched_window(impute_window_pos_list, window_center, 'A')
		match_to_A = len(match_to_A_list)
		solution_list = []
		
		if match_to_A > 0:
			match_to_B_list = self.get_matched_window(impute_window_pos_list, window_center, 'B')
			match_to_B = len(match_to_B_list)
			for ref_A in match_to_A_list:
				center_A = ref_A[window_center_index]
				for ref_B in match_to_B_list:
					center_B = ref_B[window_center_index]
					if impute_type == 'N':
						solution_list.append((center_A, center_B))
					elif impute_type == 'X' and center_A != center_B:
						solution_list.append((center_A, center_B))
		
		#print "solution_list new method", solution_list	
		
		return solution_list	 
		
	def impute_X_0(self, impute_window_pos_list, window_center):
		impute_window_pos_list.sort()
		total_haplotype_number = len(self.ref_title_info.split()) - 2
		#print "total_haplotype_number", total_haplotype_number
		ref_list = []
		ref_x_list = []
		
		window_center_index = impute_window_pos_list.index(window_center)
		#print "window_center_index", window_center_index
		
		"""
		check A first, if > 0, then check B
		"""
		
		for i in range(2, total_haplotype_number):
			temp_line = ""
			#temp_line_x = ""
			for pos in impute_window_pos_list:
				temp_line += self.ref_dict[pos][i]
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
		#print hap_xn_A
		#print hap_xn_B
		
		
		match_to_A = 0
		match_to_B = 0
		#print hap_xn_A[window_center_index], hap_xn_B[window_center_index]
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
		solution_list = []
		
		if match_to_A > 0 and match_to_B > 0:
			match_to_A_index_list = [index for index, ref in enumerate(ref_x_list) if ref == hap_xn_A]
			match_to_B_index_list = [index for index, ref in enumerate(ref_x_list) if ref == hap_xn_B]
			for index_A in match_to_A_index_list:
				center_A = ref_list[index_A][window_center_index]
				for index_B in match_to_B_index_list:
					center_B = ref_list[index_B][window_center_index]
					if center_A != center_B:
						solution_list.append((center_A, center_B))
		#print "solution_list", solution_list
		"""
		for ref in ref_list:
			print ref,
			if ref == hap_xn_A:
				print ref, "match"
				#break
		"""
			
		#print match_to_A, match_to_B	
		return solution_list
	
	def run(self):
			
		self.load_seed_geno_ref()	
		#elapse_time = time.time() - start_time
		#print "load_seed_geno_ref time is: " + str(format(elapse_time, "0.3f")) + "s"

		self.merge_pos()
		#elapse_time = time.time() - elapse_time
		#print "merge_pos time is: " + str(format(elapse_time, "0.3f")) + "s"
		self.hg_merge()
		#elapse_time = time.time() - elapse_time
		#print "hg_merge time is: " + str(format(elapse_time, "0.3f")) + "s"
		self.get_maf()
		#elapse_time = time.time() - elapse_time
		#print "get_maf time is: " + str(format(elapse_time, "0.3f")) + "s"
		
		maf_num_list = self.maf_dict.keys()
		#print maf_num_list
		
		maf_num_list.sort()
		maf_num_list.reverse()
		print maf_num_list
		self.combine_msg(186)
		for maf_num in maf_num_list:
			
			# improve add new pos to the current list
			#self.combine_msg(maf_num)
			if maf_num >= 0 and maf_num != 186:
				print "current maf_num:", maf_num
				self.pos_to_impute.extend(self.maf_dict[maf_num])
				#print len(self.pos_to_impute)
				pre_imputed_size = len(self.h_xn_imputed_dict)
				#print "pre_imputed_size", pre_imputed_size
				self.to_impute_window('X')
				self.to_impute_window('N')
				print "newly_imputed_size", len(self.h_xn_imputed_dict) - pre_imputed_size
		"""
		self.combine_msg(180)
		#elapse_time = time.time() - elapse_time
		#print "combine_msg time is: " + str(format(elapse_time, "0.3f")) + "s"
		start_time = time.time()
		self.to_impute_window('X')
		self.to_impute_window('N')
		elapse_time = time.time() - start_time
		print "***********************to_impute_window time is: " + str(format(elapse_time, "0.3f")) + "s"
		"""
		self.output_dict("h_xn_new.txt", self.h_xn_imputed_dict)
		seed_std_compare("h_xn_new.txt", self.chr_name)
		
		self.output_dict("h_xn_1.txt", self.h_xn_dict)
		seed_std_compare("h_xn_1.txt", self.chr_name)
		hifiAccuCheck("h_xn_1.txt", self.chr_name)
		

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











