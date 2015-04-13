#!/usr/bin/python

# ######################################################################################
# Guoxing Fu Jan 28, 2015
# Laser project, to impute haplotype from trio data
# 1. impute child haplotype from parent haplotype
# 2. impute parent haplotype from child haplotype
# 3, 4. impute all haplotype from populaton data
#######################################################################################

import os
import glob
import subprocess
import random
import operator
import time
import sys
import copy
import shutil
import fnmatch
from optparse import OptionParser

class parameters:
	def __init__(self):
		self.person_dict = {}
		self.rsID_dict = {}
		self.pos_list = {}
		self.father_list = []
		self.mather_list = []
		self.parent_id_list = []
		self.children_list = []
		self.children_dict = {}

		self.fragment_dict = {}
		self.fragment_list = []

		self.children_hap_file = "child_hap_laser_1.txt"
		self.parent_hap_file = "parent_hap_laser_2.txt"

		self.ref_file_name = ""
		self.ori_ref_title = ""
		self.ori_ref_dict = {}
		self.hifi_file = "/hifi_fu_laser.laser"


class persons:
	def __init__(self):
		self.ID = ""
		self.father = ""
		self.mather = ""
		self.children = {}
		self.genotype_dict = {}
		self.haplotype = {}
		self.hetero_pos_list = []
		self.std_hap_dict = {}

		self.fragment_dict = {}
		self.common_fragment_dict = {}
		self.common_fragment_pos_dict = {}


class fragments:
	def __init__(self):
		self.ID = ""
		self.start = ""
		self.end = ""
		self.length = ""


def list_to_line(list):
	line = ""
	for a in list:
		line += str(a).strip() + "\t"
	return line.strip()


def sort_dict_by_key(input_dict):
	sorted_list = []
	sorted_list = [x for x in input_dict.iteritems()]
	sorted_list.sort(key=lambda x: x[0])
	return sorted_list


def load_pedi(pedi_name):
	with open(pedi_name, "r") as pedi_file:
		for line in pedi_file:
			if line.startswith("PersonID"):
				title_info = list_to_line(line.strip().split())
			else:
				try:
					elements = line.strip().split()
					person = persons()
					person.ID = elements[0]
					person.father = elements[1]
					person.mather = elements[2]
					if person.father != "N/A" and person.father in parameter.person_dict and person.ID not in \
							parameter.person_dict[person.father].children:
						parameter.person_dict[person.father].children[person.ID] = ""
					if person.mather != "N/A" and person.mather in parameter.person_dict and person.ID not in \
							parameter.person_dict[person.mather].children:
						parameter.person_dict[person.mather].children[person.ID] = ""
					if person.ID not in parameter.person_dict:
						parameter.person_dict[person.ID] = person
				except:
					print "error in ", line, pedi_name


def load_geno(geno_name):
	with open(geno_name, "r") as geno_file:
		for line in geno_file:
			if line.startswith("rs#"):
				ID_list = line.strip().split()[2:]
			else:
				try:
					elements = line.strip().split()
					rsID = elements[0]
					position = int(elements[1])
					parameter.rsID_dict[position] = rsID
					genotype = elements[2:]
					for index, ID in enumerate(ID_list):
						parameter.person_dict[ID].genotype_dict[position] = genotype[index]
						parameter.person_dict[ID].std_hap_dict[position] = (genotype[index][0], genotype[index][1])
				except:
					print "error in geno", line, pedi_name


def parents_to_children():
	for ID in parameter.person_dict:
		person = parameter.person_dict[ID]
		if person.father != "N/A" and person.mather != "N/A":

			pos_list = person.genotype_dict.keys()
			pos_list.sort()
			for pos in pos_list:
				c_geno = person.genotype_dict[pos]
				f_geno = parameter.person_dict[person.father].genotype_dict[pos]
				m_geno = parameter.person_dict[person.mather].genotype_dict[pos]

				if True:
					f_set = set(f_geno)
					m_set = set(m_geno)
					c_set = set(c_geno)

					if len(c_set) == 1 and c_geno != "NN":
						if (c_geno[0] == f_geno[0] or c_geno[0] == f_geno[1]) and (
										c_geno[0] == m_geno[0] or c_geno[0] == m_geno[1]):
							parameter.person_dict[ID].haplotype[pos] = (c_geno[0], c_geno[1])
					else:
						if f_geno != "NN":
							if m_geno != "NN":
								if c_geno != "NN":
									if len(f_set) == 1:
										if len(m_set) == 1:
											parameter.person_dict[ID].haplotype[pos] = (f_geno[0], m_geno[1])
										elif len(m_set) != 1:
											cf_hap = f_geno[0]
											cm_hap = c_geno[1] if c_geno[0] == f_geno[0] else c_geno[0]
											parameter.person_dict[ID].haplotype[pos] = (cf_hap, cm_hap)
									elif len(f_set) != 1:
										if len(m_set) == 1:
											cf_hap = c_geno[1] if c_geno[0] == m_geno[0] else c_geno[0]
											cm_hap = m_geno[0]
											parameter.person_dict[ID].haplotype[pos] = (cf_hap, cm_hap)
										elif len(m_set) != 1:
											# all hetero, cannot determine
											parameter.person_dict[ID].haplotype[pos] = ("X", "X")

								elif c_geno == "NN":
									if len(f_set) == 1:
										if len(m_set) == 1:
											parameter.person_dict[ID].haplotype[pos] = (f_geno[0], m_geno[1])
										elif len(m_set) != 1:
											parameter.person_dict[ID].haplotype[pos] = (f_geno[0], "N")
									elif len(f_set) != 1:
										if len(m_set) == 1:
											parameter.person_dict[ID].haplotype[pos] = ("N", m_geno[0])
										elif len(m_set) != 1:
											parameter.person_dict[ID].haplotype[pos] = ("N", "N")
								else:
									print "error 1", pos

							elif m_geno == "NN":
								if c_geno != "NN":
									if len(f_set) == 1:
										cf_hap = f_geno[0]
										cm_hap = c_geno[1] if c_geno[0] == f_geno[0] else c_geno[0]
										parameter.person_dict[ID].haplotype[pos] = (cf_hap, cm_hap)
									elif len(f_set) != 1:
										parameter.person_dict[ID].haplotype[pos] = ("N", "N")
								elif c_geno == "NN":
									if len(f_set) == 1:
										parameter.person_dict[ID].haplotype[pos] = (f_geno[0], "N")
									elif len(f_set) != 1:
										parameter.person_dict[ID].haplotype[pos] = ("N", "N")
							else:
								print "error 2", pos

						elif f_geno == "NN":
							if m_geno != "NN":
								if c_geno != "NN":
									if len(m_set) == 1:
										cf_hap = c_geno[1] if c_geno[0] == m_geno[0] else c_geno[0]
										cm_hap = m_geno[0]
										parameter.person_dict[ID].haplotype[pos] = (cf_hap, cm_hap)
									elif len(m_set) != 1:
										parameter.person_dict[ID].haplotype[pos] = ("N", "N")
								if c_geno == "NN":
									if len(m_set) == 1:
										parameter.person_dict[ID].haplotype[pos] = ("N", m_geno[0])
									elif len(m_set) != 1:
										parameter.person_dict[ID].haplotype[pos] = ("N", "N")
							if m_geno == "NN":
								parameter.person_dict[ID].haplotype[pos] = ("N", "N")
						else:
							print "3", pos


def prepare_id_list():
	for ID in parameter.person_dict.keys():
		person = parameter.person_dict[ID]
		if len(person.children) > 0:
			parameter.children_list.extend(list(person.children.keys()))

	parameter.children_list = list(set(parameter.children_list))
	parameter.children_list.sort()

	for ID in parameter.person_dict.keys():
		person = parameter.person_dict[ID]
		if person.father != "N/A":
			parameter.father_list.append(person.father)
		if person.mather != "N/A":
			parameter.mather_list.append(person.mather)

	parameter.father_list = list(set(parameter.father_list))
	parameter.father_list.sort()

	parameter.mather_list = list(set(parameter.mather_list))
	parameter.mather_list.sort()

	parameter.parent_id_list.extend(parameter.father_list)
	parameter.parent_id_list.extend(parameter.mather_list)
	parameter.parent_id_list.sort()

	parameter.pos_list = parameter.rsID_dict.keys()
	parameter.pos_list.sort()


def find_common_fragment(fragment_list_1, fragment_list_2):
	common_fragment_list = []

	for f1 in fragment_list_1:
		for f2 in fragment_list_2:
			if f1.start == f2.start and f1.end == f2.end:
				not_in_common = True
				for cf in common_fragment_list:
					if f1.start == cf.start and f1.end == cf.end:
						not_in_common = False
				if not_in_common:
					common_fragment_list.append(f1)
	return common_fragment_list


def sort_fragment(p_id):
	children_list = parameter.person_dict[p_id].children.keys()
	children_list.sort()

	fragment_startpos_dict = {}

	for child_id in children_list:
		list = parameter.person_dict[p_id].fragment_dict[child_id]
		for fragment in list:
			start_pos = fragment.start
			end_pos = fragment.end
			if start_pos not in fragment_startpos_dict:
				fragment_startpos_dict[int(start_pos)] = []
			fragment_startpos_dict[int(start_pos)].append(fragment)

	fragment_startpos_sorted_list = sort_dict_by_key(fragment_startpos_dict)
	return fragment_startpos_sorted_list


def children_to_parents(p_id, p_code):
	children_list = parameter.person_dict[p_id].children.keys()
	children_list.sort()

	for pos in parameter.pos_list:
		f_geno = parameter.person_dict[p_id].genotype_dict[pos]
		if f_geno != "NN" and f_geno[0] != f_geno[1]:
			parameter.person_dict[p_id].hetero_pos_list.append(pos)

	parameter.person_dict[p_id].hetero_pos_list.sort()
	#print "hetero pos_list size", len(parameter.person_dict[p_id].hetero_pos_list)

	compare_child_hap(p_id, p_code, children_list[0], children_list[1])
	compare_child_hap(p_id, p_code, children_list[0], children_list[2])
	compare_child_hap(p_id, p_code, children_list[1], children_list[2])

	child_1_2 = find_common_fragment(parameter.person_dict[p_id].fragment_dict[children_list[0]],
	                                 parameter.person_dict[p_id].fragment_dict[children_list[1]])
	child_1_2_3 = find_common_fragment(child_1_2, parameter.person_dict[p_id].fragment_dict[children_list[2]])


	fragment_sort_dict = {}
	for fragment in child_1_2_3:
		#print fragment.ID, fragment.start, fragment.end
		start_pos = fragment.start
		end_pos = fragment.end
		if start_pos not in fragment_sort_dict:
			fragment_sort_dict[int(start_pos)] = []
		fragment_sort_dict[int(start_pos)].append(fragment)

		for pos in parameter.person_dict[p_id].hetero_pos_list:
			if start_pos <= pos <= end_pos:
				parameter.person_dict[p_id].common_fragment_pos_dict[pos] = 0

	#print "fragment hetero size", p_id, len(parameter.person_dict[p_id].common_fragment_pos_dict)
	fragment_sorted_list = sort_dict_by_key(fragment_sort_dict)

	#output_fragment(p_id, fragment_sorted_list)
	#print len(parameter.person_dict[p_id].hetero_pos_list)
	for f1 in child_1_2_3:
		for pos in parameter.person_dict[p_id].hetero_pos_list:
			if int(f1.start) <= pos <= int(f1.end):
				parameter.person_dict[p_id].common_fragment_dict[pos] = 0
	#print len(parameter.person_dict[p_id].common_fragment_dict)

	temp_parent_hap_A = {}
	temp_parent_hap_B = {}

	fragment_startpos_sorted_list = sort_fragment(p_id)

	status = ""
	for data in fragment_startpos_sorted_list:
		start_pos = data[0]
		fragment_list = data[1]

		#print start_pos
		for fragment in fragment_list:
			#print "******", fragment.ID, fragment.start, fragment.end
			status, temp_parent_hap_A = update_parent_hap(p_id, fragment, temp_parent_hap_A, p_code)
			if status == "same_to_A":
				pass
			elif status == "same_to_B":
				status, temp_parent_hap_B = update_parent_hap(p_id, fragment, temp_parent_hap_B, p_code)
			elif status == "mix":
				pass

	count = 0
	for pos in parameter.pos_list:
		f_geno = parameter.person_dict[p_id].genotype_dict[pos]
		if f_geno[0] == f_geno[1]:
			temp_parent_hap_A[pos] = f_geno[0]
			temp_parent_hap_B[pos] = f_geno[0]
		else:
			if pos in temp_parent_hap_A and pos in temp_parent_hap_B:
				if temp_parent_hap_A[pos] == "X" and temp_parent_hap_B[pos] != "X":
					temp_parent_hap_A[pos] = f_geno[0] if temp_parent_hap_B[pos] == f_geno[1] else f_geno[1]
				if temp_parent_hap_B[pos] == "X" and temp_parent_hap_A[pos] != "X":
					temp_parent_hap_B[pos] = f_geno[0] if temp_parent_hap_A[pos] == f_geno[1] else f_geno[1]

			elif pos in temp_parent_hap_A and pos not in temp_parent_hap_B:
				if temp_parent_hap_A[pos] == "X":
					temp_parent_hap_B[pos] = "X"
				else:
					temp_parent_hap_B[pos] = f_geno[0] if temp_parent_hap_A[pos] == f_geno[1] else f_geno[1]

			elif pos not in temp_parent_hap_A and pos in temp_parent_hap_B:
				if temp_parent_hap_B[pos] == "X":
					temp_parent_hap_A[pos] = "X"
				else:
					temp_parent_hap_A[pos] = f_geno[0] if temp_parent_hap_B[pos] == f_geno[1] else f_geno[1]

			elif pos not in temp_parent_hap_A and pos not in temp_parent_hap_B:
				temp_parent_hap_A[pos] = "N"
				temp_parent_hap_B[pos] = "N"
				count += 1

	for pos in parameter.pos_list:
		if pos not in parameter.person_dict[p_id].haplotype:
			parameter.person_dict[p_id].haplotype[pos] = ["", ""]

		if pos in temp_parent_hap_A:
			parameter.person_dict[p_id].haplotype[pos][0] = temp_parent_hap_A[pos]
		else:
			parameter.person_dict[p_id].haplotype[pos][0] = "N"
		if pos in temp_parent_hap_B:
			parameter.person_dict[p_id].haplotype[pos][1] = temp_parent_hap_B[pos]
		else:
			parameter.person_dict[p_id].haplotype[pos][1] = "N"
	"""
	with open(p_id + "_hap.txt", "w") as x_file:
		print >> x_file, "rs#", "pos", parameter.person_dict[p_id].ID + "_A", parameter.person_dict[p_id].ID + "_B"
		for pos in parameter.pos_list:
			if pos in temp_parent_hap_A and pos in temp_parent_hap_B \
					and temp_parent_hap_A[pos] != "N" and temp_parent_hap_B[pos] != "N" \
					and temp_parent_hap_A[pos] != "X" and temp_parent_hap_B[pos] != "X":
				#if pos not in parameter.person_dict[p_id].common_fragment_dict:
				print >> x_file, parameter.rsID_dict[pos], pos, temp_parent_hap_A[pos], temp_parent_hap_B[pos]
	"""

def output_fragment(p_id, fragment_sorted_list):
	with open(p_id + "_fragment.txt", "w") as f_file:
		for data in fragment_sorted_list:
			pos = data[0]
			fragment_list = data[1]
			for f1 in fragment_list:
				print >> f_file, p_id, f1.start, f1.end


def update_parent_hap(f_id, fragment, temp_parent_hap_A, p_code):
	if len(temp_parent_hap_A) == 0:
		# for temp_parent_hap_B
		for pos in parameter.pos_list:
			if pos >= fragment.start and pos <= fragment.end:
				temp_parent_hap_A[pos] = parameter.person_dict[fragment.ID].haplotype[pos][p_code]
		return "same_to_A", temp_parent_hap_A
	else:
		current_start, current_end = get_parent_hap_end(temp_parent_hap_A)

		overlap_start = 0
		overlap_end = 0

		if fragment.end <= current_start or fragment.start >= current_end:
			return "no_overlap", temp_parent_hap_A
		else:
			if fragment.start < current_start and fragment.end < current_end:
				overlap_start = current_start
				overlap_end = fragment.end
			elif fragment.start > current_start and fragment.end < current_end:
				overlap_start = fragment.start
				overlap_end = fragment.end
			else:
				overlap_start = fragment.start
				overlap_end = current_end

		same = 0
		not_same = 0
		for pos in parameter.person_dict[f_id].hetero_pos_list:
			if pos >= overlap_start and pos <= overlap_end:
				if parameter.person_dict[fragment.ID].haplotype[pos][p_code] != "X" and temp_parent_hap_A[pos] != "X" \
						and parameter.person_dict[fragment.ID].haplotype[pos][p_code] != "N" and temp_parent_hap_A[
					pos] != "N":
					#print "gxxx", pos, parameter.person_dict[fragment.ID].haplotype[pos][0], temp_parent_hap[pos]
					if parameter.person_dict[fragment.ID].haplotype[pos][p_code] == temp_parent_hap_A[pos]:
						same += 1
					else:
						not_same += 1

		percentage = float(same) / (same + not_same + 1)
		if percentage >= 0.99:
			for pos in parameter.person_dict[f_id].hetero_pos_list:
				if pos > fragment.start and pos < fragment.end and pos not in temp_parent_hap_A:
					temp_parent_hap_A[pos] = parameter.person_dict[fragment.ID].haplotype[pos][p_code]
			return "same_to_A", temp_parent_hap_A
		elif percentage <= 0.01:
			return "same_to_B", temp_parent_hap_A
		else:
			return "mix", temp_parent_hap_A


def get_parent_hap_end(temp_parent_hap):
	pos_list = temp_parent_hap.keys()
	pos_list.sort()
	return pos_list[0], pos_list[-1]


def compare_child_hap(p_id, p_code, child_ID_1, child_ID_2):
	child_fragment_dict = {}
	child_fragment_dict[child_ID_1] = []
	child_fragment_dict[child_ID_2] = []

	child_1 = parameter.person_dict[child_ID_1]
	child_2 = parameter.person_dict[child_ID_2]

	child_fragment_dict[child_ID_1].append(parameter.pos_list[0])
	child_fragment_dict[child_ID_2].append(parameter.pos_list[0])

	compare = ""

	for pos in parameter.pos_list:
		pos = int(pos)
		f_geno = parameter.person_dict[p_id].genotype_dict[pos]

		if f_geno[0] != "N" and f_geno[1] != "N" and f_geno[0] != f_geno[1] and f_geno[0] != "X" and f_geno[1] != "X":

			cf_hap_1 = child_1.haplotype[pos][p_code]
			cf_hap_2 = child_2.haplotype[pos][p_code]

			if cf_hap_1 != "X" and cf_hap_2 != "X" and cf_hap_1 != "N" and cf_hap_2 != "N":

				if compare == "":
					compare = "same" if cf_hap_1 == cf_hap_2 else "different"

				elif compare == "same":
					if cf_hap_1 != cf_hap_2:
						child_fragment_dict[child_ID_1].append(pos)
						child_fragment_dict[child_ID_2].append(pos)
						compare = "different"
				elif compare == "different":
					if cf_hap_1 == cf_hap_2:
						child_fragment_dict[child_ID_1].append(pos)
						child_fragment_dict[child_ID_2].append(pos)
						compare = "same"

	child_fragment_dict[child_ID_1].append(parameter.pos_list[-1])
	child_fragment_dict[child_ID_2].append(parameter.pos_list[-1])

	for id in child_ID_1, child_ID_2:
		for index, pos in enumerate(child_fragment_dict[id]):
			if index < len(child_fragment_dict[id]) - 1:
				fragment = fragments()
				fragment.ID = id
				fragment.start = int(child_fragment_dict[id][index])
				fragment.end = int(child_fragment_dict[id][index + 1])
				fragment.length = parameter.pos_list.index(fragment.end) - parameter.pos_list.index(fragment.start)

				if id not in parameter.person_dict[p_id].fragment_dict:
					parameter.person_dict[p_id].fragment_dict[id] = []
				parameter.person_dict[p_id].fragment_dict[id].append(fragment)


def output_hap_std_geno():
	id_list = parameter.person_dict.keys()
	id_list.sort()
	for id in id_list:
		with open(id + "_std.txt", "w") as std_file:
			print >> std_file, "rs#", "pos", parameter.person_dict[id].ID + "_A", parameter.person_dict[
				                                                                        id].ID + "_B"
			for pos in parameter.pos_list:
				print >> std_file, parameter.rsID_dict[pos], pos,
				geno = parameter.person_dict[id].genotype_dict[pos]
				print >> std_file, geno[0], geno[1]

		with open(id + "_geno.txt", "w") as geno_file:
			print >> geno_file, "rs#", "pos", parameter.person_dict[id].ID
			for pos in parameter.pos_list:
				print >> geno_file, parameter.rsID_dict[pos], pos,
				geno = parameter.person_dict[id].genotype_dict[pos]
				print >> geno_file, geno


def output_hap_std(id):
	with open(id + "_std.txt", "w") as std_file:
		print >> std_file, "rs#", "pos", parameter.person_dict[id].ID + "_A", parameter.person_dict[id].ID + "_B"
		for pos in parameter.pos_list:
			print >> std_file, parameter.rsID_dict[pos], pos,
			geno = parameter.person_dict[id].genotype_dict[pos]
			print >> std_file, geno[0], geno[1]


def output_l3_ref(p_id):

	with open(p_id + "_ref.txt", "w") as c_hap_file:
		print >> c_hap_file, "rsID", "pos",
		for id in sorted(parameter.person_dict[p_id].children.keys()):
			print >> c_hap_file, id + "_F", id + "_M",
		print >> c_hap_file, ""

		for pos in parameter.pos_list:
			keep_pos = True
			for child_id in parameter.person_dict[p_id].children.keys():
				if parameter.person_dict[child_id].haplotype[pos][0] == "X" or \
					parameter.person_dict[child_id].haplotype[pos][1] == "X" or \
					parameter.person_dict[child_id].haplotype[pos][0] == "N" or \
					parameter.person_dict[child_id].haplotype[pos][1] == "N":
					keep_pos = False
					break
			if keep_pos:
				print >> c_hap_file, parameter.rsID_dict[pos], pos,
				for child_id in sorted(parameter.person_dict[p_id].children.keys()):
					print >> c_hap_file, parameter.person_dict[child_id].haplotype[pos][0], \
						parameter.person_dict[child_id].haplotype[pos][1],
				print >> c_hap_file, ""

def output_child_hap():

	with open(parameter.children_hap_file, "w") as c_hap_file:
		print >> c_hap_file, "rs_ID", "pos",
		for id in parameter.children_list:
			print >> c_hap_file, id + "_F", id + "_M",
		print >> c_hap_file, ""

		for pos in parameter.pos_list:
			print >> c_hap_file, parameter.rsID_dict[pos], pos,
			for child_id in parameter.children_list:
				print >> c_hap_file, parameter.person_dict[child_id].haplotype[pos][0], \
					parameter.person_dict[child_id].haplotype[pos][1],
			print >> c_hap_file, ""

def output_each_child_hap():

	for child_id in parameter.children_list:
		with open(child_id + "_hap.txt", "w") as hap_file:
				print >> hap_file, "rs#", "pos", parameter.person_dict[id].ID
				for pos in parameter.pos_list:
					child_hap = parameter.person_dict[child_id].haplotype[pos]
					if child_hap[0] != "X" and child_hap[1] != "X" and child_hap[0] != "N" and child_hap[1] != "N":
						print >> hap_file, parameter.rsID_dict[pos], pos,
						print >> hap_file, child_hap[0], child_hap[1]

def output_parent_hap():
	with open(parameter.parent_hap_file, "w") as parent_hap:
		print >> parent_hap, "rs#", "pos",
		for id in parameter.parent_id_list:
			print >> parent_hap, parameter.person_dict[id].ID + "_A", parameter.person_dict[id].ID + "_B",
		print >> parent_hap, ""

		for pos in parameter.pos_list:
			print >> parent_hap, parameter.rsID_dict[pos], pos,
			for id in parameter.parent_id_list:
				print >> parent_hap, parameter.person_dict[id].haplotype[pos][0], \
					parameter.person_dict[id].haplotype[pos][1],
			print >> parent_hap, ""

def output_seed():
	id_list = parameter.person_dict.keys()
	id_list.sort()
	for id in id_list:
		with open(id + "_seed.txt", "w") as seed_file:
			print >> seed_file, "rs#", "pos", parameter.person_dict[id].ID + "_A"
			for pos in parameter.pos_list:
				hap = parameter.person_dict[id].haplotype[pos]
				if hap[0] != "N" and hap[0] != "X" and hap[0] != "":
					print >> seed_file, parameter.rsID_dict[pos], pos, hap[0]

def genome_laser(pedi_name, geno_name):
	#pedi_name = "pedi.txt"
	#geno_name = "geno.txt"

	load_pedi(pedi_name)
	load_geno(geno_name)

	prepare_id_list()

	print "Laser I"
	start_time = time.time()
	parents_to_children()
	print "Laser I run time: ", round(time.time() - start_time, 2), "s"

	output_child_hap()
	#output_hap_std_geno()

	print "Laser II"
	start_time = time.time()

	for f_id in parameter.father_list:
		children_to_parents(f_id, 0)
	for m_id in parameter.mather_list:
		children_to_parents(m_id, 1)
	print "Laser II run time: ", round(time.time() - start_time, 2), "s"

	output_parent_hap()

	#output_seed()

	laser_three()
	laser_four()

def laser_three():
	print "Laser III"
	current_path = os.getcwd() + "/"
	os.system("mkdir -p " + current_path + "laser_3")
	os.chdir(current_path + "laser_3")

	for id in parameter.parent_id_list:
		os.system("mkdir -p " + current_path + "laser_3/" + id)
		os.chdir(id)
		#os.system("cp " + current_path + id + "_std.txt ./")
		output_l3_ref(id)

		refMerger_laser2(id, parameter)

		os.system("rm " + id + "_ref.txt")
		
		#start_time = time.time()
		hifi_run_code = subprocess_execute(parameter.hifi_file)
		#print "Laser III run time: ", id, round(time.time() - start_time, 2), "s"
		if hifi_run_code != 0:
			print "Laser 3 is done", id
		else:
			print "Laser 3 is running", id
			pass

		os.chdir(current_path + "laser_3")
	os.chdir(current_path)


def laser_four():
	print "Laser IV"
	current_path = os.getcwd() + "/"
	os.system("mkdir -p " + current_path + "laser_4")
	os.chdir(current_path + "laser_4")

	for id in parameter.person_dict.keys():
		os.system("mkdir -p " + current_path + "laser_4/" + id)
		os.chdir(id)
		#os.system("cp " + current_path + id + "_std.txt ./")
		refMerger(id, parameter)

		#start_time = time.time()
		hifi_run_code = subprocess_execute(parameter.hifi_file)
		#print "Laser IV run time: ", id, round(time.time() - start_time, 2), "s"

		if hifi_run_code != 0:
			print "Laser 4 is done", id
		else:
			print "Laser 4 is running", id
			pass
		os.chdir(current_path + "laser_4")
	os.chdir(current_path)


def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-p", "--pedi", type="string", dest="pedi_name", help="Input pedi name", default="null")
	parser.add_option("-g", "--geno", type="string", dest="geno_file", help="Input file name", default="null")
	parser.add_option("-r", "--reference", type="string", dest="ref_file", help="Input reference name", default="null")
	(options, args) = parser.parse_args()

	return options


####################
# #def for refMerge #
####################


class refs:
	def __init__(self):
		self.seed_title = ""
		self.seed_dict = {}
		self.geno_title = ""
		self.geno_dict = {}
		self.ref_title = ""
		self.ref_dict = {}


def load_hap_ref_data_single(ref_file_name):
	raw_data_format_ref = "list"
	ref_title_info, hap_ref_dict = load_raw_data(ref_file_name, raw_data_format_ref)
	return ref_title_info, hap_ref_dict

def removeN(hifi_std_dict):
	temp_dict = {}
	for position, elements in hifi_std_dict.iteritems():
		if elements[2].strip() != "N" and elements[3].strip() != "N":
			temp_dict[position] = hifi_std_dict[position]
	return temp_dict

def load_raw_data(file_name, raw_data_format="list"):
	title_info = ""
	data = {}
	with open(file_name, "r") as fp:
		for line in fp:
			if line.startswith("rsID"):
				title_info = list_to_line(line.strip().split())
			else:
				elements = line.strip().split()
				try:
					if raw_data_format == "list":
						data[int(elements[1])] = elements
					elif raw_data_format == "string":
						data[int(elements[1])] = line.strip()
				except:
					#print "error in ", line, file_name
					pass
	return title_info, data

def compare_geno_ref(geno_dict, hap_ref_dict):

	geno_x_dict = {}
	geno_n_dict = {}
	geno_ref_not_consistent = {}
	ref_homo_dict = {}

	# to remove homo snps in ref
	for position in hap_ref_dict.keys():
		#print position
		alleles = hap_ref_dict[position]
		#print "gxxx", alleles
		alleles = alleles[2:]
		unique_alleles = list(set(alleles))
		n_alleles = len(unique_alleles)
		if n_alleles == 1:
			ref_homo_dict[position] = list_to_line(unique_alleles)

	# to remove snps that are conflict in ref and geno
	for position, snp in geno_dict.iteritems():
		if position in hap_ref_dict:
			exist_in_ref = False
			geno_A = snp[0]
			geno_B = snp[1]
			alleles = hap_ref_dict[position]
			alleles = alleles[2:]
			unique_alleles = list(set(alleles))

			n_alleles = len(unique_alleles)
			if n_alleles == 0 or n_alleles > 2:
				pass
			else:
				try:
					if geno_A == geno_B:
						if n_alleles == 2:
							if geno_A == unique_alleles[0] or geno_A == unique_alleles[1]:
								exist_in_ref = True
						if n_alleles == 1:
							if geno_A == unique_alleles[0]:
								exist_in_ref = True
					else:
						if n_alleles == 2:
							if (geno_A == unique_alleles[0] and geno_B == unique_alleles[1]) or (geno_A == unique_alleles[1] and geno_B == unique_alleles[0]):
								exist_in_ref = True
							if n_alleles == 1:
								pass
					if not exist_in_ref:
						if geno_A == 'X':
							geno_x_dict[position] = list_to_line(snp)
						if geno_A == 'N':
							geno_n_dict[position] = list_to_line(snp)
						else:
							geno_ref_not_consistent[position] = list_to_line(snp)
				except:
					print position, unique_alleles, n_alleles
					pass

	return ref_homo_dict, geno_ref_not_consistent, geno_n_dict

def dict_substract(large_dict, small_dict):
	return {index: value for index, value in large_dict.iteritems() if index not in small_dict}

def ref_preprocess(geno_dict, hap_ref_dict):
	ref_homo_dict, geno_ref_not_consistent, geno_n_dict = compare_geno_ref(geno_dict, hap_ref_dict)
	hap_ref_dict = dict_substract(hap_ref_dict, geno_ref_not_consistent)
	hap_ref_dict = dict_substract(hap_ref_dict, geno_n_dict)
	#hap_ref_dict = dict_substract(hap_ref_dict, ref_homo_dict)
	return hap_ref_dict

def group_seed(seed_dict, geno_dict):
	seed_homo_dict = {}
	seed_hetero_dict = {}
	for position, snp in seed_dict.iteritems():
		if position in geno_dict:  # pos in seed may not be in geno
			geno_allele = geno_dict[position]
			if geno_allele[0] == geno_allele[1]:
				seed_homo_dict[position] = snp
			else:
				seed_hetero_dict[position] = snp
		else:
			seed_hetero_dict[position] = snp
	return (seed_homo_dict, seed_hetero_dict)


def output_hap(file_name, parameter, seed_dict):
	with open(file_name, "w") as output:
		print >> output, ref.seed_title
		for pos in sorted(seed_dict.keys()):
			print >> output, parameter.rsID_dict[pos], pos, seed_dict[pos][0]

def output_geno(file_name, parameter, geno_dict):
	with open(file_name, "w") as output:
		print >> output, ref.geno_title
		for pos in sorted(geno_dict.keys()):
			print >> output, parameter.rsID_dict[pos], pos, geno_dict[pos][0]+geno_dict[pos][1]

def output_ref(file_name, ref_dict):
	with open(file_name, "w") as output:
		print >> output, ref.ref_title
		for pos in sorted(ref_dict.keys()):
			print >> output, ref_dict[pos].strip()

def make_hifi_files(ref, parameter):

	seed_dict = ref.seed_dict
	geno_dict = ref.geno_dict
	hap_ref_dict = ref.ref_dict

	hifi_seed_dict = {}
	hifi_geno_dict = {}
	hifi_ref_dict = {}

	last_seed_position = sort_dict_by_key(seed_dict)[-1][0]

	hap_ref_sorted_list = sort_dict_by_key(hap_ref_dict)

	for snp in hap_ref_sorted_list:
		position = snp[0]
		if position <= int(last_seed_position):
			hifi_ref_dict[position] = list_to_line(hap_ref_dict[position])
			if position in seed_dict:
				hifi_seed_dict[position] = seed_dict[position][0]
			elif position in geno_dict and geno_dict[position][0] == geno_dict[position][1]:
				hifi_seed_dict[position] = geno_dict[position][0]
			if position in geno_dict:
				hifi_geno_dict[position] = geno_dict[position]

	output_hap("haplotype.txt", parameter, hifi_seed_dict)
	output_geno("genotype.txt", parameter, hifi_geno_dict)
	output_ref("refHaplos.txt", hifi_ref_dict)


def find_index(a, b):
	index_list = []
	for str_b in b:
		for id in a:
			if id in str_b:
				index_list.append(b.index(str_b))
	return index_list


def remove_element(index_list, list):
	for index in sorted(index_list, reverse=True):
		del list[index]


def remove_hap_in_ref(id, parameter):
	ref_title = copy.deepcopy(parameter.ori_ref_title)
	ref_dict = copy.deepcopy(parameter.ori_ref_dict)

	id_list = []
	id_list.append(id)
	if parameter.person_dict[id].father != "N/A":
		id_list.append(parameter.person_dict[id].father)
	if parameter.person_dict[id].mather != "N/A":
		id_list.append(parameter.person_dict[id].mather)
	#print id, id_list

	index_list = find_index(id_list, ref_title)
	remove_element(index_list, ref_title)

	for pos in ref_dict.keys():
		remove_element(index_list, ref_dict[pos])

	return ref_title, ref_dict


def refMerger(id, parameter):

	person = parameter.person_dict[id]

	global ref
	ref = refs()
	#print id
	ref.ref_title, ref.ref_dict = remove_hap_in_ref(id, parameter)

	ref.seed_title = "rsID pos " + id + "_A"

	for pos in person.haplotype.keys():
		hap = person.haplotype[pos]
		if hap[0] != "N" and hap[0] != "X" and hap[0] != "":
			ref.seed_dict[pos] = hap[0]

	ref.geno_title = "rsID pos " + id
	ref.geno_dict = person.genotype_dict

	ref.ref_dict = ref_preprocess(ref.geno_dict, ref.ref_dict)

	make_hifi_files(ref, parameter)


def refMerger_laser2(id, parameter):

	person = parameter.person_dict[id]

	global ref
	ref = refs()

	ref.ref_title, ref.ref_dict = load_hap_ref_data_single(id + "_ref.txt")

	ref.seed_title = "rsID pos " + id + "_A"

	for pos in person.haplotype.keys():
		hap = person.haplotype[pos]
		if hap[0] != "N" and hap[0] != "X" and hap[0] != "":
			ref.seed_dict[pos] = hap[0]

	ref.geno_title = "rsID pos " + id
	ref.geno_dict = person.genotype_dict

	ref.ref_dict = ref_preprocess(ref.geno_dict, ref.ref_dict)

	make_hifi_files(ref, parameter)


def subprocess_execute(command, time_out=300):
	c = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	# print c.communicate()
	t = 0
	while t <= time_out and c.poll() is None:
		time.sleep(1)  # (comment 1)
		t += 1
		returncode = 1
	if c.poll() is None:
		c.terminate()
		returncode = 0  # (comment 2)
	else:
		eturncode = c.poll()
		# print "c.poll()", c.poll()
	return returncode


if __name__ == '__main__':
	options = get_args()
	pedi_name = options.pedi_name
	geno_name = options.geno_file
	global raw_data_format
	raw_data_format = "list"

	global parameter
	parameter = parameters()
	parameter.ref_file_name = options.ref_file

	#pedi_name = "pedi.txt"
	#geno_name = "geno.txt"
	#parameter.ref_file_name = "hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_ceu.phased"

	parameter.ori_ref_title, parameter.ori_ref_dict = load_raw_data(parameter.ref_file_name)

	start_time = time.time()

	current_path = os.getcwd()
	parameter.hifi_file = current_path + parameter.hifi_file
	print parameter.hifi_file
	genome_laser(pedi_name, geno_name)

	print "elapsed_time is: ", round(time.time() - start_time, 2), "s"
