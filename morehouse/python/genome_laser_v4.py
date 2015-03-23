#!/usr/bin/python
# ######################################################################################
# Guoxing Fu Jan 28, 2015
# Laser project, to impute haplotype from trio data
# Mar. 15 2015, Add hifi to the package
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy, shutil, fnmatch
from optparse import OptionParser
from hifiAccuCheck_v2 import *
from refMerg_laser import *

class parameters:
	def __init__(self):
		self.person_dict = {}
		self.rsID_dict = {}
		self.pos_list = {}
		self.father_list = []
		self.mather_list = []
		self.children_list = []
		self.children_dict = {}

		self.fragment_dict = {}
		self.fragment_list = []

		self.children_hap_file = "child_hap.txt"

		self.ref_file_name = ""
		self.ori_ref_title = ""
		self.ori_ref_dict = {}


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
		#list = parameter.fragment_dict[child_id]
		for fragment in list:
			#print fragment.ID, fragment.start, fragment.end
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
	#print children_list

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

	fragment__sorted_list = sort_dict_by_key(fragment_sort_dict)

	with open(p_id + "_fragment.txt", "w") as f_file:
		for data in fragment__sorted_list:
			pos = data[0]
			fragment_list = data[1]
			for f1 in fragment_list:
				print >> f_file, p_id, f1.start, f1.end

	#print len(parameter.person_dict[p_id].hetero_pos_list)
	#with open(p_id + "_fragment.txt", "w") as f_file:
	for f1 in child_1_2_3:
		#print >> f_file, p_id, f1.start, f1.end, f1.length
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

	with open(p_id + "_hap.txt", "w") as x_file:
		print >> x_file, "rs#", "pos", parameter.person_dict[p_id].ID + "_A", parameter.person_dict[p_id].ID + "_B"
		for pos in parameter.pos_list:
			if pos in temp_parent_hap_A and pos in temp_parent_hap_B \
					and temp_parent_hap_A[pos] != "N" and temp_parent_hap_B[pos] != "N" \
					and temp_parent_hap_A[pos] != "X" and temp_parent_hap_B[pos] != "X":
				#if pos not in parameter.person_dict[p_id].common_fragment_dict:
				print >> x_file, parameter.rsID_dict[pos], pos, temp_parent_hap_A[pos], temp_parent_hap_B[pos]


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
					#print "xxxx", pos, parameter.person_dict[fragment.ID].haplotype[pos][0], temp_parent_hap[pos]
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

"""
def output_fragment():
	p_id = "NA07347"
	child_id_1 = "NA11881_NA07347_1"
	child_id_2 = "NA11881_NA07347_2"
	child_id_3 = "NA11881_NA07347_3"

	with open("NA11881_NA07347_1_0_frag.txt", "w") as c_frag:
		print >> c_frag, "pos", "f_geno_0", "f_geno_1", "c_hap_0"
		fragment = parameter.person_dict[p_id].fragment_dict[child_id_1][3]
		for pos in parameter.person_dict[p_id].hetero_pos_list:
			if fragment.start <= pos <= fragment.end and parameter.person_dict[child_id_1].haplotype[pos][0] != "X":
				f_geno = parameter.person_dict[p_id].genotype_dict[pos]
				print >> c_frag, pos, f_geno, \
					parameter.person_dict[child_id_1].genotype_dict[pos], \
					parameter.person_dict[child_id_1].haplotype[pos][0], \
					parameter.person_dict[child_id_2].genotype_dict[pos], \
					parameter.person_dict[child_id_2].haplotype[pos][0], \
					parameter.person_dict[child_id_3].genotype_dict[pos], \
					parameter.person_dict[child_id_3].haplotype[pos][0]
"""

def output_hap_std_geno():
	#parent_id_list = []
	#parent_id_list.extend(parameter.father_list)
	#parent_id_list.extend(parameter.mather_list)
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

def output_parent_hap():
	parent_id_list = []
	parent_id_list.extend(parameter.father_list)
	parent_id_list.extend(parameter.mather_list)
	parent_id_list.sort()

	with open("parent_hap.txt", "w") as parent_hap:
		print >> parent_hap, "rs#", "pos",
		for id in parent_id_list:
			print >> parent_hap, parameter.person_dict[id].ID + "_A", parameter.person_dict[id].ID + "_B",
		print >> parent_hap, ""

		for pos in parameter.pos_list:

			print >> parent_hap, parameter.rsID_dict[pos], pos,
			for id in parent_id_list:
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
				print >> seed_file, parameter.rsID_dict[pos], pos,
				hap = parameter.person_dict[id].haplotype[pos]
				if hap[0] != "N" and hap[0] != "X" and hap[0] != "":
					print >> seed_file, hap[0]


def genome_laser(pedi_name, geno_name):
	pedi_name = "pedi.txt"
	geno_name = "geno.txt"

	load_pedi(pedi_name)
	load_geno(geno_name)

	prepare_id_list()
	parents_to_children()
	output_child_hap()

	output_hap_std_geno()

	for f_id in parameter.father_list:
		children_to_parents(f_id, 0)
		#print f_id
		hifiAccuCheck_file(f_id + "_hap.txt", f_id + "_std.txt")

	for m_id in parameter.mather_list:
		children_to_parents(m_id, 1)
		#print m_id
		hifiAccuCheck_file(m_id + "_hap.txt", m_id + "_std.txt")

	output_parent_hap()

	output_seed()

	refMerger("NA12763", parameter)



def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-p", "--pedi", type="string", dest="pedi_name", help="Input pedi name", default="null")
	parser.add_option("-g", "--geno", type="string", dest="geno_file", help="Input file name", default="null")
	parser.add_option("-r", "--reference", type="string", dest="ref_file", help="Input reference name", default="null")
	(options, args) = parser.parse_args()

	return options



##################
#refMerg
###########

"""
def load_hap_ref_data_single(ref_file_name):
	raw_data_format_ref = "list"
	ref_title_info, hap_ref_dict = load_raw_data(ref_file_name, raw_data_format_ref)
	return ref_title_info, hap_ref_dict

def compare_geno_ref(geno_dict, hap_ref_dict):

	geno_x_dict = {}
	geno_n_dict = {}
	geno_ref_not_consistent = {}
	ref_homo_dict = {}

	# to remove homo snps in ref
	for position in hap_ref_dict.keys():
		alleles = hap_ref_dict[position]
		alleles = alleles[2:]
		unique_alleles = list(set(alleles))
		n_alleles = len(unique_alleles)
		if n_alleles == 1:
			ref_homo_dict[position] = list_to_line(unique_alleles)

	# to remove snps that are conflict in ref and geno
	for position, snp in geno_dict.iteritems():
		if position in hap_ref_dict:
			exist_in_ref = False
			geno_A = snp[2][0]
			geno_B = snp[2][1]
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
							if n_alleles == 1:	# hetero_geno, homo_ref
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

	print "geno_x_dict: ", len(geno_x_dict)
	print "geno_n_dict: ", len(geno_n_dict)
	print "geno_ref_not_consistent: ", len(geno_ref_not_consistent)
	print "ref_homo", len(ref_homo_dict)
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
			geno_allele = geno_dict[position][2]
			if geno_allele[0] == geno_allele[1]:
				seed_homo_dict[position] = snp
			else:
				seed_hetero_dict[position] = snp
		else:
			seed_hetero_dict[position] = snp
	return (seed_homo_dict, seed_hetero_dict)

def output_files(file_name, title_info, dict):
	outpuf_file = open(currentPath + file_name, "w")    # for hifi
	print >> outpuf_file, title_info
	sorted_list = sort_dict_by_key(dict)
	for element in sorted_list:
		print >> outpuf_file, element[1]
	outpuf_file.close()

def make_hifi_files(seed_dict, geno_dict, hap_ref_dict):
	common_snp_number = 0
	hifi_seed_dict = {}
	hifi_geno_dict = {}
	hifi_ref_dict = {}


	last_seed_position = sort_dict_by_key(seed_dict)[-1][0]
	#print "last_haplotype_position", last_seed_position
	hap_ref_sorted_list = sort_dict_by_key(hap_ref_dict)

	seed_homo_dict, seed_hetero_dict = group_seed(seed_dict, geno_dict)

	hap_ref_sorted_list = [x for x in hap_ref_sorted_list if x[0] not in seed_hetero_dict]


	hap_ref_size = len(hap_ref_dict)
	for i in range(int(remPercent*hap_ref_size)):
		if len(hap_ref_sorted_list) > 1:
			random_index = random.randrange(0, (len(hap_ref_sorted_list)-1))
			position = hap_ref_sorted_list[random_index][0]
			if position in hap_ref_dict:
				del hap_ref_dict[int(position)]
				del hap_ref_sorted_list[random_index]

	hap_ref_sorted_list = sort_dict_by_key(hap_ref_dict)

	for snp in hap_ref_sorted_list:
		position = snp[0]
		if position <= int(last_seed_position):
			hifi_ref_dict[position] = list_to_line(hap_ref_dict[position])
			if position in seed_dict:
				hifi_seed_dict[position] = list_to_line(seed_dict[position])
			elif position in geno_dict and geno_dict[position][2][0] == geno_dict[position][2][1]:
				hifi_seed_dict[position] = geno_dict[position][0] + " " + geno_dict[position][1] + " " + geno_dict[position][2][0]
			if position in geno_dict:
				hifi_geno_dict[position] = list_to_line(geno_dict[position])


	output_files("haplotype.txt", seed_title_info, hifi_seed_dict)
	output_files("genotype.txt", geno_title_info, hifi_geno_dict)
	output_files("refHaplos.txt", ref_title_info, hifi_ref_dict)

	seed_not_in_ref = dict_substract(seed_dict, hap_ref_dict)
	output_files("seed_not_in_ref.txt", seed_title_info, seed_not_in_ref)


def refMerger(haplotype_file, chr_name, remPercent):
	global seed_title_info
	global seed_dict
	global geno_title_info
	global geno_dict
	global hap_ref_dict
	global ref_title_info
	global raw_data_format_ref

	raw_data_format_ref = "list"

	ref_title_info, hap_ref_dict = load_hap_ref_data_single(chr_name)

	geno_title_info, geno_dict = load_raw_data(genotype_file, raw_data_format_ref)
	#print "genotype_file", genotype_file
	#print "total_geno_number: ", len(geno_dict)

	seed_file = haplotype_file
	seed_title_info, seed_dict = load_raw_data(seed_file, raw_data_format_ref)
	print "seed_file is :", seed_file

	print "hap_ref_dict before process", len(hap_ref_dict)
	hap_ref_dict = ref_preprocess(geno_dict, hap_ref_dict)
	print "hap_ref_dict after process", len(hap_ref_dict)

	make_hifi_files(remPercent)
"""
####################
# #def for HIFI# #
####################


def load_raw_data(hap_file):
	title = ""
	data = {}
	with open(hap_file, "r") as fp:
		for line in fp:
			if line != "":
				elements = line.strip().split()
				if line.startswith("rs_ID") or line.startswith("rs#") or line.startswith("rsID"):
					title = elements
				else:
					# print elements[1]
					data[int(elements[1])] = elements
	return title, data


def find_index(a, b):
	index_list = []
	for str_b in b:
		if a in str_b:
			# print b.index(str_b)
			index_list.append(b.index(str_b))
	return index_list


def remove_sample_from_ref(line, rm_list):  # Deleting multiple elements from a list by index
	# new_line = [i for j, i in enumerate(line) if i not in rm_list]
	new_line = line
	# print "length", len(line)
	for index in sorted(rm_list, reverse=True):
		# print index
		del new_line[index]
	return new_line


def make_child_HIFI_files(sample_id, raw_ref_title, raw_ref_data):
	sample_id_list = []
	id_rm_list = []
	# define the person temp
	person_temp = parameter.person_dict[sample_id]

	# make the ref
	# ref_title = copy.copy(raw_ref_title)
	# ref_data = copy.copy(raw_ref_data)
	ref_title = copy.deepcopy(raw_ref_title)
	ref_data = copy.deepcopy(raw_ref_data)
	sample_id_list.append(person_temp.father)
	sample_id_list.append(person_temp.mather)
	for item in sample_id_list:
		index = find_index(item, ref_title)
		for num in index:
			id_rm_list.append(num)
	# print "F&M", id_rm_list
	output_ref = open("_ref", "w")
	print >> output_ref, "rsID" + "\t" + list_to_line(remove_sample_from_ref(ref_title, id_rm_list)[1:])
	for pos in sorted(ref_data.keys()):
		# print pos
		print >> output_ref, list_to_line(remove_sample_from_ref(ref_data[pos], id_rm_list))
	output_ref.close()

	# make the geno
	output_geno = open("_geno", "w")
	print >> output_geno, "rsID" + "\t" + "pos#" + "\t" + sample_id
	for pos in parameter.pos_list:
		print >> output_geno, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.genotype_dict[pos]
	output_geno.close()

	# make the haplo
	output_haplo = open("_haplo", "w")
	print >> output_haplo, "rsID" + "\t" + "pos#" + "\t" + sample_id + "_A"
	for pos in parameter.pos_list:
		if pos in person_temp.haplotype:
			if person_temp.haplotype[pos][0] != "N" and person_temp.haplotype[pos][0] != "X":
				print >> output_haplo, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.haplotype[pos][0]
	output_haplo.close()

	# make the haplo after laser
	output_haplo = open("_haplo_laser", "w")
	print >> output_haplo, "rsID" + "\t" + "pos#" + "\t" + sample_id + "_A"
	for pos in parameter.pos_list:
		if pos in person_temp.haplotype:
			if person_temp.haplotype[pos][0] != "N" and person_temp.haplotype[pos][0] != "X" or \
									person_temp.haplotype[pos][1] != "N" and person_temp.haplotype[pos][1] != "X":
				print >> output_haplo, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.haplotype[pos][
					0] + "\t" + person_temp.haplotype[pos][1]
	output_haplo.close()

	# make the std_haplo
	output_std_haplo = open("_std_haplo", "w")
	print >> output_std_haplo, "rsID" + "\t" + "pos#" + "\t" + sample_id + "_A" + "\t" + sample_id + "_B"
	for pos in parameter.pos_list:
		print >> output_std_haplo, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.genotype_dict[pos][
			0] + "\t" + person_temp.genotype_dict[pos][1]
	output_std_haplo.close()

	print sample_id, "HIFI files ready"


def make_parent_HIFI_files(sample_id, raw_ref_title, raw_ref_data):
	# make the ref
	# ref_title = copy.copy(raw_ref_title)
	# ref_data = copy.copy(raw_ref_data)
	ref_title = copy.deepcopy(raw_ref_title)
	ref_data = copy.deepcopy(raw_ref_data)
	id_rm_list = find_index(sample_id, ref_title)
	# print id_rm_list
	output_ref = open("_ref", "w")
	print >> output_ref, "rsID" + "\t" + list_to_line(remove_sample_from_ref(ref_title, id_rm_list)[1:])
	for pos in sorted(ref_data.keys()):
		print >> output_ref, list_to_line(remove_sample_from_ref(ref_data[pos], id_rm_list))
	output_ref.close()

	# define the person temp
	person_temp = parameter.person_dict[sample_id]

	# make the geno
	output_geno = open("_geno", "w")
	print >> output_geno, "rsID" + "\t" + "pos#" + "\t" + sample_id
	for pos in parameter.pos_list:
		print >> output_geno, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.genotype_dict[pos]
	output_geno.close()

	# make the haplo
	output_haplo = open("_haplo", "w")
	print >> output_haplo, "rsID" + "\t" + "pos#" + "\t" + sample_id + "_A"
	for pos in parameter.pos_list:
		if pos in person_temp.haplotype:
			if person_temp.haplotype[pos][0] != "N" and person_temp.haplotype[pos][0] != "X":
				print >> output_haplo, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.haplotype[pos][0]
	output_haplo.close()

	# make the haplo after laser
	output_haplo = open("_haplo_laser", "w")
	print >> output_haplo, "rsID" + "\t" + "pos#" + "\t" + sample_id + "_A"
	for pos in parameter.pos_list:
		if pos in person_temp.haplotype:
			if person_temp.haplotype[pos][0] != "N" and person_temp.haplotype[pos][0] != "X" or \
									person_temp.haplotype[pos][1] != "N" and person_temp.haplotype[pos][1] != "X":
				print >> output_haplo, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.haplotype[pos][
					0] + "\t" + person_temp.haplotype[pos][1]
	output_haplo.close()

	# make the std_haplo
	output_std_haplo = open("_std_haplo", "w")
	print >> output_std_haplo, "rsID" + "\t" + "pos#" + "\t" + sample_id + "_A" + "\t" + sample_id + "_B"
	for pos in parameter.pos_list:
		print >> output_std_haplo, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.genotype_dict[pos][
			0] + "\t" + person_temp.genotype_dict[pos][1]
	output_std_haplo.close()

	print sample_id, "HIFI files ready"


def key_compare(k1, k2):
	temp = {}
	for a in k1.keys():
		if a in k2:
			temp[a] = 0
			# print "comon key is ", a
	return temp


def cut_by_common_append():
	ref_title, ref_data = load_raw_data("refHaplos_after_merge")
	geno_title, geno_data = load_raw_data("genotype_after_merge")
	haplo_title, haplo_data = load_raw_data("haplotype_after_merge")
	com_id1 = key_compare(ref_data, geno_data)
	com_id2 = key_compare(com_id1, haplo_data)
	mini_p = min(com_id2.keys())
	max_p = max(com_id2.keys())
	# print mini_p
	# print type(mini_p)
	# print max_p
	# print type(max_p)
	geno_output = open("genotype.txt", "w")
	haplo_output = open("haplotype.txt", "w")
	ref_output = open("refHaplos.txt", "w")
	ref_list_ori = sorted(ref_data, key=int)
	geno_list_ori = sorted(geno_data, key=int)
	haplo_list_ori = sorted(haplo_data, key=int)
	ref_list = []
	geno_list = []
	haplo_list = []
	# print haplo_list
	for i in ref_list_ori:
		if int(mini_p) <= int(i) <= int(max_p):
			ref_list.append(i)
	for j in geno_list_ori:
		if int(mini_p) <= int(j) <= int(max_p):
			geno_list.append(j)
	for q in haplo_list_ori:
		if int(mini_p) <= int(q) <= int(max_p):
			haplo_list.append(q)
	# print "new list"
	# print "ref_list[0]", ref_list[0]
	# print "ref_list[-1]", ref_list[-1]
	# print "geno_list[0]", geno_list[0]
	# print "geno_list[-1]", geno_list[-1]
	# print "haplo_list[0]", haplo_list[0]
	# print "haplo_list[-1]", haplo_list[-1]
	print >> geno_output, "rsID" + "\t" + "phys_position" + "\t" + list_to_line(geno_title[2:])
	print >> haplo_output, "rsID" + "\t" + "phys_position" + "\t" + list_to_line(haplo_title[2:])
	print >> ref_output, "rsID" + "\t" + "phys_position"
	for k in sorted(ref_list):
		print >> ref_output, list_to_line(ref_data[k])
	for d in sorted(geno_list):
		print >> geno_output, list_to_line(geno_data[d])
	for c in sorted(haplo_list):
		print >> haplo_output, list_to_line(haplo_data[c])
	geno_output.close()
	haplo_output.close()
	ref_output.close()


def subprocess_execute(command, time_out=600):
	c = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	# print c.communicate()
	t = 0
	while t <= time_out and c.poll() is None:
		time.sleep(1)  # (comment 1)
		t += 1
		returncode = 1
		# print t
	if c.poll() is None:
		c.terminate()
		returncode = 0  # (comment 2)
	else:
		eturncode = c.poll()
		# print "c.poll()", c.poll()
	return returncode


def check_hapfile_run_HIFI_child():
	os.chdir(current_path)
	if os.path.isfile("child_hap.txt") and os.path.isfile("parent_hap.txt"):
		print "prepare for HIFI running..."
		# raw_ref_title, raw_ref_data = load_raw_data(ref_name)
		sample_id_list = parameter.children_list
		# processing samples data start
		for sample_id in set(sample_id_list):
			print sample_id, "processing..."

			if os.path.exists(sample_id):  # remove the exist folder and create a new one
				shutil.rmtree(sample_id)
			os.makedirs(sample_id)

			os.chdir(sample_id)
			# print sample_id
			# ref_title = copy.copy(raw_ref_title)
			# ref_data = copy.copy(raw_ref_data)
			make_child_HIFI_files(sample_id, raw_ref_title, raw_ref_data)
			print "refMerging..."
			print subprocess.Popen("python " + refMerger + " -i _haplo -n _geno -r _ref", shell=True,
			                       stdout=subprocess.PIPE).stdout.read()
			os.rename("refHaplos.txt", "refHaplos_after_merge")
			os.rename("genotype.txt", "genotype_after_merge")
			os.rename("haplotype.txt", "haplotype_after_merge")
			cut_by_common_append()  # format required by hifi
			os.system("unix2dos refHaplos.txt")
			os.system("unix2dos genotype.txt")
			os.system("unix2dos haplotype.txt")
			hifi_starttime = time.time()
			hifi_run_code = subprocess_execute(hifi_file)
			print "hifi run time with", sample_id, round((time.time() - hifi_starttime), 6), "s"
			# print hifi_run_code
			if hifi_run_code != 0:
				print "==done=="
				print ""
			else:
				print "hifi has an issue", sample_id
				pass
			os.chdir(current_path)
			# processing samples data end
	else:
		print "Warning! child_hap.txt or parent_hap.txt file missing..."


def check_hapfile_run_HIFI_parent():
	os.chdir(current_path)
	if os.path.isfile("child_hap.txt") and os.path.isfile("parent_hap.txt"):
		print "prepare for HIFI running..."
		sample_id_list = parameter.person_dict.keys()
		# processing samples data start
		for sample_id in set(sample_id_list):
			if sample_id not in parameter.children_list:
				print sample_id, "processing..."

				if os.path.exists(sample_id):  # remove the exist folder and create a new one
					shutil.rmtree(sample_id)
				os.makedirs(sample_id)

				os.chdir(sample_id)
				# print sample_id
				# ref_title = copy.copy(raw_ref_title)
				# ref_data = copy.copy(raw_ref_data)
				make_parent_HIFI_files(sample_id, raw_ref_title, raw_ref_data)
				print "refMerging..."
				print subprocess.Popen("python " + refMerger + " -i _haplo -n _geno -r _ref", shell=True,
				                       stdout=subprocess.PIPE).stdout.read()
				os.rename("refHaplos.txt", "refHaplos_after_merge")
				os.rename("genotype.txt", "genotype_after_merge")
				os.rename("haplotype.txt", "haplotype_after_merge")
				cut_by_common_append()  # format required by hifi
				os.system("unix2dos refHaplos.txt")
				os.system("unix2dos genotype.txt")
				os.system("unix2dos haplotype.txt")
				hifi_starttime = time.time()
				hifi_run_code = subprocess_execute(hifi_file)
				print "hifi run time with", sample_id, round((time.time() - hifi_starttime), 6), "s"
				# print hifi_run_code
				if hifi_run_code != 0:
					print "==done=="
					print ""
				else:
					print "hifi has an issue", sample_id
					pass
				os.chdir(current_path)
				# processing samples data end
	else:
		print "Warning! child_hap.txt or parent_hap.txt file missing..."


if __name__ == '__main__':
	options = get_args()
	pedi_name = options.pedi_name
	geno_name = options.geno_file
	ref_name = options.ref_file
	#global person_dict
	#person_dict = {}

	global parameter
	parameter = parameters()

	parameter.ref_file_name = ref_name

	parameter.ref_file_name = "hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_ceu.phased"

	parameter.ori_ref_title, parameter.ori_ref_dict = load_raw_data(parameter.ref_file_name)

	pedi_name = "pedi.txt"
	geno_name = "geno.txt"

	start_time = time.time()
	print "Laser I II ..."
	genome_laser(pedi_name, geno_name)

	###################
	# # call HIFI # #
	###################
	print "Laser I II done: ", round(time.time() - start_time, 2), "s"
	current_path = os.getcwd()
	output_path = current_path + "/output_file"
	ref_folder = current_path + "/ref_pool"
	hifi_file = current_path + "/scripts/hifi_fu_ref_for_HIFILOCA.ref"
	refMerger = current_path + "/scripts/refMerger_v5_wli_filltheend_remove_extra.py"
	acc_check_file = current_path + "/scripts/hifiAccuCheck_v3_pos_num.py"
	# #HIFI processing
	#raw_ref_title, raw_ref_data = load_raw_data(ref_name)
	#check_hapfile_run_HIFI_child()
	#check_hapfile_run_HIFI_parent()
	print "elapsed_time is: ", round(time.time() - start_time, 2), "s"
