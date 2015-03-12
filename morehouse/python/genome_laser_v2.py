#!/usr/bin/python
#######################################################################################
# Guoxing Fu Jan 28, 2015
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
#from tools import *
from hifiAccuCheck_v2 import *

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

class persons:
	def __init__(self):
		self.ID = ""
		self.father = ""
		self.mather = ""
		self.children = {}
		self.genotype_dict = {}
		self.haplotype = {}
		self.hetero_pos_list = []

		self.fragment_dict = {}


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
					if person.father != "N/A" and person.father in parameter.person_dict and person.ID not in parameter.person_dict[person.father].children:
						parameter.person_dict[person.father].children[person.ID] = ""
					if person.mather != "N/A" and person.mather in parameter.person_dict and person.ID not in parameter.person_dict[person.mather].children:
						parameter.person_dict[person.mather].children[person.ID] = ""
					if person.ID not in person_dict:
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
				except:
					print "error in geno", line, pedi_name
					#pass

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
						if (c_geno[0] == f_geno[0] or c_geno[0] == f_geno[1]) and (c_geno[0] == m_geno[0] or c_geno[0] == m_geno[1]):
							parameter.person_dict[ID].haplotype[pos] = (c_geno[0], c_geno[1])
						"""
						# keep this part. for non-NN c_geno, has discrepancy with f_geno or m_geno. uesful for future.
						else:
							if f_geno == "NN" and m_geno == "NN":
								parameter.person_dict[ID].haplotype[pos] = (c_geno[0], c_geno[1])
							elif f_geno == "NN" and m_geno != "NN" and (c_geno[0] == m_geno[0] or c_geno[0] == m_geno[1]):
								parameter.person_dict[ID].haplotype[pos] = (c_geno[0], c_geno[1])
							elif f_geno != "NN" and m_geno == "NN" and (c_geno[0] == f_geno[0] or c_geno[0] == f_geno[1]):
								parameter.person_dict[ID].haplotype[pos] = (c_geno[0], c_geno[1])
							else:
								#pass
								print "child geno not in parent geno", ID, pos, f_geno, m_geno, c_geno
								#sys.exit(1)
						"""
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
					try:
						#if len(c_set) != 1:
						#print f_geno, m_geno, c_geno, parameter.person_dict[ID].haplotype[pos][0], parameter.person_dict[ID].haplotype[pos][1]
						pass
					except:
						pass
						#print "error 4", pos, ID, f_geno, m_geno, c_geno
						#print "child geno not in parent geno", ID, pos, f_geno, m_geno, c_geno
						#sys.exit(1)

def output_child_hap():
	pos_list = parameter.rsID_dict.keys()
	pos_list.sort()

	with open(parameter.children_hap_file, "w") as c_hap_file:
		print >> c_hap_file, "rs_ID", "pos",
		for id in parameter.children_list:
			print >> c_hap_file, id + "_F", id + "_M",
		print >> c_hap_file, ""

		for pos in pos_list:
			print >> c_hap_file, parameter.rsID_dict[pos], pos,
			for child_id in parameter.children_list:
				print >> c_hap_file, parameter.person_dict[child_id].haplotype[pos][0], parameter.person_dict[child_id].haplotype[pos][1],
			print >> c_hap_file, ""

def prepare_id_list():
	# children list
	for ID in parameter.person_dict.keys():
		person = parameter.person_dict[ID]
		if len(person.children) > 0:
			parameter.children_list.extend(list(person.children.keys()))
			#print parameter.children_list

	parameter.children_list = list(set(parameter.children_list))
	parameter.children_list.sort()
	#print parameter.children_list

	# parent list
	for ID in parameter.person_dict.keys():
		person = parameter.person_dict[ID]
		if person.father != "N/A":
			parameter.father_list.append(person.father)
		if person.mather != "N/A":
			parameter.mather_list.append(person.mather)

	parameter.father_list = list(set(parameter.father_list))
	parameter.father_list.sort()
	#print parameter.father_list

	parameter.mather_list = list(set(parameter.mather_list))
	parameter.mather_list.sort()
	#print parameter.mather_list

	#parameter.pos_list = [int(x) for x in parameter.rsID_dict.keys()]
	parameter.pos_list = parameter.rsID_dict.keys()
	parameter.pos_list.sort()

def sort_fragment(p_id):

	children_list = parameter.person_dict[p_id].children.keys()
	children_list.sort()
	#print children_list

	fragment_startpos_dict = {}
	"""
	# keep for remove duplicates
	total_fragment = 0
	unique_fragment_total = 0

	for child_id in children_list:
		unique_fragment = []
		list = parameter.fragment_dict[child_id]
		for fragment in list:
			total_fragment += 1
			in_unique_list = False
			#print fragment.ID, fragment.start, fragment.end
			start_pos = fragment.start
			end_pos = fragment.end
			for uni_frangment in unique_fragment:
				if uni_frangment.start_pos == start_pos and uni_frangment.end_pos == end_pos:
					in_unique_list = True
					break

			if not in_unique_list:
				unique_fragment_total += 1
				if start_pos not in fragment_startpos_dict:
					fragment_startpos_dict[int(start_pos)] = []
				fragment_startpos_dict[int(start_pos)].append(fragment)

	print "total_fragment", total_fragment
	print "unique_fragment_total", unique_fragment_total

	fragment_startpos_sorted_list = sort_dict_by_key(fragment_startpos_dict)
	"""
	# To sort fragment_list
	for child_id in children_list:
		#list = parameter.person_dict[p_id].fragment_dict[child_id]
		list = parameter.fragment_dict[child_id]
		for fragment in list:
			#print fragment.ID, fragment.start, fragment.end
			start_pos = fragment.start
			end_pos = fragment.end
			if start_pos not in fragment_startpos_dict:
				fragment_startpos_dict[int(start_pos)] = []
			fragment_startpos_dict[int(start_pos)].append(fragment)

	fragment_startpos_sorted_list = sort_dict_by_key(fragment_startpos_dict)

	"""
	for data in fragment_startpos_sorted_list:
		start_pos = data[0]
		fragment_list = data[1]

		print start_pos,
		for fragment in fragment_list:
			print fragment.ID, fragment.start
	"""

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
	print "hetero pos_list size", len(parameter.person_dict[p_id].hetero_pos_list)

	compare_child_hap(p_id, p_code, children_list[0], children_list[1])
	compare_child_hap(p_id, p_code, children_list[0], children_list[2])
	compare_child_hap(p_id, p_code, children_list[1], children_list[2])

	temp_parent_hap_A = {}
	temp_parent_hap_B = {}

	fragment_startpos_sorted_list = sort_fragment(p_id)
	#group_parent_hap(p_id, fragment_startpos_sorted_list)

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

	#print "temp_parent_hap_A", len(temp_parent_hap_A)
	#print "temp_parent_hap_B", len(temp_parent_hap_B)


	count = 0
	for pos in parameter.pos_list:
		f_geno = parameter.person_dict[p_id].genotype_dict[pos]
		if f_geno[0] == f_geno[1]:
			#if pos not in temp_parent_hap_A:
			temp_parent_hap_A[pos] = f_geno[0]
			#if pos not in temp_parent_hap_B:
			temp_parent_hap_B[pos] = f_geno[0]
		else:
			if pos in temp_parent_hap_A and pos in temp_parent_hap_B:
				if temp_parent_hap_A[pos] == "X" and temp_parent_hap_B[pos] != "X":
					#print "xxxA000", temp_parent_hap_A[pos], temp_parent_hap_B[pos], f_geno[0], f_geno[1]
					temp_parent_hap_A[pos] = f_geno[0] if temp_parent_hap_B[pos] == f_geno[1] else f_geno[1]
					#print "xxxA111", temp_parent_hap_A[pos],temp_parent_hap_B[pos]
				if temp_parent_hap_B[pos] == "X" and temp_parent_hap_A[pos] != "X":
					temp_parent_hap_B[pos] = f_geno[0] if temp_parent_hap_A[pos] == f_geno[1] else f_geno[1]
					#print "xxx", temp_parent_hap_A[pos], temp_parent_hap_B[pos]

			elif pos in temp_parent_hap_A and pos not in temp_parent_hap_B:
				temp_parent_hap_B[pos] = f_geno[0] if temp_parent_hap_A[pos] == f_geno[1] else f_geno[1]

			elif pos not in temp_parent_hap_A and pos in temp_parent_hap_B:
				temp_parent_hap_A[pos] = f_geno[0] if temp_parent_hap_B[pos] == f_geno[1] else f_geno[1]

			elif pos not in temp_parent_hap_A and pos not in temp_parent_hap_B:
				temp_parent_hap_A[pos] = "N"
				temp_parent_hap_B[pos] = "N"
				count += 1


	#print "count N", count


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


	# to output data
	with open(p_id + "_hap.txt", "w") as x_file:
		print >> x_file, "rs#", "pos", parameter.person_dict[p_id].ID + "_A", parameter.person_dict[p_id].ID + "_B"
		for pos in parameter.pos_list:
			print >> x_file, parameter.rsID_dict[pos], pos,
			f_geno = parameter.person_dict[p_id].genotype_dict[pos]
			#if f_geno[0] != f_geno[1]:
			if True:
				if pos in temp_parent_hap_A:
					print >> x_file, temp_parent_hap_A[pos],
				else:
					print >> x_file, "N",
				if pos in temp_parent_hap_B:
					print >> x_file, temp_parent_hap_B[pos]
				else:
					print >> x_file, "N"

	#print "temp_parent_hap_A final", len(temp_parent_hap_A)
	#print "temp_parent_hap_B final", len(temp_parent_hap_B)


def output_hap_std():
	parent_id_list = []
	parent_id_list.extend(parameter.father_list)
	parent_id_list.extend(parameter.mather_list)
	parent_id_list.sort()
	for p_id in parent_id_list:
		with open(p_id + "_std.txt", "w") as std_file:
			print >> std_file, "rs#", "pos", parameter.person_dict[p_id].ID + "_A", parameter.person_dict[p_id].ID + "_B"
			for pos in parameter.pos_list:
				print >> std_file, parameter.rsID_dict[pos], pos,
				f_geno = parameter.person_dict[p_id].genotype_dict[pos]
				print >> std_file, f_geno[0], f_geno[1]


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
				print >> parent_hap, parameter.person_dict[id].haplotype[pos][0], parameter.person_dict[id].haplotype[pos][1],
			print >> parent_hap, ""


def group_parent_hap(p_id, fragment_startpos_sorted_list):

	group_hap_A = {}
	group_hap_B = {}
	fragment_group_A_list = []
	fragment_group_B_list = []
	"""
	children_list = parameter.person_dict[p_id].children.keys()
	children_list.sort()
	size = len(parameter.fragment_dict[children_list[0]] + parameter.fragment_dict[children_list[1]] + parameter.fragment_dict[children_list[2]])

	print "ori_f_list_size", size
	"""

	fragment_sorted_list = []
	#fragment_startpos_sorted_list = sort_fragment(p_id)
	status = ""
	for data in fragment_startpos_sorted_list:
		start_pos = data[0]
		fragment_list = data[1]
		fragment_sorted_list.extend(fragment_list)
	print "fragment_sorted_list", len(fragment_sorted_list)
	for fragment in fragment_sorted_list:
		#print fragment.ID, fragment.start, fragment.end, fragment.length
		pass


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

		#print "current_start, current_end, overlap_start, overlap_end", current_start, current_end, overlap_start, overlap_end

		same = 0
		not_same = 0
		#f_id = "1NA19702"
		for pos in parameter.person_dict[f_id].hetero_pos_list:
			if pos >= overlap_start and pos <= overlap_end:
				if parameter.person_dict[fragment.ID].haplotype[pos][p_code] != "X" and temp_parent_hap_A[pos] != "X" \
						and parameter.person_dict[fragment.ID].haplotype[pos][p_code] != "N" and temp_parent_hap_A[pos] != "N":
					#print "xxxx", pos, parameter.person_dict[fragment.ID].haplotype[pos][0], temp_parent_hap[pos]
					if parameter.person_dict[fragment.ID].haplotype[pos][p_code] == temp_parent_hap_A[pos]:
						same += 1
					else:
						not_same += 1

		#print "before", len(temp_parent_hap_A)
		percentage = float(same)/(same + not_same + 1)
		#print percentage
		if percentage >= 0.99:
			for pos in parameter.person_dict[f_id].hetero_pos_list:
				if pos >= fragment.start and pos <= fragment.end and pos not in temp_parent_hap_A:
					temp_parent_hap_A[pos] = parameter.person_dict[fragment.ID].haplotype[pos][p_code]
			#print "added snp", len(temp_parent_hap_A)
			#print "percentage A", percentage
			return "same_to_A", temp_parent_hap_A
		elif percentage <= 0.01:
			#print "percentage B", percentage
			return "same_to_B", temp_parent_hap_A
		else:
			#print fragment.ID, fragment.start, fragment.end, fragment.length, percentage
			#print "percentage =", percentage
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

			#parameter.person_dict[p_id].hetero_pos_list.append(pos)

			if cf_hap_1 != "X" and cf_hap_2 != "X" and cf_hap_1 != "N" and cf_hap_2 != "N":

				if compare == "":
					#print f_geno, cf_hap_1, cf_hap_2
					compare = "same" if cf_hap_1 == cf_hap_2 else "different"
					#print "initial compare status:", compare

				#print "all", pos, cf_hap_1, cf_hap_2
				elif compare == "same":
					if cf_hap_1 != cf_hap_2:
						child_fragment_dict[child_ID_1].append(pos)
						child_fragment_dict[child_ID_2].append(pos)
						compare = "different"
						#print "different",  child_ID_1, child_ID_2, pos, cf_hap_1, cf_hap_2
				elif compare == "different":
					if cf_hap_1 == cf_hap_2:
						child_fragment_dict[child_ID_1].append(pos)
						child_fragment_dict[child_ID_2].append(pos)
						compare = "same"
						#print "same", child_ID_1, child_ID_2, pos, cf_hap_1, cf_hap_2

	child_fragment_dict[child_ID_1].append(parameter.pos_list[-1])
	child_fragment_dict[child_ID_2].append(parameter.pos_list[-1])

	#print len(child_fragment_dict[child_ID_1])
	#print len(child_fragment_dict[child_ID_2])

	for id in child_ID_1, child_ID_2:
		for index, pos in enumerate(child_fragment_dict[id]):
			if index < len(child_fragment_dict[id]) - 1:
				fragment = fragments()
				fragment.ID = id
				fragment.start = int(child_fragment_dict[id][index])
				fragment.end = int(child_fragment_dict[id][index + 1])
				fragment.length = parameter.pos_list.index(fragment.end) - parameter.pos_list.index(fragment.start)
				#print fragment.ID, fragment.start, fragment.end, fragment.length, child_1.haplotype[pos][p_code], child_2.haplotype[pos][p_code]
				print child_ID_1, "vs", child_ID_2, fragment.start, fragment.end, fragment.length

				#if id not in parameter.person_dict[p_id].fragment_dict:
				#	parameter.person_dict[p_id].fragment_dict[id] = []
				#parameter.person_dict[p_id].fragment_dict[id].append(fragment)

				if id not in parameter.fragment_dict:
					parameter.fragment_dict[id] = []
				parameter.fragment_dict[id].append(fragment)

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
					parameter.person_dict[child_id_1].genotype_dict[pos],\
					parameter.person_dict[child_id_1].haplotype[pos][0], \
					parameter.person_dict[child_id_2].genotype_dict[pos],\
					parameter.person_dict[child_id_2].haplotype[pos][0], \
					parameter.person_dict[child_id_3].genotype_dict[pos],\
					parameter.person_dict[child_id_3].haplotype[pos][0]

	"""
	for child in parameter.person_dict[p_id].fragment_dict.keys():
		print child
		for fragment in parameter.person_dict[p_id].fragment_dict[child]:
			for pos in parameter.person_dict[p_id].hetero_pos_list:


			fragment.ID

	"""
	#f_geno = parameter.person_dict[p_id].genotype_dict[pos]
	#	if f_geno != "NN" and f_geno[0] != f_geno[1]:
	#		parameter.person_dict[p_id].hetero_pos_list.append(pos)


def compare_child_hap_old(p_id, p_code, child_ID_1, child_ID_2):
	"""keep it for now"""

	child_fragment_dict = {}
	child_fragment_dict[child_ID_1] = []
	child_fragment_dict[child_ID_2] = []

	fragment_list = []

	break_point_dict = {}
	break_point_list = []

	child_1 = parameter.person_dict[child_ID_1]
	child_2 = parameter.person_dict[child_ID_2]

	child_fragment_dict[child_ID_1].append(parameter.pos_list[0])
	child_fragment_dict[child_ID_2].append(parameter.pos_list[0])

	compare = ""

	for pos in parameter.pos_list:
		f_geno = parameter.person_dict[p_id].genotype_dict[pos]

		if f_geno != "NN" and f_geno[0] != f_geno[1]:

			cf_hap_1 = child_1.haplotype[pos][p_code]
			cf_hap_2 = child_2.haplotype[pos][p_code]

			parameter.person_dict[p_id].hetero_pos_list.append(pos)

			if cf_hap_1 != "X" and cf_hap_2 != "X" and cf_hap_1 != "N" and cf_hap_2 != "N":

				if compare == "":
					#print f_geno, cf_hap_1, cf_hap_2
					compare = "same" if cf_hap_1 == cf_hap_2 else "different"
					#print "initial compare status:", compare

				#print "all", pos, cf_hap_1, cf_hap_2
				elif compare == "same":
					if cf_hap_1 != cf_hap_2:
						child_fragment_dict[child_ID_1].append(pos)
						child_fragment_dict[child_ID_2].append(pos)
						compare = "different"
						#print "different", pos, cf_hap_1, cf_hap_2
				elif compare == "different":
					if cf_hap_1 == cf_hap_2:
						child_fragment_dict[child_ID_1].append(pos)
						child_fragment_dict[child_ID_2].append(pos)
						compare = "same"
						#print "same", pos, cf_hap_1, cf_hap_2

	child_fragment_dict[child_ID_1].append(parameter.pos_list[-1])
	child_fragment_dict[child_ID_2].append(parameter.pos_list[-1])

	#print child_fragment_dict[child_ID_1]
	#print child_fragment_dict[child_ID_2]

	for id in child_ID_1, child_ID_2:
		for index, pos in enumerate(child_fragment_dict[id]):
			if index < len(child_fragment_dict[id]) - 1:
				fragment = fragments()
				fragment.ID = id
				fragment.start = child_fragment_dict[id][index]
				fragment.end = child_fragment_dict[id][index + 1]
				fragment.length = parameter.pos_list.index(fragment.end) - parameter.pos_list.index(fragment.start)
				#print fragment.ID, fragment.start, fragment.end, fragment.length, child_1.haplotype[pos][p_code], child_2.haplotype[pos][p_code]
				print fragment.ID, fragment.start, fragment.end, fragment.length

				if id not in parameter.person_dict[p_id].fragment_dict:
					parameter.person_dict[p_id].fragment_dict[id] = []
				parameter.person_dict[p_id].fragment_dict[id].append(fragment)

def genome_laser(pedi_name, geno_name):

	pedi_name = "pedi.txt"
	geno_name = "geno.txt"

	load_pedi(pedi_name)
	load_geno(geno_name)

	prepare_id_list()
	parents_to_children()
	output_child_hap()

	output_hap_std()

	f_code = 0
	m_code = 1

	for f_id in parameter.father_list:
		children_to_parents(f_id, f_code)
		print f_id
		#hifiAccuCheck_file(f_id+"_hap.txt", f_id+"_std.txt")


	for m_id in parameter.mather_list:
		children_to_parents(m_id, m_code)
		#hifiAccuCheck_file(m_id+"_hap.txt", m_id+"_std.txt")
		pass

	output_parent_hap()

	#hifiAccuCheck_file("NA07347_hap.txt", "NA07347_std.txt")

	#output_fragment()

def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-p", "--pedi", type="string", dest="pedi_name", help="Input pedi name", default="null")
	parser.add_option("-g", "--geno", type="string", dest="geno_name", help="Input file name", default="null")
	(options, args) = parser.parse_args()

	return options

if __name__ == '__main__':
	options = get_args()
	pedi_name = options.pedi_name
	geno_name = options.geno_name
	global person_dict
	person_dict = {}

	global parameter
	parameter = parameters()

	start_time = time.time()

	genome_laser(pedi_name, geno_name)

	print "elapsed_time is: ", round(time.time() - start_time, 2), "s"