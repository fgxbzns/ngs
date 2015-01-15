#!/usr/bin/python
#######################################################################################
# Guoxing Fu Jan 28, 2015
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

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

class fragments:
	def __init__(self):
		self.ID = ""
		self.start = ""
		self.end = ""
		self.length = ""

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
					print "error in ", line, pedi_name
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
								print "2", pos

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
						#pass
						#print "error 4", pos, ID, f_geno, m_geno, c_geno
						print "child geno not in parent geno", ID, pos, f_geno, m_geno, c_geno
						#sys.exit(1)

def output_child_hap():
	"""
	for ID in parameter.person_dict.keys():
		person = parameter.person_dict[ID]
		if len(person.children) > 0:
			parameter.children_list.extend(list(person.children.keys()))
			#print parameter.children_list

	parameter.children_list = list(set(parameter.children_list))
	parameter.children_list.sort()
	print parameter.children_list
	"""

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

def children_to_parents_old():

	f_id = "10NA19836"

	children_list = parameter.person_dict[f_id].children.keys()
	children_list.sort()
	print children_list

	#pos_list = parameter.rsID_dict.keys()
	#pos_list.sort()

	compare = "same"
	fragment_list = []
	temp_fragment_A = {}
	temp_fragment_B = {}

	break_point_list = []

	for pos in parameter.pos_list:

		cf_hap_1 = parameter.person_dict[children_list[0]].haplotype[pos][0]
		cf_hap_2 = parameter.person_dict[children_list[1]].haplotype[pos][0]
		cf_hap_3 = parameter.person_dict[children_list[2]].haplotype[pos][0]

		if (cf_hap_1 != cf_hap_2 or cf_hap_2 != cf_hap_3) and cf_hap_1 != "X" and cf_hap_2 != "X" and cf_hap_3 != "X":
			#print pos, cf_hap_1, cf_hap_2, cf_hap_3, parameter.person_dict[f_id].genotype_dict[pos]
			pass

		if cf_hap_1 != cf_hap_2 and cf_hap_1 != "X" and cf_hap_1 != "N" and cf_hap_2 != "X" and cf_hap_2 != "N":
			#print pos
			pass

		if cf_hap_1 != cf_hap_3 and cf_hap_1 != "X" and cf_hap_1 != "N" and cf_hap_3 != "X" and cf_hap_3 != "N":
			#print pos
			pass

		if cf_hap_1 != "X" and cf_hap_2 != "X":
			if compare == "same":
				if cf_hap_1 == cf_hap_2:
					temp_fragment_A[pos] = cf_hap_1
				else:
					fragment_list.append(temp_fragment_A)
					temp_fragment_A = {}
					compare = "different"


		"""
		try:
			print parameter.person_dict[children_list[0]].haplotype[pos][0], \
				parameter.person_dict[children_list[1]].haplotype[pos][0], \
				parameter.person_dict[children_list[2]].haplotype[pos][0]
		except:
			print pos, children_list
		"""


def children_to_parents():
	f_id = "10NA19836"

	children_list = parameter.person_dict[f_id].children.keys()
	children_list.sort()
	print children_list

	for pos in parameter.pos_list:
		f_geno = parameter.person_dict[f_id].genotype_dict[pos]
		if f_geno != "NN" and f_geno[0] != f_geno[1]:
			parameter.person_dict[f_id].hetero_pos_list.append(pos)

	print "hetero pos_list size", len(parameter.person_dict[f_id].hetero_pos_list)

	#for child_id in children_list:

	compare_child_hap(children_list[0], children_list[1])
	compare_child_hap(children_list[0], children_list[2])
	compare_child_hap(children_list[1], children_list[2])

	max_length_fragment = ""
	max_length = 0

	for child_id in children_list:
		list = parameter.fragment_dict[child_id]
		for fragment in list:
			if max_length == 0:
				max_length = fragment.length
				max_length_fragment = fragment
			else:
				if fragment.length > max_length:
					max_length = fragment.length
					max_length_fragment = fragment
			#print fragment.ID, fragment.start, fragment.end, fragment.length
	#print max_length, max_length_fragment.ID

	temp_parent_hap_A = {}
	temp_parent_hap_B = {}
	temp_parent_hap_list = [temp_parent_hap_A, temp_parent_hap_B]

	for pos in parameter.pos_list:
		if pos >= max_length_fragment.start and pos < max_length_fragment.end:
			temp_parent_hap_A[pos] = parameter.person_dict[max_length_fragment.ID].haplotype[pos][0]
			#print pos, temp_parent_hap_A[pos]

	child_id = children_list[1]
	if True:
	#for child_id in children_list:
		if child_id != max_length_fragment.ID:
			list = parameter.fragment_dict[child_id]
			fragment = list[0]
			if True:
			#for fragment in list:
				#current_start, current_end = get_parent_hap_end(temp_parent_hap_A)
				update_parent_hap(fragment, temp_parent_hap_A)
				pass

def update_parent_hap(fragment, temp_parent_hap):
	current_start, current_end = get_parent_hap_end(temp_parent_hap)

	overlap_start = 0
	overlap_end = 0

	if fragment.end <= current_start or fragment.start >= current_end:
		return "no_overlap", temp_parent_hap
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

	print current_start, current_end, overlap_start, overlap_end

	same = 0
	not_same = 0
	f_id = "10NA19836"
	for pos in parameter.person_dict[f_id].hetero_pos_list:
		if pos >= overlap_start and pos <= overlap_end:
			if parameter.person_dict[fragment.ID].haplotype[pos][0] != "X" and temp_parent_hap[pos] != "X" \
					and parameter.person_dict[fragment.ID].haplotype[pos][0] != "N" and temp_parent_hap[pos] != "N":
				#print "xxxx", pos, parameter.person_dict[fragment.ID].haplotype[pos][0], temp_parent_hap[pos]
				if parameter.person_dict[fragment.ID].haplotype[pos][0] == temp_parent_hap[pos]:
					same += 1
				else:
					not_same += 1
	print "dddddddddddd", same, not_same

def get_parent_hap_end(temp_parent_hap):
	pos_list = temp_parent_hap.keys()
	pos_list.sort()
	#print pos_list[0], pos_list[-1]
	return(pos_list[0], pos_list[-1])

def compare_child_hap(child_ID_1, child_ID_2):

	f_id = "10NA19836"

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
		f_geno = parameter.person_dict[f_id].genotype_dict[pos]

		if f_geno != "NN" and f_geno[0] != f_geno[1]:

			cf_hap_1 = child_1.haplotype[pos][0]
			cf_hap_2 = child_2.haplotype[pos][0]

			parameter.person_dict[f_id].hetero_pos_list.append(pos)

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

	print child_fragment_dict[child_ID_1]
	print child_fragment_dict[child_ID_2]

	for id in child_ID_1, child_ID_2:
		for index, pos in enumerate(child_fragment_dict[id]):
			if index < len(child_fragment_dict[id]) - 1:
				fragment = fragments()
				fragment.ID = id
				fragment.start = child_fragment_dict[id][index]
				fragment.end = child_fragment_dict[id][index + 1]
				fragment.length = parameter.pos_list.index(fragment.end) - parameter.pos_list.index(fragment.start)
				#print fragment.ID, fragment.start, fragment.end, fragment.length
				if id not in parameter.fragment_dict:
					parameter.fragment_dict[id] = []
				parameter.fragment_dict[id].append(fragment)












	pass

def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-p", "--pedi", type="string", dest="pedi_name", help="Input pedi name", default="null")
	parser.add_option("-g", "--geno", type="string", dest="geno_file", help="Input file name", default="null")
	(options, args) = parser.parse_args()
	"""
	if options.hg18_name == "null" or options.hg19_name == "null" or options.del_name == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	"""
	return options

if __name__ == '__main__':
	options = get_args()
	pedi_name = options.pedi_name
	global person_dict
	person_dict = {}

	global parameter
	parameter = parameters()

	start_time = time.time()
	pedi_name = "pedi.txt"
	geno_name = "geno.txt"

	load_pedi(pedi_name)

	load_geno(geno_name)
	print len(person_dict)
	for ID in parameter.person_dict:
		person = parameter.person_dict[ID]
		#print person.ID, person.father, person.mather,
		#print person.children.keys()
	#print person_dict["1NAC1002"].genotype_dict['9935312']
	#print person_dict["1NAC1002"].genotype_dict['10014103']
	#print person_dict["1NAC1002"].genotype_dict['10065514']


	prepare_id_list()
	parents_to_children()

	#output_child_hap()

	#children_to_parents()
	#compare_child_hap()

	children_to_parents()

	print "elapsed_time is: ", round(time.time() - start_time, 2), "s"