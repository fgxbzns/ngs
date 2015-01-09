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
		self.father_list = []
		self.mather_list = []
		self.children_list = []
		self.children_dict = {}

class persons:
	def __init__(self):
		self.ID = ""
		self.father = ""
		self.mather = ""
		self.children = {}
		self.genotype_dict = {}
		self.haplotype = {}

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
					position = elements[1]
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

	for ID in parameter.person_dict.keys():
		person = parameter.person_dict[ID]
		if len(person.children) > 0:
			parameter.children_list.extend(list(person.children.keys()))
			#print parameter.children_list

	parameter.children_list = list(set(parameter.children_list))
	parameter.children_list.sort()
	print parameter.children_list

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
		print person.ID, person.father, person.mather,
		print person.children.keys()
	#print person_dict["1NAC1002"].genotype_dict['9935312']
	#print person_dict["1NAC1002"].genotype_dict['10014103']
	#print person_dict["1NAC1002"].genotype_dict['10065514']

	#parents_to_children()

	output_child_hap()

	print "elapsed_time is: ", round(time.time() - start_time, 2), "s"