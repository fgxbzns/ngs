#!/usr/bin/python
# ######################################################################################
# Guoxing Fu Jan 28, 2015
# To generate simulation data for laser genome
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser

class parameters:
	def __init__(self):
		self.person_dict = {}
		self.personID_list = []
		self.rsID_dict = {}
		self.pos_list = []
		self.father_list = []
		self.mather_list = []
		self.children_list = []
		self.children_dict = {}

		self.fragment_dict = {}
		self.fragment_list = []

		self.pedi_list = []

class pedis:
	def __init__(self):
		self.fater = ""
		self.mater = ""
		self.children_list = []

class persons:
	def __init__(self):
		self.ID = ""
		self.A_dict = {}
		self.B_dict = {}

def list_to_line(list):
	line = ""
	for a in list:
		line += str(a).strip() + "\t"
	return line.strip()

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
					# convert the position to int for sorting
					if raw_data_format == "list":
						data[int(elements[1])] = elements
					elif raw_data_format == "string":
						data[int(elements[1])] = line.strip()
				except:
					#print "error in ", line, file_name
					pass
	return (title_info, data)

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

def load_ref(ref_file_name):
	ref_title_info, hap_ref_dict = load_raw_data(ref_file_name)
	print "total_ref_number pos: ", len(hap_ref_dict)

	for index, ID in enumerate(ref_title_info.split()[2:]):
		if index % 2 == 0:
			person_ID = ID[:-2]
			#print rsID
			person = persons()
			person.ID = person_ID
			person.A_dict = {}
			person.B_dict = {}
			parameter.person_dict[person_ID] = person
			parameter.personID_list.append(person_ID)
	#print parameter.personID_list

	for pos in hap_ref_dict.keys():
		parameter.pos_list.append(pos)
		elements = hap_ref_dict[pos]
		rsID = elements[0]
		parameter.rsID_dict[pos] = rsID
		allele_list = elements[2:]
		for index, personID in enumerate(parameter.personID_list):
			parameter.person_dict[personID].A_dict[pos] = allele_list[index * 2]
			parameter.person_dict[personID].B_dict[pos] = allele_list[index * 2 + 1]

	parameter.pos_list.sort()

def get_crossover_pos():
	crossover_times = 3
	crossover_pos_list = []

	temp_pos_list = parameter.pos_list[:]
	last_pos = parameter.pos_list[-1]
	#print last_pos
	crossover_pos_list.append(temp_pos_list[0])
	del temp_pos_list[0]
	del temp_pos_list[-1]
	#print len(temp_pos_list)
	for i in range(crossover_times):
		crossover_pos_index = random.randrange(0, len(temp_pos_list))
		crossover_pos = temp_pos_list[crossover_pos_index]
		del temp_pos_list[crossover_pos_index]
		crossover_pos_list.append(crossover_pos)
	crossover_pos_list.append(last_pos)
	crossover_pos_list.sort()
	print crossover_pos_list
	return crossover_pos_list

def get_hap(person, crossover_pos_list):
	hap_dict = {}
	for pos in parameter.pos_list:
		if pos >= crossover_pos_list[0] and pos < crossover_pos_list[1]:
			hap_dict[pos] = person.A_dict[pos]
		elif pos >= crossover_pos_list[1] and pos < crossover_pos_list[2]:
			hap_dict[pos] = person.B_dict[pos]
		elif pos >= crossover_pos_list[2] and pos < crossover_pos_list[3]:
			hap_dict[pos] = person.A_dict[pos]
		elif pos >= crossover_pos_list[3] and pos <= crossover_pos_list[4]:
			hap_dict[pos] = person.B_dict[pos]
		#elif pos >= crossover_pos_list[4] and pos <= crossover_pos_list[5]:
		#	hap_dict[pos] = person.A_dict[pos]
		else:
			#print "error"
			pass
	return hap_dict

def get_hap_various_crossover(person, crossover_pos_list):
	hap_dict = {}
	crossover_pos_index = 0
	hap_side = "A"

	for pos in parameter.pos_list:
		if crossover_pos_index < len(crossover_pos_list):
			if pos >= crossover_pos_list[crossover_pos_index] and pos <= crossover_pos_list[crossover_pos_index+1]:
				if hap_side == "A":
					hap_dict[pos] = person.A_dict[pos]
				elif hap_side == "B":
					hap_dict[pos] = person.B_dict[pos]
			else:
				if hap_side == "A":
					hap_dict[pos] = person.A_dict[pos]
				elif hap_side == "B":
					hap_dict[pos] = person.B_dict[pos]
				hap_side = "A" if hap_side == "B" else hap_side == "A"
				crossover_pos_index += 1

	return hap_dict


def generate_pedi():

	f_index = random.randrange(0, len(parameter.personID_list))
	father = parameter.person_dict[parameter.personID_list[f_index]]
	del parameter.personID_list[f_index]

	m_index = random.randrange(0, len(parameter.personID_list))
	mather = parameter.person_dict[parameter.personID_list[m_index]]
	del parameter.personID_list[m_index]

	pedi = pedis()
	pedi.father = father.ID
	pedi.mather = mather.ID

	for i in range(1, 4):

		person = persons()
		person.ID = father.ID + "_" + mather.ID + "_" + str(i)
		crossover_pos_list = get_crossover_pos()
		person.A_dict = get_hap(father, crossover_pos_list)
		crossover_pos_list = get_crossover_pos()
		person.B_dict = get_hap(mather, crossover_pos_list)
		pedi.children_list.append(person)

	parameter.pedi_list.append(pedi)

def output_pedi(pedi_num):

	for i in range(pedi_num):
		generate_pedi()

	with open("pedi.txt", "w") as pedi_output:
		print >> pedi_output, "PersonID", "Father", "Mather"
		for pedi in parameter.pedi_list:
			print >> pedi_output, pedi.father, "N/A", "N/A"
			print >> pedi_output, pedi.mather, "N/A", "N/A"
			for child in pedi.children_list:
				print >> pedi_output, child.ID, pedi.father, pedi.mather

	with open("geno.txt", "w") as geno_output:
		print >> geno_output, "rs#", "pos",
		for pedi in parameter.pedi_list:
			print >> geno_output, pedi.father, pedi.mather,
			for child in pedi.children_list:
				print >> geno_output, child.ID,
		print >> geno_output, ""

		for pos in parameter.pos_list:
			print >> geno_output, parameter.rsID_dict[pos], pos,
			for pedi in parameter.pedi_list:
				print >> geno_output, parameter.person_dict[pedi.father].A_dict[pos] \
					+ parameter.person_dict[pedi.father].B_dict[pos] \
					, parameter.person_dict[pedi.mather].A_dict[pos] + parameter.person_dict[pedi.mather].B_dict[pos],
				for child in pedi.children_list:
					print >> geno_output, child.A_dict[pos] + child.B_dict[pos],
			print >> geno_output, ""

def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-r", "--ref", type="string", dest="ref_file", help="Input file name", default="null")
	parser.add_option("-n", "--num", type="int", dest="pedi_num", help="Input file name", default="1")

	(options, args) = parser.parse_args()
	return options

if __name__ == '__main__':
	options = get_args()
	global person_dict
	person_dict = {}

	global parameter
	parameter = parameters()
	ref_file_name = options.ref_file
	pedi_num = options.pedi_num
	ref_file_name = "hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_ceu.phased"

	start_time = time.time()
	load_ref(ref_file_name)
	output_pedi(pedi_num)

	print "elapsed_time is: ", round(time.time() - start_time, 2), "s"


