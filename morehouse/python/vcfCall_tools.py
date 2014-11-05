#!/usr/bin/python

# for zebra fish data
# Latest version of snpPick, use temp dict to store pos

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser
from tools import *
import sqlite as sql
import sqlite3 as lite

class parameter:
	def __init__(self):
		self.sam_file = ""
		self.sam_file_name = ""
		self.chr_name = ""
		self.snp_in_mimi = 0
		self.total_mimi = 0
		self.ref_file = ""
		self.temp_data_dict = {}
		self.temp_data_dict_counter = 0
		self.temp_data_dict_max_limit = 5000
		self.temp_data_dict_min_limit = 1000
		self.total_called_pos = 0
		self.output_file = ""

		self.db_name = ""
		self.db_base_name = ""
		self.ref_file = ""
		self.second_largest_allele_depth_cutoff = 2
		self.quality_score_threshold = 30

		self.indel_dict = {}
		self.ref_chr_seq = ""

class position_data:
	def __init__(self):
		self.pos = 0
		self.chr_name = ""
		self.ref_allele = ""
		self.A_depth = 0
		self.T_depth = 0
		self.C_depth = 0
		self.G_depth = 0

		self.A_qs = ""
		self.T_qs = ""
		self.C_qs = ""
		self.G_qs = ""

class indel:
	def __init__(self):
		self.pos = 0
		self.chr_name = ""
		self.ref_allele = ""
		self.deletion = {}
		self.insertion = {}

def get_ref_geno(parameters):
	chr_seq = ""
	input_file = open(parameters.ref_file, "r")
	for lines in input_file:
		if not lines.startswith(">"):
			chr_seq += lines.strip()
	print "total base number: ", len(chr_seq)
	return chr_seq

def variant_call_pair_end_sql(sam_file):
	"""
	the sequence pair has already been processed
	now treat the read as single end
	using sqlite
	"""

	total_reads_number = wccount(sam_file)
	percentage_of_total_file = 0

	chr_seq = get_ref_geno(chr_name)

	global table_name
	con = lite.connect(db_name)
	with con:
		cur = con.cursor()

		inputfile_sam = open(currentPath + sam_file, "r")
		sam_line_first = inputfile_sam.readline()  # the first read line in a pair
		total_reads_num = 0
		covered_snp_total_number = 0

		insert_size_lower_bond = 0
		insert_size_upper_bond = 1000

		while sam_line_first != '':
			if not sam_line_first.startswith("@"):
				current_percent = int(float(total_reads_number * percentage_of_total_file) / 100)
				if total_reads_num == current_percent:
					print "current progress: ", percentage_of_total_file
					percentage_of_total_file += 10

				total_reads_num += 1
				elements_first = sam_line_first.strip().split()
				try:
					read_ID_first = elements_first[0].strip()
					chrName_first = elements_first[2].strip()
					insert_size_first = abs(int(elements_first[8].strip()))  # insert_size for second read is negative
				except:
					print "error in first read:", sam_line_first
				if (insert_size_first >= insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):
					if True:  # this is for pair-end verification. Not needed at this moment.
						if chrName_first.startswith(chr_name):
							# if chrName_first == chr_name:
							# first read
							qName_first = elements_first[0].strip()
							flag_first = elements_first[1].strip()
							start_position_first = int(elements_first[3].strip())
							read_sequence_first = elements_first[9].strip()
							read_length_first = len(read_sequence_first)
							quality_score_sequence_first = elements_first[10].strip()

							for i in range(read_length_first):
								current_base_position = start_position_first + i
								A_depth = 0
								T_depth = 0
								C_depth = 0
								G_depth = 0

								covered_snp = read_sequence_first[i]  # ith position is the covered snp
								quality_score_symbol = quality_score_sequence_first[i]
								# update each read, or save certain amount of reads in an array and update together.
								# can be improved if needed.
								if (not covered_snp == 'N') and (
											(ord(quality_score_symbol) - 33) > parameters.quality_score_threshold):  # check quality_score
									#print ord(quality_score_symbol) - 33
									if covered_snp == "A":
										A_depth += 1
									elif covered_snp == "T":
										T_depth += 1
									elif covered_snp == "C":
										C_depth += 1
									elif covered_snp == "G":
										G_depth += 1

									cur.execute("SELECT *  from " + table_name + " where position=" + str(
										current_base_position))
									row = cur.fetchone()
									if row == None:

										try:
											inset_querry = "INSERT INTO " + table_name + \
											               " (position, chr, ref_allele, A_depth, T_depth, C_depth, G_depth ) VALUES (" + \
											               str(current_base_position) + \
											               ",'" + chrName_first + "','" + chr_seq[
												               current_base_position - 1] + "'," + str(
												A_depth) + "," + str(
												T_depth) \
											               + "," + str(C_depth) + "," + str(G_depth) + ")"

											cur.execute(inset_querry)
										except:
											print sam_line_first
											print parameters.sam_file
											print inset_querry

									else:
										A_depth += int(row[3])
										T_depth += int(row[4])
										C_depth += int(row[5])
										G_depth += int(row[6])
										update_querry = "UPDATE " + table_name + " set A_depth=" + str(A_depth) + \
										                ", T_depth=" + str(T_depth) + ", C_depth=" + str(
											C_depth) + ", G_depth=" + \
										                str(G_depth) + " where position=" + str(current_base_position)
										#print update_querry
										cur.execute(update_querry)
					else:
						print "first and second read ID do not match", read_ID_first
			sam_line_first = inputfile_sam.readline()
		inputfile_sam.close()
	return total_reads_num


def snpPick(parameters):
	sam_file_name = parameters.sam_file[:(len(parameters.sam_file) - 4)]
	print "sam file: ", parameters.sam_file_name
	total_reads_num = variant_call_pair_end(parameters)
	print "total_reads_num", total_reads_num
	print "total called pos", parameters.total_called_pos


def get_data(db_name, table_name, start_line, end_line):
	con = lite.connect(db_name)
	with con:
		cur = con.cursor()
		querry = "SELECT * FROM " + table_name + " where position>" + start_line + " and position<" + end_line
		# print querry
		cur.execute(querry)
		rows = [[str(item) for item in results] for results in cur.fetchall()]
		#print len(rows)
		return rows

def output_data(file_name, start_line, end_line):
	total_row_number = int(end_line) - int(start_line)
	# print total_row_number
	if total_row_number <= 0:
		print "error in start point and end point"
		sys.exit(0)
	percentage_of_total = 0
	current_row = 0
	#output_file = open(currentPath + file_name, "w")
	with open(currentPath + file_name, "w") as output_file:
		#output_file.write('pos\tchr\tref_allele\t A \t T \t C \t G \n')
		print >> output_file, "pos", "chr", "ref_allele", "A", "T", "C", "G"
		rows = get_data(db_name, table_name, start_line, end_line)
		for item in rows:

			current_percent = int(float(total_row_number * percentage_of_total) / 100)
			#print current_percent
			if current_row == current_percent:
				print "current progress: ", percentage_of_total, "% current row:", current_row + int(start_line)
				percentage_of_total += 10
			current_row += 1
			#output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (item[0], item[1], item[2], item[3], item[4], item[5], item[6]))
			print >> output_file, item[0], item[1], item[2], item[3], item[4], item[5], item[6]
	print "total snp number :", current_row


def get_single_pos_data(db_name, table_name, pos):
	con = lite.connect(db_name)
	with con:
		cur = con.cursor()
		querry = "SELECT * FROM " + table_name + " where position=" + pos
		cur.execute(querry)
		rows = [[str(item) for item in results] for results in cur.fetchall()]
		return rows


def output_single_pos_data(pos_file_name):
	percentage_of_total = 0
	current_row = 0
	pos_list = []
	with open(currentPath + pos_file_name, "r") as pos_file:
		for line in pos_file:
			pos = line.strip()
			if pos != "":
				pos_list.append(pos)
	print "total pos number:", len(pos_list)

	file_name = pos_file_name + "_poslist_output"
	with open(currentPath + file_name, "w") as output_file:
		print >> output_file, "pos", "chr", "ref_allele", "A", "T", "C", "G"
		for pos in pos_list:
			rows = get_single_pos_data(db_name, table_name, pos)
			if len(rows) != 0:
				for item in rows:
					print >> output_file, item[0], item[1], item[2], item[3], item[4], item[5], item[6]
			else:
				print >> output_file, pos


def output_data_filter(file_name, start_line, end_line):
	# output data with small total number, divide the total number into small pieces
	total_row_number = int(end_line) - int(start_line)
	if total_row_number <= 0:
		print "error in start point and end point"
		sys.exit(0)
	percentage_of_total = 0
	current_row = 0
	with open(currentPath + file_name, "w") as output_file:
		print >> output_file, "pos", "chr", "ref_allele", "A", "T", "C", "G"
		rows = get_data(db_name, table_name, start_line, end_line)
		for item in rows:
			temp_list = [int(x) for x in item[3:7]]
			if temp_list.count(0) < 3:
				temp_list.sort()
				current_percent = int(float(total_row_number * percentage_of_total) / 100)
				if current_row == current_percent:
					print "current progress: ", percentage_of_total, "current row:", current_row + int(start_line)
					percentage_of_total += 10
				current_row += 1
				# pos, chr, ref_base, A, T, C, G, Allele#, Major, Minor, 3rd, 4th, rs#, SNP_allele_1, SNP_allele_2, mm_snp_overlapping
				print >> output_file, item[0], item[1], item[2], item[3], item[4], item[5], item[6], temp_list[3], \
					temp_list[2], temp_list[1], temp_list[0]
	print "total snp number :", current_row


def get_snp_info(pos, first_allele_number, second_allele_number, allele_list):
	base_list = ["A", "T", "C", "G"]
	pos = int(pos)
	if pos in hap_std_dict:
		parameters.snp_in_mimi += 1
		rs_number = hap_std_dict[pos][0]
		SNP_allele_A = hap_std_dict[pos][2]
		SNP_allele_B = hap_std_dict[pos][3]
		# 1 = yes, 2 = no
		mm_snp_overlap = 2

		first_allele = base_list[allele_list.index(first_allele_number)]
		second_allele = base_list[allele_list.index(second_allele_number)]
		# print pos, first_allele

		if first_allele == SNP_allele_A and second_allele == SNP_allele_B:
			mm_snp_overlap = 1
		elif first_allele == SNP_allele_B and second_allele == SNP_allele_A:
			mm_snp_overlap = 1
		return rs_number, SNP_allele_A, SNP_allele_B, mm_snp_overlap
	else:
		return "", "", "", ""


def data_filter(start_line, end_line):
	# prepare data portion for output_filtered_data
	start_time = time.time()
	total_row_number = int(end_line) - int(start_line)
	data_list = []
	percentage_of_total = 0
	current_row = 0
	rows = get_data(db_name, table_name, str(start_line), str(end_line))
	for item in rows:
		temp_list = [int(x) for x in item[3:7]]
		allele_list = [int(x) for x in item[3:7]]
		# if temp_list.count(0) < 3: # remove homo position
		number_of_zero = temp_list.count(0)
		if number_of_zero == 2:
			# Remove the postions with 3 zeros and filter by the second_largest_allele_depth
			temp_list.sort()
			second_largest_allele_depth = temp_list[2]
			if second_largest_allele_depth >= second_largest_allele_depth_cutoff:
				current_percent = int(float(total_row_number * percentage_of_total) / 100)
				if current_row == current_percent:
					print "current progress: ", percentage_of_total, "current row:", current_row + int(start_line)
					percentage_of_total += 10
				current_row += 1
				"""
				# This part is for solid data, to compare mimi with snp
				#snp_info = get_snp_info(pos, first_allele, second_allele, hap_std_dict)
				snp_info = get_snp_info(item[0], temp_list[3], temp_list[2], allele_list)
				# pos, chr, ref_base, A, T, C, G, Allele#, 1st, 2nd, 3rd, 4th, rs#, SNP_allele_1, SNP_allele_2, mm_snp_overlapping
				data_list.append(list((
					item[0], item[1], item[2], item[3], item[4], item[5], item[6], (4-number_of_zero), temp_list[3], temp_list[2], temp_list[1],
					temp_list[0], snp_info[0], snp_info[1], snp_info[2], snp_info[3])))
				"""
				data_list.append(list((
					item[0], item[1], item[2], item[3], item[4], item[5], item[6], (4 - number_of_zero), temp_list[3],
					temp_list[2], temp_list[1],
					temp_list[0])))
	elapse_time = time.time() - start_time
	print "running time: ", round(elapse_time, 3), "s"
	return data_list


def output_filtered_data(start_line, end_line):
	file_name = db_base_name + "_" + start_line + "_" + end_line + "_filtered_2nddepth_" + str(
		second_largest_allele_depth_cutoff) + ".txt"

	start_line = int(start_line)
	end_line = int(end_line)
	total_row_number = end_line - start_line
	with open(file_name, "w") as output_file:
		print >> output_file, "pos", "chr", "ref_allele", "A", "T", "C", "G", \
			"Allele#", "1st", "2nd", "3rd", "4th", "rs#", "SNP_allele_1", "SNP_allele_2", "mm_snp_overlapping"
		if total_row_number <= 0:
			print "error in start point and end point"
			sys.exit(0)
		elif total_row_number >= 100000:
			number_of_subfile = 100
			total_number_ceilling = int(math.ceil(float(total_row_number) / 100) * 100)
			print "total_number_ceilling: ", total_number_ceilling
			num_in_each_file = total_number_ceilling / number_of_subfile
			print "total_number_ceilling in each: ", num_in_each_file
			# seed_removed_in_last_subfile = int(math.fmod(len(seed_hetero_sorted_list), seed_removed_in_each_subfile))
			for i in range(number_of_subfile):
				if i != number_of_subfile - 1:
					print "processing ", i, start_line, start_line + num_in_each_file - 1
					data_list = data_filter(start_line, start_line + num_in_each_file - 1)
					print "number of alleles in this range: ", len(data_list)
					parameters.total_mimi += len(data_list)
					for data in data_list:
						print >> output_file, " ".join(str(x) for x in data)
				else:
					print "processing ", i, start_line, end_line
					data_list = data_filter(start_line, end_line)
					print "number of alleles in this range: ", len(data_list)
					parameters.total_mimi += len(data_list)
					for data in data_list:
						print >> output_file, " ".join(str(x) for x in data)
				start_line = start_line + num_in_each_file
		else:
			print "processing ", start_line, end_line
			data_list = data_filter(start_line, end_line)
			for data in data_list:
				print >> output_file, " ".join(str(x) for x in data)

def load_mimi_data(file_name):
	data = {}
	with open(file_name, "r") as fp:
		for line in fp:
			elements = line.strip().split()
			try:
				data[int(elements[0])] = elements[:7]
			except:
				# print "error in ", line, file_name
				pass
	return data

"""
Following is function for temp_dcit
"""


def variant_call_pair_end_old(parameters):
	"""
	the sequence pair has already been processed
	now treat the read as single end
	using temp_dict
	"""

	total_reads_number = wccount(parameters.sam_file)
	percentage_of_total_file = 0

	chr_seq = get_ref_geno(parameters)

	global table_name
	if True:

		inputfile_sam = open(parameters.sam_file, "r")
		sam_line_first = inputfile_sam.readline()  # the first read line in a pair
		total_reads_num = 0

		insert_size_lower_bond = 0
		insert_size_upper_bond = 1000

		while sam_line_first != '':
			if not sam_line_first.startswith("@"):
				current_percent = int(float(total_reads_number * percentage_of_total_file) / 100)
				if total_reads_num == current_percent:
					print "current progress: ", percentage_of_total_file
					percentage_of_total_file += 10

				total_reads_num += 1
				elements_first = sam_line_first.strip().split()
				try:
					read_ID_first = elements_first[0].strip()
					chrName_first = elements_first[2].strip()
					insert_size_first = abs(int(elements_first[8].strip()))  # insert_size for second read is negative
				except:
					print "error in first read:", sam_line_first
				if (insert_size_first >= insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):
					if True:  # this is for pair-end verification. Not needed at this moment.
						if chrName_first.startswith(parameters.chr_name):
							# if chrName_first == parameters.chr_name:
							# first read
							qName_first = elements_first[0].strip()
							flag_first = elements_first[1].strip()
							start_position_first = int(elements_first[3].strip())
							read_sequence_first = elements_first[9].strip()
							read_length_first = len(read_sequence_first)
							quality_score_sequence_first = elements_first[10].strip()

							for i in range(read_length_first):
								current_base_position = start_position_first + i

								covered_snp = read_sequence_first[i]  # ith position is the covered snp
								quality_score_symbol = quality_score_sequence_first[i]
								if (covered_snp != 'N') and (
											(ord(quality_score_symbol) - 33) > parameters.quality_score_threshold):  # check quality_score
									#print ord(quality_score_symbol) - 33
									if current_base_position not in parameters.temp_data_dict:
										temp_pos = position_data()
										temp_pos.pos = current_base_position
										temp_pos.chr_name = parameters.chr_name
										temp_pos.ref_allele = chr_seq[current_base_position - 1]
										parameters.temp_data_dict[current_base_position] = temp_pos
										parameters.temp_data_dict_counter += 1
										parameters.total_called_pos += 1

									if covered_snp == "A":
										parameters.temp_data_dict[current_base_position].A_depth += 1
									elif covered_snp == "T":
										parameters.temp_data_dict[current_base_position].T_depth += 1
									elif covered_snp == "C":
										parameters.temp_data_dict[current_base_position].C_depth += 1
									elif covered_snp == "G":
										parameters.temp_data_dict[current_base_position].G_depth += 1

								# output the data when it reaches the size limit
								if parameters.temp_data_dict_counter > parameters.temp_data_dict_max_limit:
									output_temp_dict(parameters)

					else:
						print "first and second read ID do not match", read_ID_first
			sam_line_first = inputfile_sam.readline()
		inputfile_sam.close()
	if len(parameters.temp_data_dict) > 0:
		#print "last temp dict size", len(parameters.temp_data_dict)
		output_temp_dict(parameters)
	return total_reads_num

def variant_call_pair_end(parameters):
	"""
	the sequence pair has already been processed
	now treat the read as single end
	using temp_dict
	output qs or not
	"""

	total_reads_number = wccount(parameters.sam_file)
	percentage_of_total_file = 0

	chr_seq = get_ref_geno(parameters)
	total_reads_num = 0

	with open(parameters.sam_file, "r") as inputfile_sam:
		for sam_line_first in inputfile_sam:

			insert_size_lower_bond = 0
			insert_size_upper_bond = 1000

			if not sam_line_first.startswith("@"):
				current_percent = int(float(total_reads_number * percentage_of_total_file) / 100)
				if total_reads_num == current_percent:
					print "current progress: ", percentage_of_total_file
					percentage_of_total_file += 10

				total_reads_num += 1
				elements_first = sam_line_first.strip().split()
				try:
					read_ID_first = elements_first[0].strip()
					chrName_first = elements_first[2].strip()
					insert_size_first = abs(int(elements_first[8].strip()))  # insert_size for second read is negative
				except:
					print "error in first read:", sam_line_first
				if (insert_size_first >= insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):
					if True:  # this is for pair-end verification. Not needed at this moment.
						if chrName_first.startswith(parameters.chr_name):
							# if chrName_first == parameters.chr_name:
							# first read
							qName_first = elements_first[0].strip()
							flag_first = elements_first[1].strip()
							start_position_first = int(elements_first[3].strip())
							read_sequence_first = elements_first[9].strip()
							read_length_first = len(read_sequence_first)
							quality_score_sequence_first = elements_first[10].strip()

							for i in range(read_length_first):
								current_base_position = start_position_first + i

								covered_snp = read_sequence_first[i]  # ith position is the covered snp
								quality_score_symbol = quality_score_sequence_first[i]
								if (covered_snp != 'N') and (
											(ord(quality_score_symbol) - 33) > parameters.quality_score_threshold):  # check quality_score
									#print ord(quality_score_symbol) - 33
									if current_base_position not in parameters.temp_data_dict:
										temp_pos = position_data()
										temp_pos.pos = current_base_position
										temp_pos.chr_name = parameters.chr_name
										try:
											temp_pos.ref_allele = chr_seq[current_base_position - 1]
										except:
											print "chr_seq pos error, chr_length_,current_pos:", len(chr_seq), current_base_position-1
											temp_pos.ref_allele = ""
										parameters.temp_data_dict[current_base_position] = temp_pos
										parameters.temp_data_dict_counter += 1
										parameters.total_called_pos += 1

									process_covered_allele(parameters, current_base_position, covered_snp)
									#process_covered_allele_qs(parameters, current_base_position, covered_snp, quality_score_symbol)
									"""
									if covered_snp == "A":
										parameters.temp_data_dict[current_base_position].A_depth += 1
									elif covered_snp == "T":
										parameters.temp_data_dict[current_base_position].T_depth += 1
									elif covered_snp == "C":
										parameters.temp_data_dict[current_base_position].C_depth += 1
									elif covered_snp == "G":
										parameters.temp_data_dict[current_base_position].G_depth += 1
									"""

								# output the data when it reaches the size limit
								if parameters.temp_data_dict_counter > parameters.temp_data_dict_max_limit:
									output_temp_dict(parameters)

					else:
						print "first and second read ID do not match", read_ID_first
	if len(parameters.temp_data_dict) > 0:
		#print "last temp dict size", len(parameters.temp_data_dict)
		output_temp_dict(parameters)
	return total_reads_num


def process_covered_allele(parameters, current_base_position, covered_snp):
	if covered_snp == "A":
		parameters.temp_data_dict[current_base_position].A_depth += 1
	elif covered_snp == "T":
		parameters.temp_data_dict[current_base_position].T_depth += 1
	elif covered_snp == "C":
		parameters.temp_data_dict[current_base_position].C_depth += 1
	elif covered_snp == "G":
		parameters.temp_data_dict[current_base_position].G_depth += 1

def process_covered_allele_qs(parameters, current_base_position, covered_snp, quality_score_symbol):
	"""
	output qs for each allele
	:param current_base_position:
	:param covered_snp:
	:param quality_score_symbol:
	:return:
	"""

	if covered_snp == "A":
		parameters.temp_data_dict[current_base_position].A_depth += 1
		parameters.temp_data_dict[current_base_position].A_qs += quality_score_symbol
	elif covered_snp == "T":
		parameters.temp_data_dict[current_base_position].T_depth += 1
		parameters.temp_data_dict[current_base_position].T_qs += quality_score_symbol
	elif covered_snp == "C":
		parameters.temp_data_dict[current_base_position].C_depth += 1
		parameters.temp_data_dict[current_base_position].C_qs += quality_score_symbol
	elif covered_snp == "G":
		parameters.temp_data_dict[current_base_position].G_depth += 1
		parameters.temp_data_dict[current_base_position].G_qs += quality_score_symbol


def output_temp_dict(parameters):
	#current_time = time.time()
	#print "before output", len(parameters.temp_data_dict)
	sorted_pos_list = parameters.temp_data_dict.keys()
	sorted_pos_list.sort()
	output_size = len(sorted_pos_list) - parameters.temp_data_dict_min_limit \
		if len(sorted_pos_list) >= parameters.temp_data_dict_min_limit else len(sorted_pos_list)

	if len(parameters.temp_data_dict) < parameters.temp_data_dict_max_limit:
		output_size = len(sorted_pos_list)
		#print sorted_pos_list[-1]

	for i in range(output_size):
		temp_pos_data = parameters.temp_data_dict[sorted_pos_list[i]]
		print >> parameters.output_file, temp_pos_data.pos, temp_pos_data.chr_name, temp_pos_data.ref_allele, \
				temp_pos_data.A_depth, temp_pos_data.T_depth, temp_pos_data.C_depth, temp_pos_data.G_depth,
		if temp_pos_data.A_depth > 0:
			print >> parameters.output_file, "A_qs:", temp_pos_data.A_qs,
		if temp_pos_data.T_depth > 0:
			print >> parameters.output_file, "T_qs:", temp_pos_data.T_qs,
		if temp_pos_data.C_depth > 0:
			print >> parameters.output_file, "C_qs:", temp_pos_data.C_qs,
		if temp_pos_data.G_depth > 0:
			print >> parameters.output_file, "G_qs:", temp_pos_data.G_qs,
		print >> parameters.output_file, ""

		del parameters.temp_data_dict[sorted_pos_list[i]]
	parameters.temp_data_dict_counter = 0
	#print "after output", len(parameters.temp_data_dict)
	#print "take time: ", round(time.time() - current_time, 3), "s"


def data_filter_txt(line, parameters):
	"""
	prepare data portion for output_filtered_data_txt
	"""
	item = line.strip().split()
	temp_list = [int(x) for x in item[3:7]]
	allele_list = [int(x) for x in item[3:7]]
	# if temp_list.count(0) < 3: # remove homo position
	number_of_zero = temp_list.count(0)
	if number_of_zero == 2:
		# Remove the postions with 3 zeros and filter by the second_largest_allele_depth
		temp_list.sort()
		second_largest_allele_depth = temp_list[2]
		if second_largest_allele_depth >= parameters.second_largest_allele_depth_cutoff:
			line = list_to_line(list((
					item[0], item[1], item[2], item[3], item[4], item[5], item[6], (4 - number_of_zero), temp_list[3],
					temp_list[2], temp_list[1],
					temp_list[0])))
			#print line
			return line
	return ""

def output_filtered_data_txt_all(parameters):
	"""
	to filter data from txt file, all position
	"""
	output_file_name = parameters.db_base_name + "_2nd_" + str(parameters.second_largest_allele_depth_cutoff) + ".txt"
	print "output in progress: ", parameters.db_name
	print "second largest allele: ", parameters.second_largest_allele_depth_cutoff
	with open(output_file_name, "w") as output_file:
		print >> output_file, "pos", "chr", "ref_allele", "A", "T", "C", "G", \
			"Allele#", "1st", "2nd", "3rd", "4th", "rs#", "SNP_allele_1", "SNP_allele_2", "mm_snp_overlapping"
		with open(parameters.db_name, "r") as input_file:
			for line in input_file:
				line = data_filter_txt(line, parameters)
				if line != "":
					print >> output_file, line


def output_filtered_data_txt(start_line, end_line, parameters):
	"""
	to filter data from txt file
	"""
	output_file_name = parameters.db_base_name + "_" + start_line + "_" + end_line + ".txt"

	start_line = int(start_line)
	end_line = int(end_line)
	total_row_number = end_line - start_line
	with open(output_file_name, "w") as output_file:
		print >> output_file, "pos", "chr", "ref_allele", "A", "T", "C", "G", \
			"Allele#", "1st", "2nd", "3rd", "4th", "rs#", "SNP_allele_1", "SNP_allele_2", "mm_snp_overlapping"
		if total_row_number < 0:
			print "end point smaller than start point"
			sys.exit(0)
		else:
			with open(parameters.db_name, "r") as input_file:
				for line in input_file:
					line = data_filter_txt(line)
					if line != "":
						print >> output_file, line

def output_indel(parameters):
	"""
	the sequence pair has already been processed
	now treat the read as single end
	using temp_dict
	to get the indel info
	"""

	total_reads_number = wccount(parameters.sam_file)
	percentage_of_total_file = 0

	#chr_seq = get_ref_geno(parameters)
	parameters.ref_chr_seq = get_ref_geno(parameters)
	total_reads_num = 0

	with open(parameters.sam_file, "r") as inputfile_sam:
		for sam_line_first in inputfile_sam:

			insert_size_lower_bond = 0
			insert_size_upper_bond = 1000

			if not sam_line_first.startswith("@"):
				current_percent = int(float(total_reads_number * percentage_of_total_file) / 100)
				if total_reads_num == current_percent:
					print "current progress: ", percentage_of_total_file
					percentage_of_total_file += 10

				total_reads_num += 1
				elements_first = sam_line_first.strip().split()
				try:
					read_ID_first = elements_first[0].strip()
					chrName_first = elements_first[2].strip()
					insert_size_first = abs(int(elements_first[8].strip()))  # insert_size for second read is negative
				except:
					print "error in first read:", sam_line_first
				if (insert_size_first >= insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):
					if True:  # this is for pair-end verification. Not needed at this moment.
						if chrName_first.startswith(parameters.chr_name):
							# first read
							qName_first = elements_first[0].strip()
							flag_first = elements_first[1].strip()
							start_position_first = int(elements_first[3].strip())
							read_sequence_first = elements_first[9].strip()
							read_length_first = len(read_sequence_first)
							quality_score_sequence_first = elements_first[10].strip()

							cigar = elements_first[5].strip()

							if is_indel(cigar):
								#print sam_line_first
								#print cigar

								cigar_process(chrName_first, start_position_first, read_sequence_first, cigar, parameters)
								total_reads_num += 1

	print "del info *******"
	for pos in parameters.indel_dict.keys():
		indel = parameters.indel_dict[pos]
		print pos, indel.chr_name, indel.ref_allele
		for seq in indel.deletion.keys():
			print seq, indel.deletion[seq]

	print "insertion info *******"
	for pos in parameters.indel_dict.keys():
		indel = parameters.indel_dict[pos]
		if len(indel.insertion) > 0:
			print pos, indel.chr_name, indel.ref_allele
			for seq in indel.insertion.keys():
				print seq, indel.insertion[seq]

	print "total_del_num", len(parameters.indel_dict)
	print "total_insertion_num", len([x for x in parameters.indel_dict.keys() if len(parameters.indel_dict[x].insertion) > 0])

	print "total indel number", total_reads_num
	return total_reads_num

def is_indel(cigar):
	num_indel = cigar.count("I") + cigar.count("D")
	num_N = cigar.count("N")
	is_indel = True if num_indel > 0 and num_N == 0 else False
	return is_indel

def cigar_process(chr_name, start_pos, read_seq, cigar, parameters):

	current_index = 0
	indel_length = ""
	indel_type = ""
	for info in cigar:
		try:
			int(info)
			indel_length += info
		except:
			indel_type = info
			indel_length = int(indel_length)
			if indel_type == "M":
				current_index += indel_length
			elif indel_type == "D":	# deletion, get ref seq
				deletion_pos = start_pos + current_index
				deletion_seq = parameters.ref_chr_seq[deletion_pos: deletion_pos + indel_length]
				"""
				print "deletion_pos, deletion_seq", deletion_pos, deletion_seq
				print "start_pos, current_index", start_pos, current_index
				print "read", read_seq
				print "orir", parameters.ref_chr_seq[start_pos - 1 : start_pos + len(read_seq) - 1]
				"""
				if deletion_pos not in parameters.indel_dict:
					temp_indel = indel()
					temp_indel.pos = deletion_pos
					temp_indel.chr_name = chr_name
					temp_indel.deletion[deletion_seq] = 1
					parameters.indel_dict[deletion_pos] = temp_indel
				else:
					if deletion_seq not in parameters.indel_dict[deletion_pos].deletion:
						parameters.indel_dict[deletion_pos].deletion[deletion_seq] = 1
					else:
						parameters.indel_dict[deletion_pos].deletion[deletion_seq] += 1

			elif indel_type == "I": # insertion, get inserted seq
				insertion_pos = start_pos + current_index
				insertion_seq = read_seq[current_index : current_index + indel_length]
				current_index += indel_length
				"""
				print "insertion_pos", insertion_seq
				print "start_pos, current_index", start_pos, current_index
				print "read", read_seq
				print "orir", parameters.ref_chr_seq[start_pos - 1 : start_pos + len(read_seq) - 1]
				"""

				if insertion_pos not in parameters.indel_dict:
					temp_indel = indel()
					temp_indel.pos = insertion_pos
					temp_indel.chr_name = chr_name
					temp_indel.insertion[insertion_seq] = 1
					parameters.indel_dict[insertion_pos] = temp_indel
				else:
					if insertion_seq not in parameters.indel_dict[insertion_pos].insertion:
						parameters.indel_dict[insertion_pos].insertion[insertion_seq] = 1
					else:
						parameters.indel_dict[insertion_pos].insertion[insertion_seq] += 1
			else:
				print indel_length, indel_type
			indel_length = ""
			indel_type = ""






def get_args():
	desc = "variation call"
	usage = "snpPick_fish -s sam_file -c chr -m update -d db_name -q qscore \n" \
	        "snpPick_fish -c chr -m output -b startLine -e endLine -d db_name -q qscore"
	parser = OptionParser(usage=usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile", help="Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode", help="mode", default="null")
	parser.add_option("-b", "--startLine", type="string", dest="startLine", help="start line", default="null")
	parser.add_option("-e", "--endLine", type="string", dest="endLine", help="end line", default="null")
	parser.add_option("-d", "--dbname", type="string", dest="dbname", help="db name", default="null")
	parser.add_option("-q", "--qscore", type="int", dest="qscore", help="qscore", default="40")
	parser.add_option("-p", "--pos", type="string", dest="posList", help="Input position list Name", default="null")
	(options, args) = parser.parse_args()
	if options.mode == "null" or options.chrName == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__ == '__main__':
	start_time = time.time()

	global db_name
	global db_base_name
	global ref_file
	global second_largest_allele_depth_cutoff

	options = get_args()

	chr_name = options.chrName
	db_name = options.dbname
	db_base_name = db_name
	print "db_base_name: ", db_base_name
	mode = options.mode

	parameters = parameter()
	parameters.chr_name = chr_name
	parameters.second_largest_allele_depth_cutoff = 2

	parameters.quality_score_threshold = options.qscore
	#quality_score_threshold = 13
	parameters.quality_score_threshold = 30
	print "quality_score_threshold: ", parameters.quality_score_threshold

	# gx
	# db_name = "/home/guoxing/disk2/lima/mimi_snpPick_db/" + db_name + ".db"    # for hg19 chrX mimi data
	#db_name = "/home/guoxing/disk2/lima/mimi_solid/mimi_solid_snpPick_db/" + db_name + ".db"    # for solid mimi
	#db_name = "/home/guoxing/disk2/lima/yang/mimi_yang_snpPick_db/" + db_name + ".db"    # for yang mimi

	# for zebra fish mimi data
	ref_path = "/home/guoxing/disk2/wli/ref_genome/"
	ref_file = ref_path + "danRer7_" + chr_name + ".fa"
	print ref_file

	"""
	# load hap_std_dict
	global hap_std_dict
	file_path = "/home/guoxing/disk2/solid/common_files/"
	hap_std_file = file_path + "ASW_"+chr_name+"_child_hap.txt"
	#hap_std_dict = load_raw_data(hap_std_file)[1]
	#print "hap_std_dict", len(hap_std_dict)
	"""

	global hap_std_dict
	hap_std_dict = {}

	global table_name
	table_name = "mimi"

	if (mode == "update"):
		sam_file = options.samFile
		attribute = "position INT PRIMARY KEY, chr TEXT, ref_allele TEXT, 	\
		A_depth INT, T_depth INT, C_depth INT, G_depth INT"
		#attribute = "position INT PRIMARY KEY, chr TEXT, geno_allele TEXT, total_depth INT, 	\
		#A_depth INT, T_depth INT, C_depth INT, G_depth INT, max_allele TEXT, max_allele_number INT, max_allele_percentage FLOAT"
		sql.creat_table(db_name, table_name, attribute)

		#output_data("chr1.output")
		parameters.sam_file = sam_file
		snpPick(sam_file)
	elif (mode == "output"):
		start_line = options.startLine
		end_line = options.endLine
		file_name = db_base_name + "_" + start_line + "_" + end_line + ".txt"
		#print file_name
		output_data(file_name, start_line, end_line)
	elif (mode == "filter"):
		start_line = options.startLine
		end_line = options.endLine
		file_name = db_name + "_" + start_line + "_" + end_line + "_filtered.txt"
		#print file_name
		output_data_filter(file_name, start_line, end_line)
		output_filtered_data(start_line, end_line)
	elif (mode == "mf"):
		start_line = options.startLine
		end_line = options.endLine
		output_filtered_data(start_line, end_line)
	elif (mode == "pos_list"):
		posList = options.posList
		output_single_pos_data(posList)

	elapse_time = time.time() - start_time
	print "run time: ", round(elapse_time, 3), "s"
	#print "snp_in_mimi", parameters.snp_in_mimi
	#print "total_mimi", parameters.total_mimi
	parameters.output_file.close()


