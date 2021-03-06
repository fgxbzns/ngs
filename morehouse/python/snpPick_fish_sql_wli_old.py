#!/usr/bin/python

# for zebra fish data

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser
from tools import *
import sqlite as sql
import sqlite3 as lite

# A from Father, B from Mother
class snps:
	def __init__(self):
		self.rsID = ""
		self.position = 0
		self.A = ""
		self.B = ""
		self.depth = 0
		self.covered_reads_list = []
		self.allele_dict = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
		self.consistence = True
		self.max_allele = ""
		self.max_allele_number = 0
		self.max_allele_percentage = 0


# class to store reads from sam file
class reads:
	def __init__(self):
		self.qName = ""
		self.flag = ""
		self.chrName = ""
		self.start_position = 0
		self.read_sequence = ""
		self.quality_score_sequence = ""
		self.read_length = 0
		self.covered_snp = ""


def get_ref_geno(chr_name):
	chr_seq = ""

	#print ref_file
	input_file = open(ref_file, "r")
	for lines in input_file:
		if not lines.startswith(">"):
			chr_seq += lines.strip()
	print "total base number: ", len(chr_seq)
	return chr_seq


def variant_call_pair_end(sam_file):
	"""the sequence pair has already been processed
	now treat the read as single end """

	total_reads_number = wccount(sam_file)
	percentage_of_total_file = 0

	chr_seq = get_ref_geno(chr_name)
	#print chr_seq[19986799]

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
					insert_size_first = abs(int(elements_first[8].strip()))  #  insert_size for second read is negative
				except:
					print "error in first read:", sam_line_first
				#print "this is a new read"	
				if (insert_size_first > insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):
					#if True:
					if True:
						if chrName_first.startswith(chr_name):
							#if chrName_first == chr_name:
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
								if (not covered_snp == 'N') and (
									(ord(quality_score_symbol) - 33) > quality_score_threshold):  # check quality_score
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
										inset_querry = "INSERT INTO " + table_name + \
										               " (position, chr, ref_allele, A_depth, T_depth, C_depth, G_depth ) VALUES (" + \
										               str(current_base_position) + \
										               ",'" + chrName_first + "','" + chr_seq[
											               current_base_position - 1] + "'," + str(A_depth) + "," + str(
											T_depth) \
										               + "," + str(C_depth) + "," + str(G_depth) + ")"
										#print inset_querry
										cur.execute(inset_querry)
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


def snpPick(sam_file):
	global sam_file_name
	sam_file_name = sam_file[:(len(sam_file) - 4)]
	print "sam file: ", sam_file_name
	total_reads_num = variant_call_pair_end(sam_file)
	print "total_reads_num", total_reads_num


def get_data(db_name, table_name, start_line, end_line):
	con = lite.connect(db_name)
	with con:
		cur = con.cursor()
		querry = "SELECT * FROM " + table_name + " where position>" + start_line + " and position<" + end_line
		#print querry
		cur.execute(querry)
		rows = [[str(item) for item in results] for results in cur.fetchall()]
		#print len(rows)
		return rows


def output_data(file_name, start_line, end_line):
	total_row_number = int(end_line) - int(start_line)
	#print total_row_number
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


def output_data_filter(file_name, start_line, end_line):
	# output data with small total number
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
				print >> output_file, item[0], item[1], item[2], item[3], item[4], item[5], item[6], temp_list[3], \
				temp_list[2], temp_list[1], temp_list[0]
	print "total snp number :", current_row


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
		if True:  # wli
			#if temp_list.count(0) < 3:	# lima
			temp_list.sort()
			current_percent = int(float(total_row_number * percentage_of_total) / 100)
			if current_row == current_percent:
				#print "current progress: ", percentage_of_total, "current row:", current_row + int(start_line)
				percentage_of_total += 10
			current_row += 1
			data_list.append(list((
			item[0], item[1], item[2], item[3], item[4], item[5], item[6], temp_list[3], temp_list[2], temp_list[1],
			temp_list[0])))
	elapse_time = time.time() - start_time
	print "time: ", round(elapse_time, 3), "s"
	return data_list


def output_filtered_data(start_line, end_line):
	file_name = chr_name + "_" + start_line + "_" + end_line + "_filtered.txt"
	total_snp_num = 0
	start_line = int(start_line)
	end_line = int(end_line)
	total_row_number = end_line - start_line
	with open(currentPath + file_name, "w") as output_file:
		print >> output_file, "pos", "chr", "ref_allele", "A", "T", "C", "G"
		if total_row_number <= 0:
			print "error in start point and end point"
			sys.exit(0)
		elif total_row_number >= 100000:
			number_of_subfile = 10
			total_number_ceilling = int(math.ceil(float(total_row_number) / 100) * 100)
			print "total_number_ceilling: ", total_number_ceilling
			num_in_each_file = total_number_ceilling / number_of_subfile
			print "total_number_ceilling in each: ", num_in_each_file
			for i in range(number_of_subfile):
				current_end_line = start_line + num_in_each_file - 1 if i != number_of_subfile - 1 else end_line
				print "processing ", i, start_line, current_end_line
				data_list = data_filter(start_line, current_end_line)
				print "filtered snp number: ", len(data_list)
				total_snp_num += len(data_list)
				for data in data_list:
					print >> output_file, " ".join(str(x) for x in data)
				start_line = start_line + num_in_each_file
		else:
			print "processing ", start_line, end_line
			data_list = data_filter(start_line, end_line)
			print "filtered snp number: ", len(data_list)
			total_snp_num += len(data_list)
			for data in data_list:
				print >> output_file, " ".join(str(x) for x in data)
	print "total_snp_num: ", total_snp_num


def compare_filtered_data(a_file_name, b_file_name):
	a_file = open(currentPath + a_file_name, "r")
	b_file = open(currentPath + b_file_name, "r")

	# skip the firs line
	a_line = a_file.readline()
	b_line = b_file.readline()
	a_line = a_file.readline().strip()
	b_line = b_file.readline().strip()
	total_line_kept = 0
	total_line_with_1 = 0
	output_name = a_file_name + "_combine_" + b_file_name + "_depthcutoff_" + str(depth_cutoff)
	with open(currentPath + output_name, "w") as output_file:
		while a_line != "" or b_line != "":
			if a_line != "" and b_line == "":
				#print >> output_file, " ".join(a_line.split()[:7])
				a_line = a_file.readline().strip()
			if a_line == "" and b_line != "":
				#print >> output_file, " ".join(b_line.split()[:3]), " ", " ", " ", " ", " ".join(b_line.split()[3:7])
				b_line = b_file.readline().strip()
			if a_line != "" and b_line != "":
				a_elements = a_line.split()
				b_elements = b_line.split()
				a_pos = int(a_elements[0])
				b_pos = int(b_elements[0])
				if a_pos < b_pos:
					#print >> output_file, " ".join(a_elements[:7])
					a_line = a_file.readline().strip()
				elif a_pos > b_pos:
					#print >> output_file, " ".join(a_elements[:3]), " ".join(b_elements[3:7])
					b_line = b_file.readline().strip()
				else:
					# if the index of zero is different, the snp is different


					temp_a_elements = []
					for a in a_elements:
						temp_a = "0" if a in depth_remove_list else a
						temp_a_elements.append(temp_a)
					a_elements = temp_a_elements

					temp_b_elements = []
					for b in b_elements:
						temp_b = "0" if b in depth_remove_list else b
						temp_b_elements.append(temp_b)
					b_elements = temp_b_elements

					a_index_of_zero = [i for i in range(3, 7) if a_elements[i] == "0"]
					#print a_index_of_zero
					b_index_of_zero = [i for i in range(3, 7) if b_elements[i] == "0"]
					#print b_index_of_zero

					if len(a_index_of_zero) != 4 and len(b_index_of_zero) != 4 and a_index_of_zero != b_index_of_zero:


						keep_line = False
						a_covered_base_index = [i for i in range(3, 7) if i not in a_index_of_zero]
						b_covered_base_index = [i for i in range(3, 7) if i not in b_index_of_zero]
						if len(a_covered_base_index) > len(b_covered_base_index):
							different_index = [i for i in a_covered_base_index if i not in b_covered_base_index]
							for index in different_index:
								if int(a_elements[index]) >= depth_cutoff:
									keep_line = True
						else:
							different_index = [i for i in b_covered_base_index if i not in a_covered_base_index]
							for index in different_index:
								if int(b_elements[index]) >= depth_cutoff:
									keep_line = True

						if keep_line == True:
							total_line_kept += 1
							print >> output_file, " ".join(a_elements[:3]), " ".join(a_elements[3:7]), " ".join(
								b_elements[3:7])
					else:
						total_line_with_1 += 1
					a_line = a_file.readline().strip()
					b_line = b_file.readline().strip()
				"""
				 keep the different base with number 1
				 else:
					 # if the index of zero is different, the snp is different
					 a_index_of_zero = [i for i in range(len(a_elements)) if a_elements[i] == "0"]
					 b_index_of_zero = [i for i in range(len(b_elements)) if b_elements[i] == "0"]
					 if a_index_of_zero != b_index_of_zero:
						 print >> output_file, " ".join(a_elements[:3]), " ".join(a_elements[3:7]), " ".join(b_elements[3:7])
					 a_line = a_file.readline().strip()
					b_line = b_file.readline().strip()
				"""
	print "total_line_kept: ", total_line_kept
	print "total_line_with_1: ", total_line_with_1


def get_args():
	desc = "variation call"
	usage = "snpPick_fish -s sam_file -c chr -m update -d db_name -q qscore \nsnpPick_fish -c chr -m output -b startLine -e endLine -d db_name -q qscore"
	parser = OptionParser(usage=usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile", help="Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode", help="mode", default="null")
	parser.add_option("-b", "--startLine", type="string", dest="startLine", help="start line", default="null")
	parser.add_option("-e", "--endLine", type="string", dest="endLine", help="end line", default="null")
	parser.add_option("-d", "--dbname", type="string", dest="dbname", help="db name", default="null")
	parser.add_option("-q", "--qscore", type="int", dest="qscore", help="qscore", default="30")

	parser.add_option("-x", "--afile", type="string", dest="afile", help="afile", default="null")
	parser.add_option("-y", "--bfile", type="string", dest="bfile", help="bfile", default="null")
	parser.add_option("-z", "--dcut", type="string", dest="depthcutoff", help="depth remove cutoff", default="2")

	(options, args) = parser.parse_args()
	if options.mode == "null" or options.chrName == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options


if __name__ == '__main__':
	start_time = time.time()

	options = get_args()
	mode = options.mode

	if (mode == "compare"):
		global depth_remove_list
		global depth_cutoff
		a_file_name = options.afile
		b_file_name = options.bfile
		for depth_cutoff in range(1, 5):
			print "depth_cutoff", depth_cutoff
			depth_remove_list = []
			depth_cutoff = int(options.depthcutoff)
			for i in range(1, depth_cutoff + 1):
				depth_remove_list.append(str(i))

			print "depth_remove_list", depth_remove_list

			compare_filtered_data(a_file_name, b_file_name)
	else:
		global db_name
		global ref_file

		chr_name = options.chrName
		db_name = options.dbname

		# gx
		"""
		db_name = "/home/guoxing/disk2/" + db_name + ".db"
		ref_path = "/home/guoxing/disk2/zebra_fish/ref_genome/"
		ref_file = ref_path + "danRer7_" + chr_name + ".fa"
		ref_file = ref_path + "lm_" + chr_name + ".fa"
		
		
		# lm
		db_name = "/home/lima/disk2_node3/" + db_name + ".db"
		ref_path = "/home/lima/Public/"
		ref_file = ref_path + "chrX.fa"
		"""
		# wli

		db_name = "/home/wli/nfs1_node2/" + db_name + ".db"
		ref_path = "/home/wli/nfs1_node2/ref_genome/"
		ref_file = ref_path + "danRer7_" + chr_name + ".fa"

		global quality_score_threshold
		quality_score_threshold = options.qscore
		print "quality_score_threshold: ", quality_score_threshold

		global table_name
		table_name = "zebra"

		if (mode == "update"):
			sam_file = options.samFile
			attribute = "position INT PRIMARY KEY, chr TEXT, ref_allele TEXT, 	\
			A_depth INT, T_depth INT, C_depth INT, G_depth INT"
			#attribute = "position INT PRIMARY KEY, chr TEXT, geno_allele TEXT, total_depth INT, 	\
			#A_depth INT, T_depth INT, C_depth INT, G_depth INT, max_allele TEXT, max_allele_number INT, max_allele_percentage FLOAT"
			sql.creat_table(db_name, table_name, attribute)

			#output_data("chr1.output")
			snpPick(sam_file)
		elif (mode == "output"):
			start_line = options.startLine
			end_line = options.endLine
			file_name = chr_name + "_" + start_line + "_" + end_line + ".txt"
			#print file_name
			output_data(file_name, start_line, end_line)
		elif (mode == "filter"):
			start_line = options.startLine
			end_line = options.endLine
			file_name = chr_name + "_" + start_line + "_" + end_line + "_filtered.txt"
			#print file_name
			output_data_filter(file_name, start_line, end_line)
			output_filtered_data(start_line, end_line)
		elif (mode == "mf"):
			start_line = options.startLine
			end_line = options.endLine
			output_filtered_data(start_line, end_line)

	print "run time is: ", round((time.time() - start_time), 3), "s"
	

