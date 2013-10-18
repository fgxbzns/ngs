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
		self.allele_dict = {'A':0, 'T':0, 'C':0, 'G':0}
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
	#ref_path = "/home/guoxing/disk2/zebra_fish/ref_genome/"
	#ref_file = ref_path + "danRer7_" + chr_name + ".fa"
	#ref_file = ref_path + "lm_" + chr_name + ".fa"
	
	ref_path = "/home/lima/Public/"
	ref_file = ref_path + "chrX.fa"
	
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
		sam_line_first = inputfile_sam.readline() # the first read line in a pair
		total_reads_num = 0
		covered_snp_total_number = 0
		
		insert_size_lower_bond = 0
		insert_size_upper_bond = 1000
	
		while sam_line_first!='':
			if not sam_line_first.startswith("@"):
				current_percent = int(float(total_reads_number * percentage_of_total_file)/100)
				if total_reads_num == current_percent:
					print "current progress: ", percentage_of_total_file
					percentage_of_total_file += 10
					
				total_reads_num += 1	
				elements_first = sam_line_first.strip().split()
				try:
					read_ID_first = elements_first[0].strip()
					chrName_first = elements_first[2].strip()
					insert_size_first = abs(int(elements_first[8].strip()))			#  insert_size for second read is negative
				except:
					print "error in first read:", sam_line_first
				#print "this is a new read"	
				if (insert_size_first > insert_size_lower_bond) and (insert_size_first <= insert_size_upper_bond):
				#if True:
					if True:
						#if True:
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
								current_base_position = start_position_first+i
								A_depth = 0
								T_depth = 0
								C_depth = 0
								G_depth = 0
								
								covered_snp = read_sequence_first[i]			# ith position is the covered snp
								quality_score_symbol = quality_score_sequence_first[i]
								if (not covered_snp == 'N') and ((ord(quality_score_symbol)-33) > quality_score_threshold):	# check quality_score
									if covered_snp == "A":
										A_depth += 1
									elif covered_snp == "T":
										T_depth += 1
									elif covered_snp == "C":
										C_depth += 1
									elif covered_snp == "G":
										G_depth += 1
					
									cur.execute("SELECT *  from " + table_name + " where position=" + str(current_base_position))
									row = cur.fetchone()
									if row == None:
										inset_querry = "INSERT INTO " + table_name + \
										" (position, chr, ref_allele, A_depth, T_depth, C_depth, G_depth ) VALUES (" + \
										str(current_base_position) + \
										",'" + chrName_first + "','" + chr_seq[current_base_position-1] + "'," + str(A_depth) + "," + str(T_depth) \
										 + "," + str(C_depth) + "," + str(G_depth) + ")"
										#print inset_querry
										cur.execute(inset_querry)
									else:
										A_depth += int(row[3])
										T_depth += int(row[4])
										C_depth += int(row[5])
										G_depth += int(row[6])
										update_querry = "UPDATE " + table_name + " set A_depth=" + str(A_depth) + \
										", T_depth=" + str(T_depth) + ", C_depth=" + str(C_depth) + ", G_depth=" + \
										str(G_depth) + " where position=" + str(current_base_position)
										#print update_querry
										cur.execute(update_querry)									
					else:
						print "first and second read ID do not match", read_ID_first, read_ID_second					
			sam_line_first = inputfile_sam.readline()
		inputfile_sam.close()
	return total_reads_num

def snpPick(sam_file):
	start = time.time()		
	
	global sam_file_name
	
	sam_file_name = sam_file[:(len(sam_file)-4)]
	
	print "sam file: ", sam_file_name

	total_reads_num = variant_call_pair_end(sam_file)
	print "total_reads_num", total_reads_num
	
	end = time.time()
	run_time = str(format((end - start), "0.3f"))
	print "run time is: " + run_time + "s"

def get_data(db_name, table_name, start_line, end_line):
	con = lite.connect(db_name)
	with con:    
	    cur = con.cursor()
	    querry = "SELECT * FROM " + table_name + " where position>" + start_line + " and position<"  + end_line
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
			
			current_percent = int(float(total_row_number * percentage_of_total)/100)
			#print current_percent
			if current_row == current_percent:
				print "current progress: ", percentage_of_total, "% current row:", current_row + int(start_line)
				percentage_of_total += 10			
			current_row += 1	
			#output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (item[0], item[1], item[2], item[3], item[4], item[5], item[6]))
			print >> output_file, item[0], item[1], item[2], item[3], item[4], item[5], item[6]
	print "total snp number :", current_row

def output_data_filter(file_name, start_line, end_line):
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
				current_percent = int(float(total_row_number * percentage_of_total)/100)
				if current_row == current_percent:
					print "current progress: ", percentage_of_total, "current row:", current_row + int(start_line)
					percentage_of_total += 10			
				current_row += 1
				print >> output_file, item[0], item[1], item[2], item[3], item[4], item[5], item[6], temp_list[3], temp_list[2], temp_list[1], temp_list[0] 
	print "total snp number :", current_row

def get_args():
	desc="variation call"
	usage = "snpPick_fish -s sam_file -c chr -m update \n snpPick_fish -c chr -m output -b startLine -e endLine" 
	parser = OptionParser(usage = usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode",help = "mode", default="null")
	parser.add_option("-b", "--startLine", type="string", dest="startLine",help = "start line", default="null")
	parser.add_option("-e", "--endLine", type="string", dest="endLine",help = "end line", default="null")
	parser.add_option("-d", "--dbname", type="string", dest="dbname",help = "db name", default="null")
	parser.add_option("-q", "--qscore", type="int", dest="qscore",help = "qscore", default="40")
	(options, args) = parser.parse_args()
	if options.mode == "null" or options.chrName == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options
	
if __name__=='__main__':
	start_time = time.time()
	options = get_args()
	
	chr_name = options.chrName
	mode = options.mode	
	
	global db_name
	# gx
	#db_name = "/home/guoxing/disk2/ngs_" + chr_name + ".db"
	# lm
	db_name = options.dbname
	db_name = "/home/lima/disk2_node3/" + db_name + ".db"
	
	global quality_score_threshold
	quality_score_threshold = options.qscore	

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
		file_name = "zebra_" + chr_name + "_" + start_line + "_" + end_line + ".txt"
		#print file_name
		output_data(file_name, start_line, end_line)
	elif (mode == "filter"):
		start_line = options.startLine
		end_line = options.endLine
		file_name = "zebra_" + chr_name + "_" + start_line + "_" + end_line + "_filtered.txt"
		#print file_name
		output_data_filter(file_name, start_line, end_line)
	
	elapse_time = time.time() - start_time
	print "run time: " + str(format(elapse_time, "0.3f")) + "s"
	

