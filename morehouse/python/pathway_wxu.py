#!/usr/bin/python

import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser

currentPath = os.getcwd() + '/'

def load_microarray_data(microarray_file_name):
	microarray_dict = {}
	with open(microarray_file_name, "r") as file:
		for line in file:
			if not line.startswith("SYMBOL"):
				elements = line.strip().split(",")
				#print elements[0].upper()
				microarray_dict[elements[0].upper()] = elements
	return microarray_dict

def load_pathway_data(pathway_file_name):
	pathway_dict = {}
	with open(pathway_file_name, "r") as file:
		for line in file:
			if not line.startswith("Name"):
				elements = line.strip().split(",")
				pathway_dict[elements[3]] = elements
	return pathway_dict

def load_process_data(process_file_name):
	process_dict = {}
	process_list = []
	with open(process_file_name, "r") as file:
		for line in file:
			#if not line.startswith("Aorta"):
			if not line.startswith("Heart-"):
				elements = line.strip().split(",")
				try:
					if int(elements[4]) > 1:
						process_list.append(elements[0])
						if elements[0] != "":
							process_dict[elements[0]] = elements
				except:
					print "load_process_data error", line
	return process_dict, process_list

def data_process():
	#data_process(microarray_file_name, pathway_file_name, process_file_name)
	microarray_file_name = "C:\wxu_test\microarray.csv"
	pathway_file_name = "C:\wxu_test\disease_488_pathway.csv"
	process_file_name = "C:\wxu_test\disease_488.csv"
	
	microarray_dict = load_microarray_data(microarray_file_name)
	print "microarray_dict", len(microarray_dict)
	pathway_dict = load_pathway_data(pathway_file_name)
	print "pathway_dict", len(pathway_dict)
	print pathway_dict.keys()
	process_dict, process_list = load_process_data(process_file_name)
	print "process_dict", len(process_dict)
	#print process_dict
	
	combined_data_dict = {}
	gene_not_match_microarray_pathway_studio = 0
	for cell_process in process_dict.keys():
		gene_dict = {}
		if cell_process in pathway_dict:
			for gene in pathway_dict[cell_process][4].strip().split(";"):
				if gene in microarray_dict:
					gene_dict[gene] = [float(x) for x in microarray_dict[gene][1:]]
				else:
					#print gene
					gene_not_match_microarray_pathway_studio += 1
			#calculate_mean_sd(gene_dict)
			combined_data_dict[cell_process] = calculate_mean_sd(gene_dict)
		else:
			print cell_process
	print "combined_data_dict", len(combined_data_dict)	
	print "gene_not_match_microarray_pathway_studio", gene_not_match_microarray_pathway_studio
	
	with open("C:\wxu_test\output.txt", "w") as output:
		for process in process_list:
			if process != "" and process in pathway_dict:
				print >> output, process + ":" + str(format(combined_data_dict[process][0], "0.3f")) + "+-" + str(format(combined_data_dict[process][1], "0.3f"))
			else:
				print >> output, ""
				
def calculate_mean_sd(gene_dict):
	data_list = []
	for gene in gene_dict.keys():
		if sum(gene_dict[gene][:3]) != 0:
			data_list.append(sum(gene_dict[gene][3:])/sum(gene_dict[gene][:3]))
		else:
			#print "control is zero"
			data_list.append(0.)
	mean = sum(data_list)/float(len(data_list))
	#print mean

	variance = map(lambda x: (x - mean)**2, data_list)
	#print len(variance)
	variance_mean = sum(variance)/float(len(variance))
	standard_deviation = math.sqrt(variance_mean)
	#print "standard_deviation", standard_deviation
	return mean, standard_deviation
	
def get_args():
	desc="combine vcf file"
	usage = "file_combine -c chr_name" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
	parser.add_option("-a", "--afile", type="string", dest="a_file_name",help = "Input file Name", default="null")
	parser.add_option("-b", "--bfile", type="string", dest="b_file_name",help = "Input file Name", default="null")
	(options, args) = parser.parse_args()
	
	if options.chrName == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	
	return options

if __name__=='__main__':
	start = time.time()	
	data_process()
	end = time.time()
	run_time = str(format((end - start), "0.3f"))
	print "run time is: " + run_time + "s"
	
