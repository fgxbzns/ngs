#!/usr/bin/python
#######################################################################################
# Guoxing Fu July 7, 2013
# for rna-seq data
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy
from optparse import OptionParser
from tools import *

class parameters:
	def __init__(self):
		self.vcf_dict = {}


def covert_line(line):
	#position_key = 0
	gtf_line = ""
	elements = line.strip().split(',')
	try:

		gene_id = elements[10][:-3]
		chr_name = "chr" + gene_id[2]
		#transcript_id = elements[10][-1:]
		transcript_id = elements[10]
		source_name = "ncbi"
		seq_type = "exon"
		start = elements[-2]
		end = elements[-1]

		gene_name = elements[11]
		gene_product = elements[12]
		go_number = elements[8]
		gene_function = elements[9]

		strand = "+" if int(start) <= int(end) else "-"
		#position_key = int(start) if int(start) >= int(end) else int(end)
		if strand == "-":
			start = elements[-1]
			end = elements[-2]

		gtf_line = chr_name + "\t" + source_name + "\t" + seq_type + "\t" \
			+ start + "\t" + end + "\t.\t" + strand + "\t.\t" + "gene_id \"" + gene_id \
			+ "\"; " + "transcript_id \"" + transcript_id \
			+ "\"; " + "gene_name \"" + gene_name \
			+ "\"; " + "gene_product \"" + gene_product  \
			+ "\"; " + "GO_number \"" + go_number + "\"; " \
			+ "gene_function \"" + gene_function + "\"; "
			# for sugar beet only. save t_id as gene name
			#+ "gene_name_ori \"" + gene_name + "\"; "

	except:
		print "error in ", line
		sys.exit(1)

	return chr_name, gtf_line

def load_vcf(vcf_file):

	with open(vcf_file, "r") as input_file:
		for line in input_file:
			if not line.startswith("#"):
				#print line
				elements = line.strip().split()
				#print type(elements)
				try:
					chr_name = elements[0]
					#print chr_name
					pos = elements[1]
					#print pos
					if chr_name not in parameter.vcf_dict:
						parameter.vcf_dict[chr_name] = {}
					if pos not in parameter.vcf_dict[chr_name]:
						parameter.vcf_dict[chr_name][pos] = elements
				except:
					#print "error in ", line.strip()
					pass


def filter_line(line):
	csv_line = ""
	elements = line.strip().split(',')
	depth_info = elements[-1][1:-1].split()
	chr_name = elements[0]
	pos = elements[1]

	#if True:
	try:
		kgoct_all_af = elements[13][1:-1] if elements[13] != "." else elements[13]
		if kgoct_all_af == ".":
			csv_line = line.replace(elements[-1], "").strip() + depth_info[0] + "," \
			           + depth_info[1] + "," + depth_info[2]
			if pos in parameter.vcf_dict[chr_name]:
				depth_info = parameter.vcf_dict[chr_name][pos][-1].split(":")[1].split(",")
				ref_depth = depth_info[0]
				alt_depth = depth_info[1]
				csv_line += "," + ref_depth + "," + alt_depth

		elif float(kgoct_all_af) <= 0.3:
			csv_line = line.replace(elements[-1], "").replace(elements[13], "\"" + kgoct_all_af + "\"").strip() + \
			           depth_info[0] + "," + depth_info[1] + "," + depth_info[2]
			#print csv_line
			if pos in parameter.vcf_dict[chr_name]:
				depth_info = parameter.vcf_dict[chr_name][pos][-1].split(":")[1].split(",")
				ref_depth = depth_info[0]
				alt_depth = depth_info[1]
				csv_line += "," + ref_depth + "," + alt_depth
	except:
	 	#print "error in filter", line.strip()
		pass
	return csv_line

def convert_file(file):
	with open(file[:-4] + "_filtered.csv", "w") as output_file:
		with open(file, "r") as input_file:
			head_line = "Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,GeneDetail.refGene,ExonicFunc." \
				            "refGene,AAChange.refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all," \
				            "1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,SIFT_score,SIFT_pred," \
				            "Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred," \
				            "LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score," \
				            "MutationAssessor_pred,FATHMM_score,FATHMM_pred,RadialSVM_score,RadialSVM_pred,LR_score," \
				            "LR_pred,VEST3_score,CADD_raw,CADD_phred,GERP++_RS,phyloP46way_placental,phyloP100way_vertebrate," \
				            "SiPhy_29way_logOdds,clinvar_20140211,cosmic70,nci60,Genotype,QUAL,Total_depth,Ref_depth,Alt_depth"
			print >> output_file, head_line
			for line in input_file:
				csv_line = filter_line(line)
				if csv_line != "":
					print >> output_file, csv_line

def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-v", "--vcffile", type="string", dest="vcf_name", help="Input vcf file name", default="null")
	parser.add_option("-a", "--annfile", type="string", dest="ann_name", help="Input annovar file name", default="null")
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
	vcf_name = options.vcf_name
	ann_name = options.ann_name

	global parameter
	parameter = parameters()

	#file = "10.csv"

	start_time = time.time()
	load_vcf(vcf_name)
	convert_file(ann_name)
	#os.system("sort -k 1,1 -k 4,4n " + file + ".gtf > sorted_ref.gtf")
	elapsed_time = time.time() - start_time
	print "elapsed_time is: ", round(elapsed_time, 2), "s"
	
	
	
	
	
	
	