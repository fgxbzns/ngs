#!/usr/bin/python

# calculate the seed percentage and hifi error in a centain range.


import os, glob, subprocess, random, operator, time
from optparse import OptionParser
from tools import *

class data_class():
	def __init__(self):
		self.chr_name = ""

		self.seed_file_name = ""
		self.seed_file_prefix = ""
		self.seed_title_info = ""
		self.seed_dict = {}
		self.seed_error_list = []

		self.hifi_result_file = ""
		self.hifi_result_dict = {}
		self.hifi_error_list = []

		self.seed_hifi_error_dis_list = []

		self.hap_std_dict = {}

		self.hap_std_file_name = ""

		self.partition_size = 2000

		self.bucket_size = 50000
		self.bucket_number = 0
		self.bucket_dict = {}

def process_html(fbo_file_name):
	output_file_name = fbo_file_name[:-4] + "_fixed.txt"

	with open(output_file_name, "w") as output_file:
		with open(fbo_file_name, "w") as fbo_file:
			for line in fbo_file:
				if line.startswith("<"):


	pass



def calculate_distance_abs():
	data.seed_error_list.sort()
	data.hifi_error_list.sort()

	for hifi_error_pos in data.hifi_error_list:
		data.seed_hifi_error_dis_list.append(min([abs(x - hifi_error_pos) for x in data.seed_error_list]))

def calculate_distance():
	data.seed_error_list.sort()
	data.hifi_error_list.sort()

	for hifi_error_pos in data.hifi_error_list:
		temp_distance_list = [abs(hifi_error_pos - x) for x in data.seed_error_list]
		index_closest_seed_error = temp_distance_list.index(min(temp_distance_list))
		data.seed_hifi_error_dis_list.append(hifi_error_pos - data.seed_error_list[index_closest_seed_error])


def fill_bucket():
	data.bucket_number = math.ceil(float(max([abs(x) for x in data.seed_hifi_error_dis_list])/data.bucket_size))
	print data.bucket_number

	data.bucket_dict = {x: 0 for x in range((-2-int(data.bucket_number)), int(data.bucket_number)+2)}
	#print data.bucket_dict.keys()
	for distance in data.seed_hifi_error_dis_list:
		data.bucket_dict[distance / data.bucket_size] += 1
		"""
		if distance == 0:
			data.bucket_dict[0] += 1
		elif distance > 0:
			data.bucket_dict[(distance / data.bucket_size) + 1] += 1
		else:
			data.bucket_dict[(distance / data.bucket_size) - 1] += 1
		"""

	output_dict("bucket_" + str(data.bucket_size) + ".txt", data.bucket_dict)

def output_dict(output_file_name, dict):
	with open(output_file_name, "w") as output_file:
		for key in dict.keys():
			if dict[key] == 0:
				print >> output_file, key
			else:
				print >> output_file, key, dict[key]

def get_args():
	desc = ""
	usage = ""
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-c", "--chr", type="string", dest="chrName", help="Input chr Name", default="null")
	parser.add_option("-i", "--seed", type="string", dest="hifiSeed", help="Input seed file Name", default="null")
	parser.add_option("-r", "--result", type="string", dest="hifiResult", help="Input result file Name", default="null")
	(options, args) = parser.parse_args()
	"""
	if options.chrName == "null" or options.hifiSeed == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	"""
	return options


if __name__ == '__main__':
	options = get_args()

	start_time = time.time()
	global data
	data = data_class()
	data.chr_name = options.chrName
	data.seed_file_name = options.hifiSeed

	data.hifi_result_file = options.hifiResult

	seed_distribution()
	print len(data.seed_error_list)
	print len(data.hifi_error_list)

	calculate_distance()
	print len(data.seed_hifi_error_dis_list)
	#print data.seed_hifi_error_dis_list
	print get_average(data.seed_hifi_error_dis_list)

	fill_bucket()


	print "run time is: ", round((time.time() - start_time), 3), "s"
