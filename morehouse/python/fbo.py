#!/usr/bin/python

# calculate the seed percentage and hifi error in a centain range.


import os, glob, subprocess, operator, time, sys
from optparse import OptionParser

class data_class():
	def __init__(self):
		self.researved_word = ["PRESOL", "MOD", "AWARD", "SRCSGT", "SNOTE"]

def process_html_old(fbo_file):
	"""
	old methods, do not output multiple lines in DESC
	:param fbo_file:
	:return:
	"""
	output_file_name = fbo_file[:-4] + "_fixed.txt"
	temp_desc = ""

	with open(output_file_name, "w") as output_file:
		with open(fbo_file, "r") as fbo_file:
			for line in fbo_file:
				line = line.strip()
				if line.startswith("<"):
					if line.startswith("</"):
						print >> output_file, line
					else:
						label = line[1:line.find(">")]
						if label in data.researved_word:
							print >> output_file, line
						else:
							line += "<" + label + ">"
							print >> output_file, line
				else:
					print >> output_file, line

def process_html(fbo_file_name):
	"""
	method for processing multiple lines in DESC
	:param fbo_file_name:
	:return:
	"""
	remline_file_name = fbo_file_name[:-5] + "_remline.txt"
	with open(remline_file_name, "w") as output_file:
		with open(fbo_file_name, "r") as fbo_file:
			for line in fbo_file:
				if line.strip() != "":
					print >> output_file, line.strip()

	output_file_name = fbo_file_name[:-5] + "_fixed.txt"

	fbo_file = open(remline_file_name, "r")
	with open(output_file_name, "w") as output_file:
		line = fbo_file.readline().strip()
		while line != "":
			if line.startswith("<"):
				if line.startswith("</"):
					# closing tag
					print >> output_file, line
					line = fbo_file.readline().strip()
				else:
					label = line[1:line.find(">")]
					if label in data.researved_word:
						print >> output_file, line
						line = fbo_file.readline().strip()
					elif label == "DESC" or label == "AWARDEE":
						next_line = fbo_file.readline().strip()
						if next_line.startswith("<"):
							# only one line for DESC
							line += "<" + label + ">"
							print >> output_file, line

							label = next_line[1:next_line.find(">")]
							next_line += "<" + label + ">"
							print >> output_file, next_line
							line = fbo_file.readline().strip()

						else:
							# multiple lines for DESC
							print >> output_file, line
							line = next_line.strip()
							while not line.startswith("<") and line.strip() != "":
								# loop over all the lines for DESC
								next_line = fbo_file.readline().strip()
								if next_line.startswith("<"):
									line += "<" + label + ">"
									print >> output_file, line
									line = next_line.strip()
								else:
									print >> output_file, line
									line = next_line.strip()
					else:
						# for all other tags
						line += "<" + label + ">"
						print >> output_file, line
						line = fbo_file.readline().strip()
			else:
				print "please check ", line
				sys.exit(1)
	fbo_file.close()


def get_args():
	desc = ""
	usage = "fbo -i fbo.html"
	parser = OptionParser(usage=usage, description=desc)
	parser.add_option("-i", "--fbo", type="string", dest="fbo_file", help="Input fbo file Name", default="null")
	(options, args) = parser.parse_args()

	if options.fbo_file == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)

	return options


if __name__ == '__main__':
	options = get_args()

	start_time = time.time()
	global data
	data = data_class()
	fbo_file = options.fbo_file

	process_html(fbo_file)

	print "run time is: ", round((time.time() - start_time), 3), "s"
