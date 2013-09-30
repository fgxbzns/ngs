#!/usr/bin/python
#######################################################################################
# hifi v2
#######################################################################################

import os, glob, subprocess, random, operator, time, sys
from optparse import OptionParser

from tools import *
from calculate_maf import calculate_maf

from pandas import Series, DataFrame
import pandas as pd

class hifi:
	
	def __init__(self):
		self.seed_file = "haplotype.txt"
		self.seed_file_name = ""
		self.seed_title_info = ""
		self.seed_dict = {}
		self.seed_homo_dict = {}
		self.seed_hetero_dict = {}
		self.seed_frame = ""
		
		self.geno_file_name = "genotype.txt"
		self.geno_title_info = ""
		self.geno_dict = {}
		self.geno_homo_dict = {}
		self.geno_hetero_dict = {}
		self.geno_frame = ""
		
		self.ref_file_name = "refHaplos.txt"
		self.ref_title_info = ""
		self.hap_ref_dict = {}
		self.ref_frame = ""
		
		self.hgr_pos = ""
	
	def load_raw_data(self):
		#self.seed_frame = pd.read_table(self.seed_file, index_col = 1)
		#self.geno_frame = pd.read_table(self.geno_file_name, index_col = 1)
		columns_name = ['rsID', 'position', 'base']
		self.seed_frame = pd.read_table(self.seed_file, delim_whitespace=True, names=columns_name)
		self.geno_frame = pd.read_table(self.geno_file_name, delim_whitespace=True, names=columns_name)
		self.ref_frame = pd.read_table(self.ref_file_name, delim_whitespace=True)
		ref_columns = [col for col in self.ref_frame.columns]
		ref_columns[1] = "position"
		self.ref_frame = self.ref_frame.reindex(columns=ref_columns)
				
		#print self.seed_frame.head()
		#print self.ref_frame.head()
		#self.print_ds(self.seed_frame)
		print self.seed_frame.shape
		print self.geno_frame.shape
		print self.ref_frame.shape		
		
	def merge_pos(self):
		"""
		This step is already done in refMerger, all positions in hap
		and geno are in ref.		
		"""
		
		#seed_pos = self.seed_frame['position']
		seed_pos = self.seed_frame[self.seed_frame.columns[1]]
		#print seed_pos
		#print len(seed_pos)
		geno_pos = self.geno_frame['position']
		#print geno_pos
		#print len(geno_pos)
		ref_pos = self.ref_frame['position']
		#ref_pos.sort()
		#print ref_pos
		#print len(ref_pos)

		self.hgr_pos = pd.concat([seed_pos, geno_pos, ref_pos]).drop_duplicates().copy()
		self.hgr_pos.sort()
		print len(self.hgr_pos)
		print self.hgr_pos
	
	def hg_merge(self):
		
		#seed_pos = [pos for pos in self.seed_frame['position']]
		#geno_pos = [pos for pos in self.geno_frame['position']]
		
		#self.hgr_pos = [pos for pos in self.ref_frame['position']]
		
		seed_pos = self.seed_frame['position']
		geno_pos = self.geno_frame['position']
		self.hgr_pos = self.ref_frame['position']
		
		seed_pos_frame = pd.DataFrame(self.seed_frame['base'].copy(), index=self.seed_frame['position'])
		print seed_pos_frame.head()
		hap_xn_A = []
		hap_xn_B = []
		a = 0
		b = 0
		for pos in self.hgr_pos:
			if pos in seed_pos:
				a += 1
				#print type(pos)
			else:
				b += 1
		print a, b
			
		
		"""
		print self.seed_frame.head()
		b = []
		for pos in self.seed_frame['position']:
			b.append(pos)
		#b = self.seed_frame['position'].copy()
		c = self.seed_frame.drop('position', axis=1)
		print c.head()
		a = c.reindex(index=b)
		print a.head()
		
		f = pd.DataFrame(self.seed_frame, index = b)
		print f.head()
		#self.print_ds(a)
		
		pass
		"""
	
	 	
	def print_ds(self, data_structure):
		if isinstance(data_structure, pd.core.frame.DataFrame):
			for index, row in data_structure.iterrows():
				print index,
				for data in row:
					print data,
				print ""
		elif isinstance(data_structure, pd.core.series.Series):
			for data in data_structure:
				print data
		"""
		for row_index in df.index:
				row = df.ix[row_index]
				for column in row.index:	
					print row.ix[column],
				print ""
		"""
	
	
	def run(self):
		self.load_raw_data()	
		#self.merge_pos()
		self.hg_merge()

def get_args():
	desc="Compare seed and std hap, to check purity of seed"

	usage = "seed_std_compare -i seed_file -c chr#" 
	parser = OptionParser(usage = usage, description=desc) 
	parser.add_option("-c", "--chr", type="string", dest="chrName",help = "Input chr Name", default="null")
	parser.add_option("-i", "--seed", type="string", dest="hifiSeed",help = "Input seed file Name", default="null")
	(options, args) = parser.parse_args()
	if options.chrName == "null" or options.hifiSeed == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options

if __name__=='__main__':
	"""
	options = get_args()
	chr_name = options.chrName
	seed_input_file = options.hifiSeed	
	seed_std_compare(seed_input_file, chr_name)
	"""
	hifi = hifi()
	hifi.run()

















"""
delim_whitespace: Parse whitespace-delimited (spaces or tabs) file (much faster than using a regular expression)

geno_frame = geno_frame.set_index(['position'], append=True)
print seed_frame.shape[1] returns tuple of the data frame shape

put df in list and concat them. remove duplicates
c = [a, b]
pd.concat(c).drop_duplicates()

seed_frame = pd.read_table(self.seed_file, sep = ' ', index_col = 1)


	 	for item, frame in geno_frame.iteritems():
	 		print item
	 		print frame

	for row_index, row in geno_frame.iterrows():
	 		#print '%s\t%s' % (row_index, row)
	 		print row_index, row
	 	for i, row_index in enumerate(seed_frame.index):
	 		print i
	 		print row_index
	 	print seed_frame.shape[1]
	 	#print seed_frame
	 	#print geno_frame
	 	a = seed_frame[0:2]
	 	b = geno_frame[1:4]
	 	print a
	 	print b
	 	print len(b)
	 	#print geno_frame.ix[29].index
	 	c = [a, b]
	 	print pd.concat(c)
	 	print pd.concat(c).drop_duplicates()
	 	
	 	
	 	hg_frame = pd.merge(seed_frame, geno_frame, how="outer")
	 	print hg_frame.columns
	 	for row_index in hg_frame.index:
	 		row = hg_frame.ix[row_index]

	 		for column in row.index:
	 			pass
	 			#print row.ix[column],
	 		#print ""

"""