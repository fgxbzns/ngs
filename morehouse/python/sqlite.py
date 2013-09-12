#!/usr/bin/python

# this is for pre-processing the sam file before snpPick
"""
retrieved data from sqlite has a "u" in front of the output
[(u'arun',), (u'kiran',), (u'kulku',), (u'sugu',)]
The u prefix indicates that the string is an Unicode string
If you're sure that the string is a basic ASCII string, 
you can transform the unicode string into a non-unicode one with a simple str().
"""



import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser
from tools import *

import sqlite3 as lite
import sys


def check_version(db_name):
	"""
	con = None
	try:
	    con = lite.connect(db_name)
	    cur = con.cursor()    
	    cur.execute('SELECT SQLITE_VERSION()')
	    data = cur.fetchone()	    
	    print "SQLite version: %s" % data                	    
	except lite.Error, e:	    
	    print "Error %s:" % e.args[0]
	    sys.exit(1)	    
	finally:
	    if con:
	        con.close()
	"""  
	con = lite.connect(db_name)
	with con:
	    cur = con.cursor()    
	    cur.execute('SELECT SQLITE_VERSION()')    
	    data = cur.fetchone()
	    print "SQLite version: %s" % data 
	

def creat_table(db_name, table_name, attribute):
	try:
	    con = lite.connect(db_name)
	    cur = con.cursor()    
	    cur.execute("DROP TABLE IF EXISTS " + table_name)
	    cur.execute("CREATE TABLE "+table_name+" ("+attribute+")")
	    #con.commit()    
	except lite.Error, e:
	    if con:
	        con.rollback()
	    print "Error %s:" % e.args[0]
	    sys.exit(1)
	finally:
	    if con:
	        con.close() 

def execute_querry(con, querry):
	#con = lite.connect(db_name)
	with con:    
	    cur = con.cursor()    
	    cur.execute(querry)

def execute_querry_0(db_name, querry):
	con = lite.connect(db_name)
	with con:    
	    cur = con.cursor()    
	    cur.execute(querry)
	    #rows = cur.fetchall()
	    #rows = [[str(item) for item in results] for results in cur.fetchall()]
	    #for row in rows:
	    #   print row

def write_data(db_name, table_name, value):
	"""
	try:
	    con = lite.connect(db_name , isolation_level=None)
	    cur = con.cursor()
	    cur.execute("INSERT INTO "+table_name+" VALUES (" + value + ")")    
	    #con.execute("INSERT INTO zebra_fish (position, chr, read_ID) \
      #VALUES (1001, 'chr', '1234')");
	    con.commit()    
	except lite.Error, e:
	    if con:
	        con.rollback()
	    print "Error %s:" % e.args[0]
	    sys.exit(1)
	finally:
	    if con:
	        con.close()
	"""
	
	con = lite.connect(db_name)
	with con:
	    cur = con.cursor()    
	    cur.execute("INSERT INTO "+table_name+" VALUES (" + value + ")")    

def retrive_data(db_name, table_name):
	con = lite.connect(db_name)
	with con:    
	    cur = con.cursor()    
	    cur.execute("SELECT * FROM " + table_name)
	    #rows = cur.fetchall()
	    rows = [[str(item) for item in results] for results in cur.fetchall()]
	    for row in rows:
	        print row

def select_data(db_name, table_name, querry):
	con = lite.connect(db_name)
	with con:
		cur = con.cursor()
		cur.execute("SELECT *  from " + table_name + " where position=1002")
		#rows = cur.fetchall()
		rows = [[str(item) for item in results] for results in cur.fetchall()]
		#for row in rows:
		#	print row

def check_existing_data(db_name, table_name, pos):
	con = lite.connect(db_name)
	with con:
		cur = con.cursor()
		cur.execute("SELECT *  from " + table_name + " where position=" + str(pos))
		rows = cur.fetchall()
		if len(rows) != 0:
			return True
		else:
			return False
	        
def get_args():
	desc="variation call"
	usage = "snpPick_fish -s sam_file" 
	parser = OptionParser(usage = usage)
	parser.add_option("-s", "--sam", type="string", dest="samFile",help = "Input File Name", default="null")
	parser.add_option("-m", "--mode", type="string", dest="mode",help = "", default="null")
	(options, args) = parser.parse_args()
	if options.samFile == "null":
		print "parameters missing..."
		print usage
		sys.exit(1)
	return options


	

def sam_process():
	options = get_args()
	sam_file = options.samFile
	mode = options.mode
	sam_file_name = sam_file[:(len(sam_file)-4)]
	
	if mode == "sp":
		pair_end_indel(sam_file)
	elif mode == "mp":
		pair_end_indel_multiple()
	elif mode == "sep":
		seperate_by_chr(sam_file)
	
	
if __name__=='__main__':
	
	#start = time.time()
	
	global db_name
	db_name = "/home/guoxing/disk2/ngs.db"
	
	check_version(db_name)
	
	table_name = "test"
	attribute = "position INT PRIMARY KEY, chr TEXT, read_ID TEXT"
	creat_table(db_name, table_name, attribute)
	
	value = "1001, 'chr1', '1234'"
	write_data(db_name, table_name, value)
	value = "1002, 'chr2', '12355555545'"
	write_data(db_name, table_name, value)
	
	retrive_data(db_name, table_name)
	
	select_data(db_name, table_name, "a")
	
	#check_existing_data(db_name, table_name, 1002)
	"""
	db_name = "/home/guoxing/disk2/ngs.db"
	table_name = "zebra_fish1"
	attribute = "position INT PRIMARY KEY, read_ID TEXT, chr TEXT, geno_allele TEXT, total_depth INT, 	\
	A_depth INT, T_depth INT, C_depth INT, G_depth INT, max_allele TEXT, max_allele_number INT, max_allele_percentage FLOAT"
	creat_table(db_name, table_name, attribute)
	
	querry = "INSERT INTO zebra_fish1 (position, read_ID, chr, A_depth,T_depth, C_depth, G_depth ) VALUES \
	(48079216,'ILLUMINA-8ACF53:9:FC66462AAXX:1:13:7935:4947','chr1',0,0,1,0)"

	execute_querry(db_name, querry)
	"""










 
