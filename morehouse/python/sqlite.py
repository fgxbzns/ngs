#!/usr/bin/python

# this is for pre-processing the sam file before snpPick

import os, glob, subprocess, random, operator, time, sys, math
from optparse import OptionParser
from tools import *

import sqlite3 as lite
import sys


def check_version():
	con = None
	try:
	    con = lite.connect('ngs.db')
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
	con = lite.connect('ngs.db')
	with con:
	    cur = con.cursor()    
	    cur.execute('SELECT SQLITE_VERSION()')    
	    data = cur.fetchone()
	    print "SQLite version: %s" % data 
	"""

def creat_table(table_name, attribute):
	try:
	    con = lite.connect('ngs.db')
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

def write_data(table_name, value):
	try:
	    con = lite.connect('ngs.db' , isolation_level=None)
	    cur = con.cursor()
	    cur.execute("INSERT INTO "+table_name+" VALUES ("+value+")")    
	    #con.commit()    
	except lite.Error, e:
	    if con:
	        con.rollback()
	    print "Error %s:" % e.args[0]
	    sys.exit(1)
	finally:
	    if con:
	        con.close() 

def retrive_data(table_name):
	con = lite.connect('ngs.db')
	with con:    
	    cur = con.cursor()    
	    cur.execute("SELECT * FROM " + table_name)
	    rows = cur.fetchall()
	    for row in rows:
	        print row
	        
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
	check_version()
	
	table_name = "zebra_fish"
	attribute = "position INT PRIMARY KEY, chr TEXT, read_ID TEXT"
	creat_table(table_name, attribute)
	
	value = "1001, chr2, 1234"
	#write_data(table_name, value)
	retrive_data(table_name)


















 
