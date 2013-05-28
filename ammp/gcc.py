#!/usr/bin/python
import os, glob, subprocess, random, operator
from optparse import OptionParser

homePath=os.getenv("HOME")

currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-f", "--file", type="string", dest="fileName",help = "Input File Name", default="null")
(options, args) = parser.parse_args()

fileName = options.fileName

os.system('gcc -g '+currentPath + fileName+'.c -o '+fileName)
os.system('chmod +x '+currentPath + fileName)
