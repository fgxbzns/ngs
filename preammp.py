#!/usr/bin/python
import os, glob, subprocess, random
from optparse import OptionParser

homePath=os.getenv("HOME")
projectPath = homePath+'/project/dadh/ammp/'
ligandPath = homePath+'/project/dadh/ammp/IAA/'
ammpPath = homePath+'/tool/sadgeom/ammp'
atomstuna = homePath+'/tool/preammp/atoms.tuna'

templatePath = '/home/guoxing/tool/preammp/new_charge/'
#templatePath = homePath+'/tool/preammp/DAA/'
#templatePath = homePath+'/tool/preammp/IAA/'
preammpPath = homePath+'/tool/preammp/presp4'
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-f", "--file", type="string", dest="fileName",help = "Input File Name")
(options, args) = parser.parse_args()

ligandName = options.fileName

preammpFile = open (currentPath+'preammp.ammp','w')
preammpFile.write(atomstuna + '\n')
preammpFile.write(templatePath + '\n')
preammpFile.write(currentPath + ligandName + '.pdb' + '\n')
preammpFile.write(currentPath + ligandName + '.ammp' +'\n')
preammpFile.close()

os.system(preammpPath + ' < preammp.ammp')
os.system('rm '+currentPath + 'preammp.ammp')

