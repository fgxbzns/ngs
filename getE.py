#!/usr/bin/python
import os, glob, subprocess, random, operator
from optparse import OptionParser

homePath=os.getenv("HOME")
#projectPath = homePath+'/project/dadh/ammp/'
#ligandPath = homePath+'/project/dadh/ammp/IAA/'
#ammpPath = homePath+'/tool/sadgeom/ammp'
#preammpPath = homePath+'/tool/preammp/preammp.py'

currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1" 
parser = OptionParser(usage = usage) 
parser.add_option("-f", "--file", type="string", dest="fileName",help = "Input File Name", default="null")
(options, args) = parser.parse_args()
ligandName = options.fileName

class energy:
	def __init__(me, ligand, energy):
		me.ligand = ligand
		me.energy = energy

def sortDict(adict):
    items = adict.items()
    items.sort()
    print [value for key, value in items]


def getEnergy(ligandName):
	f = open(currentPath+'dadh_'+ligandName + '.log', "r")
	for line in f:
		if 'Vnonbon total external' in line :
			x = float(line.strip().split()[3])
			return x
			
if(ligandName != "null"):
	print ligandName+" "+str(getEnergy(ligandName))
else:
	energyList=[]	
	for infile in glob.glob(os.path.join(currentPath,'*.log')):
		ligandName = infile[(infile.find('.')-3):infile.find('.')]
		#print 'current file is '+ ligandName 
		energyList.append(energy(ligandName, getEnergy(ligandName)))
		energyList.sort(key=operator.attrgetter('energy'))
	for a in energyList:
		print a.ligand + ' ' + str(a.energy)

		
		



