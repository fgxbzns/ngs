#!/usr/bin/python
import os, glob, subprocess, random
from optparse import OptionParser

homePath=os.getenv("HOME")
projectPath = homePath+'/project/dadh/ammp/'
ligandPath = homePath+'/project/dadh/ammp/dadhIAA/'
ammpPath = homePath+'/tool/sadgeom/ammp'
currentPath = os.getcwd() + '/'

# Reading options
usage = "usage: %prog [options] arg1 arg2" 
parser = OptionParser(usage = usage) 
parser.add_option("-f", "--file", type="string", dest="fileName",help = "Input File Name")
parser.add_option("-c", type ="string",dest = "ligandCharge", help = "Charge of the ligand", default='0')
(options, args) = parser.parse_args()

ligandName = options.fileName
ligandCharge = options.ligandCharge

# Prepare ammp file
ammpFile = open (currentPath+ligandName+'_minimize.ammp','w')
ammpFile.write('output dadh_'+ligandName+'.log;'+'\n')
ammpFile.write('echo off;'+'\n')
ammpFile.write('read dadh.ammp;'+'\n')
ammpFile.write('read fad.ammp;'+'\n')
ammpFile.write('read '+ligandName+'.ammp;'+'\n')
ammpFile.write('\n')

ammpFile.write('echo;'+'\n')
ammpFile.write('seti ires 995;'+'\n')
ammpFile.write('loopi end: ires 1376;'+'\n')
ammpFile.write('mov jres ires;'+'\n')
ammpFile.write('add jres 1;'+'\n')
ammpFile.write('read peplink.amp;'+'\n')
ammpFile.write('end:;'+'\n')
ammpFile.write('\n')

ammpFile.write('nzinactive 0 10000000;'+'\n')
ammpFile.write('abuild 15;'+'\n')
ammpFile.write('\n')

ammpFile.write('momadd 150000 150099;'+'\n')
ammpFile.write('mom 0 15;'+'\n')
ammpFile.write('momadd 160000 160099;'+'\n')
ammpFile.write('mom '+ligandCharge+' 15;'+'\n')
ammpFile.write('\n')

ammpFile.write('setf mmbox 10.;'+'\n')
ammpFile.write('setf mxdq 0.75;'+'\n')
ammpFile.write('\n')

ammpFile.write('use none angle;'+'\n')
ammpFile.write('cngdel 200;'+'\n')
ammpFile.write('use bond;'+'\n')
ammpFile.write('cngdel 200;'+'\n')
ammpFile.write('use hybrid;'+'\n')
ammpFile.write('cngdel 200;'+'\n')
ammpFile.write('\n')

ammpFile.write('use none angle;'+'\n')
ammpFile.write('cngdel 200;'+'\n')
ammpFile.write('use bond;'+'\n')
ammpFile.write('cngdel 200;'+'\n')
ammpFile.write('use hybrid;'+'\n')
ammpFile.write('cngdel 200;'+'\n')
ammpFile.write('monitor;'+'\n')
ammpFile.write('\n')

ammpFile.write('active 0 100000000;'+'\n')
ammpFile.write('use none angle bond torsion hybrid nonbon;'+'\n')
ammpFile.write('cngdel 500;'+'\n')
ammpFile.write('use none angle bond torsion hybrid nonbon;'+'\n')
ammpFile.write('cngdel 500;'+'\n')
ammpFile.write('monitor;'+'\n')
ammpFile.write('\n')

#ammpFile.write('use none nonbon;'+'\n')
#ammpFile.write('analyze 150000 150099;'+'\n')

ammpFile.write('use none nonbon;'+'\n')
ammpFile.write('analyze 160000 160099;'+'\n')
ammpFile.write('monitor;'+'\n')
ammpFile.write('\n')

ammpFile.write('echo off; output dadh_'+ligandName+'.ammp;'+'\n')
ammpFile.write('dump atom bond angle hybrid restrain torsion tether ;'+'\n')
ammpFile.write('close;'+'\n')
ammpFile.write('output dadh_'+ligandName+'.pdb; dump pdb; close;'+'\n')
ammpFile.write('exit;'+'\n')
ammpFile.close()

ammp = ammpPath+' < '+ligandName+'_minimize.ammp'
ammpProcess = subprocess.Popen(ammp, shell=True)
ammpProcess.wait()

os.system('rm '+currentPath + ligandName+'_minimize.ammp')
