#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import re
import subprocess as sp
import os

print("python version is {0}".format(sys.version))

if len(sys.argv) == 2:
    dataDir = sys.argv[1]
else : print("Please enter root directory of data as parameters")

def main():

    files = walkDirForMutData(dataDir)

    if len(files) == 0 : print("No files found. Only TSV and CSV is currently supported. Bro, did you even LIFT?")

    print(os.getuid())

    for aFile in files :

        baseName = os.path.basename(aFile)
        provider = os.path.basename(os.path.dirname(os.path.dirname(aFile)))
	
	runLog = "{}/{}.log".format(os.path.dirname(aFile),baseName)	
	errorLog =  "{}/{}.err".format(os.path.dirname(aFile),baseName) 

    	print(provider)
    	jobName = "{0}-{1}".format(provider,baseName)

    	cmd = 'bsub -u afollette@ebi.ac.uk -J {0} -n 4 -M 15000 -R "rusage[mem=15000]" -e {1} -o {2}  python2.7 /nfs/nobackup/spot/mouseinformatics/pdx/pdx_Annotater/AnnotationMerger.py {3}'.format(jobName, errorLog, runLog, aFile) 
    	print(cmd)
    	sp.call(cmd,shell=True)

def makeFileIfNoneExistent(file):
	if not os.path.exists(runLog):
        	os.makedirs(runLog)

def walkDirForMutData(dataDir):

    matches = []
    for dirpath, dirnames, filenames in os.walk(dataDir):
        for x in dirnames + filenames:
            filepath = os.path.join(dirpath, x)
            if isMutationData(filepath):
                print("regex match:" + filepath)
                matches.append(filepath)
    return matches

def isMutationData(filepath):
    return re.match(".+/mut/.+(tsv|csv)$", filepath)

main()
