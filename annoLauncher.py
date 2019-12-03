#!/usr/local/bin/python
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

    for file in files :

        baseName = os.path.basename(file)
        provider = os.path.basename(os.path.dirname(os.path.dirname(file)))

        jobName = "{0}-{1}".format(provider,baseName)
        jobScript = '/nfs/nobackup/spot/mouseinformatics/pdx/omicAnno/annoJob.sh'

        print('bsub bash' + jobScript + jobName + ' ' + file)

        cmd = ['bsub ','bash ', jobScript,' ', jobName,' ', file] 
        sp.call(cmd)


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
