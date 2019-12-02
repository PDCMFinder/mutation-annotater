import csv
import operator
import re
import sys


def sortFun(x):

    decMatch = "^[0-9]{1,2}$"
    isDec = re.match(decMatch, x)

    if isDec:
        return int(x)
    elif len(x) == 1:
        return ord(x)
    else :
        return 0

def sortLocation(x):
    return int(x) if x != 'pos' else 0

def sort(aVCFfile, saveTo):
    with open(aVCFfile, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        secondKey = sorted(reader, key=lambda x: sortLocation(x[1]))
        sortedFile = sorted(secondKey, key=lambda x: sortFun(x[0]))

    with open(saveTo, 'w') as outFile:
        writer = csv.writer(outFile, delimiter='\t')
        for row in sortedFile:
            writer.writerow(row)
