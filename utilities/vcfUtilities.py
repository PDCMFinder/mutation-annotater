import csv
import re
import pandas as pa

def formatChromo(givenChromo):
    incorrectChrFormat = "(?i)^([0-9]{1,2}|[xym]{1}|mt|un)$"
    isMatch = re.match(incorrectChrFormat, givenChromo)
    if isMatch:
        chromo = "chr" + isMatch.group(1)
    else:
        chromo = givenChromo
    return chromo

def sortInPlace(aVCFfile):
    with open(aVCFfile, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        secondKey = sorted(reader, key=lambda x: sortLocation(x[1]))
        sortedFile = sorted(secondKey, key=lambda x: sortChromo(x[0]))
    with open(aVCFfile, 'w') as outFile:
        writer = csv.writer(outFile, delimiter='\t')
        for row in sortedFile:
            writer.writerow(row)

def sortChromo(x):
    decMatch = "^[0-9]{1,2}$"
    isDec = re.match(decMatch, x)

    if isDec:
        return int(x)
    elif len(x) == 1:
        return ord(x)
    else :
        return 0

def sortLocation(x):
    firstTenD = "^[0-9]{1,10}"
    loc = re.search(firstTenD, x)
    return int(loc.group(0)) if x != 'pos' and loc != None else 0

def dropDuplicates(vcfFile):
    vcfDf = pa.read_csv(vcfFile, sep='\t', keep_default_na=False, na_values=[''], dtype=str)
    vcfDf.drop_duplicates(inplace=True)
    vcfDf.to_csv(vcfFile, sep='\t', index=False, na_rep='')