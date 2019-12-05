import csv
import re



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
    return int(loc) if x != 'pos' and loc != None else 0

def sort(aVCFfile, saveTo):
    with open(aVCFfile, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        secondKey = sorted(reader, key=lambda x: sortLocation(x[1]))
        sortedFile = sorted(secondKey, key=lambda x: sortChromo(x[0]))

    with open(saveTo, 'w') as outFile:
        writer = csv.writer(outFile, delimiter='\t')
        for row in sortedFile:
            writer.writerow(row)
