#!/usr/local/bin/python
# -*- coding: utf-8 -*-
import pandas as pa
import IOutilities
import re
import sys


def run(annoRows, parentDirectory):
    
    print("*****versions*************")
    print(pa.__version__)
    print("pythong version is {0}".format(sys.version))

    EMBLrows = filterRowsByDB(annoRows, "EMBL")
    NCBIrows = filterRowsByDB(annoRows, "NCBI")

    EMBLcanon = getCanon(EMBLrows)
    NCBIcanon = getCanon(NCBIrows)

    EMBLcanonCount = len(EMBLrows)
    NCBIcanonCount = len(NCBIrows)

    if not rowsExistForBoth(EMBLcanon,NCBIcanon):
        if len(EMBLrows) == 1 and len(NCBIrows) == 1:
            rowsReadyToBuild = annoRows
        else:
            filteredRows = filterToCompleteData(annoRows,EMBLrows,NCBIrows,parentDirectory)
            rowsReadyToBuild = selectColumnsByCriteria(filteredRows)
    elif EMBLcanonCount == 1 and NCBIcanonCount == 1:
        rowsReadyToBuild = pa.concat([EMBLcanon, NCBIcanon])
    elif (EMBLcanonCount > 0 and NCBIcanonCount > 0):
        filteredRows = pa.concat([EMBLcanon,NCBIcanon])
        rowsReadyToBuild = selectColumnsByCriteria(filteredRows)
    return rowsReadyToBuild

def filterToCompleteData(annoRows,EMBLrows, NCBIrows, parentDirectory):

    EMBLcanon = getCanon(EMBLrows)
    NCBIcanon = getCanon(NCBIrows)

    EMBLcanonCount = len(EMBLcanon)
    NCBIcanonCount = len(NCBIcanon)

    firstFilter = pa.DataFrame
    filteredRows = pa.DataFrame()

    if rowsExistInExclusiveOr(EMBLcanon, NCBIcanon):
        if EMBLcanonCount == 0 and NCBIcanonCount > 0:
            firstFilter = concatIfSecondDataFrameHasrows(EMBLrows, NCBIcanon)
        elif EMBLcanonCount > 0 and NCBIcanonCount == 0:
            firstFilter = concatIfSecondDataFrameHasrows(EMBLcanon, NCBIrows)
        filteredRows = firstFilter
    elif rowsExistInEither(EMBLrows, NCBIrows):
        filteredRows = annoRows
    else:
        logDroppedPoint(parentDirectory)

    return filteredRows

def selectColumnsByCriteria(filteredRows):

    EMBLrows = filterRowsByDB(filteredRows, "EMBL")
    NCBIrows = filterRowsByDB(filteredRows, "NCBI")

    if rowsExistInExclusiveOr(EMBLrows,NCBIrows):
        if len(EMBLrows) > 0: selectedRow = EMBLrows
        else : selectedRow = NCBIrows
        selectedAnnoRows = selectAnnotationByCriteria(selectedRow)
    elif rowsExistForBoth(EMBLrows, NCBIrows):
        selectedAnnoRows = selectAnnotationByMatch(EMBLrows,NCBIrows)
    else :
        selectedAnnoRows = pa.DataFrame()

    return selectedAnnoRows

def selectAnnotationByMatch(EMBLrows,NCBIrows):

    fixedNCBIrows = NCBIrows.reset_index(drop=True)
    fixedEMBLrows = EMBLrows.reset_index(drop=True)

    scoreSeriesCondensed = fixedEMBLrows.apply(lambda x: returnTopMatchingScore(x, fixedNCBIrows), axis=1)
    scoreSeries = pa.concat(scoreSeriesCondensed.tolist(), ignore_index=True)
    highestScoredEMBL = scoreSeries[scoreSeries['Score'] == scoreSeries['Score'].max()]

    if len(highestScoredEMBL) != 2:
        assert("ROWS HAVE THE SAME SCORING")

    return highestScoredEMBL.drop(['Score'],axis=1)


def returnTopMatchingScore(row, NCBIrows):

    extras = row.loc['Extra']

    symbolRe = "SYMBOL=[A-Za-z0-9]{0,15}"
    biotypeRe = "BIOTYPE=[A-Za-z_]{0,50}"
    impactRe = "IMPACT=(HIGH|MODERATE|MODIFIER|LOW)"
    canonicalRe = "CANONICAL=(YES|NO)"

    symbol = "NOT-FOUND"
    biotype = "NOT-FOUND"
    impact = "NOT-FOUND"
    isCanonical = "NOT-FOUND"

    symbolMatch = re.search(symbolRe, extras)
    if symbolMatch: symbol = symbolMatch.group(0)
    biotypeMatch = re.search(biotypeRe,extras)
    if biotypeMatch: biotype = biotypeMatch.group(0)
    impactMatch = re.search(impactRe,extras)
    if impactMatch : impact = impactMatch.group(0)
    isCanonicalMatch = re.search(canonicalRe,extras)
    if isCanonicalMatch: isCanonical = isCanonicalMatch.group(0)

    condensedRows[:] = NCBIrows.apply(lambda x: calculateMatchingScore(x,row, symbol, biotype, impact, isCanonical), axis=1)
    scoredRows = pa.concat(condensedRows.tolist(), ignore_index=True)

    highestScore = scoredRows[scoredRows['Score'] == scoredRows['Score'].max()]
    return highestScore

def calculateMatchingScore(x,row, symbol, isCanonical, biotype, impact):

    x['Score'] = 0

    if len(x) > 1:
        extras = x.loc['Extra']

        if re.search(symbol,extras): x.loc['Score'] += 1000
        if re.search(isCanonical,extras): x.loc['Score'] += 100
        if re.search(biotype,extras): x.loc['Score'] += 10
        if re.search(impact,extras): x.loc['Score'] += 1

    row['Score'] = x['Score']

    ncbi = pa.DataFrame(x).transpose()
    embl = pa.DataFrame(row).transpose()

    return pa.concat([embl,ncbi])

def selectAnnotationByCriteria(row):

    biotypePC = "BIOTYPE=protein_coding"
    highImpact = "IMPACT=HIGH"
    modifier = "IMPACT=MODIFIER"
    moderateImpact = "IMPACT=MODERATE"

    row['Score'] = 0

    row.loc[row['Extra'].str.contains(biotypePC), 'Score'] += 100
    row.loc[row['Extra'].str.contains(highImpact), 'Score'] += 10
    row.loc[row['Extra'].str.contains(moderateImpact), 'Score'] += 2
    row.loc[row['Extra'].str.contains(modifier), 'Score'] += 1

    return row.loc[row['Score'] == row['Score'].max()].drop('Score',axis=1).iloc[[0]]


def filterRowsByDB(annoRows, emblOrNCBI):
    dbRows = pa.DataFrame()
    if not (len(annoRows) == 0):
        if emblOrNCBI == "EMBL":
            dbRows = annoRows[annoRows['Gene'].str.contains("^ENS") & annoRows['Feature'].str.contains("ENS")]
        elif emblOrNCBI == "NCBI":
            dbRows = annoRows[~annoRows['Gene'].str.contains("^ENS") & ~annoRows['Feature'].str.contains("^ENS")]
    return dbRows

def getCanon(annoRows):
    return annoRows[annoRows['Extra'].str.contains("CANONICAL=YES")] if len(annoRows) > 0 else pa.DataFrame()

def matchRowsBetweenDB(firstFilter):
    return firstFilter

def rowsExistInExclusiveOr(EMBLrows, NCBIrows):
    return (len(EMBLrows) > 0 and not len(NCBIrows) > 0) or (not len(EMBLrows) > 0 and len(NCBIrows) > 0)

def rowsExistInEither(EMBLrows, NCBIrows):
    return len(EMBLrows) > 0 or len(NCBIrows) > 0

def rowsExistForBoth(EMBLrows, NCBIrows):
    return len(EMBLrows) > 0 and len(NCBIrows) > 0

def concatIfSecondDataFrameHasrows(row,rowToCheck):
    if len(rowToCheck) > 0:
        concatRows = pa.concat([row,rowToCheck])
    else :
        concatRows = row
    return concatRows

def logDroppedPoint(parentDirectory):
    message = "No matching annotations found. May be an error in genomic coordinates"
    IOutilities.logMessage(parentDirectory, message)
