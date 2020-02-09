#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import pandas as pa
import IOutilities
import re


def run(annoRows, fullFileName, parentDirectory):

    fileBaseName = fullFileName[:-4]

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
            filteredRows = filterToCompleteData(annoRows,EMBLrows,NCBIrows,fileBaseName, parentDirectory)
            rowsReadyToBuild = selectColumnsByCriteria(filteredRows)
    elif EMBLcanonCount == 1 and NCBIcanonCount == 1:
        rowsReadyToBuild = pa.concat([EMBLcanon, NCBIcanon])
    elif (EMBLcanonCount > 0 and NCBIcanonCount > 0):
        filteredRows = pa.concat([EMBLcanon,NCBIcanon])
        rowsReadyToBuild = selectColumnsByCriteria(filteredRows).reset_index(drop=True)
    return rowsReadyToBuild

def filterToCompleteData(annoRows,EMBLrows, NCBIrows, fileBaseName, parentDirectory):

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
        logDroppedPoint(fileBaseName, parentDirectory)

    return filteredRows

def selectColumnsByCriteria(filteredRows):

    EMBLrows = filterRowsByDB(filteredRows, "EMBL")
    NCBIrows = filterRowsByDB(filteredRows, "NCBI")

    if rowsExistInExclusiveOr(EMBLrows,NCBIrows):
        if len(EMBLrows) > 0: selectedRow = EMBLrows
        else : selectedRow = NCBIrows
        selectedAnnoRows = selectAnnotationByCriteria(selectedRow)
    elif rowsExistForBoth(EMBLrows, NCBIrows):
        selectedAnnoRows = selectAnnotationHighestScoringMatch(EMBLrows, NCBIrows)
    else :
        selectedAnnoRows = pa.DataFrame()

    return selectedAnnoRows

def selectAnnotationHighestScoringMatch(EMBLrows, NCBIrows):

    fixedNCBIrows = NCBIrows.reset_index(drop=True)
    fixedEMBLrows = EMBLrows.reset_index(drop=True)

    scoreSeries = returnTopMatchingPairScore(fixedEMBLrows, fixedNCBIrows)
    droppedScore = scoreSeries.drop(['Score'],axis=1)

    return droppedScore.iloc[0:2]


def returnTopMatchingPairScore(EMBLrows, NCBIrows):

    allPairedScoredRows = pa.DataFrame()

    for index, singleEmblRow in EMBLrows.iterrows():
        extraVariables = extractDataFromExtrasColumn(singleEmblRow)
        transposedEmbleRow = pa.DataFrame(singleEmblRow).transpose().reset_index(drop=True)
        pairedScoredRows = calculateScoreOfMatchesBetweenNcbiRowsAndOneEmblRow(NCBIrows, transposedEmbleRow, extraVariables)
        allPairedScoredRows = pa.concat([pairedScoredRows,allPairedScoredRows])

    highestScore = allPairedScoredRows[(allPairedScoredRows['Score'] == allPairedScoredRows['Score'].max())]
    return highestScore

def extractDataFromExtrasColumn(singleEmblRow):

    extras = singleEmblRow.loc['Extra']

    symbolRe = "SYMBOL=[A-Za-z0-9]{0,15}"
    biotypeRe = "BIOTYPE=[A-Za-z_]{0,50}"
    impactRe = "IMPACT=(HIGH|MODERATE|MODIFIER|LOW)"
    canonicalRe = "CANONICAL=(YES|NO)"

    symbol = "NOT-FOUND"
    biotype = "NOT-FOUND"
    impact = "NOT-FOUND"
    isCanonical = "NOT-FOUND"

    symbolMatch = re.search(symbolRe, extras)
    biotypeMatch = re.search(biotypeRe, extras)
    impactMatch = re.search(impactRe, extras)
    isCanonicalMatch = re.search(canonicalRe, extras)

    if symbolMatch: symbol = symbolMatch.group(0)
    if biotypeMatch: biotype = biotypeMatch.group(0)
    if impactMatch: impact = impactMatch.group(0)
    if isCanonicalMatch: isCanonical = isCanonicalMatch.group(0)

    return [symbol,biotype,impact,isCanonical]


def calculateScoreOfMatchesBetweenNcbiRowsAndOneEmblRow(NCBIrows, EmbleRow, extraVariables):

   concatEmbl = pa.DataFrame()
   concatNcbi = pa.DataFrame()
   ncbi = pa.DataFrame()
   embl = pa.DataFrame()

   for index,NCBIrow in NCBIrows.iterrows():

       transposedNcbiRow = pa.DataFrame(NCBIrow).transpose().reset_index(drop=True)
       transposedNcbiRow.loc[:,'Score'] = 0

       if len(transposedNcbiRow) == 1:
            extras = transposedNcbiRow.at[0,'Extra']

            if re.search(extraVariables[0],extras): transposedNcbiRow.at[0, 'Score'] += 1000
            if re.search(extraVariables[1],extras): transposedNcbiRow.at[0,'Score'] += 100
            if re.search(extraVariables[2],extras): transposedNcbiRow.at[0,'Score'] += 10
            if re.search(extraVariables[3],extras): transposedNcbiRow.at[0,'Score'] += 1

       EmbleRow.at[0, 'Score'] = transposedNcbiRow.at[0, 'Score']

       concatNcbi=pa.concat([ncbi,pa.DataFrame(transposedNcbiRow)])
       concatEmbl=pa.concat([embl, EmbleRow])

       ncbi = concatNcbi
       embl = concatEmbl


   return pa.concat([concatEmbl,concatNcbi])

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
    uniqDbRows = pa.DataFrame()
    if not (len(annoRows) == 0):
        if emblOrNCBI == "EMBL":
            dbRows = annoRows[annoRows['Gene'].str.contains("^ENS") & annoRows['Feature'].str.contains("ENS")]
            uniqDbRows = dbRows.drop_duplicates()
        elif emblOrNCBI == "NCBI":
            dbRows = annoRows[~annoRows['Gene'].str.contains("^ENS") & ~annoRows['Feature'].str.contains("^ENS")]
            uniqDbRows = dbRows.drop_duplicates()
    return uniqDbRows

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

def logDroppedPoint(baseName,parentDirectory):
    message = "No matching annotations found. May be an error in genomic coordinates"
    IOutilities.logMessage(parentDirectory, baseName, message)
