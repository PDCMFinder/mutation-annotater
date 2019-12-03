#!/usr/local/bin/python
# -*- coding: utf-8 -*-
import filter

import csv
import json
import os
import sys
import time

import IOutilities
import re
import pandas as pa

mergedPointsMissed = 0


if len(sys.argv) > 1:
    TSVfilePath = sys.argv[1]
    parentDirectory = os.path.dirname(sys.argv[1])
else :
    sys.stderr.write("Warning: Merger is being ran without file input. This should only be used for testing")

def run():

    print("AnnoMergers python version {0}".format(sys.version))

    mergeRowsAndWrite()
    IOutilities.saveToExcel(TSVfilePath)


def isColumnHeader(line):
    return ((line[0] == '#') and (line[1] != '#'))


def isCanonical(line):
    return (line.find("CANONICAL=YES") != -1)


def mergeRowsAndWrite():
    with open(TSVfilePath + ".hmz", 'w+') as finalTemplate, \
            open(TSVfilePath + ".ANN", 'r') as annoFile, \
            open(TSVfilePath, 'r') as tsvFile:

        print('***********************************************')
        print(tsvFile)

        outFileWriter = csv.writer(finalTemplate, delimiter="\t")
        TSVreader = csv.DictReader(tsvFile, delimiter="\t")
        annoReader = pa.read_csv(annoFile, delimiter='\t', error_bad_lines=False, header=97)

        message = "Mergings original data : {0} /n and annotated data : {1} at {2}".format(TSVfilePath, annoFile, time.ctime())
        IOutilities.logMessage(parentDirectory, message)
        print(message)

        headers = buildHeaders()
        outFileWriter.writerow(headers)

        rowNum = 0
        rowAdded = 0

        for row in TSVreader:
            rowNum += 1

            if rowIsValidForMerge(row):
                mergedRow = mergeRows(row, annoReader)
                if not len(mergedRow) < 26:
                    outFileWriter.writerow(mergedRow)
                    rowAdded += 1
                else:
                    message = ("Info: Dropping row for incorrect template size. RowNum {0} - Len {1}".format(rowNum,
                                                                                                             len(
                                                                                                                 mergedRow)))
                    IOutilities.logMessage(parentDirectory, message)
            else:
                message2 = ("Info: row {0} is broken or legacy".format(rowNum)) 
                print(row)
                IOutilities.logMessage(parentDirectory, message2)

        masterLog = "{0}/masterLog".format(os.path.dirname(os.path.dirname(parentDirectory)))
        message = "{0} The completed file file {1} has {2} data points (including header)".format(time.ctime(),finalTemplate, rowAdded)
        IOutilities.masterlogMessage(masterLog, message)


def rowIsValidForMerge(row):
    return rowIsHg38(row) and getFromRow(row, "chromosome") and getFromRow(row, "seq_start_position")


def rowIsHg38(row):
    hg38Regex = "(?i)(hg38|grch38|38)"
    return re.match(hg38Regex, getFromRow(row, "genome_assembly"))


def mergeRows(row, annoReader):
    annoRows = compareKeysOfFileAndReturnMatchingRows(row, annoReader)
    if (len(annoRows) == 0):
        builtRow = pa.Series()
    else:
        twoMatchingRows = filter.run(annoRows, parentDirectory)
        builtRow = buildFinalTemplate(twoMatchingRows, row)
    return builtRow


def compareKeysOfFileAndReturnMatchingRows(row, annoReader):
    chrStartPosKey = formatChrPosKey(row)

    resultdf = annoReader[annoReader['Location'].str.contains(chrStartPosKey)].drop_duplicates()


    if len(resultdf) == 0:
        logMissedPosition(row, chrStartPosKey)
    return resultdf


def formatChrPosKey(row):
    chr = IOutilities.formatChromo(getFromRow(row, "chromosome"))
    seqStart = getFromRow(row, "seq_start_position")
    ref = getFromRow(row, "ref_allele")
    alt = getFromRow(row, "alt_allele")

    if len(ref) > 0 and len(alt) > 0 and (ref[0] == alt[0]):
        adjustedSeq = (str)((int)(seqStart) + 1)
        if len(ref) == 1:
            chrPosKey = "{0}-{2}".format(chr,seqStart, adjustedSeq)
        else :
            chrPosKey = "{0}:{1}".format(chr, adjustedSeq)
    else :
        chrPosKey ="{0}:{1}".format(chr, seqStart)

    return chrPosKey


def extraColumnToJSON(extra):
    semiToComma = re.sub(';', '","', extra)
    equalsToColon = re.sub("=", '":"', semiToComma)
    addCurlyToStart = re.sub('(?m)^', '{"', equalsToColon)
    JSONstr = re.sub('(?m)$', '"}', addCurlyToStart)

    return json.loads(JSONstr) if extra != "" else pa.Series()


def buildHeaders():
    return ["model_id", "sample_id", "sample_origin", "host_strain_nomenclature", "passage", "symbol", "biotype",
            "coding_sequence_change",
            "variant_class", "codon_change", "amino_acid_change", "consequence", "functional_prediction", "read_depth",
            "allele_frequency",
            "chromosome", "seq_start_position", "ref_allele", "alt_allele", "ucsc_gene_id", "ncbi_gene_id",
            "ncbi_transcript_id", "ensemble_gene_id",
            "ensemble_transcript_id", "variation_id", "genome_assembly", "platform"]


def buildFinalTemplate(twoMatchingRows, row):
    annoRow = twoMatchingRows.iloc[0]
    NCBIrow = twoMatchingRows.iloc[1]

    extra = getFromRow(annoRow, 'Extra')
    extraAnno = extraColumnToJSON(extra)

    return [getFromRow(row, 'Model_ID'), getFromRow(row, 'Sample_ID'), getFromRow(row, 'sample_origin'),
            getFromRow(row, 'host strain nomenclature'),
            getFromRow(row, 'Passage'), getFromRow(extraAnno, 'SYMBOL'), getFromRow(extraAnno, 'BIOTYPE'),
            parseHGSVc(getFromRow(extraAnno, 'HGVSc')), getFromRow(extraAnno, 'VARIANT_CLASS'),
            getFromRow(annoRow, 'Codons'),
            buildAminoAcidChange(getFromRow(annoRow, 'Amino_acids'), getFromRow(annoRow, 'Protein_position')),
            getFromRow(annoRow, 'Consequence'),
            parseFunctionalPredictions(getFromRow(extraAnno, 'PolyPhen'), getFromRow(extraAnno, 'SIFT')),
            getFromRow(row, 'read_depth'), getFromRow(row, 'Allele_frequency'), getFromRow(row, 'chromosome'),
            getFromRow(row, 'seq_start_position'),
            getFromRow(row, 'ref_allele'), getFromRow(row, 'alt_allele'), getFromRow(row, 'ucsc_gene_id'),
            getFromRow(NCBIrow, 'Gene'),
            getFromRow(NCBIrow, 'Feature'), getFromRow(annoRow, 'Gene'), getFromRow(annoRow, 'Feature'),
            getFromRow(annoRow, 'Existing_variation'),
            getFromRow(row, 'genome_assembly'), getFromRow(row, 'Platform')]


def getFromRow(row, attributeID):
    returnStr = ""
    if row.get(attributeID) and (type(row.get(attributeID)) == str):
        returnStr = row.get(attributeID)

    return returnStr


def parseHGSVc(HGSV):
    regexToRemoveAccession = "(?m)c\\.(.+$)"
    hgsvMatch = re.findall(regexToRemoveAccession, HGSV)
    return hgsvMatch[0] if len(hgsvMatch) > 0 else ""


def buildAminoAcidChange(aminoAcids, protienPosition):
    return aminoAcids[0] + protienPosition + aminoAcids[2] if (
            len(aminoAcids) == 3 and protienPosition.isdecimal()) else ""


def parseFunctionalPredictions(polyphen, sift):
    return "PolyPhen: {0} | SIFT: {1}".format(polyphen, sift) if (polyphen and sift) else ""


def logMissedPosition(row, chrStartPosKey):
    global mergedPointsMissed
    mergedPointsMissed += 1

    message = "Total dropped: {0} could not find {1} in annotations. Ref: {2} Alt {3} for sample {4}".format(
        mergedPointsMissed,
        chrStartPosKey,
        row["ref_allele"],
        row["alt_allele"], row['Sample_ID'])
    IOutilities.logMessage(os.path.dirname(TSVfilePath), message)

def filterrun(annoRows, parentDirectory):
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

    condensedRows = NCBIrows.apply(lambda x: calculateMatchingScore(x,row, symbol, biotype, impact, isCanonical), axis=1)
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

if len(sys.argv) > 1:
    run()
