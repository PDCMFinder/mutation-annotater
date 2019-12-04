#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import filter
import IOutilities

import csv
import json
import os
import sys
import time

import re
import pandas as pa

mergedPointsMissed = 0

if len(sys.argv) > 1:
    TSVfilePath = sys.argv[1]
    parentDirectory = os.path.dirname(sys.argv[1])
else:
    sys.stderr.write("Warning: Merger is being ran without file input. This should only be used for testing")


def run():
    mergeRowsAndWrite()
    IOutilities.saveToExcel(TSVfilePath)


def isColumnHeader(line):
    return ((line[0] == '#') and (line[1] != '#'))


def isCanonical(line):
    return (line.find("CANONICAL=YES") != -1)


def mergeRowsAndWrite():
    with open(TSVfilePath + ".hmz", 'w') as finalTemplate, \
            open(TSVfilePath + ".ANN", 'r') as annoFile, \
            open(TSVfilePath, 'r') as tsvFile:

        outFileWriter = csv.writer(finalTemplate, delimiter="\t")
        if TSVfilePath.endswith(".tsv"):
            reader = csv.DictReader(tsvFile, delimiter="\t")
        elif TSVfilePath.endswith(".csv"):
            reader = csv.DictReader(tsvFile, delimiter=",")
        annoReader = pa.read_csv(annoFile, delimiter='\t', error_bad_lines=False, header=97)

        message = "Merging original data : {0} /n and annotated data : {1} at {2}".format(TSVfilePath, annoFile,
                                                                                          time.ctime())
        IOutilities.logMessage(parentDirectory, message)

        headers = buildHeaders()
        outFileWriter.writerow(headers)

        rowNum = 0
        rowAdded = 0

        for row in reader:
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
                    print(row)
                    IOutilities.logMessage(parentDirectory, message)
            else:

                message2 = ("Info: row {0} is broken or legacy".format(rowNum))
                print(row)
                IOutilities.logMessage(parentDirectory, message2)

        masterLog = "{0}/masterLog".format(os.path.dirname(os.path.dirname(parentDirectory)))
        message = "{0} The completed file file {1} has {2} data points (including header)".format(time.ctime(),
                                                                                                  finalTemplate,
                                                                                                  rowAdded)
        IOutilities.masterlogMessage(masterLog, message)


def rowIsValidForMerge(row):
    return rowIsHg38(row) and getFromRow(row, "chromosome") and getFromRow(row, "seq_start_position")


def rowIsHg38(row):
    hg38Regex = "(?i)(hg38|GRCh38|38)"
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
            chrPosKey = "{0}:{1}-{2}".format(chr, seqStart, adjustedSeq)
        else:
            chrPosKey = "{0}:{1}".format(chr, adjustedSeq)
    else:
        chrPosKey = "{0}:{1}".format(chr, seqStart)

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
    inputLen = len(twoMatchingRows)
    annoRow = pa.DataFrame()
    emblGeneColumnName = 'Gene'
    emblFeatureColumnName = 'Feature'

    if inputLen > 0:

        if inputLen > 1:
            annoRow = twoMatchingRows.iloc[0]
            NCBIrow = twoMatchingRows.iloc[1]
        elif inputLen == 1:
            annoRow = twoMatchingRows.iloc[0]
            reGene = getFromRow(twoMatchingRows, 'Gene')
            reFeature = getFromRow(twoMatchingRows, 'Feature')

            if re.search("^ENS",reGene) and (re.search("^ENS",reFeature)):
                NCBIrow = pa.DataFrame()
            else:
                emblGeneColumnName = 'Gene'
                emblFeatureColumnName = 'Feature'

        extra = getFromRow(annoRow, 'Extra')
        extraAnno = extraColumnToJSON(extra)

        builtRow = [getFromRow(row, 'Model_ID'), getFromRow(row, 'Sample_ID'), getFromRow(row, 'sample_origin'),
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
                    getFromRow(NCBIrow, 'Feature'), getFromRow(annoRow, emblGeneColumnName),
                    getFromRow(annoRow, emblFeatureColumnName),
                    getFromRow(annoRow, 'Existing_variation'),
                    getFromRow(row, 'genome_assembly'), getFromRow(row, 'Platform')]

    else:
        builtRow = list()

    return builtRow


def getFromRow(row, attributeID):
    returnStr = ""
    if row.get(attributeID) and ((type(row.get(attributeID)) == str) or type(row.get(attributeID)) == unicode):
        returnStr = row.get(attributeID)

    return returnStr


def parseHGSVc(HGSV):
    regexToRemoveAccession = "(?m)c\\.(.+$)"
    hgsvMatch = re.findall(regexToRemoveAccession, HGSV)
    return hgsvMatch[0] if len(hgsvMatch) > 0 else ""


def buildAminoAcidChange(aminoAcids, protienPosition):
    return aminoAcids[0] + protienPosition + aminoAcids[2] if (
            len(aminoAcids) == 3 and protienPosition.isdigit()) else ""


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


if len(sys.argv) > 1:
    run()
