#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import csv
import json
import os
import sys
import time
import logging
import re
import pandas as pa

if len(sys.argv) > 1:
    tsvFilePath = sys.argv[1]
    tsvFileName = os.path.basename(tsvFilePath)
    parentDirectory = os.path.dirname(sys.argv[1])
    provider = os.path.dirname(parentDirectory)
    Updog = os.path.dirname(provider)

    logging.basicConfig(filename='{}.log'.format(tsvFilePath), filemode='a+', level=logging.DEBUG)
    logging.info(" Beginning to compile annotation with provider data ")
else:
    sys.stderr.write(" Warning: Merger is being ran without file input. This should only be used for testing ")


def run():
    mergeRowsAndWrite()
    logging.info("Merge complete")

def mergeRowsAndWrite():
    with open(tsvFilePath + ".hmz", 'w') as finalTemplate, \
            open(tsvFilePath + ".ANN", 'r') as annoFile, \
            open(tsvFilePath, 'r') as tsvFile:
        outFileWriter = csv.writer(finalTemplate, delimiter="\t")
        tsvReader = getTsvReader(tsvFile)
        annoReader = pa.read_csv(annoFile, delimiter='\t', error_bad_lines=False)
        headers = buildHeaders()
        outFileWriter.writerow(headers)
        logBeginningOfMerge(tsvFilePath)
        iterateThroughRowsAndMerge(tsvReader,annoReader, outFileWriter)

def getTsvReader(tsvFile):
    if tsvFilePath.endswith(".csv"):
        tsvReader = csv.DictReader(tsvFile, delimiter=",")
    else:
        tsvReader = csv.DictReader(tsvFile, delimiter="\t")
    return tsvReader

def logBeginningOfMerge(tsvFilePath):
    message = ("Merging original data :"
               " {0} /n and annotated data : {1} at {2}".format(tsvFilePath, tsvFilePath + "ANN", time.ctime()))
    logging.info(message)

def iterateThroughRowsAndMerge(reader, annoReader, outFileWriter):
    rowNum = 0
    rowAdded = 0
    for row in reader:
        rowNum += 1
        if rowIsValidForMerge(row):
            mergedRow = mergeRows(row, annoReader, rowNum)
            if len(mergedRow) != 0:
                outFileWriter.writerow(mergedRow)
                rowAdded += 1
            else:
                message = ("Info: Row Number {0}: "
                           "Dropping row for being invalid (size or column headers) "
                           "or missing match in annotations (Chromosome position error). "
                           "Row: {0} - Length {1} - data {2}".format(rowNum, len(mergedRow), mergedRow.to_string()))
                logging.warning(message)
        else:
            message2 = ("Info: Row Number {0}: is broken or legacy".format(rowNum))
            logging.warning(message2)
            logging.warning(row.to_string())
    message3 = ("{0} The completed file file {1} has {2}"
                " data points out of {3}".format(time.ctime(), tsvFileName + ".ANN", rowAdded, rowNum))
    logging.info(message3)

def rowIsValidForMerge(row):
    return rowIsHg38(row) and getFromRow(row, "chromosome") and getFromRow(row, "seq_start_position")


def rowIsHg38(row):
    hg38Regex = "(?i)(hg38|GRCh38|38)"
    return re.match(hg38Regex, getFromRow(row, "genome_assembly"))


def mergeRows(row, annoReader, rowNum):
    annoRow = returnMatchingRows(row, annoReader)
    if len(annoRow) == 0:
        builtRow = pa.Series()
    else:
        builtRow = buildRow(annoRow, row, rowNum)
    return builtRow


def returnMatchingRows(row, annoReader):
    annotationKey = createAnnotationKey(row)
    resultdf = annoReader[annoReader['id'] == annotationKey]
    if len(resultdf) == 0:
        logMissedPosition(row, annotationKey)
    return resultdf

def createAnnotationKey(row):
    return "{}_{}_{}_{}".format(formatChromo(row["chromosome"]),
                                row["seq_start_position"], row["ref_allele"], row["alt_allele"])

def formatChromo(givenChromo):
    formattedChromo = givenChromo
    incorrectChrFormat = "(?i)^([0-9]{1,2}|[xym]{1}|mt|un)$"
    isMatch = re.match(incorrectChrFormat, givenChromo)
    if isMatch:
        formattedChromo = "chr" + isMatch.group(1)
    return formattedChromo

def buildHeaders():
    return ["model_id", "sample_id", "sample_origin", "host_strain_nomenclature", "passage", "symbol", "biotype",
            "coding_sequence_change",
            "variant_class", "codon_change", "amino_acid_change", "consequence", "functional_prediction", "read_depth",
            "allele_frequency",
            "chromosome", "seq_start_position", "ref_allele", "alt_allele", "ucsc_gene_id", "ncbi_gene_id",
            "ncbi_transcript_id", "ensembl_gene_id",
            "ensembl_transcript_id", "variation_id", "genome_assembly", "platform"]

def buildRow(row, annoRow):
    return [getEitherFromRow(row, 'model_id', 'Model_ID'), getEitherFromRow(row, 'sample_id', 'Sample_ID'),
            getFromRow(row, 'sample_origin'),
            getEitherFromRow(row, 'host_strain_nomenclature', 'host strain nomenclature'),
            getEitherFromRow(row, 'passage', 'Passage'), getFromRow(annoRow, 'SYMBOL', ),
            getFromRow(annoRow, 'BIOTYPE'), parseHGSVc(getFromRow(annoRow, 'HGVSc')),
            getFromRow(annoRow, 'VARIANT_CLASS'), getFromRow(annoRow, 'Codons'),
            buildAminoAcidChange(getFromRow(annoRow, 'Amino_acids'),
                                             getFromRow(annoRow, 'Protein_position')),
            getFromRow(annoRow, 'Consequence'),
            parseFunctionalPredictions(getFromRow(annoRow, 'PolyPhen'), getFromRow(annoRow, 'SIFT')),
            getFromRow(row, 'read_depth'), getEitherFromRow(row, 'Allele_frequency', 'allele_frequency'),
            getFromRow(row, 'chromosome'),
            getFromRow(annoRow, 'pos'),
            getFromRow(annoRow, 'ref'), getFromRow(annoRow, 'alt'), getFromRow(row, 'ucsc_gene_id'),
                        "",
                        "", getFromRow(annoRow, 'Gene'),
            getFromRow(annoRow, 'Feature'),
            getFromRow(annoRow, 'Existing_variation'),
            getFromRow(row, 'genome_assembly'), getEitherFromRow(row, 'platform', 'Platform')]

def getEitherFromRow(row, attributeId, alternativeId):
    returnStr = getFromRow(row, attributeId)
    if not returnStr or returnStr == "":
        returnStr = getFromRow(row, alternativeId)
    return returnStr

def getFromRow(row, attributeID):
    returnStr = ""
    attribute = str(row.get(attributeID))
    attributeIsStrOrUnicode = (type(attribute) == str or type(attributeID) == unicode)
    if attribute and attributeIsStrOrUnicode:
        returnStr = row.get(attributeID)
    return returnStr

def parseHGSVc(HGSV):
    regexToRemoveAccession = "(?m)c\\.(.+$)"
    hgsvMatch = re.findall(regexToRemoveAccession, str(HGSV))
    return hgsvMatch[0] if len(hgsvMatch) > 0 else ""


def buildAminoAcidChange(aminoAcids, protienPosition):
    return aminoAcids[0] + protienPosition + aminoAcids[2] if (aminoAcids and protienPosition and
                                                               len(
                                                                   aminoAcids) == 3) else ""

def parseFunctionalPredictions(polyphen, sift):
    return "PolyPhen: {0} | SIFT: {1}".format(polyphen, sift) if (polyphen and sift) else ""


def isColumnHeader(line):
    return ((line[0] == '#') and (line[1] != '#'))


def logMissedPosition(row, chrStartPosKey):
    global mergedPointsMissed
    mergedPointsMissed += 1

    message = "Total dropped: {0} could not find {1} in annotations:".format(mergedPointsMissed,chrStartPosKey)
    logging.warning(message)
    logging.warning(row.to_string())

if len(sys.argv) > 1:
    run()
