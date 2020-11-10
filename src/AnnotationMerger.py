#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import csv
import os
import sys
import time
import logging
import re
import pandas as pa

if len(sys.argv) > 1:
    tsvFilePath = sys.argv[1]
    tsvFileName = os.path.basename(tsvFilePath)
    annotationFilePath = "{}.ANN".format(tsvFilePath)
    outFilePath = "{}.hmz".format(tsvFilePath)
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
    with open(outFilePath, 'w+') as outFile, \
            open(tsvFilePath, 'r') as tsvFile:
        outFileWriter = csv.writer(outFile, delimiter="\t")
        tsvReader = csv.DictReader(tsvFile, delimiter="\t")
        headers = buildHeaders()
        outFileWriter.writerow(headers)
        logBeginningOfMerge()
        annoReader = pa.read_csv(annotationFilePath, delimiter='\t', error_bad_lines=False)

        iterateThroughRowsAndMerge(tsvReader, annoReader, outFileWriter)


def logBeginningOfMerge():
    message = ("Merging original data :"
               " {0} /n and annotated data : {1} at {2}".format(tsvFilePath, annotationFilePath, time.ctime()))
    logging.info(message)


def iterateThroughRowsAndMerge(reader, annoReader, outFileWriter):
    rowAdded = 0
    index = 0
    for row in reader:
        rowAdded = validateMergeAndWrite(index, row, annoReader, outFileWriter, rowAdded)
        index += 1
    message3 = ("{0} The completed file file {1} has {2}"
                " data points out of {3} attempted".format(time.ctime(), tsvFileName + ".ANN", rowAdded, index))
    logging.info(message3)


def validateMergeAndWrite(index, row, annoReader, outFileWriter, rowAdded):
    if checkRowForValidity(index, row):
        mergedRow = mergeRows(row, annoReader)
        if len(mergedRow) != 0:
            outFileWriter.writerow(mergedRow)
            rowAdded += 1
        else:
            message = ("Info: Row Dropped at index {0}: "
                       "Dropping row for being invalid or failed to annotated"
                       "data {1}".format(index, mergedRow.to_string()))
            logging.warning(message)
    return rowAdded


def checkRowForValidity(index, row):
    isValid = False
    assembly = getFromRow(row, "genome_assembly")
    chromosome = getFromRow(row, "chromosome")
    pos = getFromRow(row, "chromosome")
    if strIsHg38(assembly) and chromosome and pos:
        isValid = True
    else:
     message2 = ("Info: Row Index {} is missing essential value or legacy. Assembly: {} chromosome: {} pos: {} "
                 .format(index,assembly, chromosome,pos))
     logging.warning(message2)
    return isValid

def strIsHg38(assembly):
    hg38Regex = "(?i)(hg38|GRCh38|38)"
    return re.match(hg38Regex, assembly)


def mergeRows(row, annoReader):
    annoRow = returnMatchingRows(row, annoReader)
    if len(annoRow) == 0:
        builtRow = pa.Series()
    else:
        builtRow = buildRow(row, annoRow)
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
    attribute = ""
    if type(row) is pa.DataFrame:
        rawAttribute = row.get(attributeID).iloc[0]
    else:
        rawAttribute = row.get(attributeID)
    if bool(rawAttribute) and str(rawAttribute) != 'nan':
        attribute = str(rawAttribute)
    return attribute


def parseHGSVc(HGSV):
    regexToRemoveAccession = "(?m)c\\.(.+$)"
    hgsvMatch = re.findall(regexToRemoveAccession, str(HGSV))
    return hgsvMatch[0] if len(hgsvMatch) > 0 else ""


def buildAminoAcidChange(aminoAcids, protienPosition):
    aminoAcidChange = ""
    if (bool(aminoAcids) and bool(protienPosition)) and len(aminoAcids) == 3:
        aminoAcidChange = aminoAcids[0] + protienPosition + aminoAcids[2]
    return aminoAcidChange


def parseFunctionalPredictions(polyphen, sift):
    functionalPredictions = ""
    if bool(polyphen) and bool(sift):
        functionalPredictions = "PolyPhen: {0} | Sift: {1}".format(polyphen, sift)
    return functionalPredictions


def isColumnHeader(line):
    return ((line[0] == '#') and (line[1] != '#'))


def logMissedPosition(row, chrStartPosKey):
    global mergedPointsMissed
    mergedPointsMissed += 1

    message = "Total dropped: {0} could not find {1} in annotations:".format(mergedPointsMissed, chrStartPosKey)
    logging.warning(message)
    logging.warning(row.to_string())


if len(sys.argv) > 1:
    run()
