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

import AnnotationFilter
import IOutilities

mergedPointsMissed = 0

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

def isColumnHeader(line):
    return ((line[0] == '#') and (line[1] != '#'))


def isCanonical(line):
    return (line.find("CANONICAL=YES") != -1)


def mergeRowsAndWrite():
    with open(tsvFilePath + ".hmz", 'w') as finalTemplate, \
            open(tsvFilePath + ".ANN", 'r') as annoFile, \
            open(tsvFilePath, 'r') as tsvFile:


        outFileWriter = csv.writer(finalTemplate, delimiter="\t")
        if tsvFilePath.endswith(".tsv"):
            reader = csv.DictReader(tsvFile, delimiter="\t")
        elif tsvFilePath.endswith(".csv"):
            reader = csv.DictReader(tsvFile, delimiter=",")

        annoReader = pa.read_csv(annoFile, delimiter='\t', error_bad_lines=False, header=97)

        message = "Merging original data : {0} /n and annotated data : {1} at {2}".format(tsvFilePath, annoFile,
                                                                                          time.ctime())
        logging.info(message)

        headers = buildHeaders()
        outFileWriter.writerow(headers)

        rowNum = 0
        rowAdded = 0

        for row in reader:
            rowNum += 1

            print(".")

            if rowIsValidForMerge(row):
                mergedRow = mergeRows(row, annoReader)
                if not len(mergedRow) < 26:
                    outFileWriter.writerow(mergedRow)
                    rowAdded += 1
                else:
                    message = ("Info: Dropping row for being invalid (size or column headers) or missing match in annotations (Chromosome position error). RowNum {0} - Len {1}".format(rowNum,
                                                                                                             len(
                                                                                                                 mergedRow)))
                    print(row)
                    logging.warn(message)
            else:

                message2 = ("Info: row {0} is broken or legacy".format(rowNum))
                print(row)
                logging.warn(message2)

        message3 = "{0} The completed file file {1} has {2} data points (including header)".format(time.ctime(),
                                                                                                  finalTemplate,
                                                                                                  rowAdded)
        logging.info(message3)


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
        twoMatchingRows = AnnotationFilter.run(annoRows, tsvFileName, parentDirectory)
        builtRow = buildFinalTemplate(twoMatchingRows, row)
    return builtRow


def compareKeysOfFileAndReturnMatchingRows(row, annoReader):
    chrStartPosKey = formatChrPosKey(row)

    resultdf = annoReader[annoReader['Location'].str.contains(chrStartPosKey)].drop_duplicates()

    if len(resultdf) == 0:
        logMissedPosition(row, chrStartPosKey)
    return resultdf


def formatChrPosKey(row):
    formatedchr = IOutilities.formatChromo(getFromRow(row, "chromosome"))
    seqStart = getFromRow(row, "seq_start_position")
    ref = getFromRow(row, "ref_allele")
    alt = getFromRow(row, "alt_allele")

    if len(ref) > 0 and len(alt) > 0 :
        adjustedSeq = (str)((int)(seqStart) + 1)
        if(ref[0] == alt[0]):
            if len(ref) == 1:
                chrPosKey = "{0}:{1}-{2}".format(formatedchr, seqStart, adjustedSeq)
            else:
                chrPosKey = "{0}:{1}".format(formatedchr, adjustedSeq)
        elif(ref[0] == '-'):
            chrPosKey = "{0}:{1}-{2}".format(formatedchr, seqStart, adjustedSeq)
            logging.warn( "Attemping to adjust for improper insertion format {}".format(chrPosKey))
        else :
            chrPosKey = "{0}:{1}".format(formatedchr, seqStart)
    else:
        chrPosKey = "{0}:{1}".format(formatedchr, seqStart)

    return chrPosKey


def extraColumnToJSON(extra):
    if extra:
        semiToComma = re.sub(';', '","', extra)
        equalsToColon = re.sub("=", '":"', semiToComma)
        addCurlyToStart = re.sub('(?m)^', '{"', equalsToColon)
        JSONstr = re.sub('(?m)$', '"}', addCurlyToStart)
        extraJson = json.loads(JSONstr) if extra != "" else pa.Series()
    else :
        extraJson = pa.Series()
    return extraJson

def buildHeaders():
    return ["model_id", "sample_id", "sample_origin", "host_strain_nomenclature", "passage", "symbol", "biotype",
            "coding_sequence_change",
            "variant_class", "codon_change", "amino_acid_change", "consequence", "functional_prediction", "read_depth",
            "allele_frequency",
            "chromosome", "seq_start_position", "ref_allele", "alt_allele", "ucsc_gene_id", "ncbi_gene_id",
            "ncbi_transcript_id", "ensemble_gene_id",
            "ensemble_transcript_id", "variation_id", "genome_assembly", "platform"]


def parseFilteredRows(twoMatchingRows):
    annoRow = pa.DataFrame()
    NCBIrow = pa.DataFrame()
    if len(twoMatchingRows) > 1:
        annoRow = twoMatchingRows.iloc[0]
        NCBIrow = twoMatchingRows.iloc[1]
    elif len(twoMatchingRows) == 1:
        annoRow = twoMatchingRows.iloc[0]
    return [annoRow, NCBIrow]


def isEnsemblData(row):
    geneId = getFromRow(row, 'Gene')
    transcriptId = getFromRow(row, 'Feature')
    return re.match("ENS", geneId) and re.match("ENS", transcriptId)

def buildFinalTemplate(twoMatchingRows, row):
    NCBIrow = pa.DataFrame()
    builtRow = []
    if len(twoMatchingRows) > 0:
        parsedRows = parseFilteredRows(twoMatchingRows)
        if parsedRows[0].size > 0 and isEnsemblData(parsedRows[0]):
            annoRow = parsedRows[0]
            if parsedRows[1].size > 0 and not isEnsemblData(parsedRows[1]):
                NCBIrow = parsedRows[1]

            extra = getFromRow(annoRow, 'Extra')
            if extra == None:
                extra = ""
            extraAnno = extraColumnToJSON(extra)

            builtRow = [getFromRow(row, 'model_id'), getFromRow(row, 'sample_id'), getFromRow(row, 'sample_origin'),
                        getFromRow(row, 'host_strain_nomenclature'),
                        getFromRow(row, 'passage'), getFromRow(extraAnno, 'SYMBOL'), getFromRow(extraAnno, 'BIOTYPE'),
                        parseHGSVc(getFromRow(extraAnno, 'HGVSc')), getFromRow(extraAnno, 'VARIANT_CLASS'),
                        getFromRow(annoRow, 'Codons'),
                        buildAminoAcidChange(getFromRow(annoRow, 'Amino_acids'), getFromRow(annoRow, 'Protein_position')),
                        getFromRow(annoRow, 'Consequence'),
                        parseFunctionalPredictions(getFromRow(extraAnno, 'PolyPhen'), getFromRow(extraAnno, 'SIFT')),
                        getFromRow(row, 'read_depth'), getFromRow(row, 'Allele_frequency'), getFromRow(row, 'chromosome'),
                        getFromRow(row, 'seq_start_position'),
                        getFromRow(row, 'ref_allele'), getFromRow(row, 'alt_allele'), getFromRow(row, 'ucsc_gene_id'),
                        getFromRow(NCBIrow, 'Gene'),
                        getFromRow(NCBIrow, 'Feature'), getFromRow(annoRow, 'Gene'),
                        getFromRow(annoRow, 'Feature'),
                        getFromRow(annoRow, 'Existing_variation'),
                        getFromRow(row, 'genome_assembly'), getFromRow(row, 'platform')]

    else:
        logging.info("No annotations found for row with values : {}".format(row.values()))
        builtRow = list()

    return builtRow


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
    logging.warn(message)


if len(sys.argv) > 1:
    run()
