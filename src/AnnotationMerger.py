#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import csv
import glob
import os
import sys
import time
import logging
import re
import pandas as pa


class AnnotationMerger:

    def __init__(self, mutTarget):
        self.tsvFilePath = mutTarget
        self.tsvFileName = os.path.basename(self.tsvFilePath)
        self.annotationFilePath = "{}.ANN".format(self.tsvFilePath)
        self.outFilePath = "{}.hmz".format(self.tsvFilePath)
        self.parentDirectory = os.path.dirname(sys.argv[1])
        self.provider = os.path.dirname(self.parentDirectory)
        self.Updog = os.path.dirname(self.provider)

    def run(self):
        self.mergeRowsAndWrite()
        logging.info("Merge complete")


    def mergeRowsAndWrite(self):
        with open(self.outFilePath, 'w+') as outFile, \
                open(self.tsvFilePath, 'r') as tsvFile:
            outFileWriter = csv.writer(outFile, delimiter="\t")
            tsvReader = csv.DictReader(tsvFile, delimiter="\t")
            headers = self.buildHeaders()
            outFileWriter.writerow(headers)
            self.logBeginningOfMerge()
            annoReader = pa.read_csv(self.annotationFilePath, delimiter='\t', error_bad_lines=False)

            self.iterateThroughRowsAndMerge(tsvReader, annoReader, outFileWriter)


    def logBeginningOfMerge(self):
        message = ("Merging original data :"
                   " {0} /n and annotated data : {1} at {2}".format(self.tsvFilePath, self.annotationFilePath, time.ctime()))
        logging.info(message)


    def iterateThroughRowsAndMerge(self, reader, annoReader, outFileWriter):
        rowAdded = 0
        index = 0
        for row in reader:
            rowAdded = self.validateMergeAndWrite(index, row, annoReader, outFileWriter, rowAdded)
            index += 1
        message3 = ("{0} The completed file file {1} has {2}"
                    " data points out of {3} attempted".format(time.ctime(), self.tsvFileName + ".ANN", rowAdded, index))
        logging.info(message3)


    def validateMergeAndWrite(self, index, row, annoReader, outFileWriter, rowAdded):
        if self.checkRowForValidity(index, row):
            mergedRow = self.mergeRows(row, annoReader)
            if len(mergedRow) != 0:
                outFileWriter.writerow(mergedRow)
                rowAdded += 1
            else:
                message = ("Info: Row Dropped at index {0}: "
                           "Dropping row for being invalid or failed to annotated"
                           "data {1}".format(index, mergedRow.to_string()))
                logging.warning(message)
        return rowAdded


    def checkRowForValidity(self, index, row):
        isValid = False
        chromosome = self.getFromRow(row, "chromosome")
        pos = self.getFromRow(row, "chromosome")
        if chromosome and pos:
            isValid = True
        else:
         message2 = ("Info: Row Index {} is missing essential value or legacy. chromosome: {} pos: {} "
                     .format(index,chromosome,pos))
         logging.warning(message2)
        return isValid

    def mergeRows(self, row, annoReader):
        annoRow = self.returnMatchingRows(row, annoReader)
        if len(annoRow) == 0:
            builtRow = pa.Series()
        else:
            builtRow = self.buildRow(row, annoRow)
        return builtRow


    def returnMatchingRows(self, row, annoReader):
        annotationKey = self.createAnnotationKey(row)
        resultdf = annoReader[annoReader['id'] == annotationKey]
        if len(resultdf) == 0:
            self.logMissedPosition(row, annotationKey)
        return resultdf


    def createAnnotationKey(self, row):
        return "{}_{}_{}_{}".format(self.formatChromo(row["chromosome"]),
                                    row["seq_start_position"], row["ref_allele"], row["alt_allele"])


    def formatChromo(self, givenChromo):
        formattedChromo = givenChromo
        incorrectChrFormat = "(?i)^([0-9]{1,2}|[xym]{1}|mt|un)$"
        isMatch = re.match(incorrectChrFormat, givenChromo)
        if isMatch:
            formattedChromo = "chr" + isMatch.group(1)
        return formattedChromo


    def buildHeaders(self):
        return ["sample_id","symbol", "biotype",
                "coding_sequence_change",
                "variant_class", "codon_change", "amino_acid_change", "consequence", "functional_prediction", "read_depth",
                "allele_frequency",
                "chromosome", "seq_start_position", "ref_allele", "alt_allele", "ucsc_gene_id", "ncbi_gene_id",
                "ncbi_transcript_id", "ensembl_gene_id",
                "ensembl_transcript_id", "variation_id", "platform_id"]


    def buildRow(self, row, annoRow):
        return [
                self.getFromRow(row, 'sample_id'),
                self.getFromRow(annoRow, 'SYMBOL', ),
                self.getFromRow(annoRow, 'BIOTYPE'),
                self.parseHGSVc(self.getFromRow(annoRow, 'HGVSc')),
                self.getFromRow(annoRow, 'VARIANT_CLASS'),
                self.getFromRow(annoRow, 'Codons'),
                self.buildAminoAcidChange(self.getFromRow(annoRow, 'Amino_acids'),
                self.getFromRow(annoRow, 'Protein_position')),
                self.getFromRow(annoRow, 'Consequence'),
                self.parseFunctionalPredictions(self.getFromRow(annoRow, 'PolyPhen'),
                self.getFromRow(annoRow, 'SIFT')),
                self.getFromRow(row, 'read_depth'),
                self.getFromRow(row, 'Allele_frequency'),
                self.getFromRow(row, 'chromosome'),
                self.getFromRow(annoRow, 'pos'),
                self.getFromRow(annoRow, 'ref'),
                self.getFromRow(annoRow, 'alt'),
                self.getFromRow(row, 'ucsc_gene_id'),
                "",
                "",
                self.getFromRow(annoRow, 'Gene'),
                self.getFromRow(annoRow, 'Feature'),
                self.getFromRow(annoRow, 'Existing_variation'),
                self.getFromRow(row, 'platform_id')]


    def getFromRow(self, row, attributeID):
        attribute = ""
        if type(row) is pa.DataFrame:
            rawAttribute = row.get(attributeID).iloc[0]
        else:
            rawAttribute = row.get(attributeID)
        if bool(rawAttribute) and str(rawAttribute) != 'nan':
            attribute = str(rawAttribute)
        return attribute


    def parseHGSVc(self, HGSV):
        regexToRemoveAccession = "(?m)c\\.(.+$)"
        hgsvMatch = re.findall(regexToRemoveAccession, str(HGSV))
        return hgsvMatch[0] if len(hgsvMatch) > 0 else ""


    def buildAminoAcidChange(self, aminoAcids, protienPosition):
        aminoAcidChange = ""
        if (bool(aminoAcids) and bool(protienPosition)) and len(aminoAcids) == 3:
            aminoAcidChange = aminoAcids[0] + protienPosition + aminoAcids[2]
        return aminoAcidChange


    def parseFunctionalPredictions(self, polyphen, sift):
        functionalPredictions = ""
        if bool(polyphen) and bool(sift):
            functionalPredictions = "PolyPhen: {0} | Sift: {1}".format(polyphen, sift)
        return functionalPredictions


    def isColumnHeader(self, line):
        return ((line[0] == '#') and (line[1] != '#'))


    def logMissedPosition(self, row, chrStartPosKey):
        global mergedPointsMissed
        mergedPointsMissed += 1

        message = "Total dropped: {0} could not find {1} in annotations:".format(mergedPointsMissed, chrStartPosKey)
        logging.warning(message)
        logging.warning(row.to_string())




def cmdline_runner():
    if len(sys.argv) > 1:
        mutTarget = sys.argv[1]
        if os.path.isfile(mutTarget):
            logging.basicConfig(filename='{}.log'.format(mutTarget), filemode='a+', level=logging.DEBUG)
            logging.info("Starting merge of annotations")
            AnnotationMerger(mutTarget).run()
        elif os.path.isdir(mutTarget):
           globForTsv = os.path.join(mutTarget, "*tsv")
           for mutFile in glob.iglob(globForTsv):
               mutfile_path = os.path.join(mutTarget,mutFile)
               print(mutfile_path)
               AnnotationMerger(mutfile_path).run()
    else:
        logging.info("Please pass the absolute path of the file to annotate")
cmdline_runner()