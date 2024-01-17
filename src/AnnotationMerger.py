#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import csv
import glob
import os
import sys
import time
import logging
import re
import pandas as pd
from math import isnan


class AnnotationMerger:
    global mergedPointsMissed
    mergedPointsMissed = 0

    def __init__(self, mutTarget, run_type, local):
        self.parentDirectory = os.path.dirname(mutTarget)
        self.annotationFilePath = "{}.ANN".format(os.path.join(self.parentDirectory, 'mut/annotations/merged'))
        self.provider = os.path.dirname(self.parentDirectory)
        self.Updog = os.path.dirname(self.provider)
        self.run_type = run_type
        self.local = local
        self.read_annotation_file()

    def run(self, mutTarget):
        self.tsvFilePath = mutTarget
        self.tsvFileName = os.path.basename(self.tsvFilePath)
        self.outFilePath = "{}.hmz".format(self.tsvFilePath)
        self.mergeRowsAndWrite()
        #logging.info("Merge complete")

    def read_annotation_file(self):
        start = time.time()
        self.annoReader = pd.read_csv(self.annotationFilePath, delimiter='\t', low_memory=False)
        self.annoReader.columns = self.annoReader.columns.str.lower()
        self.annoReader = self.annoReader.apply(lambda x: self.generate_annotation_columns(x), axis=1)
        mapper = {'codons': 'codon_change', 'pos': 'seq_start_position', 'ref': 'ref_allele',
                  'alt': 'alt_allele', 'existing_variation': 'variation_id'}
        self.annoReader.rename(columns=mapper, inplace=True)
        cols = ['id', 'symbol', 'biotype', 'coding_sequence_change', 'variant_class', 'codon_change',
                'amino_acid_change', 'consequence', 'functional_prediction', 'strand', 'seq_start_position',
                'ref_allele', 'alt_allele', "ncbi_gene_id", "ncbi_transcript_id", "ensembl_gene_id",
                'ensembl_transcript_id', "variation_id"]
        self.annoReader = self.annoReader[cols]
        logging.info("{0}: Annotation file processed in {1} mins!".format(time.ctime(), round((time.time()-start)/60)))

    def mergeRowsAndWrite(self):
        self.iterateThroughRowsAndMerge()#tsvReader, outFileWriter)

    def generate_annotation_columns(self, row):
        row['amino_acid_change'] = self.buildAminoAcidChange(row['amino_acids'], row['protein_position'])
        row['functional_prediction'] = self.parseFunctionalPredictions(row['polyphen'], row['sift'])
        row['coding_sequence_change'] = self.parseHGSVc(row['hgvsc'])
        ncbi_gene_id = ''
        ncbi_transcript_id = ''
        ensembl_gene_id = ''
        ensembl_transcript_id = ''
        source = row['source']
        if isinstance(source, str) and pd.notna(source):
            if str(source).lower().__contains__('ensembl'):
                ensembl_gene_id = row['gene']
                ensembl_transcript_id = row['feature']
            elif str(source).lower().__contains__('refseq'):
                ncbi_gene_id = row['gene']
                ncbi_transcript_id = row['feature']
        row['ncbi_transcript_id'] = ncbi_transcript_id
        row['ncbi_gene_id'] = ncbi_gene_id
        row['ensembl_gene_id'] = ensembl_gene_id
        row['ensembl_transcript_id'] = ensembl_transcript_id
        strand = ''
        if not isnan(row['strand']):
            strand = int(row['strand'])
        row['strand'] = strand
        row['chromosome'] = row['#chrom'].replace('chr', '')
        return row

    def iterateThroughRowsAndMerge(self):#, reader, outFileWriter):
        mut_raw, out_cols, mut_size = self.process_raw_data()
        annotated = mut_raw.merge(self.annoReader, left_on='annotation_key',
                                               right_on='id', how='left', indicator=True)
        rows_without_match = annotated[annotated['_merge'] == 'left_only']
        annotated = annotated[annotated['_merge'] == 'both']
        annotated[out_cols].to_csv(self.outFilePath, sep='\t', index=False)
        indices_to_drop = rows_without_match.index
        for index in indices_to_drop:
            message = ("Info: Row Dropped at index {0}: "
                   "Dropping row for being invalid or failed to annotated"
                   "data {1}, {2}".format(index, rows_without_match.at[index, 'sample_id'],
                                          rows_without_match.at[index, 'annotation_key']))
            logging.warning(message)
        message3 = ("{0}: Annotated file {1} has {2}"
                    " data points out of {3} attempted".format(time.ctime(), self.tsvFileName + ".ANN", annotated.shape[0],
                                                               mut_size))
        logging.info(message3)

    def process_raw_data(self):
        mut_raw = pd.read_csv(self.tsvFilePath, sep='\t', low_memory=False)
        out_cols = mut_raw.columns
        mut_size = mut_raw.shape[0]
        indices_to_drop = mut_raw[mut_raw[['chromosome', 'seq_start_position']].isna().any(axis=1)].index
        for index in indices_to_drop:
            chromosome = mut_raw.at[index, 'chromosome']
            pos = mut_raw.at[index, 'seq_start_position']
            message2 = ("Info: Row Index {} is missing essential value or legacy. chromosome: {} pos: {} "
                        .format(index, chromosome, pos))
            logging.warning(message2)
        mut_raw = mut_raw.dropna(subset=['chromosome', 'seq_start_position'])
        mut_raw['chromosome'] = mut_raw['chromosome'].astype(str)
        mut_raw['annotation_key'] = mut_raw.fillna('').apply(self.createAnnotationKey, axis=1)
        mut_raw['chromosome'] = mut_raw['chromosome'].str.replace('chr', '')
        annotations = ['sample_id', 'annotation_key', 'platform_id', 'ucsc_gene_id', 'read_depth', 'allele_frequency',
                       'chromosome']
        return mut_raw[annotations], out_cols, mut_size

    def createAnnotationKey(self, row):
        if self.run_type == 'hgvs':
            annotation_key = "{}:c.{}".format(row["ncbi_transcript_id"], row["coding_sequence_change"])
        else:
            annotation_key = "{}_{}_{}_{}".format(self.formatChromo(row["chromosome"]),
                                                  row["seq_start_position"], row["ref_allele"],
                                                    row["alt_allele"])
        return annotation_key

    def formatChromo(self, givenChromo):
        formattedChromo = givenChromo
        incorrectChrFormat = "(?i)^([0-9]{1,2}|[xym]{1}|mt|un)$"
        isMatch = re.match(incorrectChrFormat, str(givenChromo))
        if isMatch:
            formattedChromo = "chr" + isMatch.group(1)
        return formattedChromo


    def parseHGSVc(self, HGSV):
        regexToRemoveAccession = "(?m)(c|n)\\.(.+$)"
        hgsvMatch = re.findall(regexToRemoveAccession, str(HGSV))
        return hgsvMatch[0][1] if len(hgsvMatch) > 0 else ""

    def buildAminoAcidChange(self, aminoAcids, protienPosition):
        aminoAcidChange = ""
        if (bool(aminoAcids) and bool(protienPosition)) and isinstance(aminoAcids, str) and pd.notna(aminoAcids) and len(aminoAcids) == 3:
            try:
                protienPosition = str(int(protienPosition))
            except:
                protienPosition = str(protienPosition)
            aminoAcidChange = aminoAcids[0] + protienPosition + aminoAcids[2]
        return aminoAcidChange

    def parseFunctionalPredictions(self, polyphen, sift):
        functionalPredictions = ""
        if bool(polyphen) and bool(sift) and isinstance(polyphen, str) and pd.notna(polyphen) and isinstance(sift, str) and pd.notna(sift):
            functionalPredictions = "PolyPhen: {0} | Sift: {1}".format(polyphen, sift)
        return functionalPredictions


def cmdline_runner():
    if len(sys.argv) > 1:
        mutTarget = sys.argv[1]
        run_type = sys.argv[2]
        local = sys.argv[3]
        mutDir = os.path.dirname(mutTarget)
        merger = AnnotationMerger(mutDir, run_type, local)
        if os.path.isfile(mutTarget):
            logging.basicConfig(filename='{}.log'.format(mutTarget), filemode='a+', level=logging.DEBUG)
            logging.info("Starting merge of annotations")
            #AnnotationMerger(mutTarget, run_type, local).run()
            merger.run(mutTarget)
        elif os.path.isdir(mutTarget):
            globForTsv = os.path.join(mutTarget, "*tsv")
            for mutFile in glob.iglob(globForTsv):
                mutfile_path = os.path.join(mutTarget, mutFile)
                print(mutfile_path)
                merger.run(mutfile_path)
                #AnnotationMerger(mutfile_path, run_type, local).run()
    else:
        logging.info("Please pass the absolute path of the file to annotate")
#cmdline_runner()
