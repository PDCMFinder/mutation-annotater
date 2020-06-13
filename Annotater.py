#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import os
import csv
import subprocess as sp
import sys
import re
import logging
import pandas as pa

import IOutilities
import vcfSorter

if len(sys.argv) > 1:
    file = sys.argv[1]
    fileName = os.path.basename(file)
    parentDirectoryPath = os.path.dirname(file)
    provider = os.path.dirname(parentDirectoryPath)
    Updog = os.path.dirname(provider)
    vcfFilePath = file + '.vcf'
    logging.basicConfig(filename='{}.log'.format(file), filemode='a+', level=logging.DEBUG)
    logging.info(" Starting annotation pipleline ")

else :
    logging.info("Please pass the absolute path of the file to annotate")

def run():
    formatToVCFAndSave(file)
    annotateVCF(vcfFilePath,file)
    logging.info("Annotating is complete")

def formatToVCFAndSave(filePath):
    tsvOrCsvFile = open(filePath)
    vcfFile = open(vcfFilePath, "w+")
    IOutilities.chmodFile(vcfFilePath)

    if filePath.endswith(".csv"):
        reader = csv.DictReader(tsvOrCsvFile, delimiter=",")
    else :
        reader = csv.DictReader(tsvOrCsvFile, delimiter="\t")
        if not filePath.endswith(".tsv"):
            logging.info("File {}  is not suffixed as tsv... procceding anyways".format(filePath))

    print("Writing {0} to VCF".format(filePath))
    vcfFile.write("#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfo\n")
    rowCount = 0
    try:
        for row in reader:
            rowCount += 1
            proccessRowToVCF(row, vcfFile)
    except EOFError:
        logging.info("End of file at {}".format(rowCount))

    IOutilities.flushCloseFile(tsvOrCsvFile)
    IOutilities.flushCloseFile(vcfFile)
    vcfSorter.sort(vcfFilePath, vcfFilePath)
    dropDuplicates(vcfFilePath)

    message = "The file {0} has {1} data points (including header)".format(filePath, rowCount)
    logging.info(message)

def dropDuplicates(vcfFile):
    with open(vcfFile, 'r') as vcf:
        vcfDf = pa.read_csv(vcf, sep='\t')
        vcfDf.drop_duplicates(inplace=True)
        vcfDf.to_csv(vcfFile, sep='\t', index=False)

def proccessRowToVCF(row, vcfFile) :
    hg38RE = "(?i)(hg38|grch38|38)"
    if not bool(re.match(hg38RE, row["genome_assembly"])):
        logging.warning("Warning found legacy data : {0}".format(row.items()))
    elif (anyGenomicCoordinateAreMissing(row)) :
        logging.info(
            "Row has incomplete data : {0} in file {1} caused by missing chro,seq start, ref or alt allele data"
                .format(row.items(), vcfFilePath))
    elif allGenomicDataIsMissing(row):
        raise EOFError
    else:
        formatRowToVCFAndWrite(row, vcfFile)


def formatRowToVCFAndWrite(row, vcfFile) :
    chromo = IOutilities.formatChromo(row["chromosome"])
    alleles = formatImproperInserions(row["ref_allele"],row["alt_allele"])
    posId =  createPosId(row)
    vcfRow = "{0}\t{1}\t{2}\t{3}\t{4}\t.\t.\t.\t\n".format(chromo, posId, row["seq_start_position"],
                                                         alleles[0], alleles[1])
    vcfFile.write(vcfRow)

def formatImproperInserions(refAllele, altAllele) :
    if not refAllele[0] == " " and refAllele[0] == "-":
        formatedRefAllele = "A"
        formatedAltAllele = "A{}".format(altAllele)
    elif '-' in refAllele:
        logging.info("Ref allele {} not supported".format(refAllele))
    else:
        formatedRefAllele = refAllele
        formatedAltAllele = altAllele
    return [formatedRefAllele,formatedAltAllele]

def createPosId(row):
    return "{}_{}_{}_{}".format(row["chromosome"], row["seq_start_position"],row["ref_allele"],row["alt_allele"])

def anyGenomicCoordinateAreMissing (row) :
    return not row["chromosome"] or not row["seq_start_position"] or not row["ref_allele"] or not row["alt_allele"]

def allGenomicDataIsMissing (row) :
    return not row["chromosome"] and not row["seq_start_position"] and not row["ref_allele"] and not row["alt_allele"]

def annotateVCF(vcfFile, targetFile):
    fastaDir = "/nfs/nobackup/spot/mouseinformatics/pdx/vepDBs/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    alleleDB = "/nfs/nobackup/spot/mouseinformatics/pdx/vepDBs/homo_sapiens_vep_98_GRCh38"
    singularityVepImage = "./ensembl-vep.simg"

    if not os.path.isfile(fastaDir) and not os.path.isfile(alleleDB):
        fastaDir = "/home/afollette/vepWD/db/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        alleleDB = "/home/afollette/vepWD/db/homo_sapiens_vep_98_GRCh38"

    if not os.path.isfile(singularityVepImage):
        singularityVepImage = "/homes/afollette/ensembl-vep.simg"

    vepIn = vcfFile
    vepWarningFile = targetFile + ".vepWarnings"
    vepOut = targetFile + ".ANN"

    vepCMD = """vep -e -q -check_existing  -symbol -polyphen -sift -merged -use_transcript_ref —hgvs —hgvsg —variant_class \
    -canonical -fork 4 -format vcf -force -offline -no_stats --warning_file {0} -vcf\
     -cache -dir_cache {1} -fasta {2} -i {3} -o {4} 2>> {5}.log""".format(vepWarningFile,alleleDB,fastaDir,vepIn, vepOut,fileName)
    logging.debug("singularity exec {0} {1}".format(singularityVepImage, vepCMD))
    returnSignal = sp.call(
        "singularity exec {0} {1}".format(singularityVepImage, vepCMD), shell=True)
    if(returnSignal != 0):
        raise Exception("Vep returned a non-zero exit code {}".format(returnSignal))

if len(sys.argv) > 1:
    run()
