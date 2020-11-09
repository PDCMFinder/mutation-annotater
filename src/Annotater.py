#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import os
import csv
import subprocess as sp
from multiprocessing import cpu_count
import sys
import re
import logging
import pandas as ps
import yaml

if len(sys.argv) > 1:
    mutFile = sys.argv[1]
    fileName = os.path.basename(mutFile)
    parentDirectoryPath = os.path.dirname(mutFile)
    vcfFilePath = mutFile + '.vcf'
    ensemblFilePath = mutFile + '.ensembl'
    logging.basicConfig(filename='{}.log'.format(mutFile), filemode='a+', level=logging.DEBUG)
    logging.info(" Starting annotation pipleline ")

else:
    logging.info("Please pass the absolute path of the file to annotate")

def run():
    formatToVcfOrEnsemblAndSave()
    processFiles()
    annotateFile(vcfFilePath, "vcf")
    annotateFile(ensemblFilePath, "ensembl")
    mergeResultAnnos(vcfFilePath, ensemblFilePath)
    logging.info("Annotating is complete")


def formatToVcfOrEnsemblAndSave():
    with open(vcfFilePath, "w") as vcfFile, \
            open(ensemblFilePath, "w") as ensembl:
        reader = ps.read_csv(mutFile, sep='\t', dtype=str)
        reader['chromsome'] = reader.apply(lambda x: formatChromo(x.chromosome), axis=1)
        logging.info("Writing {0} to VCF".format(mutFile))
        vcfFile.write("#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfo\n")
        ensembl.write("#chrom\tpos\tend\tref/alt\tstrand\tid\n")
        for index, row in reader.iterrows():
            rowDict = row.to_dict()
            if isRowValidForProcessing(rowDict):
                if rowIsEnsembl(rowDict):
                    formatRowToEnsemblAndWrite(rowDict, ensembl)
                else:
                    formatRowToVCFAndWrite(rowDict, vcfFile)
        message = "The file {0} has {1} data points (including header)".format(mutFile, (index + 1))
        logging.info(message)


def rowIsEnsembl(row):
    refAllele = row["ref_allele"]
    altAllele = row["alt_allele"]
    return bool(re.search('-', refAllele)) or bool(re.search('-', altAllele))


def isRowValidForProcessing(row):
    isValid = True
    hg38RE = "(?i)(hg38|grch38|38)"
    if not bool(re.match(hg38RE, row["genome_assembly"])):
        logging.warning("Warning found legacy data : {0}".format(row.items()))
        isValid = False
    elif anyGenomicCoordinateAreMissing(row):
        logging.info(
            "Row has incomplete data : {0} in file {1} caused by missing chro,seq start, ref or alt allele data"
                .format(row.items(), vcfFilePath))
        isValid = False
    return isValid


def processFiles():
    logging.info("removing duplicates in VCF/Ensembl files")
    sortInPlace(vcfFilePath)
    sortInPlace(ensemblFilePath)
    dropDuplicates(vcfFilePath)
    dropDuplicates(ensemblFilePath)


def sortInPlace(aVCFfile):
    with open(aVCFfile, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        secondKey = sorted(reader, key=lambda x: sortLocation(x[1]))
        sortedFile = sorted(secondKey, key=lambda x: sortChromo(x[0]))
    with open(aVCFfile, 'w') as outFile:
        writer = csv.writer(outFile, delimiter='\t')
        for row in sortedFile:
            writer.writerow(row)

def sortChromo(x):
    decMatch = "^[0-9]{1,2}$"
    isDec = re.match(decMatch, x)

    if isDec:
        return int(x)
    elif len(x) == 1:
        return ord(x)

    return 0

def sortLocation(x):
    firstTenD = "^[0-9]{1,10}"
    loc = re.search(firstTenD, x)
    return int(loc.group(0)) if x != 'pos' and loc != None else 0

def dropDuplicates(vcfFilePath):
    vcfDf = ps.read_csv(vcfFilePath, sep='\t', keep_default_na=False, na_values=[''], dtype=str)
    vcfDf.drop_duplicates(inplace=True)
    vcfDf.to_csv(vcfFilePath, sep='\t', index=False, na_rep='')


def formatRowToVCFAndWrite(row, vcfFile):
    vcfRow = "{0}\t{1}\t{2}\t{3}\t{4}\t.\t.\t.\n".format(row['chromosome'], row["seq_start_position"], createPosId(row),
                                                         row["ref_allele"], row["alt_allele"])
    vcfFile.write(vcfRow)


def formatRowToEnsemblAndWrite(row, ensemblFile):
    startPos, endPos = resolveEnsemblEndPos(row)
    ensemblRow = "{0}\t{1}\t{2}\t{3}/{4}\t+\t{5}\n".format(row['chromosome'], startPos, endPos,
                                                           row["ref_allele"], row["alt_allele"], createPosId(row))
    ensemblFile.write(ensemblRow)


def resolveEnsemblEndPos(row):
    startPos = row['seq_start_position']
    endPos = startPos
    if bool(re.search("-", row["ref_allele"])):
        # insertion rule: Start - 1 = end coordinate
        endPos = int(startPos) - 1
    elif bool(re.search("-", row["alt_allele"])):
        # Deletion rule: Endpos = startPos + ( len(refAllele) - 1 )
        endPos = int(startPos) + len(row["ref_allele"]) - 1
    return startPos, endPos


def createPosId(row):
    return "{}_{}_{}_{}".format(formatChromo(row["chromosome"]),
                row["seq_start_position"],row["ref_allele"], row["alt_allele"])

def anyGenomicCoordinateAreMissing(row):
    return not row["chromosome"] or not row["seq_start_position"] or not row["ref_allele"] or not row["alt_allele"]


def allGenomicDataIsMissing(row):
    return not row["chromosome"] and not row["seq_start_position"] and not row["ref_allele"] and not row["alt_allele"]


def formatChromo(givenChromo):
    formattedChromo = givenChromo
    incorrectChrFormat = "(?i)^([0-9]{1,2}|[xym]{1}|mt|un)$"
    isMatch = re.match(incorrectChrFormat, givenChromo)
    if isMatch:
        formattedChromo = "chr" + isMatch.group(1)
    return formattedChromo


def annotateFile(vepIn, format):
    fastaDir, alleleDB, singularityVepImage, vepArgumentsList = getVepConfigurations()
    vepArguments = " ".join(vepArgumentsList)
    vepWarningFile = vepIn + ".vepWarnings"
    vepOut = vepIn + ".ANN"
    threads = cpu_count() * 2
    vepCMD = """vep {0} --format {1} --fork={2} --warning_file {3} -cache -dir_cache {4} -fasta {5} -i {6} -o {7} 2>> {8}.log""" \
        .format(vepArguments, format, threads, vepWarningFile, alleleDB, fastaDir, vepIn, vepOut, mutFile)
    logging.info("singularity exec {0} {1}".format(singularityVepImage, vepCMD))
    returnSignal = sp.call(
        "singularity exec {0} {1}".format(singularityVepImage, vepCMD), shell=True)
    if (returnSignal != 0):
        raise Exception("Vep returned a non-zero exit code {}".format(returnSignal))


def getVepConfigurations():
    with open('./config.yaml', 'r') as config:
        configDirs = yaml.safe_load(config)
        fastaDir = configDirs.get("fastaDir")
        alleleDB = configDirs.get("alleleDB")
        singularityVepImage = configDirs.get("vepSingularityImage")
        vepArguments = configDirs.get("vepArguments")
        if not os.path.exists(fastaDir):
            raise IOError("Fasta database does not exist at {}".format(fastaDir))
        if not os.path.exists(alleleDB):
            raise IOError("vep data base does not exist at {}".format(alleleDB))
        if not os.path.exists(singularityVepImage):
            raise IOError("singularity vep image does not exist at {}".format(singularityVepImage))
    return fastaDir, alleleDB, singularityVepImage, vepArguments

def mergeResultAnnos(vcfPath, ensemblPath):
       vcfAnnos = vcfPath + ".ANN"
       ensemblAnnos = ensemblPath + ".ANN"
       mergedAnnos = mutFile + ".ANN"
       with open(vcfAnnos, 'r') as vcfFile:
           vcfFile.readline()
           vcfFile.readline()
           matchGroup = re.search("Format:(.+$)", vcfFile.readline())
           infoColumnsHeaders = matchGroup.group(1).split("|")

           vcfDf = ps.read_csv(vcfAnnos, sep='\t', header=3)
           headers = vcfDf.columns
           mergedAnnosDf = (ps
                .read_csv(ensemblAnnos, sep='\t', header=3, names=headers)
                .append(vcfDf, ignore_index=True)
           )

           infoColumns = mergedAnnosDf['info'].str.split("|").tolist()
           infoColumnsDf = ps.DataFrame(infoColumns, columns=infoColumnsHeaders)
           mergedAnnosDf.join(infoColumnsDf).to_csv(mergedAnnos, sep='\t', index=False)
           os.remove(vcfAnnos)
           os.remove(ensemblAnnos)


if len(sys.argv) > 1:
    run()
