#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

#BSUB -j $1_annotater_$(date)
#BSUB --mail-user=afollette@ebi.ac.uk
#BSUB  -B -N
#BSUB -e /homes/afollette/$i.err.%j
#BSUB -o /homes/afollette/$i.out.%j
#BSUB -M 10000
#BSUB -n 4



import IOutilities
import vcfSorter

import os
import csv
import subprocess as sp
import sys

file = sys.argv[1]

def run():

    global parentDirectoryPath
    parentDirectoryPath = os.path.dirname(file)
    global vcfFilePath
    vcfFilePath = file + '.vcf'

    formatToVCFAndSave(file)
    annotateVCF(vcfFilePath,file)
    sp.call(("bash mergeWrapper.sh {0}".format(file)), shell=True)



def formatToVCFAndSave(filePath):

    tsvOrCsvFile = open(filePath)
    vcfFile = open(vcfFilePath, "w+")
    IOutilities.chmodFile(vcfFilePath)

    if filePath.endswith(".tsv"):
        reader = csv.DictReader(tsvOrCsvFile, delimiter="\t")
    elif filePath.endswith(".csv"):
        reader = csv.DictReader(tsvOrCsvFile, delimiter=",")

    assert(reader != None)

    vcfFile.write("#chrom\tpos\tid\tref\talt\tqual\tfilter\tinfo\n")

    rowCount = 0

    for row in reader:
        rowCount += 1
        if attemptToWriteRowToVCFisNotSuccessful(row,vcfFile) :
            break

    IOutilities.flushCloseFile(tsvOrCsvFile)
    IOutilities.flushCloseFile(vcfFile)

    vcfSorter.sort(vcfFilePath, vcfFilePath)
    #uniqVCF(vcfFilePath,vcfFilePath)
    os.path.dirname(os.path.dirname(filePath))

    masterLog = "{0}/masterLog".format(os.path.dirname(os.path.dirname(os.path.dirname(filePath))))
    message = "The file {0} has {1} data points (including header)".format(filePath, rowCount)
    IOutilities.masterlogMessage(masterLog,message)

def uniqVCF(vcfFilePath,vcfOutFile):

    cmds = ['uniq ',vcfFilePath,' > ', vcfFilePath]
    print(cmds)
    sp.calls(cmds,shell=False)

def attemptToWriteRowToVCFisNotSuccessful (row, vcfFile) :

    isEOF_orError = False

    if (row["genome_assembly"] == "Hg38") or (row["genome_assembly"] == "GRCh38"):
        if genomeDataIsMissing(row):
            IOutilities.logMessage(parentDirectoryPath, "Row has incomplete data : {0} in file {1} caused by missing chro,seq start, ref or alt allele data".format(row.items(), vcfFilePath))
        elif allGenomicDataIsMissing(row):
            isEOF_orError = True
        else:
            formatRowToVCFAndWrite(row, vcfFile)
    else:
        IOutilities.logMessage(parentDirectoryPath, "Warning found legacy data : {0}".format(row.items()))
    return isEOF_orError

def formatRowToVCFAndWrite (row, vcfFile) :

    chromo = IOutilities.formatChromo(row["chromosome"])

    vcfRow = "{0}\t{1}\t.\t{2}\t{3}\t.\t.\t.\t\n".format(chromo, row["seq_start_position"],
                                                         row["ref_allele"], row["alt_allele"])
    vcfFile.write(vcfRow)

def genomeDataIsMissing (row) :
    return not row["chromosome"] or not row["seq_start_position"] or not row["ref_allele"] or not row["alt_allele"]

def allGenomicDataIsMissing (row) :
    return not row["chromosome"] and not row["seq_start_position"] and not row["ref_allele"] and not row["alt_allele"]

def annotateVCF(vcfFile, file):

    vcfFileName = os.path.basename(vcfFilePath)
    provider = os.path.dirname(parentDirectoryPath)
    providerName = os.path.basename(provider)

    fastaDir = "/nfs/nobackup/spot/mouseinformatics/pdx/vepDBs/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    alleleDB = "/nfs/nobackup/spot/mouseinformatics/pdx/vepDBs/homo_sapiens_vep_98_GRCh38"

    vepIn = vcfFile
    vepOut = file + ".ANN"
    annoFilename = file + ".ANN"

    open(annoFilename, "w+")
    IOutilities.chmodFile(annoFilename)

    vepCMD = """vep -e -q -check_existing  -symbol -polyphen -sift -merged --use_transcript_ref —hgvs —hgvsg —variant_class -canonical -fork 4 -format vcf -force -offline -no_stats -cache -dir_cache {0} -fasta {1} -i {2} -o {3}""".format(alleleDB,fastaDir,vepIn, vepOut)

    #print("singularity exec ensembl-vep.simg {0}".format(vepCMD))

    sp.call(
        "singularity exec ensembl-vep.simg {0}".format(vepCMD), shell=True)


run()
