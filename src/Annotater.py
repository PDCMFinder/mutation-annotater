#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import glob
import os
import csv
import subprocess as sp
from multiprocessing import cpu_count
import sys
import re
import logging
import pandas as pd
import yaml
from os.path import join


class Annotater:
    def __init__(self, mutTarget, run_type, local, configDir='./config.yaml'):
        self.mutTarget = mutTarget
        self.configDir = configDir
        self.run_type = run_type
        self.local = local
        vcf_cols = ["#chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"]
        ensembl_cols = ["#chrom", "pos", "end", "ref/alt", "strand", "id"]
        if os.path.isfile(self.mutTarget):
            target = os.path.dirname(self.mutTarget)
        else:
            target = mutTarget
        self.hgvsFilePath = os.path.join(target, 'merged.hgvs')
        self.vcfFilePath = os.path.join(target, 'merged.vcf')
        self.ensemblFilePath = os.path.join(target, 'merged.ensembl')
        self.annFilePath = os.path.join(target, 'merged.ANN')
        self.vcfDf = pd.DataFrame(columns=vcf_cols)
        self.ensemblDf = pd.DataFrame(columns=ensembl_cols)

    def run(self, mutFile):
        if mutFile != '':
            self.mutTarget = join(self.mutTarget, mutFile)
        self.fileName = os.path.basename(self.mutTarget)
        self.parentDirectoryPath = os.path.dirname(self.mutTarget)
        logging.basicConfig(filename='{}.log'.format(join(self.parentDirectoryPath, 'annotater')), filemode='a+', level=logging.DEBUG)
        logging.info("Starting merge of annotations using " + str(cpu_count())+" CPUs")
        if self.run_type == 'vcf':
            self.process_VCForEnsembl()
        elif self.run_type == 'hgvs':
            self.process_hgvs()

    def annotate(self):
        if self.run_type == 'hgvs':
            self.annotateFile(self.hgvsFilePath, 'hgvs')
            self.process_hgvs_annotations()
            os.remove(self.hgvsFilePath)
        elif self.run_type == 'vcf':
            self.annotateFile(self.vcfFilePath, "vcf")
            self.annotateFile(self.ensemblFilePath, "ensembl")
            self.mergeResultAnnos(self.vcfFilePath, self.ensemblFilePath)
        logging.info("Annotating is complete")


    def process_hgvs(self):
        self.formatHGVSFiles()
        self.dropDuplicates(self.hgvsFilePath)
        #'singularity exec -B /Users/tushar/PycharmProjects/PDCMDataAggregator/output,/Users/tushar/pdx/update-data,/Users/tushar/pdx/mutation-annotater:/Users/tushar/pdx/mutation-annotater:rw pdx-liftover_vep_release98.3.sif vep --verbose --vcf --force --check_existing --check_ref --refseq --use_given_ref --symbol --polyphen s --sift s --hgvs --variant_class --no_stats --pick --pick_order biotype,canonical,appris,tsl,ccds,rank,length,mane --fork=16 --warning_file  /Users/tushar/PycharmProjects/PDCMDataAggregator/output/CMP/pdcm_format/mut/CMP_mut.tsv.ensembl.vepWarnings --cache --dir_cache /Users/tushar/pdx/mutation-annotater/vepDB/homo_sapiens_vep_98_GRCh38 --fasta /Users/tushar/pdx/mutation-annotater/vepDB/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa --format hgvs -i /Users/tushar/PycharmProjects/PDCMDataAggregator/output/CMP/pdcm_format/mut/WGS/CMP_mut_SIDS01895.tsv.refseq -o /Users/tushar/PycharmProjects/PDCMDataAggregator/output/CMP/pdcm_format/mut/WGS/CMP_mut_SIDS01895.tsv.ANN 2>> /Users/tushar/PycharmProjects/PDCMDataAggregator/output/CMP/pdcm_format/mut/WGS/CMP_mut_SIDS01895.tsv.log'
        #self.annotateFile(self.hgvsFilePath, 'hgvs')
        #self.process_hgvs_annotations()
        #os.remove(self.hgvsFilePath)

    def process_hgvs_annotations(self):
        mergedAnnos = self.mutTarget + ".ANN"
        vcfAnnos = self.hgvsFilePath+".ANN"
        with open(vcfAnnos, 'r') as vcfFile:
            vcfFile.readline()
            vcfFile.readline()
            matchGroup = re.search("Format:(.+$)", vcfFile.readline())
            infoColumnsHeaders = matchGroup.group(1).split("|")

            vcfDf = pd.read_csv(vcfAnnos, sep='\t', header=3)
            headers = vcfDf.columns
            infoColumns = vcfDf['INFO'].str.split("|").tolist()
            infoColumnsDf = pd.DataFrame(infoColumns, columns=infoColumnsHeaders)

            vcfDf.join(infoColumnsDf).rename(columns=str.lower).to_csv(mergedAnnos, sep='\t', index=False)
            os.remove(vcfAnnos)

    def formatHGVSFiles(self):
        with open(self.hgvsFilePath, 'w') as hgvsFile:
            reader = pd.read_csv(self.mutTarget, sep='\t', dtype=str, engine='python').fillna('')
            logging.info('Writing {0} to HGVS'. format(self.mutTarget))
            hgvsFile.write('#HGVSIdentifier\n')
            for index, row in reader.iterrows():
                if row['ncbi_transcript_id'] != '' and row['coding_sequence_change'] != '':
                    out_row = row['ncbi_transcript_id']+':c.'+row['coding_sequence_change'] +'\n'
                    hgvsFile.write(out_row)
            message = "The file {0} has {1} data points (including header)".format(self.mutTarget, (index + 1))
            logging.info(message)

    def process_VCForEnsembl(self):
        self.formatToVcfOrEnsemblAndSave()
        #self.processFiles()
        #self.annotateFile(self.vcfFilePath, "vcf")
        #self.annotateFile(self.ensemblFilePath, "ensembl")
        #self.mergeResultAnnos(self.vcfFilePath, self.ensemblFilePath)

    def formatToVcfOrEnsemblAndSave(self):
        reader = pd.read_csv(self.mutTarget, sep='\t', dtype=str, engine='python')
        reader = self.isRowValidForProcessing(reader)
        reader['chromosome'] = reader.apply(lambda x: self.formatChromo(x.chromosome), axis=1)
        reader['id'] = reader.apply(
            lambda x: "{}_{}_{}_{}".format(self.formatChromo(x["chromosome"]), x["seq_start_position"],
                                           x["ref_allele"], x["alt_allele"]), axis=1)

        logging.info("Writing {0} to VCF".format(self.mutTarget))
        self.generate_ensembl_file(reader)
        self.generate_vcf_file(reader)

        #for index, row in reader.iterrows():
        #    rowDict = row.to_dict()
        #    if self.isRowValidForProcessing(rowDict):
        #        if self.rowIsEnsembl(rowDict):
        #            self.formatRowToEnsemblAndWrite(rowDict, ensembl)
        #        else:
        #            self.formatRowToVCFAndWrite(rowDict, vcfFile)

        message = "The file {0} has {1} data points (including header)".format(self.mutTarget, (reader.shape[0] + 1))
        logging.info(message)

    def isRowValidForProcessing(self, df):
        indices_before = df.index
        df = df.dropna(subset=['chromosome', 'seq_start_position', 'ref_allele', 'alt_allele'])
        indices_after = df.index
        dropped_indices = indices_before.difference(indices_after)
        for index in dropped_indices:
            logging.info(
                "Row has incomplete data : Row {0} in file {1} has missing chr, seq start, ref or alt allele data"
                .format(index, self.mutTarget))
        return df

    def processFiles(self):
        logging.info("removing duplicates in VCF/Ensembl files")
        self.vcfDf = sort_vcf_ensembl_df(self.vcfDf)
        self.ensemblDf = sort_vcf_ensembl_df(self.ensemblDf)
        self.vcfDf.to_csv(self.vcfFilePath, sep='\t', index=False)
        self.ensemblDf.to_csv(self.ensemblFilePath, sep='\t', index=False)
        #self.sortInPlace(self.vcfFilePath)
        #self.sortInPlace(self.ensemblFilePath)

    def sortInPlace(self, aVCFfile):
        with open(aVCFfile, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            secondKey = sorted(reader, key=lambda x: self.sortLocation(x[1]))
            sortedFile = sorted(secondKey, key=lambda x: self.sortChromo(x[0]))
        with open(aVCFfile, 'w') as outFile:
            writer = csv.writer(outFile, delimiter='\t')
            for row in sortedFile:
                writer.writerow(row)

    def sortChromo(self, x):
        decMatch = "^[0-9]{1,2}$"
        isDec = re.match(decMatch, x)
        if isDec:
            return int(x)
        elif len(x) == 1:
            return ord(x)
        return 0

    def sortLocation(self, x):
        firstTenD = "^[0-9]{1,10}"
        loc = re.search(firstTenD, x)
        return int(loc.group(0)) if x != 'pos' and loc != None else 0

    def dropDuplicates(self, vcfFilePath):
        vcfDf = pd.read_csv(vcfFilePath, sep='\t', keep_default_na=False, na_values=[''], dtype=str)
        vcfDf.drop_duplicates(inplace=True)
        vcfDf.to_csv(vcfFilePath, sep='\t', index=False, na_rep='')

    def formatRowToVCFAndWrite(self, row, vcfFile):
        vcfRow = "{0}\t{1}\t{2}\t{3}\t{4}\t.\t.\t.\n".format(row['chromosome'], row["seq_start_position"],
                                                             self.createPosId(row),
                                                             row["ref_allele"], row["alt_allele"])
        vcfFile.write(vcfRow)

    def generate_vcf_file(self, df):
        temp = df[~df['ref_allele'].str.contains('-') & ~df['alt_allele'].str.contains('-')]
        if len(temp) > 0:
            temp = temp[['chromosome', 'strand', 'seq_start_position', 'ref_allele', 'alt_allele', 'id']]
            mapper = {'chromosome': '#chrom', 'seq_start_position': 'pos',
                      'ref_allele': 'ref', 'alt_allele': 'alt'}
            temp.rename(columns=mapper, inplace=True)
            temp['qual'] = ''
            temp['filter'] = ''
            temp['info'] = ''
            #temp['strand'] = ['1' if r == "" else str(int(r)) for r in temp['strand']]
            cols = ["#chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"]
            self.vcfDf = pd.concat([self.vcfDf, temp[cols]])
            self.vcfDf = self.vcfDf.drop_duplicates().reset_index(drop=True)

    def generate_ensembl_file(self, df):
        temp = df[df['ref_allele'].str.contains('-') | df['alt_allele'].str.contains('-')]
        if len(temp)>0:
            temp = temp[['chromosome', 'strand', 'seq_start_position', 'ref_allele', 'alt_allele', 'id']]
            temp['strand'] = ['1' if r == "" else str(int(r)) for r in temp['strand'].fillna('')]
            temp['pos'] = temp['seq_start_position']
            temp['end'] = temp.apply(lambda x: int(x['pos']) - 1 if bool(re.search("-", x["ref_allele"])) else int(x['pos']) + len(x["ref_allele"]) - 1 if bool(re.search("-", x["alt_allele"])) else "" , axis=1)
            temp['ref/alt'] = temp['ref_allele'].astype(str) + '/' + temp['alt_allele'].astype(str)
            temp['#chrom'] = temp['chromosome']
            cols = ['#chrom', 'pos', 'end', 'ref/alt', 'strand', 'id']
            self.ensemblDf = pd.concat([self.ensemblDf, temp[cols]])
            self.ensemblDf = self.ensemblDf.drop_duplicates().reset_index(drop=True)

    def formatRowToEnsemblAndWrite(self, row, ensemblFile):
        startPos, endPos = self.resolveEnsemblEndPos(row)
        strand = '+' if row['strand'] == "" else row["strand"]
        ensemblRow = "{0}\t{1}\t{2}\t{3}/{4}\t{5}\t{6}\n".format(row['chromosome'], startPos, endPos,
                                                                 row["ref_allele"], row["alt_allele"], strand,
                                                                 self.createPosId(row))
        ensemblFile.write(ensemblRow)

    def resolveEnsemblEndPos(self, row):
        startPos = row['seq_start_position']
        endPos = startPos
        if bool(re.search("-", row["ref_allele"])):
            # insertion rule: Start - 1 = end coordinate
            endPos = int(startPos) - 1
        elif bool(re.search("-", row["alt_allele"])):
            # Deletion rule: Endpos = startPos + ( len(refAllele) - 1 )
            endPos = int(startPos) + len(row["ref_allele"]) - 1
        return startPos, endPos

    def createPosId(self, row):
        return "{}_{}_{}_{}".format(self.formatChromo(row["chromosome"]),
                                    row["seq_start_position"], row["ref_allele"], row["alt_allele"])

    def anyGenomicCoordinateAreMissing(self, row):
        return not row["chromosome"] or not row["seq_start_position"] or not row["ref_allele"] or not row["alt_allele"]

    def allGenomicDataIsMissing(self, row):
        return not row["chromosome"] and not row["seq_start_position"] and not row["ref_allele"] and not row[
            "alt_allele"]

    def formatChromo(self, givenChromo):
        formattedChromo = givenChromo
        incorrectChrFormat = "(?i)^([0-9]{1,2}|[xym]{1}|mt|un)$"
        isMatch = re.match(incorrectChrFormat, givenChromo)
        if isMatch:
            formattedChromo = "chr" + isMatch.group(1)
        return formattedChromo

    def annotateFile(self, vepIn, format):
        fastaDir, alleleDB, singularityVepImage, vepArgumentsList, mutationAnnotator, dataPath = self.getVepConfigurations()
        vepArguments = " ".join(vepArgumentsList)
        mutationAnnotator = dataPath +":"+ dataPath +","+ mutationAnnotator + ":" + mutationAnnotator + ":rw"
        vepWarningFile = vepIn + ".vepWarnings"
        vepOut = vepIn + ".ANN"
        threads = cpu_count() * 2
        vepCMD = """vep {0} --format {1} --fork={2} --warning_file {3} --cache --dir_cache {4} --fasta {5} -i {6} -o {7} 2>> {8}.log""" \
            .format(vepArguments, format, threads, vepWarningFile, alleleDB, fastaDir, vepIn, vepOut, join(self.parentDirectoryPath, 'annotater'))
        if format != 'hgvs':
            vepCMD = vepCMD + " --offline --merged " #Get Refseq and Ensembl lookup databases
        elif format == 'hgvs':
            vepCMD = vepCMD + " --refseq "
        print(vepCMD)
        if not self.local:
            logging.info("singularity exec -B {0} {1} {2}".format(mutationAnnotator, singularityVepImage, vepCMD))
            returnSignal = sp.call(
                "singularity exec -B {0} {1} {2}".format(mutationAnnotator, singularityVepImage, vepCMD), shell=True)
        else:
            logging.info("vagrant ssh -c 'singularity exec -B {0} {1} {2}'".format(mutationAnnotator, singularityVepImage, vepCMD))
            returnSignal = sp.call("cd ../vm-singularity && vagrant ssh -c 'singularity exec -B {0} {1} {2}'".format(mutationAnnotator, singularityVepImage, vepCMD), shell=True)
        if (returnSignal != 0):
            raise Exception("Vep returned a non-zero exit code {}".format(returnSignal))

    def getVepConfigurations(self):
        with open(self.configDir, 'r') as config:
            configDirs = yaml.safe_load(config)
            if not self.local:
                mutationAnnotator = configDirs.get("mutationAnnotator_codon")
                dataPath = configDirs.get("dataPath_codon")
                fastaDir = join(mutationAnnotator, configDirs.get("fastaDir_codon"))
                alleleDB = join(mutationAnnotator,configDirs.get("alleleDB_codon"))
                singularityVepImage = join(mutationAnnotator, configDirs.get("vepSingularityImage_codon"))
            else:
                mutationAnnotator = configDirs.get("mutationAnnotator")
                dataPath = configDirs.get("dataPath")
                fastaDir = join(mutationAnnotator, configDirs.get("fastaDir"))
                alleleDB = join(mutationAnnotator, configDirs.get("alleleDB"))
                singularityVepImage = join(mutationAnnotator, configDirs.get("vepSingularityImage"))

            vepArguments = configDirs.get("vepArguments")

            if not os.path.exists(dataPath):
                raise IOError("Please fix the path to the data at {}".format(dataPath))
            if not os.path.exists(mutationAnnotator):
                raise IOError("Please fix the path to the mutation-annotator at {}".format(mutationAnnotator))
            if not os.path.exists(fastaDir):
                raise IOError("Fasta database does not exist at {}".format(fastaDir))
            if not os.path.exists(alleleDB):
                raise IOError("vep data base does not exist at {}".format(alleleDB))
            if not os.path.exists(singularityVepImage):
                raise IOError("singularity vep image does not exist at {}".format(singularityVepImage))

        return fastaDir, alleleDB, singularityVepImage, vepArguments, mutationAnnotator, dataPath

    def mergeResultAnnos(self, vcfPath, ensemblPath):
        vcfAnnos = vcfPath + ".ANN"
        ensemblAnnos = ensemblPath + ".ANN"
        mergedAnnos = self.annFilePath
        with open(vcfAnnos, 'r') as vcfFile:
            vcfFile.readline()
            vcfFile.readline()
            matchGroup = re.search("Format:(.+$)", vcfFile.readline())
            infoColumnsHeaders = matchGroup.group(1).split("|")

            vcfDf = pd.read_csv(vcfAnnos, sep='\t', header=4)
            headers = vcfDf.columns
            mergedAnnosDf = pd.concat([pd.read_csv(ensemblAnnos, sep='\t', header=4, names=headers), vcfDf], ignore_index=True)
            infoColumns = mergedAnnosDf['info'].str.split("|").tolist()
            infoColumnsDf = pd.DataFrame(infoColumns, columns=infoColumnsHeaders)
            mergedAnnosDf.join(infoColumnsDf).to_csv(mergedAnnos, sep='\t', index=False)


def cmdline_runner():
    if len(sys.argv) > 1:
        mutTarget = sys.argv[1]
        run_type = sys.argv[2]
        local = sys.argv[3]
        local = local == "local"
        annotate = Annotater(mutTarget, run_type, local)
        if os.path.isfile(mutTarget):
            logging.basicConfig(filename='{}.log'.format(mutTarget), filemode='a+', level=logging.DEBUG)
            logging.info(" Starting annotation pipleline ")
            annotate.run('')
            annotate.processFiles()
            annotate.annotate()
        elif os.path.isdir(mutTarget):
            globForTsv = os.path.join(mutTarget, "*tsv")
            for mutFile in glob.iglob(globForTsv):
                #mutfile_path = os.path.join(mutTarget, mutFile)
                #print(mutfile_path)
                annotate.run(mutFile)
            annotate.processFiles()
            annotate.annotate()

    else:
        logging.info("Please pass the absolute path of the file to annotate")

def sort_vcf_ensembl_df(df):
    cols = df.columns
    df['chromosome'] = pd.Categorical(df['#chrom'], ordered=True,
                                              categories=['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY'])
    df = df.sort_values(by=['chromosome', 'pos'])
    return df[cols]

cmdline_runner()
