import time
import pandas as ps
import requests
import json
import re
import datetime
import logging
import glob
import random

logging.basicConfig(filename="../../log_results.tsv", level=logging.INFO)
logging.info("Starting data integration checks {}".format(datetime.datetime.now()))
logging.info("filename\trsid\tvariant_class\tchro\tpos\tapi_pos\tpos_results\tref_allele\tref_alt\tapi_ref_allele\tallele_result")

ncbiRootURL = 'http://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/'
ensemblApiRsRoot = 'http://rest.ensembl.org/variation/human/rs'
seqApi = 'http://rest.ensembl.org/sequence/region/human/'

updogFolders = glob.glob("/Users/tushar/pdx/pdxfinder-data/data/UPDOG/*")

def parseDBsnpJson(responseJson):
    dbsnpJson = json.load(responseJson)
    primaryShot = dbsnpJson.get("primary_snapshot_data")
    if primaryShot != None:
        alleles = primaryShot.get("placements_with_allele")[0].get("alleles")
        if alleles != None:
            allele = alleles[0].get("allele").get("spdi")
            allele_pos = allele.get("position")
            allele_ref = allele.get("deleted_sequence")
            allele_alt = allele.get("inserted_sequenc")
    return allele_pos, allele_pos, allele_ref, allele_alt


def parseEnsemblJson(responseJson):
    mappings = json.loads(responseJson).get("mappings")[0]
    if mappings != None:
        allele_start = mappings.get("start")
        allele_pos = mappings.get("location")
        allele_ref = mappings.get("ancestral_allele")
        allele_alt = mappings.get("allele_string")
    return allele_start, allele_pos, allele_ref, allele_alt


def requestJson(apiKey, database, backOff):
    response = ""
    if not backOff == 100:
        if database == 'dbsnp':
            rootApi = ncbiRootURL
            url = rootApi + apiKey
        elif database == 'ensembl':
            rootApi = ensemblApiRsRoot
            url = rootApi + apiKey + "?content-type=application/json"
        elif database == 'sequence':
            rootApi = seqApi
            url = rootApi + apiKey + "?content-type=application/json"
        print(url)
        response = requests.get(url)
        if (response.status_code != 200):
            time.sleep(backOff)
            incrBackOff = backOff + 10
            response = requestJson(apiKey, database, incrBackOff)
    return response

def parseSeq(seqJson):
    return "{} {}".format(seqJson.get("query"), seqJson.get("seq"))

def checkRsid(row, database):
    fileURI = row.loc['file']
    variationClass = row.loc['variant_class']
    variations = row.loc['variation_id']
    rsid = re.match("rs[0-9]+", variations).group(0)
    chro = row.loc['chromosome']
    pos = row.loc['seq_start_position']
    ref = row.loc['ref_allele']
    alt = row.loc['alt_allele']
    rsNum = re.sub('rs', '', rsid)
    if rsNum.isdigit():
        replyJson = requestJson(rsNum, database, 10).text
        time.sleep(.5)
        if replyJson is not None:
            if database == 'dbsnp':
                start_pos, allele_pos, allele_ref, allele_alt = parseDBsnpJson(replyJson)
            elif database == 'ensembl':
                start_pos, allele_pos, allele_ref, allele_alt = parseEnsemblJson(replyJson)


        region = "{}:{}..{}:1".format(chro, pos-3, pos+3)
        seqJson = requestJson(region, "sequence", 10).json()
        if requestJson is not None:
            surroundingSeq = parseSeq(seqJson)
            posStr = ""
            additionalStr = ""
            alleleStr = ""
            if pos == start_pos:
                posStr = 'alleles match'
            else:
                posStr = "positions do not match"
                if pos - start_pos == 1 or pos - start_pos == -1:
                    additionalStr = " off-by-one "
            logMessage = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}\t{}".format(fileURI, rsid, variationClass, chro, pos,
                                                                        allele_pos, posStr, additionalStr,
                                                                        ref,alt, allele_ref, allele_alt, surroundingSeq)
            print(logMessage)
            logging.info(logMessage)

def analyzeMutFileVariants(mutFile, sampleSize):
    if mutFile:
        df = ps.read_csv(mutFile, sep='\t', ) \
            .fillna(" ")
        rsDf = df[df.variation_id.str.contains("rs[0-9]{0,10}")]
        database = "ensembl"
        snvDf = rsDf[rsDf.variant_class == "SNV"]
        if len(snvDf) >= sampleSize:
            (snvDf
             .sample(sampleSize)
             .assign(file=mutFile)
             .apply(lambda x: checkRsid(x, database), axis=1)
             )
        insertionDf = rsDf[rsDf.variant_class == "insertion"]
        if len(insertionDf) >= sampleSize:
            (insertionDf
             .sample(sampleSize)
             .assign(file=mutFile)
             .apply(lambda x: checkRsid(x, database), axis=1)
             )
        deletionDf = rsDf[rsDf.variant_class == "deletion"]
        if len(deletionDf) >= sampleSize:
            (deletionDf
             .sample(sampleSize)
             .assign(file=mutFile)
             .apply(lambda x: checkRsid(x, database), axis=1)
             )
    else:
        print("No mut files found for {}".format(i))

def parseAllSamples():
    updogFolders = glob.glob("/home/afollette/Finder_Data_Repositories/pdxfinder-data/data/UPDOG/*")
    for i in updogFolders:
        mutFilesGlob = "{}/mut/*_mut*tsv".format(i)
        mutFiles = glob.glob(mutFilesGlob)
        randomMutFile = random.sample(mutFiles, 1)
        analyzeMutFileVariants(randomMutFile,2)


def parseParticularSamples():
    curieLc = glob.glob("/home/afollette/Finder_Data_Repositories/pdxfinder-data/data/UPDOG/Curie-LC/mut/Curie-LC_mut.tsv")[0]
    lih = glob.glob("/home/afollette/Finder_Data_Repositories/pdxfinder-data/data/UPDOG/LIH/mut/LIH_mut.tsv")[0]
    crl_rnaseq = glob.glob("/home/afollette/Finder_Data_Repositories/pdxfinder-data/data/UPDOG/CRL/mut/rnaseq/*_mut*.tsv")
    crl_wes = glob.glob("/home/afollette/Finder_Data_Repositories/pdxfinder-data/data/UPDOG/CRL/mut/wes/*_mut*tsv")
    randomCRL_rnaseq = random.sample(crl_rnaseq, 1)[0]
    randomCRL_Wes = random.sample(crl_wes, 1)[0]

    analyzeMutFileVariants(curieLc, 1)
    analyzeMutFileVariants(lih, 1)
    #analyzeMutFileVariants(randomCRL_rnaseq, 5)
    #analyzeMutFileVariants(randomCRL_Wes, 5)


parseParticularSamples()