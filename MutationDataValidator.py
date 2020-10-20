import time
import pandas as ps
import requests
import json
import re
import datetime
import logging
import glob
import random

logging.basicConfig(filename="log_results.tsv",level=logging.INFO)
logging.info("Starting data integration checks {}".format(datetime.datetime.now()))
logging.info("filename\trsid\tvariant_class\tchro\tpos\tapi_pos\tpos_results\tref_allele\tref_alt\tapi_ref_allele\tallele_result")

api_rootURL = 'https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/'
updogFolders = glob.glob("/home/afollette/Finder_Data_Repositories/pdxfinder-data/data/UPDOG/*")
variationFilter = 'deletion'

def checkRsid(row):
    backOff = 10
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
        url = api_rootURL + rsNum
        print(url)
        response = requests.get(url)
        if (response.status_code != 200):
            time.sleep(backOff)
            backOff += 10
        else:
            replyJson = json.loads(response.text)
            if replyJson != None:
                primaryShot = replyJson.get("primary_snapshot_data")
                if primaryShot != None:
                    alleles = primaryShot.get("placements_with_allele")[0].get("alleles")
                    if alleles != None:
                        allele = alleles[0].get("allele").get("spdi")
                        allele_pos = allele.get("position")
                        allele_ref_seq = allele.get("deleted_sequence")
                        posStr = ""
                        additionalStr = ""
                        alleleStr = ""
                        if pos == allele_pos:
                            posStr = 'alleles match'
                        else:
                            posStr = "positions do not match"
                            if pos - allele_pos == 1 or pos - allele_pos == -1:
                                additionalStr = " off-by-one "
                        if allele_ref_seq == ref:
                            alleleStr = "Ref sequences match"
                        else:
                            alleleStr = "Ref's do not match on first lookup allele"

                        logMessage = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}".format(fileURI, rsid, variationClass, chro, pos,
                                                                                    allele_pos, posStr, additionalStr,
                                                                                    ref,alt, allele_ref_seq,alleleStr)
                        print(logMessage)
                        logging.info(logMessage)

sampleDf = ps.DataFrame()

for i in updogFolders:
    mutFilesGlob = "{}/mut/*_mut*tsv".format(i)
    mutFiles = glob.glob(mutFilesGlob)
    if mutFiles:
        dataFile = random.sample(mutFiles,1)
        df = ps.read_csv(dataFile[0], sep='\t')
        rsDf = df[df.variation_id.str.contains("rs[0-9]{0,10}")]
        sampleSize = 10
        snvDf = rsDf[rsDf.variant_class == "SNV"]
        if len(rsDf) >= snvDf:
            (snvDf
             .sample(sampleSize)
             .assign(file=i)
             .apply(lambda x: checkRsid(x), axis=1)
             )
        insertionDf = rsDf[rsDf.variant_class == "insertion"]
        if len(insertionDf) >= sampleSize:
            (insertionDf
             .sample(sampleSize)
             .assign(file=i)
             .apply(lambda x: checkRsid(x), axis=1)
             )
        deletionDf = rsDf[rsDf.variant_class == "deletion"]
        if len(deletionDf) >= sampleSize:
            (deletionDf
             .sample(sampleSize)
             .assign(file=i)
             .apply(lambda x: checkRsid(x), axis=1)
             )
    else:
        print("No mut files found for {}".format(i))