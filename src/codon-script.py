#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import os
import sys
import logging
import time
from os.path import join, exists, isdir, isfile
from Annotater import Annotater
from AnnotationMerger import AnnotationMerger
from shutil import rmtree

def get_dirs(path):
    return [f for f in os.listdir(path) if isdir(join(path, f))]

def get_files(path):
    return [join(path, f) for f in os.listdir(path) if isfile(join(path, f)) and f.endswith('.tsv')]

def generate_mutTarget(path):
    file_list = get_files(path)
    if len(file_list) == 0:
        dirs = get_dirs(path)
        for dir in dirs:
            temp = get_files(join(path, dir))
            file_list = file_list + temp
    return file_list

def remove_files(path):
    if exists(path):
        os.remove(path)
        logging.info("File removed: " + path)
    else:
        logging.info("File does not exist: " + path)

if len(sys.argv) > 1:
    target = sys.argv[1]
    run_type = sys.argv[2]
    local = sys.argv[3]
    log_location = sys.argv[4]
    local = local == "local"
    logging.basicConfig(filename='{}.log'.format(log_location), filemode='a+', level=logging.DEBUG)
    logging.info("{0}: Starting annotations".format(time.ctime()))
    #skip_provider = ["BROD", "CCIA", "CHOP", "CMP", "CRL", "CSHL", "CUIMC", "Curie-BC", "Curie-LC", "GCCRI"]
    skip_provider = []
    if exists(target):
        for provider in sorted(get_dirs(target)):
            start = time.time()
            provider_path = join(target, provider)
            mut_path = join(provider_path, 'mut')
            if exists(mut_path) and provider not in skip_provider:
                files = generate_mutTarget(mut_path)
                annotater = Annotater(mut_path, run_type, local)
                for mutTarget in files:
                    if os.path.isfile(mutTarget):
                        #logging.info("Annotating file: " + mutTarget)
                        annotater.run(mutTarget)
                annotater.processFiles()
                mid = time.time() - start
                logging.info("{0}: Create VEP input for {1} took {2} mins!".format(time.ctime(), provider, round(mid/60)))
                annotater.annotate()
                mid = time.time() - mid
                logging.info("{0}: Annotations for {1} took {2} mins!".format(time.ctime(), provider, round(mid/60)))
                annotater = AnnotationMerger(mut_path, run_type, local)
                for mutTarget in files:
                    if os.path.isfile(mutTarget):
                        annotater.run(mutTarget)
                        os.remove(mutTarget)
                        os.rename(mutTarget + '.hmz', mutTarget)
                    else:
                        logging.info("Not a file: " + mutTarget)
                mid = time.time() - mid
                logging.info("{0}: Merge of annotations for {1} took {2} mins!".format(time.ctime(), provider, round(mid/60)))
                rmfs = ['.vcf', '.vcf.ANN', '.ensembl', '.ensembl.ANN', '.ensembl.vepWarnings', '.ANN']
                rmtree(join(mut_path, 'annotations'))
                end = round((time.time() - start)/60)
                logging.info("{0}: Annotations for {1} took {2} mins!".format(time.ctime(), provider, end))
            else:
                logging.info("Please pass the absolute path of the file to annotate")