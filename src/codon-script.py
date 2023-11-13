#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import os
import subprocess as sp
import sys
import logging
from os.path import join, exists
from Annotater import Annotater
from AnnotationMerger import AnnotationMerger

def get_dirs(path):
    return [f for f in os.listdir(path) if os.isdir(join(path, f))]

def get_files(path):
    return [join(path, f) for f in os.listdir(path) if os.isfile(join(path, f)) and f.endswith('.tsv')]

def generate_mutTarget(path):
    file_list = get_files(path)
    if len(file_list) == 0:
        dirs = get_dirs(path)
        for dir in dirs:
            temp = get_files(join(path, dir))
            file_list = file_list + temp
    return file_list


if len(sys.argv) > 1:
    logging.basicConfig(filename='{}.log'.format(), filemode='a+', level=logging.DEBUG)
    logging.info("Starting annotations")
    target = sys.argv[1]
    run_type = sys.argv[2]
    local = sys.argv[3]
    local = local == "local"

    if exists(target):
        for provider in get_dirs(target):
            provider_path = join(target, provider)
            mut_path = join(provider_path, 'mut')
            if exists(mut_path):
                files = generate_mutTarget(mut_path)
                for mutTarget in files:
                    if os.path.isfile(mutTarget):
                        logging.info("Annotating file: " + mutTarget)
                        Annotater(mutTarget, run_type, local).run()
                        logging.info("Starting merge of annotations")
                        AnnotationMerger(mutTarget, run_type, local).run()
                        logging.info("Annotations complete")
                        logging.info(sp.call("tail -n 2 "+mutTarget+".log"))
                    else:
                        logging.info("Not a file: "+ mutTarget)
            else:
                logging.info("Please pass the absolute path of the file to annotate")