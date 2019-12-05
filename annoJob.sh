#!/bin/bash

#BSUB --mail-user=afollette@ebi.ac.uk
#BSUB  -B -N
#BSUB -n 4

jobName=$1

bsub -J $jobName -n 4 -M 10000 -R "rusage[mem=10000]" -e /homes/afollette/$i.err.%j -o /homes/afolllette/$i.out.%j python2.7 Annotater.py $2 



