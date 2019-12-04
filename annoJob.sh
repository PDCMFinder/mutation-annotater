#!/bin/bash
#BSUB -J $1_annotater_$(date)
#BSUB --mail-user=afollette@ebi.ac.uk
#BSUB  -B -N
#BSUB -e /homes/afollette/$i.err.%j
#BSUB -o /homes/afollette/$i.out.%j
#BSUB -M 10000
#BSUB -n 4

jobName=$1

bsub -J $jobName python2.7 Annotater.py $2 



