#BSUB -j $1_annotater_$(date)
#BSUB --mail-user=afollette@ebi.ac.uk
#BSUB  -B -N
#BSUB -e /homes/afollette/$i.err.%j
#BSUB -o /homes/afollette/$i.out.%j
#BSUB -M 10000
#BSUB -n 4

echo $(id)
echo $(groups) 

echo $1
echo $2

source /nfs/nobackup/spot/mouseinformatics/pdx/omicAnno/venv/bin/activate
python3 Annotater.py $2 



