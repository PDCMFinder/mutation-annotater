# mutation-annotater
Workflow:

![Annotation pipeline workflow](annotation%20pipeline.png "Annotation pipeline workflow")

annotation pipeline for PDXFinder mut templates using VEP.

The pipeline was written in python 2.7. It has not been tested in higher versions.

The pipeline requires Singularity to run the containerized VEP. After installing 
Singularity to the system run the install.sh script. The install.sh script will:
    
* pdxfinder/pdx-liftover:vep_release98.3, this is the containerized instance of VEP  
* Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz, a VEP running dependency
* homo_sapiens_merged_vep_98_GRCh38.gz, another vep running dependency
* homo_sapiens_refseq_vep_98_GRCh38.gz, another vep running dependency


Once these are installed the variables in the `config.yaml` must be set.

These are:
- fastaDir, the location of the Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa file
- alleleDB, location of the homo_sapiens_merged_vep_98_GRCh38 folder
- vepSingularityImage, the location of the Singularity image

The rest of `config.yaml` contains the arguments to pass to the VEP pipeline. The important arguments have comment 
annotations. Changing an option can completely change the pipelines behavior and results in unwanted behavior (that's
difficult to debug). So please only change if you know what you are doing. All the options can be found here:
https://grch37.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic

To run the annotater, run the Annotater.py with the target mut.tsv as the argument. Several files will be created with 
various suffixes. The following are the suffixes:

* ANNO, the annotations
* log, the log of the process. Check this after the run.
* VCF and ENSEMBL, these files are the various formats the pipeline uses to convert data.

After these files are generated run the AnnotaterMerger.py with same mut.tsv argument as before. This will compile the
annotations with the original data and output the final .hmz file.

Check lines counts of the original file against the hmz to see if any data was lost. Data lost only happens when 
formatting errors occur. 100% data conversion is completely feasible.

The run time for this pipeline can be up to several hours if the data set is sufficiently large.
