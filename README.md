# Dominant antibody light chain sequence from bulk RNAseq data 
  
This repository contains scripts to identify the dominant light chain antibody cell sequence from plasma cell dyscrasias. MiXCR ([Bolotin et al 2017](https://www.nature.com/articles/nmeth.3364),  [Commercial version of MiXCR](https://docs.milaboratories.com/)) is the primary sequence identification tool used.  
  
The core of the pipeline is in the script "pipeline.qsub" that takes a text file of provided dbGaP accession numbers and performs the initial processing and quality control before running MiXCR. MiXCR will identify sequences associated with IGL, IGK, IGH gene regions, perform VDJ gene annotations, collapse closely related seqeunces into clonotypes, and assign a "count" associated with each clonotype.
  
After running the MiXCR core pipeline, a second pipeline can be run "pipeline_summary_qc.qsub" that will aggregate MiXCR results and identify what, if any, the dominant clonotype is. This pipeline can also be used to produce summary figures.
    
## To run main processing steps from downloading from samples through running MiXCR, use the following script:
### pipeline.qsub ###
This script is written to be run on a linux computing cluster, such as Boston University's shared computing cluster, as an array job. To run on a different computing environment, some tweaks may need to be made to environmental variables syntax. See end of read me   
At the top of the script is a list of variables in "Global variables" that contain information such as file locations that need to be provided/updated (including a text file that contains a list of dbGaP accession numbers, one sample per line.)   
This script will download samples from dbGaP, downsamples file, trims off any remaining illumina adapters (adapter sequence file location needs to be specified in "adapt" in global variables), runs FastQC (on full file, downsampled file, and trimmed file), renames files to be provided to MiXCR, and then runs MiXCR.  
Notes:  
- Original file names will be used in FastQC results. All MiXCR results will be used following the new names.
- Sequencing adapters to be trimmed off should be added to "pipeline/next_tru2_tru3.fa". This file currently contains a list of commonly found adapters.  
- "new_names/rename2.py" should be run first (or your own program) to create a key that will be used to convert dbGaP accessions into a new naming scheme.
    
    
## To run the main post MiXCR processing steps that run quality control as well as identify the dominant plasma cell sequence, use the following script:
### pipeline_summary_qc.qsub ###
At the top of the script is a list of variables in "Global variables" that contain information such as file locations that need to be provided/updated. In this section there are a list of boolean global variables that control what post MiXCR analysis steps should be run.     
Notes:     
- Some steps require previous steps to be run  
- "do_jobaccounting" pulls out information related to how long each sample took to process. This MUST be set to "false" if not run on a Sun Grid Engine (or similar) computing cluster. 
- Prior to running this scrip for the first time, download IMGT reference fasta files, and run imgtFASTAtoCSV.py to create reference csv files of the leader scripts and VDJC genes. 
- Troubleshooting notes:
    - "blastpair_germline.py" and "blastpair_contig_germline.py" needs to import "blastpair.py". Make sure the location specified at the top of the scripts are accurate.
   
    
## To run script to identify gene coverage based on pseudo-IMGT residue numbering, run:
### sequence_coverage/sequence_coverage_IMGT.R ###
Notes:  
- IMGT HighV-Quest MUST be run first.
- File locations need to be updated at the top of the scripts
   
## To randomly generate new sample names to be used when running MiXCR, run:
### new_names/rename2.py ###
Notes:  
- To make this reproducible, a random number seed must be set prior to using random.shuffle. It is recommended to NOT save the seed used anywhere public.
- If you want to leave names unchanged, replace the new names in the generated csv file with the original names.
   
## To produce final tables, adjust variables at tops of scripts and run (This will eventually be added to pipeline/pipeline_summary.qc.qsub):
### pipeline/final_table.py ###
    
    
## This program requires the following tools to be installed/available:  
Requirements for "pipeline.qsub":   
- python3/3.8.10
- sratoolkit/2.11.1
- fastqc/0.11.7
- seqtk/1.3
- trimmomatic/0.39
- mixcr/3.0.13
   
Requirements for "pipeline_summary_qc.qsub":   
- python3/3.8.10 (main python version)
    - Packages: 
        - [Bio](https://biopython.org/docs/1.75/api/index.html)
        - pandas, numpy, tkinter
        - Base python: os, glob, sys, time, re, distutils, io, math, random
    - blastn executable of blast+/2.12.0 needs to be installed. 
        - Location needs to be specified in "blastpair.py" in "blastplus" variable
- python3/3.7.9 
    - This version of python is required for multiqc
- multiqc/1.10.1
- R/4.0.5
    - Libraries: 
        - Biostrings, msa, stringr
        - epitools, nnet
        - tidyverse, ggpubr, ggplot2, GGally, janitor, naniar
- CAP3 executable
    - ([Huang et al 1999](https://pubmed.ncbi.nlm.nih.gov/10508846/))
    - Location should be specified in the variable "cap3_prog" in "pipeline_summary_qc.qsub"  
   
Additional requirements not listed above:  
- R/4.0.5
    - Libraries:
        - stringi, forcats


