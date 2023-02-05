# Dominant antibody light chain sequence from bulk RNAseq data 
  
This repository contains scripts to identify the dominant light chain antibody cell sequence from plasma cell dyscrasias. MiXCR ([Bolotin et al 2017](https://www.nature.com/articles/nmeth.3364),  [Commercial version of MiXCR](https://docs.milaboratories.com/)) is the primary sequence identification tool used.  
  
The core of the pipeline is in the script "pipeline.qsub" that takes a text file of provided dbGaP accession numbers and performs the initial processing and quality control before running MiXCR. MiXCR will identify sequences associated with IGL, IGK, IGH gene regions, perform VDJ gene annotations, collapse closely related seqeunces into clonotypes, and assign a "count" associated with each clonotype.
  
After running the MiXCR core pipeline, a second pipeline can be run "pipeline_summary_qc.qsub" that will aggregate MiXCR results and identify what, if any, the dominant clonotype is. This pipeline can also be used to produce summary figures.
    
## To run main processing steps from downloading from samples through running MiXCR, use the following script:
### pipeline.qsub ###
This script is written to be run on a linux computing cluster, such as Boston University's shared computing cluster, as an array job. To run on a different computing environment, some tweaks may need to be made to environmental variables syntax.   
At the top of the script is a list of variables in "Global variables" that contain information such as file locations that need to be provided/updated (including a text file that contains a list of dbGaP accession numbers, one sample per line.)   
This script will download samples from dbGaP, downsamples file, trims off any remaining illumina adapters (adapter sequence file location needs to be specified in "adapt" in global variables), runs FastQC (on full file, downsampled file, and trimmed file), renames files, and then runs MiXCR. MiXCR will be run on each sample individually (in an array job), and output the results separately for IGK, IGL, and IGH separately. To combine these results and continue processing, see "pipeline_summary_qc.qsub"  
Notes:  
- This pipeline starts by downloading files from dbgap. Place the location of the ngc key in the global variable "ngc_key"
    - The version of sratoolkit used appended a suffix to all filenames (e.g. "_dbGaP-28845.sra"). See "downloaded_name" variable and set up the expected name appropriately.
- Original file names will be used in FastQC results. All MiXCR results will use the new names.
- Sequencing adapters to be trimmed off should be added to "pipeline/next_tru2_tru3.fa". This file currently contains a list of commonly found adapters.  
- "new_names/rename2.py" should be run first (or your own program) to create a key that will be used to convert dbGaP accessions into a new naming scheme.
- This pipeline WILL rename SRA accession numbers to a new number  
    - "pipeline/lookupname.py" expects this key of name changes to be located in "all_newids.csv" (public repo only, see "pipeline/lookupname.py" to check and change if necessary)
    - To create your own dictionary of name changes, create a csv that has the original names in the column "Run" and the new names in the column "newname". If you wish to skip renaming files, just make these two columns identical.
    
    
## To run the main post MiXCR processing steps that aggregate results, run quality control, and identify the dominant plasma cell sequence, use the following script:
### pipeline_summary_qc.qsub ###
At the top of the script is a list of variables in "Global variables" that contain information such as file locations that need to be provided/updated. In this section there are a list of boolean global variables that control what post MiXCR analysis steps should be run.     
This is the script that will aggregate the MiXCR results and come up with summary statistics.   
This script calls "pipeline/categorize.R" that will split the samples into four categories regarding light chain sequence.   
- Category 1) samples that have a single unambiguous top clone. (Default: >=95% of light chain sequences)  
- Category 2) samples that have a single top clone (default >95% of LC sequence) if the top two clones were sufficiently similar to be collapsed together. Currently, collapse of the top two sequences happens if the top two clones have all of the following: matching CDR3 sequences, matched chain usage, and >=200bp alignment overlap (default) with 100% identity (default) and 100% coverage in the overlapped region.
- Category 3) samples that have a single major clone, but at a lower threshold. (Default: top clones >50% of LC sequences and the second clone <=5% of LC sequences.)  
- Category 4) samples that do not fall into any of the above categories. Samples in this category may be truly polyclonal, but it is more likely a result of technical issues. (e.g there is a minor mutation between the top two clones that prevented their collapse together, poor plasmablast enrichment / contamination, etc.)
The thresholds used for each of these four categories can be modified at the top of "pipeline/categorize.R".    
    
Notes:        
- Some steps require previous steps to be run. It is recommended to set all steps to true except "do_jobaccounting" and "do_standardoutput".  
- "do_jobaccounting" pulls out information related to how long each sample took to process. This MUST be set to "false" if not run on a Sun Grid Engine (or similar) computing cluster. 
- Prior to running this scrip for the first time, download IMGT reference fasta files, and run "imgtFASTAtoCSV.py" to create reference csv files of the leader scripts and VDJC genes. 
- Troubleshooting notes:
    - "blastpair_germline.py" and "blastpair_contig_germline.py" needs to import "blastpair.py". Make sure the location specified at the top of the scripts are accurate.
   
    
## To run script to identify gene coverage based on pseudo-IMGT residue numbering, run:  
### sequence_coverage/sequence_coverage_IMGT.R ###  
In order to assess sequence similarity and come up with a more universal residue index system, reference files from IMGT must be produced (see "imgt_setup.imgtFASTAtoCSV.py" section). 
Notes:  
- IMGT HighV-Quest MUST be run first.
- File locations need to be updated at the top of the scripts

## For the pipeline to run properly, reference IMGT sequences must be provided and the following run:  
### imgt_setup/imgtFASTAtoCSV.py ###  
This script will take provided FASTA reference files (with gaps) downloaded from [IMGT.org](https://www.imgt.org/) and convert it into a csv to be used for the other scripts in the pipeline. FASTA files should be downloaded from IMGT and stored in a single directory representing all V, D, J, C sequences for IGL, IGK, and IGL chains for the species. Store leadersequences in a separate directory. Reference sequences can be found at [IMGT's reference directory](https://www.imgt.org/vquest/refseqh.html).  
Both the directory for leaders sequences, and all V, D, J, and C sequences need to be specified at the top of the script (see "file_directory" variable). Run the script once to generate the reference csv for leader sequences (global variable "do_leader=True") and a separate time to generate the reference csv for the VDJC sequences (global variable "do_leader=False").  
- Program expects files to be named "IG*.fasta". For VDJC sequences, typical naming of fasta files would be "IGLV.fasta", "IGLJ.fasta", "IGLC.fasta", "IGHV.fasta", "IGHD.fasta", etc. For leader sequences, typical naming of leader sequences would be "IGL_leader.fasta" and "IGK_leader.fasta".

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
    - This version of python is required for multiqc only.  
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


