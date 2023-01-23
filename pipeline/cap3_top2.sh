#!/bin/bash -l

# does CAP3 assembly over all fasta files in directory:
# (location of CAP3 assembly program) (directory where fasta files are)

# Directory with fasta files:
# in_dir=sra/800MiXCR/SUMMARY/fasta_top2_play20220502_1/
in_dir=$2

# Files:
files=${in_dir}*.fasta

echo "Number of files: "${#files[@]}

# CAP3 program:
CAP3=$1
# CAP3=tools/cap3/CAP3/cap3

# Loop through and run cap3 on all:
for f in $files
do 
    echo "Process $f file..."
    ${CAP3} ${f} > ${f}.out
done