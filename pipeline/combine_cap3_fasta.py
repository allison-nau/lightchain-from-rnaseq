#!/usr/bin/env python3
"""
Combines CAP3 contig fasta results into a master fasta file. (In preparation to submit to IGMT)
"""

# Import modules:
import os
import sys
import glob

# Get command line arguments:
# (directory fasta files are in) (output fasta to store results in)
FASTA_DIR = sys.argv[1]
OUT_FILE = sys.argv[2]


# --------------------------------------------------------------------------------------------------------------------
def read_fasta(filepath: str):
    """
    Reads in fasta file and returns string. Fasta filename is adding to the beginning of each sequence header.
    :param filepath: fasta file location.
    :return: string representing results of fasta file.
    """
    # Initialize string to store contents:
    built_string = ""
    # Grab relevant filename header:
    header = str(str(filepath.split(sep="/")[-1]).split(sep=".")[0])
    with open(filepath, "r") as myfile:
        for line in myfile:
            if line[0] == ">":
                new_line = ">" + header + "_" + line[1:]
                built_string = built_string + new_line
            else:
                built_string = built_string + line
    return built_string


# --------------------------------------------------------------------------------------------------------------------
# Get list of fasta files:
fastafile_list = glob.glob(f"{FASTA_DIR}/*.fasta.cap.contigs")

# Initialize master string representing fasta contents:
master_string = ""
# For each file:
for file in fastafile_list:
    contents = read_fasta(file)
    master_string = master_string + contents

# Save contents:
with open(OUT_FILE, "w") as outfile:
    outfile.write(master_string)
