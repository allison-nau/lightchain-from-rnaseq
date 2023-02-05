#!/usr/bin/env python3
"""
Program converts IMGT fasta files into csv reference files. 
Update file directories to point to reference FASTA files, and run twice with "do_leader=True" to create reference leader sequences and "do_leader=False" to create the main VDJC reference file.
"""

# Import package
import pandas as pd
import os
import glob
import numpy as np
import sys
from Bio import SeqIO

# Globals:
# do_leader = True
do_leader = False

if not do_leader:
    # Gene V D J and C:
    file_directory = "20220402_IMGT/"  # Include trailing slash
else:
    # V leader sequences:
    file_directory = "20220418_IMGT_leadersequences/"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def fasta_to_df(fastafile):
    """
    Reads in fasta file and returns a pandas dataframe, with sequence header split into multiple columns.
    :param fastafile: fasta file name (str).
    :return: pandas dataframe.
    """
    # Create new file name:
    csvfile = fastafile.replace("fasta", "csv")
    # Read in fasta file:
    records = SeqIO.parse(fastafile, "fasta")
    rows = []
    for record in records:
        rows.append([record.id, record.description, str(record.seq)])
    # Convert to pandas dataframe:
    df = pd.DataFrame(rows, columns=["id", "description", "seq"])
    # Convert ID into multiple columns:
    # Example Fasta header:
    # >Z73653|IGLV1-36*01|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |296+45=341| | |
    df[["Accession", "Allele", "Species", "Functional",
        "h01", "h02", "h03"]] = df["description"].str.split("|", n=6, expand=True)
    # Split Allele Group
    df[["AlleleGroup", "SpecificAllele"]] = df["Allele"].str.split("*", n=1, expand=True)
    # Get Locus:
    df["Locus"] = df["Allele"].str[0:3]
    # V vs D vs J:
    df["GeneSegment"] = df["Allele"].str[4]
    # Allele Group no chain info:
    df["AlleleGroupNoChain"] = df["AlleleGroup"].str[3:]
    # AlleleGroup Number only:
    df["AlleleGroupNum"] = df["AlleleGroup"].str[4:]
    # Remove gap from sequence:
    df["SeqNoGaps"] = df["seq"].str.replace(pat="\.", repl="", regex=True)
    # df["SeqNoGaps"] = df["seq"].apply(lambda x: x.str.replace(".", repl=""))
    return df, csvfile


# Get list of filenames:
filenames = glob.glob(f"{file_directory}IG*.fasta")
print(f"Files read in:\n{filenames}\n")
# Convert FASTA files to CSV
for f in filenames:
    df, csvfile = fasta_to_df(f)
    # Save to csv:
    df.to_csv(csvfile, index=False)
    # Combine to master dataframe:
    if "large_df" in locals():
        large_df = pd.concat([large_df, df], axis=0)
    else:
        large_df = df

# Save combined dataframe:
if not do_leader:
    large_df.to_csv(file_directory + "IMGTref.csv", index=False)
else:
    large_df.to_csv(file_directory + "IMGTleadersequences.csv", index=False)
