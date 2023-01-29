#!/usr/bin/env python3
"""
Compares two sequences and calculated the coverage and identity. Identity is determined only for the coverage region.

This version focuses on combining with germline reference.

Biopython documentation:
https://biopython.org/docs/1.75/api/index.html
# Tutorial:
http://biopython.org/DIST/docs/tutorial/Tutorial.html
Forum with example:
https://www.biostars.org/p/42687/
# Blast options and default:
https://www.arabidopsis.org/Blast/BLASToptions.jsp

To see options for blastn, in command line on the SCC run:
module load blast+/2.12.0
blastn -help
"""

# Import modules:
import os

# Load pairwise blast files:
cwd = os.getcwd()
os.chdir("pipeline/")
from blastpair import *
os.chdir(cwd)

# Get command line arguments:
# (directory mixcr sample fasta files are in) (directory germline files are in) (output csv to store results in)
FASTA_DIR = sys.argv[1]
GERM_DIR = sys.argv[2]
OUT_FILE = sys.argv[3]


def run_pairs(query_fasta, subject_fasta):
    """
    Runs BLASTN for sequences in "query_fasta" against sequences in "subject_fasta".
    :param query_fasta: path to query fasta file.
    :param subject_fasta: path to sample subject fasta file.
    :return: dataframe of blast results.
    """
    # Initialize dataframe placeholder:
    df = None
    df_temp = None
    # Temporary fasta filenames:
    temp_query = "temp1_forblast.fasta"
    temp_subj = "temp2_forblast.fasta"
    # Read in original fasta files:
    fasta_query_org = list(SeqIO.parse(query_fasta, "fasta"))
    fasta_subj_org = list(SeqIO.parse(subject_fasta, "fasta"))
    # Number of sequences:
    query_count = len(fasta_query_org)
    subj_count = len(fasta_subj_org)
    # Do BLAST:
    # Compare all pairs:
    for q in range(query_count):
        for s in range(subj_count):
            # Create temporary fasta files:
            SeqIO.write(fasta_query_org[q], temp_query, "fasta")
            SeqIO.write(fasta_subj_org[s], temp_subj, "fasta")
            # Run paired BLAST:
            df_temp = run_pair_blast(fasta1=temp_query, fasta2=temp_subj)
            # Save if alignment count > 0:
            if df_temp.loc[0, 'alignment_count'] > 0:
                if "df" not in locals():
                    df = df_temp
                else:
                    df = pd.concat([df, df_temp])
    # Save whatever is last if there were no alignments at all:
    if "df" not in locals():
        df = df_temp
    return df


def run_all_files(indir, germdir, outfile):
    """
    Runs BLASTN for all fasta files in "indir".
    :param indir: directory where fasta files are located.
    :param outfile: csv file to save results to.
    """
    # Initialize list of to store resulting pandas dataframes:
    df_list = []
    # Get list of files:
    # TODO remove  list_of_samples = glob.glob(f'{indir}/*.fasta')
    list_of_samples = glob.glob(f'{germdir}/*.fasta')
    # Iterate through list and do pairwise blast:
    for germline_file in list_of_samples:
        # Original filename:
        fasta_file = indir + germline_file.split(sep="/")[-1]
        fasta_file = fasta_file.replace("_germline.fasta", ".fasta")
        df_list.append(run_pairs(query_fasta=fasta_file, subject_fasta=germline_file))
    # Concatenate all pandas dataframes:
    df = pd.concat(df_list, axis=0)
    # Create seperate columns for sample, clone id, clone rank, subsequence, etc for both query and subject
    # Query:
    temp_query = df["query"].str.split("_", expand=True)
    df["meta_query_name"] = temp_query.iloc[:, 0]
    df["meta_query_chain"] = temp_query.iloc[:, 1]
    df["meta_query_cID"] = temp_query.iloc[:, 2].str.split("-", expand=True).iloc[:, 1]
    df["meta_query_cRank"] = temp_query.iloc[:, 3].str.split("-", expand=True).iloc[:, 1]
    df["meta_query_cCount"] = temp_query.iloc[:, 4].str.split("-", expand=True).iloc[:, 1]
    df["meta_query_V"] = temp_query.iloc[:, 5].str.split("-", expand=True).iloc[:, 1]
    df["meta_query_J"] = temp_query.iloc[:, 6].str.split("-", expand=True).iloc[:, 1]
    df["meta_query_C"] = temp_query.iloc[:, 7].str.split("-", expand=True).iloc[:, 1]
    df["meta_query_subseq"] = temp_query.iloc[:, 8]
    # TODO remove! temporarily save intermediate file:
    # TODO remove df.to_csv("blast_temp_intermediate.csv")
    # Subject:
    df["meta_subject_name"] = df["subject"]
    # Save dataframe:
    df.to_csv(outfile)


def main(fasta_dir, germ_dir, out_file):
    """
    Runs all pairwise blast for each fasta file in "fasta_dir".
    :param fasta_dir: directory where are sample fasta files are located.
    :param germ_dir: directory where sample specific germline fasta files are located.
    :param out_file: file name to save csv combined output to.
    """
    # Run all files:
    run_all_files(indir=fasta_dir, germdir=germ_dir, outfile=out_file)


if __name__ == '__main__':
    main(fasta_dir=FASTA_DIR, germ_dir=GERM_DIR, out_file=OUT_FILE)
