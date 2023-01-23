#!/usr/bin/env python3
"""
Compares two sequences and calculated the coverage and identity. Identity is determined only for the coverage region.

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


# Import packages
# May need to install & load biopython first
from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import pairwise2  
from Bio.Align import MultipleSeqAlignment 
from Bio.Align import AlignInfo 
import pandas as pd
import sys
import glob


# blast executable installed:
blastplus = "/share/pkg.7/blast+/2.12.0/install/bin/blastn"
# blastprotein = "/share/pkg.7/blast+/2.12.0/install/bin/blastp"

# Get command line arguments:
# (directory fasta files are in) (output csv to store results in)
FASTA_DIR = sys.argv[1]
OUT_FILE = sys.argv[2]
# FASTA_DIR = "test_directory/"
# OUT_FILE = "sra/800MiXCR/SUMMARY/blast_results.csv"

# To redirect printing to a different file:
# STDOUT = "sra/800MiXCR/SUMMARY/blast_print.txt"
# Clear standard output text file of contents:
# open(STDOUT, "w").close()
# Redirect standard output to text file:
# sys.stdout = open(STDOUT, "wt")


def get_sequence(record_id, all_sequences):
    """
    Gets the sequence for a specific record in a fasta file.
    :param record_id: ID of desired sequence (str).
    :param all_sequences: list of all sequence records.
    :return: string of desired sequence
    """
    # Initialize variables:
    desired_string = ""
    # Look up sequence id:
    for sequence in all_sequences:
        if sequence.id == record_id:
            desired_string = str(sequence.seq)
    return desired_string


def get_all_sequences_from_fasta(fasta_file):
    """
    Gets all sequences from fasta file.
    :param fasta_file: fasta file path (str).
    :return: list of all sequence records.
    """
    # Initialize variables:
    sequences = []
    # Get fasta sequences:
    all_sequences = SeqIO.parse(open(fasta_file), "fasta")
    for sequence in all_sequences:
        sequences.append(sequence)
    return sequences


def run_pair_blast(fasta1, fasta2):
    """
    Computes identity and coverage for two sequences.
    :param fasta1: fasta file containing first sequence to be compared (query).
    :param fasta2: fasta file containing second sequence to be compared (subject).
    :return: pandas dataframe containing alignment information, including the following information:
        ["alignment_title", "query", "subject",
          "alignment_count", "query_length",
          "subject_length", "score", "e_value",
          "alignment_length", "query_start", "subject_start",
          "query_end", "subject_end",
          "identities", "identity_percent",
          "extension_available_from_query", "extension_available_from_subject",
          "query_coverage_percent", "subject_coverage_percent", "gaps",
          "query_seq", "match_seq", "subject_seq",
          "alignment_string", "assembled_string", "original_subject_seq", "original_query_seq"]
    """
    # Column headers:
    column_headers = ["alignment_title", "query", "subject",
                      "alignment_count", "query_length",
                      "subject_length", "score", "e_value",
                      "alignment_length", "query_start", "subject_start",
                      "query_end", "subject_end",
                      "identities", "identity_percent",
                      "extension_available_from_query", "extension_available_from_subject",
                      "query_coverage_percent", "subject_coverage_percent", "gaps",
                      "query_seq", "match_seq", "subject_seq",
                      "alignment_string", "assembled_string", "original_subject_seq", "original_query_seq"] 
    # Initialize list of dictionaries to store results:
    results = []
    # Get all sequence records out of fasta files:
    all_sequences1 = get_all_sequences_from_fasta(fasta_file=fasta1)
    all_sequences2 = get_all_sequences_from_fasta(fasta_file=fasta2)
    # Run BLAST and parse the output as XML
    output = NcbiblastnCommandline(cmd=blastplus, query=fasta1, subject=fasta2, outfmt=5)()[0]
    blast_result_record = NCBIXML.read(StringIO(output))
    # Number of alignments:
    alignment_count = len(blast_result_record.alignments)
    result_entry = {"query": blast_result_record.query,
                    "query_length": blast_result_record.query_letters,
                    "alignment_count": alignment_count}
    if alignment_count == 0:
        results.append(result_entry)
    else:
        print(f"Number of alignments: {alignment_count}")
    # Print some information on the result
    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:  # High Scoring Pair
            print(f"Alignment title: {alignment.title}")
            print(f"Query: {blast_result_record.query}")
            print(f"Query length: {blast_result_record.query_letters}")
            print(f"Subject: {alignment.hit_id}")
            print(f"Subject Length: {alignment.length}")
            print(f"Score: {hsp.score}")
            print(f"e value: {hsp.expect}")
            print(f"Alignment length: {hsp.align_length}")
            print(f"Identities: {hsp.identities}")
            identity_percent = hsp.identities/hsp.align_length*100
            print(f"% Identity (identities / alignment length): {identity_percent:.2f}%")
            query_coverage = hsp.align_length/blast_result_record.query_letters * 100
            subject_coverage = hsp.align_length/alignment.length*100
            print(f"Query Coverage (alignment length / query length): {query_coverage:.2f}%")
            print(f"Subject Coverage (alignment length / subject length): {subject_coverage:.2f}%")
            print(f"Gaps: {hsp.gaps}")
            print(f"Query start: {hsp.query_start}")
            print(f"Subject start: {hsp.sbjct_start}")
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
            result_entry["alignment_title"] = alignment.title
            result_entry["subject"] = alignment.hit_id
            result_entry["subject_length"] = alignment.length
            result_entry["score"] = hsp.score
            result_entry["e_value"] = hsp.expect
            result_entry["alignment_length"] = hsp.align_length
            result_entry["query_start"] = hsp.query_start
            result_entry["subject_start"] = hsp.sbjct_start
            result_entry["query_end"] = hsp.query_start + hsp.align_length - 1
            result_entry["subject_end"] = hsp.sbjct_start + hsp.align_length - 1
            result_entry["extension_available_from_query"] = result_entry["query_length"] - hsp.align_length
            result_entry["extension_available_from_subject"] = result_entry["subject_length"] - hsp.align_length
            result_entry["identities"] = hsp.identities
            result_entry["identity_percent"] = identity_percent
            result_entry["query_coverage_percent"] = query_coverage
            result_entry["subject_coverage_percent"] = subject_coverage
            result_entry["gaps"] = hsp.gaps
            result_entry["query_seq"] = hsp.query
            result_entry["match_seq"] = hsp.match
            result_entry["subject_seq"] = hsp.sbjct
            result_entry["alignment_string"] = str(alignment)
            # Build a consensus string:
            can_build_consensus = True
            # Get relevant indices:
            # (Reminder: slice notation means that python starts -1 from end of slice, so pay attention to when using index vs BLAST notation)
            query_start_index = result_entry["query_start"] - 1
            subject_start_index = result_entry["subject_start"] - 1
            query_end_index = result_entry["query_end"] - 1
            subject_end_index = result_entry["subject_end"] - 1
            # Get full sequences:
            subject_seq = get_sequence(record_id=result_entry["subject"], all_sequences=all_sequences2)
            query_seq = get_sequence(record_id=result_entry["query"], all_sequences=all_sequences1)
            result_entry["original_subject_seq"] = subject_seq
            result_entry["original_query_seq"] = query_seq
            print(f"Full Query:\n{query_seq}")
            print(f"Full Subject:\n{subject_seq}")
            # # Get string version (strip non-alphabet characters):
            # query_seq = "".join(char for char in hsp.query if char.isalpha())
            # subject_seq = "".join(char for char in hsp.sbjct if char.isalpha())
            # Don't build consensus if one of the sequences is not at edge of alignment:
            if (query_start_index != 0) & (subject_start_index != 0):
                can_build_consensus = False
            if (result_entry["query_end"] != result_entry["query_length"]) & (result_entry["subject_end"] != result_entry["subject_length"]):
                can_build_consensus = False 
            # Grab the consensus pieces:
            # Shared middle:
            query_middle = query_seq[query_start_index:result_entry["query_end"]]
            subject_middle = subject_seq[subject_start_index:result_entry["subject_end"]]
            if query_middle != subject_middle:
                can_build_consensus = False
            print(f"Query   middle: {query_middle}")
            print(f"Subject middle: {subject_middle}")
            if query_start_index > subject_start_index:
                consensus_start = query_seq[:query_start_index]
            elif query_start_index < subject_start_index: 
                consensus_start = subject_seq[:subject_start_index] 
            elif result_entry["query_start"] == result_entry["subject_start"]:
                consensus_start = ""
            else:
                consensus_start = ""
            if result_entry["subject_end"] < result_entry["subject_length"]:
                consensus_end = subject_seq[result_entry["subject_end"]:]
            elif result_entry["query_end"] < result_entry["query_length"]:
                consensus_end = query_seq[result_entry["query_end"]:]
            else:
                consensus_end = ""
            result_entry["assembled_string"] = consensus_start + subject_middle + consensus_end
            print(f"Consensus: {result_entry['assembled_string']}")
            # If consensus string didn't meet requirements, zero out:
            if not can_build_consensus:
                result_entry["assembled_string"] = "NotAvailable"
            results.append(result_entry)
    results = pd.DataFrame(results, columns=column_headers)
    return results


def run_multiple_pairs(fastafile):
    """
    Runs BLASTN for all pairs of sequences in "fastafile". This is done in both directions (i.e. each pair is compared
    twice). (Note: self comparisons are skipped.)
    :param fastafile: path to fasta file.
    :return: dataframe of blast results.
    """
    # Initialize dataframe placeholder:
    df = None
    df_temp = None
    # Temporary fasta filenames:
    temp1 = "temp1_forblast.fasta"
    temp2 = "temp2_forblast.fasta"
    # Read in original fastafile:
    fasta_org = list(SeqIO.parse(fastafile, "fasta"))
    # Number of sequences:
    seq_count = len(fasta_org)
    # Only do BLAST if there are more than 1 sequences:
    if seq_count > 0:
        # Compare all pairs (in both directions:):
        # SKIP the self comparison
        for i in range(seq_count):
            for j in range(seq_count):
                if i != j:
                    # Create temporary fasta files:
                    SeqIO.write(fasta_org[i], temp1, "fasta")
                    SeqIO.write(fasta_org[j], temp2, "fasta")
                    # Run paired BLAST:
                    df_temp = run_pair_blast(temp1, temp2)
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


def run_all_files(indir, outfile):
    """
    Runs BLASTN for all fasta files in "indir".
    :param indir: directory where fasta files are located.
    :param outfile: csv file to save results to.
    """
    # Initialize list of to store resulting pandas dataframes:
    df_list = []
    # Get list of files:
    list_of_samples = glob.glob(f'{indir}/*.fasta')
    # Iterate through list and do pairwise blast:
    for fasta_file in list_of_samples:
        df_list.append(run_multiple_pairs(fasta_file))
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
    # Subject:
    temp_subject = df["subject"].str.split("_", expand=True)
    df["meta_subject_name"] = temp_subject.iloc[:, 0]
    df["meta_subject_chain"] = temp_subject.iloc[:, 1]
    df["meta_subject_cID"] = temp_subject.iloc[:, 2].str.split("-", expand=True).iloc[:, 1]
    df["meta_subject_cRank"] = temp_subject.iloc[:, 3].str.split("-", expand=True).iloc[:, 1]
    df["meta_subject_cCount"] = temp_subject.iloc[:, 4].str.split("-", expand=True).iloc[:, 1]
    df["meta_subject_V"] = temp_subject.iloc[:, 5].str.split("-", expand=True).iloc[:, 1]
    df["meta_subject_J"] = temp_subject.iloc[:, 6].str.split("-", expand=True).iloc[:, 1]
    df["meta_subject_C"] = temp_subject.iloc[:, 7].str.split("-", expand=True).iloc[:, 1]
    df["meta_subject_subseq"] = temp_subject.iloc[:, 8]
    # Save dataframe:
    df.to_csv(outfile)


def test_code():
    """
    Test code for calculating identity and coverage. Contains sequences that should have a lot of overlap, as well as a
    sequence that should have no overlap.
    """
    # Create some sequences to compare:
    clone0_s1 = SeqRecord(Seq("GTCCAACAGGGCCACTGGCGTCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACAGAGTTCACTCTAAACATCAGGACTCTGG"
                              "AGCCTGAAGATTTTGCAGTGTATTATTGTCAGCAGTTTGGTGTCTCACGTCCGTGGACGTTCGGCCAAGGGACCAAGGTGGAA"
                              "ATGAATCGAACTGTGGCTGCACCATCTGTCTTCATCTTCCCGCCATCTGATGAGCAGTTGAAATCTGGAACTGCCTCTGTTGT"
                              "GTGCCTGCTGAATAACTTCTATCCCAGAGAGGCCAAAGTACAGTGGAAGGTGGATAACGCCCTCCAATCGGGTAACTCCCAGGA"
                              "GAGTGTCACAGAGCAGGACAGCAAGGACA"), id="clone0_s1")
    clone0_s2 = SeqRecord(Seq("TGCTCAGTTAGGACCCAGAGGGAACCATGGAAACCCCAGCGCAGCTTCTCTTCCTCCTGCTACTCTGGCTCCCAGATACCACC"
                              "GGAGAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCGCCCTCTCCTGCAGGGCCAGTCA"
                              "GAGTCTTAGCAGCAGCTTCTTAGCCTGGTACCAGCAGAAACC"), id="clone0_s2")
    clone1_s1 = SeqRecord(Seq("AGAAGAGCTGCTCAGTTAGGACCCAGAGGGAACCATGGAAACCCCAGCGCAGCTTCTCTTCCTCCTGCTACTCTGGCTCCCAGA"
                              "TACCACCGGAGAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCGCCCTCTCCTGCAGGGC"
                              "CAGTCAGAGTCTTAGCAGCAGCTTCTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTGTGTC"
                              "CAACAGGGCCACTGGCGTCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACAGAGTTCACTCTAAACATCAGGACTCTGGAGCC"
                              "TGAAGATTTTGCAGTGTATTATTGTCAGCAGTTTGGTGTCTCACGTCCGTGGACGTTCGGCCAAGGGACCAAGGTGGAAATGAAT"
                              "CGAACTGTGGCTGCACCATCTGTCTTCATCTTCCCGCCATCTGATGAGCAGTTGAAATCTGGAACTGCCTCTGTTGTGTGCCTGCT"
                              "GAATAACTTCTATCCCAGAGAGGCCAAAGTACAGTGGAAGGTGGATAACGCCCTCCA"
                              "ATCGGGTAACTCCCAGGAGAGT"), id="clone1_s1")
    SeqIO.write(clone0_s1, "clone0_s1.fasta", "fasta")
    SeqIO.write(clone0_s2, "clone0_s2.fasta", "fasta")
    SeqIO.write(clone1_s1, "clone1_s1.fasta", "fasta")
    # Run BLAST
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("No expected similarity:")
    output0s1_0s2 = run_pair_blast("clone0_s1.fasta", "clone0_s2.fasta")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Expect much similarity:")
    output0s1_1s1 = run_pair_blast("clone0_s1.fasta", "clone1_s1.fasta")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Expect much similarity, reverse direction:")
    output1s1_0s1 = run_pair_blast("clone1_s1.fasta", "clone0_s1.fasta")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Expect much similarity, smaller chunk compared with longer:")
    output0s2_1s1 = run_pair_blast("clone0_s2.fasta", "clone1_s1.fasta")
    # Combine results into one pandas dataframe
    df = pd.concat([output0s1_0s2, output0s1_1s1, output1s1_0s1, output0s2_1s1], axis=0, ignore_index=True)
    df.to_csv("blast_result.csv")


def test_code2():
    """
    Test code for testing alignment results.
    :return: 
    """
    fasta_file = "MMRF100639_IGL.fasta"
    df = run_multiple_pairs(fasta_file)
    df.to_csv("blast_results2.csv")


def main(fasta_dir, out_file):
    """
    Runs all pairwise blast for each fasta file in "fasta_dir".
    :param fasta_dir: directory where are fasta files are located.
    :param out_file: file name to save csv combined output to.
    """
    # Run all files:
    run_all_files(indir=fasta_dir, outfile=out_file)


if __name__ == '__main__':
    # Run test code:
    # test_code()
    # test_code2()
    main(fasta_dir=FASTA_DIR, out_file=OUT_FILE)
