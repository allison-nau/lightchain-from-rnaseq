#!/usr/bin/env python3
"""
Program takes text file of partially processed MiXCR results, and creates a FASTA file of clones > N% frequency
"""

# Import package
from tkinter import E
import pandas as pd
import os
import os.path
import sys
import glob
from distutils import util # String to bool conversion


# Get command line arguments:
# (giant mixcr combined tsv file) (directory to store fasta files to) (include only top2?)
# (Boolean use shorter names, for programs that do not except long fasta headers)
MIXCR_INPUT = sys.argv[1]
OUT_DIR = sys.argv[2]
TOP2_BOOL = util.strtobool(sys.argv[3]) 
SHORTNAME_BOOL = TOP2_BOOL
# For testing code:
# MIXCR_INPUT = "sra/800MiXCR/SUMMARY/mixcr_all_clonotypes_v0.txt"
# OUT_DIR = "sra/800MiXCR/SUMMARY/fasta/"


class CreateFasta:
    """
    Takes text file of partially processed MiXCR results, and creates a FASTA file of clones > N% frequency
    """
    def __init__(self, text_file, min_frac=0.01, min_count=100, chains=("IGH", "IGK", "IGL"), top2only=TOP2_BOOL,
                 shortname=SHORTNAME_BOOL):
        """
        Reads in text file of MiXCR results and creates FASTA of results.
        :param text_file: string containing text file name.
        :param min_frac: minimum clone fraction for an entry to be saved (range 0 to 1). 0 will keep all entries
        (default 0.01).
        :param min_count: minimum count for a clone to be save (int). 0 will keep all entries (default 100).
        :param chains: tuple of chains to separate into different FASTA files.
        :param top2only: include only the top 2 ranked clones?
        """
        self.text_file = text_file
        self.min_frac = min_frac
        self.min_count = min_count
        self.chains = chains
        if top2only is not None:
            self.top2only = top2only
        else:
            self.top2only = False
        if shortname is not None:
            self.shortname = shortname 
        else:
            self.shortname = False
        self.split_me_times = 9  # How many substring columns there are
        self.data = self.read_in_csv(filename=text_file, sep="\t")
        self.data_dict = {}
        self.string_dict = {}
        self.samples = None
        self.separate_filter()
        self.create_fasta()

    def separate_filter(self):
        """
        Separates out the chains, separates the different samples, and filters by clone fraction
        """
        for c in self.chains:
            self.data_dict[c] = {}
            self.samples = self.data["sample"].unique()
            for s in self.samples:
                if not self.top2only:
                    self.data_dict[c][s] = self.data.loc[(self.data['Chain'] == c)
                                                        & (self.data['cloneFraction'] >= self.min_frac)
                                                        & (self.data['cloneCount'] >= self.min_count)
                                                        & (self.data["sample"] == s)]
                else:
                    # Filter:
                    self.data_dict[c][s] = self.data.loc[(self.data['Chain'] == c)
                                    & (self.data['cloneFraction'] >= self.min_frac)
                                    & (self.data['cloneCount'] >= self.min_count)
                                    & (self.data["sample"] == s)
                                    & (self.data["cloneRank"] <= 1)]

    @staticmethod
    def read_in_csv(filename, sep=","):
        """
        Reads in csv "filename" as a pandas dataframe.
        :param filename: string containing file name of the csv file to be imported.
        :param sep: character that separate fields within the csv, default ","
        :return: pandas dataframe containing data within filename.
        """
        temp_table = pd.read_csv(filename, sep=sep, header=0, low_memory=False)
        return temp_table

    def create_fasta(self, ):
        """
        Creates fasta files based on read in data, separated by sample and chain.
        """
        # Columns to use to describe row:
        # describe_row = ["sample", "Chain", "cloneId", "fileRow", "sup1", "sup2", "sup3", "sup4", "cloneFraction"]
        # describe_row = ["sample", "Chain", "cloneId", "cloneRank", "cloneFraction"]
        if not self.shortname:
            describe_row = ["sample", "Chain", "cloneId", "cloneRank", "cloneCount", "top_Vgene_short", "top_Jgene_short", "top_Cgene_short"]
            describe_prefix = ["", "", "cID", "cRank", "cCount", "V", "J", "C"]
        else:
            describe_row = ["sample", "Chain", "cloneId", "cloneRank"]
            describe_prefix = ["", "", "cID", "cRank"]           
        # Build the string:
        for c in self.chains:
            for sample in self.samples:
                # new_filename = self.text_file.replace("mixcr_all_clonotypes.txt", f"{c}.fasta")
                new_filename = OUT_DIR + sample + "_" + c + ".fasta"
                fasta_string = ""
                for i, row in self.data_dict[c][sample].iterrows():
                    row_name = row[describe_row[0]]
                    for d in range(1, len(describe_row)):
                        if describe_row[d] == "cloneCount" or describe_row[d] == "cloneId":
                            row_name = row_name + "_" + describe_prefix[d] + "-" + str(int(row[describe_row[d]]))
                        elif describe_row[d] == "sample" or describe_row[d] == "Chain":
                            row_name = row_name + "_" + str(row[describe_row[d]])
                        else:
                            row_name = row_name + "_" + describe_prefix[d] + "-" + str(row[describe_row[d]])
                    for s in range(self.split_me_times):
                        frag_col = f"subseq{s+1}"
                        if not pd.isna(row[frag_col]):
                            if len(row[frag_col]) > 0:
                                entry_name = ">" + row_name + "_" + frag_col + "\n"
                                fasta_string += entry_name
                                fasta_string += row[frag_col] + "\n"
                if len(fasta_string) > 0:
                    with open(new_filename, "wb") as myfile:
                        myfile.write(fasta_string.encode())


def main():
    """
    Main function.
    """
    # # Get list of filenames:
    # filenames = glob.glob(f'*mixcr*.txt')
    # print(filenames)
    # # Iterate through
    # for f in filenames:
    #     CreateFasta(text_file=f)
    print(f"MiXCR file: {MIXCR_INPUT}")
    print(f"Output directory: {OUT_DIR}")
    print(f"Use only top 2?: {TOP2_BOOL}")
    print(f"Use short names?: {SHORTNAME_BOOL}")
    CreateFasta(text_file=MIXCR_INPUT)


if __name__ == '__main__':
    main()
