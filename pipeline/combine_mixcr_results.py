#!/usr/bin/env python3
"""
Program combines top MiXCR Results. By default will combine all clonotypes

Program should by called:
combine_mixcr_results.py directory_of_mixcr_results_to_combine directory_to_save_results_to
"""

# Import package
import pandas as pd
import os
import glob
import numpy as np
import sys


# Get command line arguments:
# (directory where each samples MiXCR results are stored in) (directory to output files to)
IN_DIR = sys.argv[1]
OUT_DIR = sys.argv[2]
# For testing code outside command line:
# IN_DIR = "sra/800MiXCR/"
# OUT_DIR = "sra/800MiXCR/SUMMARY/"

class MixTools:
    """
    One object represents a set of MiXCR results for multiple chains.
    """
    def __init__(self, my_dir=None, min_frac=0.0, min_count=0, receptor="BCR", chains=("IGH", "IGK", "IGL"),
                 saveto="mixcr_top_clonotypes.txt", sup1=None, sup2=None, sup3=None, sup4=None, anysup=False,
                 dontshorten=False, filepattern=0, topcloneonly=False):
        """
        Initializes one MixTools object. Will output a csv of all the results combined.
        :param my_dir: string representing directory where results are stored. If None, current working directory will
        be used (default).
        :param min_frac: minimum clone fraction for an entry to be saved (range 0 to 1). 0 will keep all entries
        (default).
        :param min_count: minimum clone count for an entry to by saved. 0 will keep all entries (default).
        :param receptor: string of the receptor type run through MiXCR (will be in the filename. e.g. "BCR", "TCR").
        :param chains: tuple of chains to read in.
        :param saveto: csv filename to save final combined data to (string).
        :param anysup: boolean specifying if there is any supplmentary run info
        :param sup1: supplemntary sample run info. If not supplied, will try in infer from filename.
        :param sup2: supplemntary sample run info. If not supplied, will try in infer from filename.
        :param sup3: supplemntary sample run info. If not supplied, will not be inferred.
        :param sup4: supplemntary sample run info. If not supplied, will not be inferred.
        :param dontshorten: boolean specifying if sample name should be specified from supplementary sample info,
        default False.
        :param filepattern: expected to be an integer 0, 1, or 2. This specifies the filepattern of MiXCR results
        to be read in. filepattern == 0 indicates using the default MiXCR output format.
        :param topcloneonly: boolean specifying if only the top clone should be included in the final output. This is
        intended to be used when only a single chain is being looked at.
        """
        if my_dir is None:
            self.dir = os.getcwd()
        else:
            self.dir = my_dir
        self.filenames_dict = {}  # Dictionary to store filenames
        self.data_dict = {}  # Dictionary to store data that is read in
        self.data_to_keep = {}  # Dictionary to store data that is kept after filtering "min_frac"
        self.data_comb = {}  # Dictionary to store combined data. Will be changed into hierarchical dataframe later
        self.receptor = receptor
        self.anysup = anysup
        self.sup1 = sup1
        self.sup2 = sup2
        self.sup3 = sup3
        self.sup4 = sup4
        self.dontshorten = dontshorten
        self.filepattern = filepattern
        self.topcloneonly = topcloneonly
        self.chain_list = list(chains)  # Chains to read in
        self.min_frac = min_frac  # Minimum frequency of clonotype in order to keep it
        self.min_count = min_count  # Minimum count of clonotype in order to keep it
        self.saveto = saveto  # File to save results to
        # Read in data:
        self.read_in()
        # Combine data:
        self.summarise()

    def read_in(self):
        """
        Reads in datasets in self.dir.
        """
        # For each chain in chain list:
        for c in self.chain_list:
            # Get list of filenames:
            if self.filepattern == 0:
                self.filenames_dict[c] = glob.glob(f'{self.dir}/*.mix.{self.receptor}.clonotypes.{c}.txt')
            elif self.filepattern == 1:
                self.filenames_dict[c] = glob.glob(f'{self.dir}/*.mix.{self.receptor}.{c}.myExport.Clns.txt')
            elif self.filepattern == 2:
                self.filenames_dict[c] = glob.glob(f'{self.dir}/*.mix.{self.receptor}.{c}.custom.myExport.Clns.txt')
            # Read in all files associated with that chain:
            self.data_dict[c] = {}
            for f in self.filenames_dict[c]:
                self.data_dict[c][f] = self.read_in_csv(filename=f, sep="\t")

    @staticmethod
    def read_in_csv(filename, sep=","):
        """
        Reads in csv "filename" as a pandas dataframe.
        :param filename: string containing file name of the csv file to be imported.
        :param sep: character that separate fields within the csv, default ","
        :return: pandas dataframe containing data within filename.
        """
        temp_table = pd.read_csv(filename, sep=sep, header=0, low_memory=False)
        # If empty, add row of nan:
        if len(temp_table.index) == 0:
            # temp_table.append(pd.Series([np.nan]), ignore_index=True)
            temp_table.loc[0, :] = np.nan
            temp_table.loc[temp_table.index[0], "cloneId"] = -10
            temp_table.loc[temp_table.index[0], "cloneFraction"] = 1
            temp_table.loc[temp_table.index[0], "targetSequences"] = ""
        return temp_table

    def save_to_csv(self, df, filename=None, sep="\t", index=True):
        """
        Saves combined dataframe to csv.
        :param df: pandas dataframe to save.
        :param filename: filename to save dataframe to.
        :param sep: delimiter, default "\t".
        :param index: boolean specifying if row index should be saved. (Default: True).
        """
        if filename is None:
            filename = self.saveto
        # Save to csv:
        df.to_csv(filename, sep=sep, index=index)

    def summarise(self):
        """
        Combines all the read in dataframes into a master dataframe and saves it to a csv file.
        """
        # For each chain:
        for c in self.chain_list:
            # Keep only the rows that are above the minimum clone fraction & count:
            self.data_to_keep[c] = {}
            for k in self.data_dict[c].keys():
                self.data_to_keep[c][k] = self.data_dict[c][k][self.data_dict[c][k]["cloneFraction"] >= self.min_frac]
                self.data_to_keep[c][k] = self.data_dict[c][k][self.data_dict[c][k]['cloneCount'] >= self.min_count]
                if self.anysup:
                # Add columns with supplementary run info:
                    if self.sup4 is not None:
                        self.data_to_keep[c][k].insert(0, "sup4", self.sup4)
                    elif not self.dontshorten:
                        if len(k.split("_")) > 4:
                            sup1 = k.split("_")[4]
                            self.data_to_keep[c][k].insert(0, "sup4", sup1)
                    else:
                        sup1 = k.split("_")[0].split("\\")[-1]
                        self.data_to_keep[c][k].insert(0, "sup4", sup1)
                    if self.sup3 is not None:
                        self.data_to_keep[c][k].insert(0, "sup3", self.sup3)
                    else:
                        sup1 = k.split("_")[3]
                        self.data_to_keep[c][k].insert(0, "sup3", sup1)
                    if self.sup2 is None:
                        sup1 = k.split("_")[2]
                        sup1 = sup1.split(".")[0]
                    else:
                        sup1 = self.sup2
                    self.data_to_keep[c][k].insert(0, "sup2", sup1)
                    if self.sup1 is None:
                        sup1 = k.split("_")[1]
                    else:
                        sup1 = self.sup1
                    self.data_to_keep[c][k].insert(0, "sup1", sup1)
                # Add column with sample name (taken from filename):
                sample = k.split("/")[-1]
                if not self.dontshorten:
                    sample = sample.split("_")[0]
                if len(sample) == 0:
                    sample = k.split(".")[0]
                else:
                    sample = sample.split(".")[0]
                self.data_to_keep[c][k].insert(0, "sample", sample)
            # Combine dataframes that are the same chain:
            self.data_comb[c] = pd.concat(self.data_to_keep[c], names=["file", "cloneRank"], copy=False)
        # Combine dataframes from all the chains together.
        self.data_comb = pd.concat(self.data_comb, names=["Chain"], copy=False)
        # Add in length of resulting length info:
        # self.data_comb["mixcr_length"] = self.data_comb["targetSequences"].apply(len)
        # Length, not considering if there was a comma:
        self.data_comb.insert(6, "mixcr_length", self.data_comb["targetSequences"].apply(len))
        # Expand out:
        split_me_times = 9
        # new = self.data_comb["targetSequences"].str.split(",", n=split_me_times-1, expand=True)
        # Split string into list and sort by length, longest first:
        new = self.data_comb["targetSequences"].str.split(",", n=split_me_times - 1, expand=False).apply(
            lambda x: sorted(x, key=len, reverse=True))
        # Expand out those lists:
        new = pd.DataFrame(new.tolist(), index=new.index)
        new.fillna(value="", inplace=True)
        cols_to_max_over = []
        num_new_cols = len(new.columns)
        for s in range(split_me_times):
            new_col_name = f"subseq{s+1}"
            if s < num_new_cols:
                self.data_comb.insert(7+s, new_col_name+"len", new[s].apply(len))
                self.data_comb[new_col_name] = new[s]
            else:
                self.data_comb.insert(7 + s, new_col_name + "len", 0)
                self.data_comb[new_col_name] = ""
            cols_to_max_over.append(new_col_name+"len")
        # Maximum sublength:
        # self.data_comb.insert(7, "subseq1len", new[0].apply(len))
        # self.data_comb.insert(8, "subseq2len", new[1].apply(len))
        self.data_comb.insert(7+split_me_times, "maxSubLength", self.data_comb[cols_to_max_over].max(axis=1))
        # self.data_comb.insert(7, "mixcr_lengthMain", self.data_comb["targetSequences"].str.split(",")  temp.apply(len).apply(max(axis=1)))
        # Number of subsequences:
        self.data_comb["subseq_count"] = (self.data_comb[cols_to_max_over] > 0).sum(axis=1)
        # Save to csv
        self.save_to_csv(df=self.data_comb)
        # Read back in to remove file column:
        self.data_comb = self.read_in_csv(filename=self.saveto, sep="\t")
        self.data_comb.drop(labels="file", inplace=True, axis=1)
        # If only keeping top clone:
        if self.topcloneonly:
            self.data_comb = self.data_comb[self.data_comb["cloneRank"] == 0]
        # Resave:
        self.save_to_csv(df=self.data_comb, index=False)


def main(input_directory=None, output_directory=""):
    """
    Main function
    :param input_directory: directory where MiXCR results are.
    :param output_directory: directory to output results to.
    TODO: test output directory & patterns
    TODO: for the other patterns, need to change the column names for the most important columns
    """
    mymix = MixTools(my_dir=input_directory, saveto=OUT_DIR + "mixcr_all_clonotypes_v0.txt", min_frac=0, filepattern=0)
    mymix = MixTools(my_dir=input_directory, saveto=OUT_DIR + "mixcr_over_100_v0.txt", min_count=100, filepattern=0)
    # mymix = MixTools(my_dir=input_directory, saveto=OUT_DIR + "mixcr_all_clonotypes_v1.txt", min_frac=0, filepattern=1)
    # mymix = MixTools(my_dir=input_directory, saveto=OUT_DIR + "mixcr_all_clonotypes_v2.txt", min_frac=0, filepattern=2)
    # Get top clone for each chain:
    for chain in ["IGH", "IGK", "IGL"]:
        mymix = MixTools(my_dir=input_directory, saveto=OUT_DIR + "mixcr_" + chain + "_topclone_v0.txt", min_frac=0,
                         filepattern=0, topcloneonly=True)
        # mymix = MixTools(my_dir=input_directory, saveto=OUT_DIR + "mixcr_" + chain + "_topclone_v2.txt", min_frac=0,
        #                  filepattern=2, topcloneonly=True)
        # mymix = MixTools(my_dir=input_directory, saveto=OUT_DIR + "mixcr_" + chain + "_topclone_v1.txt", min_frac=0,
        #                  filepattern=1, topcloneonly=True)


if __name__ == '__main__':
    main(input_directory=IN_DIR)
