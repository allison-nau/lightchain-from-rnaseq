#!/usr/bin/env python3

"""
Creates final output tables for each sample category, split by chain.

Category 1: Ready to go, >=95% of total LC is the top clone.
Category 2: Top 2 clones together. Chain usage match & 100% identity & >=200bp alignment & >=95% LC is top2.
Category 3: Major clone present, but below threshold. >50% First clone & <=5% Second clone among all LC.
Category 4: Samples not in any of the above categories.

Make sure filenames are pointing to correct locations.
"""
# Load libraries:
import pandas
import pandas as pd
import numpy as np


# Filenames:
COMBINED_DF_FILE = "sra/800MiXCR/SUMMARY/mixcrstats_topclones_blast_combined_correctedFractions_withCategorizations.txt"
STITCHED_FILE = "sra/800MiXCR/SUMMARY/stitched_and_padded.csv"
OUT_DIR = "sra/800MiXCR/SUMMARY/"
MIXCR_OUTPUT_FILE = "sra/800MiXCR/SUMMARY/mixcr_all_clonotypes_v0.txt"


class ProduceFinalOutputs:
    """
    Class produces final output tables.
    """
    def __init__(self, combined_dff: str, mixcr_file: str, stitched_file: str, outdir: str):
        """
        Initializes on ProduceFinalOutputs object. (And produces the final outputs.)
        :param combined_dff: file path/name of combined MiXCR stats table.
        :param cat1f: filename of category 1 sample list.
        :param cat2f: filename of category 2 sample list.
        :param cat3f: filename of category 3 sample list.
        :param cat4f: filename of category 4 sample list.
        :param mixcr_file: filename of combined MiXCR output file.
        :param outdir:  filepath of directory to save results to. (Note: Needs to have a trailing "/")
        """
        self.cat_descriptions = {1: "ReadyToGo",
                                 2: "Top2Together",
                                 3: "MajorClonepresentbelowthreshold",
                                 4: "uncategorizedsamples"}
        self.cat_substring = "Category"
        self.remove_categories = ["Category4_uncategorizedsamples"]
        self.outdir = outdir
        # Read in data:
        self.combined_df = self.read_in_csv(combined_dff, sep="\t")
        # print(f"Stats dataframe size: {self.combined_df.shape}")
        self.mixcr_df = self.read_in_csv(mixcr_file, sep="\t")
        # print(f"MiXCR dataframe size: {self.mixcr_df.shape}")
        self.stitched_df = self.read_in_csv(stitched_file, sep=",")
        self.stitched_df = self.stitched_df[["sample", "stitched", "stitched_padded",
                                             "piece1_seq", "piece2_seq", "piece3_seq"]]
        # Combined mixcr_df with combined_df, keeping only top row:
        self.combined_df = pd.merge(self.mixcr_df, self.combined_df, how="right",
                                    right_on=["sample", "domClone_cloneId"],
                                    left_on=["sample", "cloneId"])
        # print(f"Combined dataframe size: {self.combined_df.shape}")
        # Combined with stitched:
        self.combined_df = pd.merge(self.combined_df, self.stitched_df, how="left", on="sample")
        print(f"Combined dataframe size, after adding stitched: {self.combined_df.shape}")
        self.add_in_heavy_chain()
        self.initial_processing()
        # # Category dataframes:
        # self.cat_dfs = {}
        # for i in range(1, 5):
        #     self.cat_dfs[i] = self.process_cat(category=i)
        #     self.save_final_table(category=i)
        # Create fasta files:
        self.create_fasta(out_file=outdir+"finaltable.fasta", stitched=False)
        self.create_fasta(out_file=outdir+"finaltable_stitched_forIMGTonly.fasta", stitched=True)
        # Save table with and without category 4
        self.combined_df.to_csv(self.outdir + "finaltable_allcategories.tsv", sep="\t", index=False)
        self.combined_df[~self.combined_df["CATEGORY"].isin(self.remove_categories)].to_csv(self.outdir + "finaltable_cat1_2_3.tsv", sep="\t", index=False)

    def add_in_heavy_chain(self):
        """
        Adds in heavy chain information to final table results.
        """
        # Grab top heavy chain for each:
        top_heavy = self.mixcr_df.loc[(self.mixcr_df["Chain"] == "IGH")
            & (self.mixcr_df["cloneRank"] == 0)]
        # Add prefix to column names:
        top_heavy = top_heavy.add_prefix("topHeavy_")
        # Rename pandas column
        top_heavy = top_heavy.rename({"topHeavy_sample": "sample"}, axis=1)
        # Combine:
        self.combined_df = pd.merge(self.combined_df, top_heavy, how="outer", on="sample")

    def initial_processing(self):
        """
        Initial table processing
        """
        # Get first column names:
        first_col_names = [col for col in self.combined_df if col.startswith("first")]
        # print(f"Column names:  {first_col_names}")
        for col in reversed(first_col_names):
            new_col = col.replace("first", "main")
            self.combined_df.insert(loc=3, column=new_col, value=self.combined_df[col])
        # Move main output sequence
        main_col_names = [col for col in self.combined_df if col.startswith("main")]
        for col in main_col_names:
            print(f"{col} datatype: {self.combined_df[col]}")
        cols_to_add_numeric = ["cloneCount", ]
        cols_to_append = [""]
        for i in range(len(self.combined_df)):
            category = self.combined_df["CATEGORY"].loc[i]
            if category == "Category2_Top2Together":
                pass
        pass

    def process_cat(self, category=1):
        """
        Processes dataframe for specified category.
        :param category: category to process (int).
        :return: processed pandas dataframe.
        """
        # Drop these columns in present, regardless:
        drop_cols = []
        # Columns to rename:
        rename_cols = {}
        # Get category dataframe:
        temp_df = self.get_cat(cat_list=self.cat_lists[category])
        if category == 1:
            temp_df = self._process_cat1(temp_df)
        elif category == 2:
            pass
        elif category == 3:
            pass
        elif category == 4:
            pass
        else:
            pass
        # Get top gene
        temp_df = self.get_top_vj(temp_df)
        return temp_df

    def get_cat(self, cat_list):
        """
        Gets rows associated with a particular category.
        :param cat_list: list of sample names.
        :return: pandas dataframe.
        """
        # Get just the applicable rows:
        temp_table = self.combined_df.loc[self.combined_df["sample"].isin(cat_list)]
        return temp_table

    @staticmethod
    def drop_cols_with_prefix(df: pandas.DataFrame, prefix: str):
        """
        Drops columns that start with a prefix from dataframe.
        :param df: pandas dataframe.
        :param prefix: prefix of column names to drop
        :return: dataframe with applicable columns dropped.
        """
        return df.loc[:, ~df.columns.str.startswith(prefix)]

    def get_top_vj(self, df: pandas.DataFrame, vcol: str="first_allVHitsWithScore", jcol: str="first_allJHitsWithScore"):
        """
        Gets top V and J gene for each row. Also adds column with a shortened/reformatted string for the top gene.
        :param df: dataframe to get top gene usage for.
        :param vcol: name of column where V gene information is.
        :param jcol: name of column where J gene information is.
        :return: pandas dataframe with new columns added.
        Note: For shortened gene names, "D" is dropped. (i.e. IGKV3D-15 and IGKV3-15 would both end up being IGKV3-15).
        """
        # Create copy of dataframe:
        df1 = df.copy()
        # Get top V and J gene:
        topv = df1[vcol].str.split("(", n=1, expand=True)
        topj = df1[jcol].str.split("(", n=1, expand=True)
        df1["top_Vgene"] = topv.loc[:, 0].values
        df1["top_Jgene"] = topj.loc[:, 0].values
        # Create a shorthand for top gene:
        # Remove *00:
        topv_short = df1[vcol].str.split("*", n=1, expand=True)
        topj_short = df1[jcol].str.split("*", n=1, expand=True)
        df1["top_Vgene_short"] = topv_short.loc[:, 0].values
        df1["top_Jgene_short"] = topj_short.loc[:, 0].values
        # Remove IGK/IGL
        df1["top_Vgene_short"] = df1["top_Vgene_short"].str[3:]
        df1["top_Jgene_short"] = df1["top_Jgene_short"].str[3:]
        # Remove "D" in Vgene:
        df1["top_Vgene_short"] = df1["top_Vgene_short"].str.replace("D", "")
        # Add 0 if numerical single digit:
        for row in range(df.shape[0]):
            # Split string:
            entry = df1.iloc[row, df1.columns.get_loc("top_Vgene_short")]
            entry = entry.split(sep="-", maxsplit=1)
            # Pad string if needed:
            if len(entry[1]) < 2:
                new_entry = entry[0] + "-0" + entry[1]
            else:
                new_entry = entry[0] + "-" + entry[1]
            df1.iloc[row, df1.columns.get_loc("top_Vgene_short")] = new_entry
        return df1

    def _process_cat1(self, df):
        """
        NOTE: This is the older wath of processing, can be ignored.
        Private method that applies formatting applicable for category 1.
        :param df: pandas dataframe to reformat.
        :return: pandas dataframe.
        """
        drop_cols = []
        rename_cols = {}
        # Drop blast columns:
        df = self.drop_cols_with_prefix(df, prefix="blast")
        return df

    def _process_cat2(self):
        # TODO
        drop_cols = []
        rename_cols = {}
        pass

    def _process_cat3(self):
        # TODO
        drop_cols = []
        rename_cols = {}
        pass

    def confirm_results(self):
        # TODO
        pass

    def save_final_table(self, category, col_domchain="domClone_Chain"):
        """
        Saves final tables as tab delimited text files. A separate table will be saved for IGK and for IGL.
        :param category: category to save.
        :param col_domchain: name of column that contains whether IGK or IGL is the dominant chain.
        """
        # Save to filenames:
        save_to = OUT_DIR + "finaltable_cat" + str(category) + "_" + self.cat_descriptions[category]
        save_to_igk = save_to + "_igk.txt"
        save_to_igl = save_to + "_igl.txt"
        # Get IGK and IGL dataframes:
        igk_df = self.cat_dfs[category].loc[self.cat_dfs[category][col_domchain] == "IGK"]
        igl_df = self.cat_dfs[category].loc[self.cat_dfs[category][col_domchain] == "IGL"]
        self.save_df(df=igk_df, save_to=save_to_igk, sep="\t")
        self.save_df(df=igl_df, save_to=save_to_igl, sep="\t")

    @staticmethod
    def read_in_list(filename, sep="\n"):
        """
        Reads in text file as a list, with each newline as a separate entry.
        :param filename: file path/name to read in (string).
        :param sep: character seperating entries in the list. (Default newline "\n").
        :return: list of strings.
        """
        with open(filename, "r") as my_file:
            content = my_file.read()
        the_list = content.split(sep=sep)
        return the_list

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

    @staticmethod
    def save_df(df: pandas.DataFrame, save_to: str, sep="\t"):
        """
        Saves dataframe.
        :param df: pandas dataframe to save.
        :param save_to: file path/name to save dataframe to.
        :param sep: character that should be used to delimit entries in row, default "\t".
        """
        df.to_csv(save_to, sep=sep, index=False)

    def create_fasta(self, out_file, stitched):
        """
        Creates fasta files based on read in data, separated by sample and chain.
        """
        # Empty context of text file:
        open(out_file, "w").close()
        # Open text file for creation (with open automatically closes file later):
        with open(out_file, "a") as myfile:
            # Iterate across rows:
            for i in range(len(self.combined_df)):
                category = self.combined_df["CATEGORY"].loc[i]
                if category not in self.remove_categories:
                    sample = self.combined_df["sample"].loc[i]
                    chain = self.combined_df["first_Chain"].loc[i]
                    sequence = self.combined_df["stitched"].loc[i]
                    if not stitched:
                        if category == "Category2_Top2Together":
                            seq1 = str(self.combined_df["piece1_seq"].loc[i])
                            seq2 = str(self.combined_df["piece2_seq"].loc[i])
                            seq3 = str(self.combined_df["piece3_seq"].loc[i])
                            if seq1 != "nan":
                                sequence = max([seq1, seq2, seq3], key=len)
                        else:
                            sequence = self.combined_df["first_subseq1"].loc[i]
                    self.combined_df["output_seq"] = sequence
                    # Build fasta string:
                    fasta_string = ">" + sample + "__" + chain + "__" + category + "\n" + sequence + "\n"
                    # Add fasta string to file:
                    myfile.write(fasta_string)
        # Check how many output_seq match first


if __name__ == '__main__':
    readin = ProduceFinalOutputs(combined_dff=COMBINED_DF_FILE, mixcr_file=MIXCR_OUTPUT_FILE,
                                 stitched_file=STITCHED_FILE, outdir=OUT_DIR)
