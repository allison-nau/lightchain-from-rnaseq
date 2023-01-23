#!/usr/local/bin/python3
"""
Creates new sample names to use. Will NOT transform the actual fasta files.
"""

import math
import random
import pandas as pd


def create_newnames(start=100000, end=200000, prefix="MMRF"):
    """
    Creates list of new names to draw from when renaming dataset.
    :param start: starting number for new filenames.
    :param end: ending number for new filenames.
    :param prefix: prefix to add to new filenames.
    :return: pandas dataframe of the new sample ID list.
    """
    # Create list of numbers:
    number_list = list(range(100000, 200000, 1))
    # Shuffle list:
    # TODO set random seed
    random.shuffle(number_list)
    # Convert to strings and add prefix:
    number_list = ["MMRF" + str(x) for x in number_list]
    # Convert number list to pandas df:
    df = pd.DataFrame(number_list, columns=["newname"])
    return df


def read_in_meta_data(filepath: str):
    """
    Reads in meta data.
    :param filepath: string of the filepath to the meta data csv.
    :return: pandas dataframe of meta data.
    """
    meta_data = pd.read_csv(filepath)
    return meta_data


def main():
    """
    Reads in meta data, generates new sample ids, and saves the mapping and a list of timepoint 1 samples.
    """
    # Read in meta data:
    meta_data_file = "SampleIDs//SraRunTable_modified_all.csv"
    meta_data = read_in_meta_data(meta_data_file)
    print(meta_data.head(10))
    # Create sample ID list:
    new_names = create_newnames()
    print(new_names.head())
    # Merge two dataframes by row index:
    combined = pd.merge(new_names, meta_data, how="outer", left_index=True, right_index=True)
    print(combined.head())
    # Filter on timepoint:
    tp1_df = combined[combined.s_Timepoint == 1]
    # Save to new file:
    new_names_file = "SampleIDs//all_newids.csv"
    tp1_file = "SampleIDs//timepoint1_newids.csv"
    combined.to_csv(new_names_file)
    tp1_df.to_csv(tp1_file)
    # Save Original sample IDs in order to be run as text file:
    sampleid_file = "SampleIDs//runs1.txt"
    tp1_df["Run"].to_csv(sampleid_file, header=False, index=False, line_terminator="\n")


# Code to run in main (test code) -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

