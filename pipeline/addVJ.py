#!/usr/bin/env python3
"""
Add VJ to MiXCR results.
"""

# Import package
import pandas as pd
import os
import glob
import numpy as np
import sys

# Get command line arguments:
# (combined mixcr results file) (file to output results to)
IN_FILE = sys.argv[1]
OUT_FILE = sys.argv[2]


def get_top_vj(df: pd.DataFrame, vcol: str = "allVHitsWithScore", jcol: str = "allJHitsWithScore",
               ccol: str="allCHitsWithScore"):
    """
    Gets top V, J, and C gene for each row. Also adds column with a shortened/reformatted string for the top gene.
    :param df: dataframe to get top gene usage for.
    :param vcol: name of column where V gene information is.
    :param jcol: name of column where J gene information is.
    :param ccol: name of column where C gene information is (i.e. Constant).
    :return: pandas dataframe with new columns added.
    Note: For shortened gene names, "D" is dropped. (i.e. IGKV3D-15 and IGKV3-15 would both end up being IGKV3-15).
    """
    # Create copy of dataframe:
    df1 = df.copy()
    # Get top V and J gene:
    topv = df1[vcol].str.split("(", n=1, expand=True)
    topj = df1[jcol].str.split("(", n=1, expand=True)
    topc = df1[ccol].str.split("(", n=1, expand=True)
    df1["top_Vgene"] = topv.loc[:, 0].values
    df1["top_Jgene"] = topj.loc[:, 0].values
    df1["top_Cgene"] = topc.loc[:, 0].values
    # Create a shorthand for top gene:
    # Remove *00:
    topv_short = df1[vcol].str.split("*", n=1, expand=True)
    topj_short = df1[jcol].str.split("*", n=1, expand=True)
    topc_short = df1[ccol].str.split("*", n=1, expand=True)
    df1["top_Vgene_short"] = topv_short.loc[:, 0].values
    df1["top_Jgene_short"] = topj_short.loc[:, 0].values
    df1["top_Cgene_short"] = topc_short.loc[:, 0].values
    # Remove IGK/IGL (don't do this for C segment)
    df1["top_Vgene_short"] = df1["top_Vgene_short"].str[3:]
    df1["top_Jgene_short"] = df1["top_Jgene_short"].str[3:]
    # Remove "D" in Vgene:
    df1["top_Vgene_short"] = df1["top_Vgene_short"].str.replace("D", "")
    # Add 0 if numerical single digit:
    for row in range(df.shape[0]):
        # Split string:
        entry = df1.iloc[row, df1.columns.get_loc("top_Vgene_short")]
        entry = entry.split(sep="-", maxsplit=1)
        if len(entry) > 1:
            # Pad string if needed:
            if len(entry[1]) < 2:
                new_entry = entry[0] + "-0" + entry[1]
            else:
                new_entry = entry[0] + "-" + entry[1]
            df1.iloc[row, df1.columns.get_loc("top_Vgene_short")] = new_entry
    return df1


def read_in_csv(filename, sep=","):
    """
    Reads in csv "filename" as a pandas dataframe.
    :param filename: string containing file name of the csv file to be imported.
    :param sep: character that separate fields within the csv, default ","
    :return: pandas dataframe containing data within filename.
    """
    temp_table = pd.read_csv(filename, sep=sep, header=0, low_memory=False)
    return temp_table


if __name__ == '__main__':
    # Read in combined MiXCR results:
    mixcr_results = read_in_csv(filename=IN_FILE, sep="\t")
    # Add VJ column:
    vj_added = get_top_vj(df=mixcr_results)
    # Save to CSV:
    vj_added.to_csv(OUT_FILE, index=False, sep="\t")
