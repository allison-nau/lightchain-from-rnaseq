"""
Creates CSV from standard output from MiXCR pipeline. 
(Used to grab IGH, IGK, IGL counts as well as other statistics that MiXCR outputs to standard out.)
"""

# Import packages:
import pandas as pd
import os.path
import sys
import re
import numpy as np

# Get command line arguments:
# (standard output text file made using grep) (text file aggregated fastqc results produced by multiqc) (directory to output files to)
stdout_file = sys.argv[1]
multiqc_file = sys.argv[2]
OUT_DIR = sys.argv[3]

# File to save isolated multqc output to:
save_to_multiqc = OUT_DIR + "multiqc_fastqc.csv"
# File to save standard output to:
save_to_l = OUT_DIR + "processed_stdout_long.csv"
# File to save wide version of dataframe standard output + MultiQC to:
save_to_w = OUT_DIR + "processed_stdout_qc_wide.csv"
# File to save wide version with no protected names:
save_to_w2 = OUT_DIR + "processed_stdout_qc_wide_NOTrestricted.csv"

# Other global variables:
# Expected string associated with old variable names:
old_name = "Download sample"
# List of qsub scripts to save the information from
# (Note: this script is intended for scripts based on "pipeline.qsub")
# (Note 2: this will break if you used a named output file)
scripts = ["pipeline.qsub"]


# -----------------------------------------------------------------------------------------------------
# Functions for handling string processing:
def matcheslimitreached(line: str):
    """
    Processes standard output lines for the following MiXCR warning:
    "WARNING: too many partial alignments detected, consider skipping assemblePartial
    (enriched library?). /maxRightMatchesLimitReached/"
    :param line: string of line to be processed.
    :return: list containing values from processed line.
    """
    # Process line:
    processed1 = line.split(":", 1)
    processed2 = processed1[1].split("/", 2)
    returned_list = [processed1[0]] + processed2[0:2]
    return returned_list


def grab_name(line: str):
    """
    Extracts old or new name from line.
    :param line: string representing a line from standard output.
    :return: string representing.
    """
    # Process line:
    processed = line.split(":", 2)
    return processed


# -----------------------------------------------------------------------------------------------------
# Additional functions:
def read_in_multiqc(filename: str):
    """
    Reads in tsv of FastQC results summarized by FastQC, and returns a dataframe.
    :param filename: filename of MultiQC results.
    :return: pandas dataframe.
    """
    # Columns to check and combine status across both forward and reverse reads:
    columns_to_check = ['basic_statistics', 'per_base_sequence_quality', 'per_sequence_quality_scores',
                        'per_base_sequence_content', 'per_sequence_gc_content', 'per_base_n_content',
                        'sequence_length_distribution', 'sequence_duplication_levels', 'overrepresented_sequences',
                        'adapter_content']
    # Prefix to add for forward and reverse:
    prefix1 = "fastqc_1_"
    prefix2 = "fastqc_2_"
    # Read in dataframe:
    df_org = pd.read_csv(filename, sep="\t")
    # Add column with sample names:
    names = df_org["Sample"].str.split("_", expand=True)
    df_org["sample_main"] = names.iloc[:, 0]
    # Add column with what direction the sample is in:
    df_org["which_in_pair"] = names.iloc[:, 2]
    # Add column with minimum and maximum length:
    lengths = df_org["Sequence length"].str.split("-", expand=True)
    df_org["length_shortest"] = lengths[0]
    df_org["length_longest"] = lengths[1]
    # Separate 1 and 2 into two dataframes:
    df1 = df_org[df_org["which_in_pair"] == "1"]
    df2 = df_org[df_org["which_in_pair"] == "2"]
    # Add fastqc_1 or fastqc_2 prefix to all columns:
    df1 = df1.add_prefix(prefix1)
    df2 = df2.add_prefix(prefix2)
    # Remove prefix from "sample_main" columns:
    df1 = df1.rename(columns={"fastqc_1_sample_main": "accession"})
    df2 = df2.rename(columns={"fastqc_2_sample_main": "accession"})
    # Make wider by putting both 1 and 2 in the same row
    df_combine = pd.merge(df1, df2, on="accession")
    # Move accession to front:
    first_col = df_combine.pop("accession")
    df_combine.insert(0, "accession", first_col)
    # Add set of columns with True False for if it failed in either direction:]
    # For each MultiQC column to check:
    for col in columns_to_check:
        # Forward and reverse columns:
        col1 = prefix1 + col
        col2 = prefix2 + col
        # New column name:
        newname = "fastqc_" + col
        # Initiate column:
        df_combine[newname] = ""
        # Check conditions:
        df_combine[newname] = np.select(condlist=[(df_combine[col1] == "pass") & (df_combine[col2] == "pass"),
                                                  (df_combine[col1] == "fail") | (df_combine[col2] == "fail"),
                                                  (df_combine[col1] == "warn") | (df_combine[col2] == "warn")],
                                        choicelist=["pass", "fail", "warn"])
    return df_combine


# -----------------------------------------------------------------------------------------------------
# Read in multiqc results:
df_multi = read_in_multiqc(filename=multiqc_file)
df_multi.to_csv(save_to_multiqc)


# -----------------------------------------------------------------------------------------------------
# Initialize list of lists to store standard output results:
data_list = [["job_task", "variable", "value"]]

# Dictionary of rules (rules should be in order of priority, only one will be applied to a given line):
rules = {
    ".+max.+MatchesLimitReached.+": matcheslimitreached,
    ".+Download sample.+": grab_name,
    ".+Before MiXCR, file name will be changed to.+": grab_name
}


# Read in file and split up:
with open(stdout_file, "r") as myfile:
    # Read in lines
    lines = myfile.readlines()
    # Iterate through and handle:
    for the_line in lines:
        # Strip leading and trailing white space:
        the_line = the_line.strip()
        # Rule to follow:
        follow_rule = None
        # Check which rule should apply to line:
        for rule, my_fxn in rules.items():
            if bool(re.match(re.compile(rule), the_line)):
                follow_rule = my_fxn
                break
        # Process line according rule:
        if follow_rule is not None:
            temp_list = follow_rule(the_line)
            # Strip white space from each item in list:
            for i in range(len(temp_list)):
                temp_list[i] = temp_list[i].strip()
            # Check if this is related to the appropriate script:
            script_used = temp_list
            if temp_list[0].split(sep=".o")[0] in scripts:
                # Append results:
                data_list.append(temp_list)


# Convert to dataframe:
df = pd.DataFrame(data_list)
# Remove duplicate rows:
df.drop_duplicates(keep="first", inplace=True)
# Write file to csv:
df.to_csv(save_to_l, index=False, header=False)

# Make first row of dataframe column names:
df = pd.DataFrame(df.values[1:], columns=df.iloc[0])
print(df.head(10))
# Create a wide version of dataframe:
df_wide = df.pivot_table(index=["job_task"], columns="variable", values="value", aggfunc="first")
df_wide.reset_index(inplace=True)
df_wide = df_wide.rename(columns={"index": "job_task"})

# Combine with multiqc
df_wide = pd.merge(df_wide, df_multi, left_on=old_name, right_on="accession", how="outer")

# Save to csv:
df_wide.to_csv(save_to_w, index=False, header=True)
# Create a wide version of dataframe without old name:
df_wide2 = df_wide.drop(labels=["accession", old_name, "fastqc_1_Sample", "fastqc_1_Filename",
                                "fastqc_2_Sample", "fastqc_2_Filename"], axis=1)
# Save to csv:
df_wide2.to_csv(save_to_w2, index=False, header=True)

# for line in data_list:
#     print(line)
