"""
Creates CSV from job accounting.


# To get accounting of jobs of last 2 days, and save on SGE cluster:
qacct -o <USERNAME> -j -d 2 > accounting.txt


NOTE: this expects file to be complete. If lines are missing, the results will be misaligned.


Example long to wide:
https://stackoverflow.com/questions/57435047/convert-long-format-text-file-from-s3-to-wide-format-csv
"""

# Import packages:
import pandas as pd
import os.path
import sys
import time

# Get command line arguments:
username = sys.argv[1]
account_file = sys.argv[2]
OUT_DIR = sys.argv[3]

# Info dict:
info_d = {}
out_file = OUT_DIR + "job_account.csv"

# Read lines from file:
with open(account_file, "r") as my_file:
    for line in my_file.readlines():
        # Split by first whitespace only:
        line_split = line.strip().split(maxsplit=1)
        line_split = [entry.strip() for entry in line_split]
        # print(line_split)
        # Load dictionary:
        if (len(line_split) > 1) & (line_split[0] != "OWNER") & \
                (line_split[0] != username) & (isinstance(line_split, list)):
            if line_split[0] not in info_d.keys():
                info_d[line_split[0]] = [line_split[1]]
            else:
                info_d[line_split[0]].append(line_split[1])

# # Iterate over dictionary:
# for key, value in info_d.items():
#     print(f"{key}: {value}")
#     print(f"Length: {len(value)}")

# Convert to dataframe:
df = pd.DataFrame(info_d)

# Checkout dataframe:
print(f"Dataframe shape from accounting.txt: {df.shape}")
# print(df.head(n=10))

# Read in previous accounting, if it exists:
if os.path.exists(out_file):
    prev_account = pd.read_csv(out_file)
    print(f"Previous account dimension: {prev_account.shape}")
    print(f"New account dimension: {df.shape}")    
    try:
        df = pd.concat([prev_account, df])
    except:
        print("Can't concatenate current timings with previous, saving as a temporary file")
        temp_file = OUT_DIR + time.strftime("%Y%m%d_%H%M") + "job_account.csv"
        print(f"Temp filename: {temp_file}")
        # Drop any columns that started with unnamed:
        df = df[df.columns.drop(list(df.filter(regex="Unnamed")))]
        df.to_csv(temp_file, index=False)
    else:
        # Drop any columns that started with unnamed:
        df = df[df.columns.drop(list(df.filter(regex="Unnamed")))]
        # Drop duplicate rows:
        df = df.drop_duplicates()
        # Save dataframe:
        df.to_csv(out_file, index=False)
else:
    # Save dataframe:
    df.to_csv(out_file, index=False)
