#!/usr/bin/python3
"""
Looks up Run ID passed in script call and returns new sample name.
"""


# Import modules:
import sys
import pandas as pd

# Get command line argument:
old_name = sys.argv[1]

# Read in sample key:
key_file = "all_newids.csv"
key = pd.read_csv(key_file, low_memory=False)

# Look up old_name:
new_name = key["newname"][key["Run"] == old_name].values[0]
# print(f"Changing name: {old_name} ==>  {new_name}")

# Send to command line:
sys.stdout.write(new_name)
sys.exit()
