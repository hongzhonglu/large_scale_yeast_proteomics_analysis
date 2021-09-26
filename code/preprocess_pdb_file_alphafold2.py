import os
import pandas as pd
import shutil

# pdb dir
dir0 = "/Users/xluhon/Documents/UP000002311_559292_YEAST/"
all_file = os.listdir(dir0)

# make new dir
os.mkdir("/Users/xluhon/Documents/" + "alphafold_pdb")

# copy files
output = "/Users/xluhon/Documents/" + "alphafold_pdb"
for x in all_file:
    print(x)
    if ".pdb.gz" in x:
        try:
            shutil.copy(dir0 + str(x), output)
        except:
            pass

# uncompress the gz files

"""
cd /Users/xluhon/Documents/alphafold_pdb
gunzip -k *.gz
"""

# further remove file in .gz format
"""
rm *.pdb.gz
"""






