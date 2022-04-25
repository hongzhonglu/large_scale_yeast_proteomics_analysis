# ftp://ftp.ebi.ac.uk/pub/databases/alphafold/v2/UP000002311_559292_YEAST_v2.tar sce
# ftp://ftp.ebi.ac.uk/pub/databases/alphafold/v2/UP000002485_284812_SCHPO_v2.tar Schizosaccharomyces pombe
# ftp://ftp.ebi.ac.uk/pub/databases/alphafold/v2/UP000005640_9606_HUMAN_v2.tar human
# ftp://ftp.ebi.ac.uk/pub/databases/alphafold/v2/UP000000625_83333_ECOLI_v2.tar e.coli
# volume calculation
"java -jar /Users/xluhon/Documents/ProteinVolume_1.3/ProteinVolume_1.3.jar /Users/xluhon/Documents/alphafold_pdb"





import os
import pandas as pd
import shutil


# for the yeast
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




###############################################################
# for the e.coli
# pdb dir
dir0 = "/Users/xluhon/Documents/UP000000625_83333_ECOLI_v2/"
all_file = os.listdir(dir0)

# make new dir
os.mkdir("/Users/xluhon/Documents/" + "alphafold_pdb_ecoli")

# copy files
output = "/Users/xluhon/Documents/" + "alphafold_pdb_ecoli"
for x in all_file:
    print(x)
    if ".pdb.gz" in x:
        try:
            shutil.copy(dir0 + str(x), output)
        except:
            pass

# uncompress the gz files

"""
cd /Users/xluhon/Documents/alphafold_pdb_ecoli
gunzip -k *.gz
"""

# further remove file in .gz format
"""
rm *.pdb.gz
"""

# calculate the volume
"java -jar /Users/xluhon/Documents/ProteinVolume_1.3/ProteinVolume_1.3.jar /Users/xluhon/Documents/alphafold_pdb_ecoli"





###############################################################
# for the SCHPO
# pdb dir
dir0 = "/Users/xluhon/Documents/UP000002485_284812_SCHPO_v2/"
all_file = os.listdir(dir0)

# make new dir
os.mkdir("/Users/xluhon/Documents/" + "alphafold_pdb_SCHPO")

# copy files
output = "/Users/xluhon/Documents/" + "alphafold_pdb_SCHPO"
for x in all_file:
    print(x)
    if ".pdb.gz" in x:
        try:
            shutil.copy(dir0 + str(x), output)
        except:
            pass

# uncompress the gz files

"""
cd /Users/xluhon/Documents/alphafold_pdb_SCHPO
gunzip -k *.gz
"""

# further remove file in .gz format
"""
rm *.pdb.gz
"""

# calculate the volume
"java -jar /Users/xluhon/Documents/ProteinVolume_1.3/ProteinVolume_1.3.jar /Users/xluhon/Documents/alphafold_pdb_SCHPO"




# on the cluster
# for the ECOLI
"""
mkdir /lustre/home/acct-clslhz/clslhz/pdb_file/UP000000625_83333_ECOLI_v2/
tar -xvf /lustre/home/acct-clslhz/clslhz/pdb_file/UP000000625_83333_ECOLI_v2.tar -C /lustre/home/acct-clslhz/clslhz/pdb_file/UP000000625_83333_ECOLI_v2/
"""

# pdb dir
dir00 = "/lustre/home/acct-clslhz/clslhz/pdb_file/"
dir0 = dir00 + "UP000000625_83333_ECOLI_v2/"
all_file = os.listdir(dir0)
# make new dir
os.mkdir(dir00 + "alphafold_pdb_ecoli")

# copy files
output = dir00 + "alphafold_pdb_ecoli"
for x in all_file:
    print(x)
    if ".pdb.gz" in x:
        try:
            shutil.copy(dir0 + str(x), output)
        except:
            pass

# uncompress the gz files

"""
dir00="/lustre/home/acct-clslhz/clslhz/pdb_file/"
cd $dir00/alphafold_pdb_ecoli
gunzip *.gz
"""

# further remove file in .gz format
"""
rm *.pdb.gz
"""

# calculate the volume
"java -jar /Users/xluhon/Documents/ProteinVolume_1.3/ProteinVolume_1.3.jar /Users/xluhon/Documents/alphafold_pdb_ecoli"

