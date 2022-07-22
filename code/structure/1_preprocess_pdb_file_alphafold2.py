# ftp://ftp.ebi.ac.uk/pub/databases/alphafold/v2/UP000002311_559292_YEAST_v2.tar sce
# ftp://ftp.ebi.ac.uk/pub/databases/alphafold/v2/UP000000559_237561_CANAL_v2.tar Schizosaccharomyces pombe
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

# calculate the volume
"java -jar /Users/xluhon/Documents/ProteinVolume_1.3/ProteinVolume_1.3.jar /Users/xluhon/Documents/alphafold_pdb"





###############################################################
# for the SCHPO
# pdb dir
dir0 = "/Users/xluhon/Documents/UP000000559_237561_CANAL_v2/"
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







##########################################################
# on the cluster
##########################################################
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
os.mkdir(dir00 + "alphafold_pdb")

# copy files
output = dir00 + "alphafold_pdb"
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
cd $dir00/alphafold_pdb
gunzip *.gz
"""

# further remove file in .gz format
"""
rm *.pdb.gz
"""



# This part is running on cluster
# split the file into 40 sub-folders
import os
import pandas as pd
sub_fold_num = 1
os.system("rm -r /lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/sub_folder_*")
dir01="/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/sub_folder_"
for x in range(0, sub_fold_num):
    print(x)
    os.mkdir(dir01 + str(x))

# move the pdb file into each sub-folder
all_file = os.listdir('/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb')



# put the proprocess pdb file in the cycling...
pdb_ok = pd.read_excel("/lustre/home/acct-clslhz/clslhz/pdb_file/ecoli_structure_size_part1.xlsx",engine='openpyxl')
pdb_list = pdb_ok['Protein'].tolist()
pdb_list = [x+ ".pdb" for x in pdb_list]
all_file = list(set(all_file)-set(pdb_list)) # 2104 file remaining


import numpy
import shutil
l = numpy.array_split(numpy.array(all_file), sub_fold_num)
file_name = [dir01 + str(i) for i in range(0, sub_fold_num)]
dir02 = '/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb/'
for row, out in zip(l,file_name):
    print(row, out)
    for x in row:
        print(x)
        shutil.copy(dir02 + str(x), out)

# generate one sh file to process all datasets in parallel
def generate_sh():
    import os
    # first update the main function
    out_sh_file = '/lustre/home/acct-clslhz/clslhz/pdb_file/volume.sh'
    parallel = 1
    pdb_dir0 = "/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/"
    pdb_dir_all = os.listdir(pdb_dir0)
    # write in the start file
    start_part = "#!/bin/bash\n" \
                 "#SBATCH --job-name=test\n" \
                 "#SBATCH --partition=cpu\n" \
                 "#SBATCH -n 80\n" \
                 "#SBATCH --ntasks-per-node=40\n" \
                 "#SBATCH --output=%j.out\n" \
                 "#SBATCH --error=%j.err\n"\
                 "#SBATCH --mail-type=end\n"\
                 "#SBATCH --mail-user=hongzhonglu@sjtu.edu.cn\n"\
                 "module load lammps/2020-cpu\n" \
                 "ulimit -s unlimited\n" \
                 "ulimit -l unlimited\n"
    newfile = open(out_sh_file, "w")
    newfile.writelines(start_part)
    # template
    volume_sh = "java -jar /lustre/home/acct-clslhz/clslhz/ProteinVolume_1.3/ProteinVolume_1.3.jar /lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/sub_folder_1"
    for i, x in enumerate(pdb_dir_all):
        print(i, x)
        volume_sh_update = volume_sh.replace("sub_folder_1", x)
        if i % parallel == 0:
            newfile.writelines(volume_sh_update)
            newfile.write(" & " + "\n")
            newfile.write("\n")
        else:
            newfile.write(volume_sh_update)
    newfile.write("wait" + "\n")
    newfile.close()

generate_sh()


# collect the dataset
import glob
import shutil
ss = glob.glob("/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/*/*.txt")
for xx in ss:
    print(xx)
    shutil.copy(xx, "/lustre/home/acct-clslhz/clslhz/pdb_file/original_result/")

# combine the initial result
read_files = glob.glob("/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/*/*.txt")
with open("/lustre/home/acct-clslhz/clslhz/pdb_file/organism_result.txt", "wb") as outfile:
    for f in read_files:
        with open(f, "rb") as infile:
            outfile.write(infile.read())







# from here, all code run in mac
# preprocess txt file
ss = open("data/organism_result.txt",).readlines()
with open("data/organism_result_new.txt", "w") as outfile:
    for i,xx in enumerate(ss):
        print(i)
        if i <=7 and "Total Volume" in xx:
            print(xx)
            outfile.write(xx)
        if i >=7 and "-F1-model_v2 " in xx:
            outfile.write(xx)


# second part
# calculate the protein size based on its protein 3D structures
import pandas as pd
input_file = "data/organism_result_new.txt"
with open(input_file) as file_in:
    lines = []
    for line in file_in:
        lines.append(line)
lines0 = [x for x in lines if "Reading hydrogens is turned on" not in x]
lines1 = lines0[7:]
p1 = []
p2 = []
p3 = []
p4 = []
p5 = []
p6 = []
for x in lines1:
    print(x)
    ss = x.split(" ")
    ss0 = [x for x in ss if x is not '']
    ss0 = [x.replace("\n", "") for x in ss0]
    ss0 = [x.replace(",", ".") for x in ss0]
    s1 = ss0[0]
    s2 = float(ss0[1])
    s3 = float(ss0[2])
    s4 = float(ss0[3])
    s5 = float(ss0[4])
    s6 = float(ss0[5])
    p1.append(s1)
    p2.append(s2)
    p3.append(s3)
    p4.append(s4)
    p5.append(s5)
    p6.append(s6)
# note the original unit for the volume is Å
# 1Å = 0.1 nm; 1Å^3 = 0.001 nm^3
volume_df = pd.DataFrame({"Protein":p1,"Total_Volume":p2,"Void_Volume":p3,"VDW_Volume":p4,"Packing Density":p5,"Time_Taken_ms":p6})

# change the unit from Å to nm
volume_df0 = volume_df.copy()
volume_df0["Total_Volume"] = volume_df["Total_Volume"]/1000
volume_df0["Void_Volume"] = volume_df["Void_Volume"]/1000
volume_df0["VDW_Volume"] = volume_df["VDW_Volume"]/1000
volume_df0.to_excel("data/ecoli_structure_size_part4.xlsx")






##########################################################
# on the cluster
##########################################################
# UP000000559_237561_CANAL_v2.tar
"""
mkdir /lustre/home/acct-clslhz/clslhz/pdb_file/UP000000559_237561_CANAL_v2/
tar -xvf /lustre/home/acct-clslhz/clslhz/pdb_file/UP000000559_237561_CANAL_v2.tar -C /lustre/home/acct-clslhz/clslhz/pdb_file/UP000000559_237561_CANAL_v2/
"""

# pdb dir
import os
import shutil
dir00 = "/lustre/home/acct-clslhz/clslhz/pdb_file/"
dir0 = dir00 + "UP000000559_237561_CANAL_v2/"
all_file = os.listdir(dir0)
# make new dir
os.system("rm -r " + dir00 + "alphafold_pdb" )
os.mkdir(dir00 + "alphafold_pdb")

# copy files
output = dir00 + "alphafold_pdb"
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
cd $dir00/alphafold_pdb
gunzip *.gz
"""

# further remove file in .gz format
"""
rm *.pdb.gz
"""



# This part is running on cluster
# split the file into 40 sub-folders
import os
import numpy
import shutil
import pandas as pd
sub_fold_num = 10
os.mkdir('/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split')
os.system("rm -r /lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/sub_folder_*")
dir01="/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/sub_folder_"
for x in range(0, sub_fold_num):
    print(x)
    os.mkdir(dir01 + str(x))

# move the pdb file into each sub-folder
all_file = os.listdir('/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb')
# put the proprocess pdb file in the cycling...
pdb_ok = pd.read_excel("/lustre/home/acct-clslhz/clslhz/pdb_file/organism_structure_size_part1.xlsx",engine='openpyxl')
pdb_list = pdb_ok['Protein'].tolist()
pdb_list = [x+ ".pdb" for x in pdb_list]
all_file = list(set(all_file)-set(pdb_list))


l = numpy.array_split(numpy.array(all_file), sub_fold_num)
file_name = [dir01 + str(i) for i in range(0, sub_fold_num)]
dir02 = '/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb/'
for row, out in zip(l,file_name):
    print(row, out)
    for x in row:
        print(x)
        shutil.copy(dir02 + str(x), out)

# generate one sh file to process all datasets in parallel
def generate_sh():
    import os
    # first update the main function
    out_sh_file = '/lustre/home/acct-clslhz/clslhz/pdb_file/volume.sh'
    parallel = 1
    pdb_dir0 = "/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/"
    pdb_dir_all = os.listdir(pdb_dir0)
    # write in the start file
    start_part = "#!/bin/bash\n" \
                 "#SBATCH --job-name=test\n" \
                 "#SBATCH --partition=cpu\n" \
                 "#SBATCH -n 80\n" \
                 "#SBATCH --ntasks-per-node=40\n" \
                 "#SBATCH --output=%j.out\n" \
                 "#SBATCH --error=%j.err\n"\
                 "#SBATCH --mail-type=end\n"\
                 "#SBATCH --mail-user=hongzhonglu@sjtu.edu.cn\n"\
                 "module load lammps/2020-cpu\n" \
                 "ulimit -s unlimited\n" \
                 "ulimit -l unlimited\n"
    newfile = open(out_sh_file, "w")
    newfile.writelines(start_part)
    # template
    volume_sh = "java -jar /lustre/home/acct-clslhz/clslhz/ProteinVolume_1.3/ProteinVolume_1.3.jar /lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/sub_folder_1"
    for i, x in enumerate(pdb_dir_all):
        print(i, x)
        volume_sh_update = volume_sh.replace("sub_folder_1", x)
        if i % parallel == 0:
            newfile.writelines(volume_sh_update)
            newfile.write(" & " + "\n")
            newfile.write("\n")
        else:
            newfile.write(volume_sh_update)
    newfile.write("wait" + "\n")
    newfile.close()

generate_sh()


# collect the dataset
import glob
import shutil
ss = glob.glob("/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/*/*.txt")
for xx in ss:
    print(xx)
    shutil.copy(xx, "/lustre/home/acct-clslhz/clslhz/pdb_file/original_result/")

# combine the initial result
read_files = glob.glob("/lustre/home/acct-clslhz/clslhz/pdb_file/alphafold_pdb_split/*/*.txt")
with open("/lustre/home/acct-clslhz/clslhz/pdb_file/organism_result.txt", "wb") as outfile:
    for f in read_files:
        with open(f, "rb") as infile:
            outfile.write(infile.read())







# from here, all code run in mac
import pandas as pd
import os

# preprocess txt file
ss = open("data/organism_result.txt",).readlines()
with open("data/organism_result_new.txt", "w") as outfile:
    for i,xx in enumerate(ss):
        print(i)
        if i <=7 and "Total Volume" in xx:
            print(xx)
            outfile.write(xx)
        if i >=7 and "-F1-model_v2 " in xx:
            outfile.write(xx)


# second part
# calculate the protein size based on its protein 3D structures
import pandas as pd
input_file = "data/organism_result_new.txt"
with open(input_file) as file_in:
    lines = []
    for line in file_in:
        lines.append(line)
lines0 = [x for x in lines if "Reading hydrogens is turned on" not in x]
lines1 = lines0[1:]
p1 = []
p2 = []
p3 = []
p4 = []
p5 = []
p6 = []
for x in lines1:
    print(x)
    ss = x.split(" ")
    ss0 = [x for x in ss if x is not '']
    ss0 = [x.replace("\n", "") for x in ss0]
    ss0 = [x.replace(",", ".") for x in ss0]
    s1 = ss0[0]
    s2 = float(ss0[1])
    s3 = float(ss0[2])
    s4 = float(ss0[3])
    s5 = float(ss0[4])
    s6 = float(ss0[5])
    p1.append(s1)
    p2.append(s2)
    p3.append(s3)
    p4.append(s4)
    p5.append(s5)
    p6.append(s6)
# note the original unit for the volume is Å
# 1Å = 0.1 nm; 1Å^3 = 0.001 nm^3
volume_df = pd.DataFrame({"Protein":p1,"Total_Volume":p2,"Void_Volume":p3,"VDW_Volume":p4,"Packing Density":p5,"Time_Taken_ms":p6})

# change the unit from Å to nm
volume_df0 = volume_df.copy()
volume_df0["Total_Volume"] = volume_df["Total_Volume"]/1000
volume_df0["Void_Volume"] = volume_df["Void_Volume"]/1000
volume_df0["VDW_Volume"] = volume_df["VDW_Volume"]/1000
volume_df0.to_excel("data/organism_structure_size_part1.xlsx")



