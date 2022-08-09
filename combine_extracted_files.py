import pandas as pd
import numpy as np
import os
import glob
import sys 

phred_converted_version = False
sort_by_position = False

save_file = "/oak/stanford/groups/zihuai/fredlu/MpraScreen/aaron_AD/general_extraction_scripts/bigger_window.csv"

if phred_converted_version:
    save_file = save_file.replace(".csv","_phred_converted.csv")

base_dir = "/oak/stanford/groups/zihuai/fredlu/MpraScreen/aaron_AD/general_extraction_scripts/"
save_folder = sys.argv[1]

df = pd.DataFrame(columns=["index"])

# Start with Roadmap
for f in glob.glob(save_folder + "**/*",recursive=True):

    if not (("roadmap" in f.split("/")[-1]) and (".csv" in f.split("/")[-1])):
        continue
    # Read csv
    a = pd.read_csv(f,sep=',')

    # Drop bad columns
    a = a.drop(columns=["#Chrom","Pos_start","Pos_end","ref","alt","chr","position","chromosome","pos"] + [x for x in df.columns if "Unnamed" in x],errors='ignore')
    df = pd.merge(df,a,on="index",how="outer")

# Regbase
regbase_df = pd.DataFrame()
for f in glob.glob(save_folder + "**/*",recursive=True):

    if not (".csv" in f.split("/")[-1]):
        continue
    if not ("regbase" in f.split("/")[-1]):
        continue
    # Read csv
    a = pd.read_csv(f,sep=',')
    # Drop bad columns
    a = a.drop(columns=["#Chrom","Pos_start","Pos_end","ref","alt","chr","position","chromosome","pos"] + [x for x in df.columns if "Unnamed" in x],errors='ignore')
    a = a.groupby("index").mean()
    regbase_df = pd.concat([regbase_df,a])
df = pd.merge(df,regbase_df,on="index",how="outer",suffixes=('','_remove'))

# Eigen
eigen_df = pd.DataFrame()
for f in glob.glob(save_folder + "**/*",recursive=True):

    if not (".csv" in f.split("/")[-1]):
        continue
    if not ("eigen" in f.split("/")[-1]):
        continue
    # Read csv
    a = pd.read_csv(f,sep=',')
    # Drop bad columns
    a = a.drop(columns=["#Chrom","Pos_start","Pos_end","ref","alt","chr","position","chromosome","pos"] + [x for x in df.columns if "Unnamed" in x],errors='ignore')
    a = a.groupby("index").mean()
    eigen_df = pd.concat([eigen_df,a])
df = pd.merge(df,eigen_df,on="index",how="outer",suffixes=('','_remove'))

# FAVOR
favor_df = pd.DataFrame()
for f in glob.glob(save_folder + "**/*",recursive=True):

    if not (".csv" in f.split("/")[-1]):
        continue
    if not ("favor" in f.split("/")[-1]):
        continue
    # Read csv
    a = pd.read_csv(f,sep=',')
    # Drop bad columns
    a = a.drop(columns=["#Chrom","Pos_start","Pos_end","ref","alt","chr","position","chromosome","pos"] + [x for x in df.columns if "Unnamed" in x],errors='ignore')
    a = a.drop_duplicates(subset=["index"],keep="first")
    favor_df = pd.concat([favor_df,a])
df = pd.merge(df,favor_df,on="index",how="outer",suffixes=('','_remove'))
# QTL

for f in glob.glob(save_folder + "**/*",recursive=True):

    if not (".csv" in f.split("/")[-1]):
        continue
    if not ("qtl" in f.split("/")[-1]):
        continue
    # Read csv
    a = pd.read_csv(f,sep=',')
    # Drop bad columns
    a = a.drop(columns=["#Chrom","Pos_start","Pos_end","ref","alt","chr","position","chromosome","pos"] + [x for x in df.columns if "Unnamed" in x],errors='ignore')

    xQTL_combos = ["Positive_" + str(x) for x in list(np.unique(a["xQTL_Tissue"]))]

    b = pd.DataFrame(columns=["index"] + xQTL_combos)

    for ix,row in a.iterrows():
        i, xQTL_Tissue, pValue = row["index"], row["xQTL_Tissue"], row["Pvalue"]
        b.loc[len(b)] = [i] + [True if str(xQTL_Tissue) == tis.split("Positive_")[1] else False for tis in xQTL_combos]
    b = b.drop_duplicates(subset=["index"],keep="first")

    df = pd.merge(df,b,on="index",how="outer",suffixes=('','_uniqueidentifier_y'))
    df= df.drop(columns = [x for x in df.columns if "uniqueidentifier_y" in x])

# Drop useless columns
df = df.drop(columns=["#Chrom","Pos_start","Pos_end","ref","alt","chr","position","chromosome","pos"] + [x for x in df.columns if "Unnamed" in x],errors='ignore')

if phred_converted_version:
    for c in [x for x in df.columns if (("phred" in x) or ("PHRED" in x))]:
        print(c)
        df[c] = df[c].apply(lambda x: 1 - 10 ** (-float(x)/10))

if sort_by_position:
	# Sort to make sure they are in order
	df["POSITIONFORSORTING"] = df["index"].apply(lambda x: int(x.split(":")[1]))
	df = df.sort_values(by = 'POSITIONFORSORTING')
	df = df.drop(columns=["POSITIONFORSORTING"])

# Save
df.to_csv(save_file,sep=',')

