#!/home/sch-win1/miniconda3/envs/csun/bin/python
import re
import numpy as np
import sys
import subprocess
import pandas as pd
import pathlib
import os
import shutil
import inspect

def check_for_only_one(files, condition):
    check = [file for file in files if condition(file)]
    if len(check) == 1: 
        return check[0]
    else:
        func = inspect.getsource(condition)
        output = re.findall('lambda .+: (.+)', func)[0]
        print(f"Ambiguous -- exactly one file must meet the condition: {output}!")
        sys.exit(1)
    
def prep(dir, member_cutoff=3):
    if dir[-1] == "/":
        dir = dir[:-1]
    
    dir_files = os.listdir(dir)

    cloanalyst_dir = check_for_only_one(dir_files, lambda filename: filename[:2] == 'CR')
    heavy_file = check_for_only_one(dir_files, lambda filename: filename[-10:] == "_heavy.fas")
    light_file = check_for_only_one(dir_files, lambda filename: filename[-10:] == "_light.fas")

    print(f"Preparing {dir} folder...")

    cloanalyst_outputs = ['HeavyChainClonogram.nwk',
            'HeavyChainIntermediates.fasta',
            'HeavyChainIntermediates.npmf',
            'HeavyChainMarkedAlignment.fasta',
            'LightChainClonogram.nwk',
            'LightChainIntermediates.fasta',
            'LightChainIntermediates.npmf',
            'LightChainMarkedAlignment.fasta']

    [shutil.move(f"{dir}/{file}", f"{dir}/{cloanalyst_dir}/{file}") for file in dir_files if file in cloanalyst_outputs]
    
    with open(f"{dir}/{heavy_file}") as f: 
        heavy_data = f.readlines()
    with open(f"{dir}/{light_file}") as f: 
        light_data = f.readlines()

    df_assignments = pd.read_csv(f"{dir}/{cloanalyst_dir}/CloneAssignments.txt", sep="\s+")

    df = pd.read_csv(f"{dir}/{cloanalyst_dir}/Clones.txt", sep="\s+")
    subset = df[df["#Members"] >= member_cutoff]
    CloneIDs = subset["CloneID"].tolist()
    
    for clone_id in CloneIDs:
        try:
            os.mkdir(f"{dir}/clone_{clone_id}")
        except:
            pass

        subset = df_assignments[df_assignments.CloneID == clone_id]
        ReadIDs = subset.ReadID.tolist()
        CDR3s = subset.CDR3.tolist()

        with open(f"{dir}/clone_{clone_id}/members.txt", "w+") as f:
            f.write("Members\n")
            for read_id in ReadIDs:
                f.write(f"{read_id}\n")

        with open(f"{dir}/clone_{clone_id}/heavy_na.fna", "w+") as f:
            for read_id in ReadIDs:
                for id, seq in zip(heavy_data[::2], heavy_data[1::2]):

                    if read_id == id.rstrip()[1:]:
                        f.write(id)
                        f.write(seq)
                        break

        with open(f"{dir}/clone_{clone_id}/light_na.fna", "w+") as f:
            for read_id in ReadIDs:
                for id, seq in zip(light_data[::2], light_data[1::2]):
                    if read_id == id.rstrip()[1:]:
                        f.write(id)
                        f.write(seq)
                        break

        with open(f"{dir}/clone_{clone_id}/heavy_cdr3.fna", "w+") as f:
             for read_id, cdr3 in zip(ReadIDs, CDR3s):
                f.write(f">{read_id}\n")
                f.write(f"{cdr3}\n")
                    
    return CloneIDs

if __name__ == '__main__':
    prep(sys.argv[1])