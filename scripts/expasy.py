#!/home/sch-win1/miniconda3/envs/csun/bin/python
import re
import numpy as np
import sys
import subprocess
import pandas as pd
import pathlib
import os
import textwrap

from igblast import *

def expasy_query(seq):
    cmd = f"curl -s -d 'dna_sequence={seq}&output_format=fasta' https://web.expasy.org/cgi-bin/translate/dna2aa.cgi"
    output = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    temp = output.stdout.split("\n")
    temp = [_.rstrip() for _ in temp]
    temp = [f"\n{_}\n" if ">" in _ else _ for _ in temp]
    temp[0] = temp[0].replace("\n", "", 1)
    temp[-1] += "\n"
    temp = "".join(temp)

    return temp


def expasy(dir, CloneIDs=None, member_cutoff=3):
    if dir[-1] == "/":
        dir = dir[:-1]

    if not CloneIDs:
        for path, subdirs, files in os.walk(dir):
            for file in files:
                if file == "Clones.txt":
                    fullpath = os.path.join(path, file)
                    df = pd.read_csv(fullpath, sep="\s+")
                    subset = df[df["#Members"] >= member_cutoff]
                    CloneIDs = subset["CloneID"].tolist()

    for clone_id in CloneIDs:
        print(f"Examining clone group {clone_id}...")
        aa = []

        with open(f"{dir}/clone_{clone_id}/heavy_na.fna", "r+") as h:
            data = h.readlines()
            ids = data[::2]
            sequences = data[1::2]

        with open(f"{dir}/clone_{clone_id}/heavy_cdr3.fna", "r+") as c:
            data = c.readlines()
            CDR3s = data[1::2]

        for id, cdr3, seq in zip(ids, CDR3s, sequences):
            
            fragment = expasy_query(cdr3)
            fragment = fragment.split("\n")[1::2]
            fragment = [x for x in fragment if "-" not in x]

            full = expasy_query(seq)
            full = full.split("\n")[1::2]

            found_match = False
            for fragment_aa in fragment:
                for count, full_aa in enumerate(full):
                    match = re.findall(f".{fragment_aa}W", full_aa)
                    if match:
                        found_match = True
                        pre, post = re.findall(f"(.+){fragment_aa}(.+)", full_aa)[0]
                        aa.append([id, pre, fragment_aa, post, count, seq])
            if not found_match:
                aa.append([id, "No W-end CDR3 pattern found."])

        df = pd.DataFrame(aa)
        maxlen = df[1].map(len).max()

        aa = []

        with open(f"{dir}/clone_{clone_id}/heavy_aa.faa", 'w+') as f:
            for index, row in df.iterrows():
                f.write(f"{row[0]}")
                if not row[1] == "No W-end CDR3 pattern found.":
                    f.write(f"{row[1]:>{maxlen}}{row[2]}{row[3]}\n".replace(" ","-"))
                else:
                    f.write(f"{row[1]}\n")

        with open(f"{dir}/clone_{clone_id}/heavy_mixed.txt", 'w+') as f:
            for index, row in df.iterrows():
                f.write(f"{row[0]}")
                if row[4] in [3,4,5]:
                    na = row[5][row[4]-3:]
                    aa = "".join(row[1:4])
                    na_wrap = textwrap.wrap(na, 60)
                    aa_wrap = textwrap.wrap(aa, 20)
                    for a, b in zip(na_wrap, aa_wrap):
                        f.write(f"{a}\n")
                        f.write(f"{''.join([f' {i} ' for i in b])}\n")
                f.write(f"\n")
if __name__ == '__main__':
    expasy(sys.argv[1])