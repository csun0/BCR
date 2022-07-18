#!/home/sch-win1/miniconda3/envs/csun/bin/python
import re
import numpy as np
import sys
import subprocess
import pandas as pd
import pathlib
import os


def expasy(seq):
    cmd = f"curl -s -d 'dna_sequence={seq}&output_format=fasta' https://web.expasy.org/cgi-bin/translate/dna2aa.cgi"
    output = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    temp = output.stdout.split("\n")
    temp = [_.rstrip() for _ in temp]
    temp = [f"\n{_}\n" if ">" in _ else _ for _ in temp]
    temp[0] = temp[0].replace("\n", "", 1)
    temp[-1] += "\n"
    temp = "".join(temp)

    return temp


def cloanalyst_tracer(dir, member_cutoff=3):
    if dir[-1] != "/":
        dir += "/"
    with open(f"{dir}SimpleMarkedUAs.fasta") as f:
        ua = f.readlines()

    df = pd.read_csv(f"{dir}Clones.txt", sep="\s+")
    subset = df[df["#Members"] >= member_cutoff]
    CloneIDs = subset["CloneID"].tolist()

    for clone_id in CloneIDs:
        try:
            os.mkdir(f"{dir}clone_{clone_id}")
        except:
            pass

        df = pd.read_csv(f"{dir}CloneAssignments.txt", sep="\s+")
        subset = df[df.CloneID == clone_id]
        ReadIDs = subset.ReadID.tolist()
        CDR3s = subset.CDR3.tolist()
        CDR3Lengths = subset["CDR3Length"].tolist()

        nt = []
        aa = []

        for id, cdr3seq in zip(ReadIDs, CDR3s):
            print(f"Examining {id}...")
            for ids, sequences in zip(ua[::6], ua[1::6]):
                if id == ids.rstrip()[1:]:
                    seq = sequences
                    nt.append(">" + id + '\n')
                    nt.append(seq)
                    break

            fragment = expasy(cdr3seq)
            fragment = fragment.split("\n")[1::2]
            fragment = [x for x in fragment if "-" not in x]

            full = expasy(seq)
            full = full.split("\n")[1::2]

            found_match = False
            for x in fragment:
                for y in full:
                    match = re.findall(f".{x}W", y)
                    if match:
                        found_match = True
                        a, b = re.findall(f"(.+){x}(.+)", y)[0]
                        aa.append([id, a, x, b])
            if not found_match:
                aa.append([id, "No W-end CDR3 pattern found."])

        df = pd.DataFrame(aa)
        maxlen = df[1].map(len).max()
        aa = []
        for index, row in df.iterrows():
            aa.append(f">{row[0]}\n")
            if "CDR3" in row[1]:
                aa.append(f"{row[1]}\n")
            else:
                aa.append(f"{row[1]:>{maxlen}}{row[2]}{row[3]}\n")

        with open(f"{dir}clone_{clone_id}/heavy_nt.fas", 'w+') as f:
            f.writelines(nt)
        with open(f"{dir}clone_{clone_id}/heavy_aa.fas", 'w+') as f:
            f.writelines(aa)

if __name__ == '__main__':
    cloanalyst_tracer(sys.argv[1])