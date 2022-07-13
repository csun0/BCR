#!/home/sch-win1/miniconda3/envs/csun/bin/python
import re
import numpy as np
import sys
import subprocess
import pandas as pd

def expasy(seq):
    cmd = f"curl -s -d 'dna_sequence={seq}&output_format=fasta' https://web.expasy.org/cgi-bin/translate/dna2aa.cgi"
    output = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    temp = output.stdout.split("\n")
    temp = [_.rstrip() for _ in temp]
    temp = [f"\n{_}\n" if ">" in _ else _ for _ in temp ]
    temp[0] = temp[0].replace("\n", "", 1)
    temp[-1] += "\n"
    temp = "".join(temp)

    return temp

def cloanalyst_tracer(path,member_cutoff=3):
    if path[-1] != "/":
        path += "/"

    df = pd.read_csv(f"{path}Clones.txt", sep="\s+")
    subset = df[df["#Members"] >= member_cutoff]
    CloneIDs = subset["CloneID"].tolist()

    with open(f"{path}SimpleMarkedUAs.fasta") as f:
        ua = f.readlines()

    for clone_id in CloneIDs:
        df = pd.read_csv(f"{path}CloneAssignments.txt", sep="\s+")
        subset = df[df.CloneID==clone_id]
        ReadIDs = subset.ReadID.tolist()
        CDR3s = subset.CDR3.tolist()
        CDR3Lengths = subset["CDR3Length"].tolist()

        for id, cdr3seq in zip(ReadIDs, CDR3s):
            for ids, sequences in zip(ua[::6], ua[1::6]):
                if id == ids.rstrip()[1:]:
                    seq = sequences
                    break

            fragment = expasy(cdr3seq)
            fragment = fragment.split("\n")[1::2]
            fragment = [x for x in fragment if "-" not in x]

            full = expasy(seq)
            full = full.split("\n")[1::2]
            for x in fragment:
                for y in full:
                    a = re.findall(f"C{x}W",y)
                    if a:
                        print(f"{id}: group {clone_id}")
                        print(seq)
                        print(y)
                        print("\n")

if __name__ == '__main__':
    cloanalyst_tracer(sys.argv[1])