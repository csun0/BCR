#!/home/csun/miniconda3/bin/python3
import os
import re
import sys

def fasta_formatter(filename):
    with open(filename, "r+") as f:
        temp = f.readlines()

    subject = re.findall("S\d",filename)[0]
    type = re.findall("(Ig.)",filename)[0]

    temp = [_.rstrip() for _ in temp]
    temp = [f"\n{_}\n" if ">" in _ else _ for _ in temp ]
    temp[0] = temp[0].replace("\n", "", 1)
    temp[-1] += "\n"

    with open(f"{subject}_{type[:-1]}{type[-1].upper()}.fas", "w+") as f:
        f.writelines(temp)

if __name__ == '__main__':
    fasta_formatter(sys.argv[1])