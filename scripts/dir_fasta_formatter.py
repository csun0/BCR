#!/home/csun/miniconda3/bin/python3
import sys
import os
from fasta_formatter import *


inputdir = sys.argv[1]

for path, subdirs, files in os.walk(inputdir):
    for name in files:
        if name[-4:] == ".fas": 
            fullpath = os.path.join(path, name)
            print(fullpath)
            fasta_formatter(fullpath)