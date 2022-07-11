#!/home/sch-win1/miniconda3/envs/csun/bin/python
import sys
import os
from txt_formatter import *


inputdir = sys.argv[1]

for path, subdirs, files in os.walk(inputdir):
    for name in files:
        if (name[-4:] == ".txt"): 
            fullpath = os.path.join(path, name)
            print(fullpath)
            if 'RecombinationSummaries.txt' in name:
                with open(fullpath, 'r+') as f:
                    temp = f.read()
                    temp = temp.replace('No D', 'No_D')
                    temp = temp.replace(' (H)', '[H]')
                    temp = temp.replace(' (L)', '[L]]')
                with open(fullpath, 'w+') as f:
                    f.write(temp)
            txt_formatter(fullpath)
