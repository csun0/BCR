#!/home/sch-win1/miniconda3/envs/csun/bin/python
import re
import numpy as np
import sys
import subprocess
import pandas as pd
import pathlib
import os
import threading

from igblast import *
from expasy import *
from prep import *

def cloanalyst_tracer(dir, member_cutoff=3):
    # input dir needs two .fas files ending in _heavy.fas and _light.fas

    CloneIDs = prep(dir)

    t1 = threading.Thread(target=expasy, args=[dir, CloneIDs])
    t2 = threading.Thread(target=igblast, args=[dir, CloneIDs])

    t1.start()
    t2.start()

    # Pass once all threads end
    t1.join()
    t2.join()

if __name__ == '__main__':
    cloanalyst_tracer(sys.argv[1])