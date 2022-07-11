#!/home/csun/miniconda3/envs/csun/bin/python
import re
import numpy
import os
import sys

def heavy_formatter(H, L, K):
    with open(H) as f:
        heavy = f.readlines()
    with open(L) as f:
        lambdas = f.readlines()
    with open(K) as f:
        kappas = f.readlines()

    
    data = heavy
    ids = data[::2]
    seqs = data[1::2]

    total_num = int(re.findall('-(\d+)', ids[-1])[0])
    id_store_h = numpy.full(total_num+1, "", object)
    seq_store_h = numpy.full(total_num+1, "", object)

    for id, seq in zip(ids, seqs):
        id_store_h[int(re.findall('-(\d+)', id)[0])] = id
        seq_store_h[int(re.findall('-(\d+)', id)[0])] = seq

    data = kappas
    ids = data[::2]
    hk = []
    for id in ids:
        index = int(re.findall('-(\d+)', id)[0])
        hk.append(id_store_h[index])
        hk.append(seq_store_h[index])

    data = lambdas
    ids = data[::2]
    hl = []
    for id in ids:
        index = int(re.findall('-(\d+)', id)[0])
        hl.append(id_store_h[index])
        hl.append(seq_store_h[index])

    root, ext = os.path.splitext(H)
    with open(root + 'K' + ext, 'w+') as f:
        f.writelines(hk)
    with open(root + 'L' + ext, 'w+') as f:
        f.writelines(hl)


if __name__ == '__main__':
    heavy_formatter(sys.argv[1], sys.argv[2], sys.argv[3])