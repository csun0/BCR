#!/home/csun/miniconda3/envs/csun/bin/python
import pandas as pd
import sys


def txt_formatter(filename):
    with open(filename) as f:
        file = f.readlines()
    if file == []:
        return
        
    data = [_.split() for _ in file]
    df = pd.DataFrame(data, columns=data[0])
    col_len = [df[_].str.len().max() for _ in df.columns]

    temp = []
    for value, length in zip(data[0], col_len):
        temp.append(f"{value:<{length+3}}")
    temp.append("\n")

    output = []
    output.append("".join(temp))
    for line in data[1:]:
        temp = []
        for value, length in zip(line, col_len):
            temp.append(f"{value:<{length+3}}")
        temp.append("\n")
        output.append("".join(temp))

    with open(f"{filename}_formatted.txt", "w+") as f:
        f.writelines(output)

if __name__ == '__main__':
    txt_formatter(sys.argv[1])