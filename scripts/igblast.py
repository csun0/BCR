#!/home/sch-win1/miniconda3/envs/csun/bin/python
import re
import numpy as np
import sys
import subprocess
import pandas as pd
import pathlib
import os



def igblast(dir, member_cutoff=3):
    if dir[-1] != "/":
        dir += "/"

    df = pd.read_csv(f"{dir}Clones.txt", sep="\s+")
    subset = df[df["#Members"] >= member_cutoff]
    CloneIDs = subset["CloneID"].tolist()

    for clone_id in CloneIDs:
        print(f"Searching clone group {clone_id}...")
        with open(f"{dir}clone_{clone_id}/igblast_output.txt", 'w+') as f:
            f.write(f"-------------------- Group {clone_id} --------------------\n")
            pwd = f"{os.getcwd()}/{dir}clone_{clone_id}/heavy_nt.fas"
            cmd = f"bin/igblastn -germline_db_V database/IGV -germline_db_J database/IGJ -germline_db_D database/IGD -organism human -query {pwd} -show_translation"
            output = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=f"{pathlib.Path.home()}/igblast/")
            sequences = ["Query=" + _ for _ in output.stdout.split("Query=")][1:]

            for chunk in sequences:
                id = re.findall("Query= ([a-zA-Z0-9-_]+)",chunk)[0]
                id_len = len(id)
                genes = re.findall("Value([\S\s]+)Domain", chunk)[0]
                genes = genes.split("\n")
                genes = [_ for _ in genes if _ != '']
                genes = [_.split() for _ in genes]
                df = pd.DataFrame(genes)
                
                v_names = []
                d_names = []
                j_names = []

                v_prob = []
                d_prob = []
                j_prob = []

                temp = df.to_numpy()
                for count, line in enumerate(temp):
                    if "V" in line[0]:
                        v_names.append(temp[count][0])
                        v_prob.append(temp[count][2])
                    elif "D" in line[0]:
                        d_names.append(temp[count][0])
                        d_prob.append(temp[count][2])
                    else:
                        j_names.append(temp[count][0])
                        j_prob.append(temp[count][2])

                first = re.findall("(.+)\*", v_names[0])
                if re.findall("(.+)\*", v_names[1]) == first: 
                    v_names[1] = v_names[1].replace(first[0], "")
                    if re.findall("(.+)\*", v_names[2]) == first:
                        v_names[2] = v_names[2].replace(first[0], "")

                v_names = re.sub('[IGHVDJ\[\]\'\",]', "", str(v_names))
                v_prob = re.sub('[IGHVDJ\[\]\'\",]', "", str(v_prob))

                first = re.findall("(.+)\*", d_names[0])
                if re.findall("(.+)\*", d_names[1]) == first: 
                    d_names[1] = d_names[1].replace(first[0], "")
                    if re.findall("(.+)\*", d_names[2]) == first:
                        d_names[2] = d_names[2].replace(first[0], "")

                d_names = re.sub('[IGHVDJ\[\]\'\",]', "", str(d_names))
                d_prob = re.sub('[IGHVDJ\[\]\'\",]', "", str(d_prob))

                first = re.findall("(.+)\*", j_names[0])
                if re.findall("(.+)\*", j_names[1]) == first: 
                    j_names[1] = j_names[1].replace(first[0], "")
                    if re.findall("(.+)\*", j_names[2]) == first:
                        j_names[2] = j_names[2].replace(first[0], "")

                j_names = re.sub('[IGHVDJ\[\]\'\",]', "", str(j_names))
                j_prob = re.sub('[IGHVDJ\[\]\'\",]', "", str(j_prob))

                

                def empty_blast(e, names):
                    present = int(re.findall("\d.+(\d)",str(e))[0])
                    a = names.split()
                    for i in range(3-present):
                        a.append("")
                    return a

                try: v1, v2, v3 = v_names.split()
                except ValueError as e: v1, v2, v3 = empty_blast(e, v_names)
                try: d1, d2, d3 = d_names.split()
                except ValueError as e: d1, d2, d3 = empty_blast(e, d_names)
                try: j1, j2, j3 = j_names.split()
                except ValueError as e: j1, j2, j3 = empty_blast(e, j_names)

                # v_len = [len(_) for _ in v_names.split()]
                # d_len = [len(_) for _ in d_names.split()]
                # j_len = [len(_) for _ in j_names.split()]
                # lens = np.array([v_len, d_len, j_len])
                lens = np.array([[v1,v2,v3], [d1,d2,d3], [j1,j2,j3]])
                nplen = np.vectorize(len)
                lens = np.amax(np.array(list(map(nplen, lens))),0)


                v_len = [len(_) for _ in v_prob.split()]
                d_len = [len(_) for _ in d_prob.split()]
                j_len = [len(_) for _ in j_prob.split()]
                lens_prob = np.array([v_len, d_len, j_len])
                lens_prob = np.amax(lens_prob, 0)

                try: v1_prob, v2_prob, v3_prob = v_prob.split()
                except ValueError as e: v1_prob, v2_prob, v3_prob = empty_blast(e, v_prob)
                try: d1_prob, d2_prob, d3_prob = d_prob.split()
                except ValueError as e: d1_prob, d2_prob, d3_prob = empty_blast(e, d_prob)
                try: j1_prob, j2_prob, j3_prob = j_prob.split()
                except ValueError as e: j1_prob, j2_prob, j3_prob = empty_blast(e, j_prob)

                f.write(f"{id}     {'Match':<{6+lens.sum()}}   E-value\n")
                f.write(f"{'IGHV':<{id_len}}  |  {v1:<{lens[0]}}  {v2:<{lens[1]}}  {v3:<{lens[2]}}  |  {v1_prob:<{lens_prob[0]}}  {v2_prob:<{lens_prob[1]}}  {v3_prob:<{lens_prob[2]}}\n")
                f.write(f"{'IGHD':<{id_len}}  |  {d1:<{lens[0]}}  {d2:<{lens[1]}}  {d3:<{lens[2]}}  |  {d1_prob:<{lens_prob[0]}}  {d2_prob:<{lens_prob[1]}}  {d3_prob:<{lens_prob[2]}}\n")
                f.write(f"{'IGHJ':<{id_len}}  |  {j1:<{lens[0]}}  {j2:<{lens[1]}}  {j3:<{lens[2]}}  |  {j1_prob:<{lens_prob[0]}}  {j2_prob:<{lens_prob[1]}}  {j3_prob:<{lens_prob[2]}}\n\n")




if __name__ == '__main__':
    igblast(sys.argv[1])