#!/home/sch-win1/miniconda3/envs/csun/bin/python
import re
import numpy as np
import sys
import subprocess
import pandas as pd
import pathlib
import os



def igblast(dir, CloneIDs=None, member_cutoff=3):
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
        print(f"Blasting clone group {clone_id}...")

        try:
            os.mkdir(f"{dir}/clone_{clone_id}/auxiliary")
        except:
            pass

        # Clear files
        with open(f"{dir}/clone_{clone_id}/igblast.report", 'w+'): pass
        
        with open(f"{dir}/clone_{clone_id}/auxiliary/igblast_full.report", "w+"): pass
        
        with open(f"{dir}/clone_{clone_id}/igblast.report", 'a+') as f:
            for type in ['heavy', 'light']:
                f.write(f"{'-'*30} Group {clone_id} ({type}) {'-'*30}\n\n")

                seq_path = f"{os.getcwd()}/{dir}/clone_{clone_id}/{type}_na.fna"
                cmd = f"bin/igblastn -germline_db_V database/IGV -germline_db_J database/IGJ -germline_db_D database/IGD -organism human -query {seq_path} -show_translation"
                output = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=f"{pathlib.Path.home()}/igblast/")
                
                sequences = ["Query=" + _ for _ in output.stdout.split("Query=")][1:]
                df_genes = pd.DataFrame()
                ids = []

                for chunk in sequences:
                    genes = re.findall("Value([\S\s]+)Domain", chunk)[0]
                    genes = genes.split("\n")
                    genes = [_ for _ in genes if _ != '']
                    genes = [_.split() for _ in genes]
                    df = pd.DataFrame(genes)

                    id = re.findall("Query= ([a-zA-Z0-9-_]+)",chunk)[0]
                    ids.append(id)
                    df['id'] = id
                    df['id_len'] = len(id)
                    df_genes = pd.concat([df_genes,df])

                df_genes['VDJ'] = [re.findall("IG.(.)",_)[0] for _ in df_genes[0].tolist()]
                df_genes['HKL'] = [re.findall("IG(.).",_)[0] for _ in df_genes[0].tolist()]
                df_genes['gene'] = [re.findall("IG..(.+)",_)[0] for _ in df_genes[0].tolist()]
                cond = [True if float(_) > 1 else False for _ in df_genes[2].tolist()]
                df_genes['gene'] = [f"({gene})" if yn else gene for yn, gene in zip(cond, df_genes['gene'].tolist())]
                df_genes['gene_len'] = [len(_) for _ in df_genes['gene'].tolist()]
                df_genes['e_len'] = [len(_) for _ in df_genes[2].tolist()]

                gene_len = df_genes['gene_len'].max()
                id_len = df_genes['id_len'].max()
                e_len = df_genes['e_len'].max()

                for id in ids:
                    df_subset = df_genes[df_genes['id'] == id]
                    hkl = df_subset['HKL'].tolist()[0]

                    f.write(f"{id:<{id_len}}     {'Gene':<{gene_len*3+4}}     E-value\n")
                    for gene in ['V', 'D', 'J']:
                        df = df_subset[df_subset['VDJ'] == gene]
                        genes = df['gene'].tolist()
                        for i in range(3 - len(genes)):
                            genes.append("--")
                        g1, g2, g3 = genes

                        e = df[2].tolist()
                        for i in range(3 - len(e)):
                            e.append("--")
                        e1, e2, e3 = e

                        f.write(f"{'IG'+hkl+gene:<{id_len}}  |  {g1:<{gene_len}}  {g2:<{gene_len}}  {g3:<{gene_len}}  |  {e1:>{e_len}}  {e2:>{e_len}}  {e3:>{e_len}}\n")
                    f.write("\n")

                with open(f"{dir}/clone_{clone_id}/auxiliary/igblast_full.report", "a+") as a:
                    a.write(output.stdout)

            f.write(f"{'-'*30}-{'-'*5}-{'-'*len(str(clone_id))}--{'-'*len(type)}--{'-'*30}\n\n")

            
            
                

if __name__ == '__main__':
    igblast(sys.argv[1])