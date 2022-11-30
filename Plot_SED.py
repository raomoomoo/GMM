import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import subprocess

centre  = []
monomer = []
sphere   = []
#Agg = {centre : 1, monomer: 0.15, sphere:1.23}
Aggs = [centre, monomer, sphere]
wavs = np.arange(550, 2500, 50).tolist()
PATH_TO_RUNS = os.path.expanduser('~/runs/Agg_models/AGG_3_wav')
for i in wavs:
    for j in Aggs:
        subprocess.run("cd ${PATH_TO_RUNS}")
        subprocess.run("cd ${i}") 

        with open("gmm01f.in") as f:
            lines = f.readlines()

            lines # ['This is the first line.\n', 'This is the second line.\n']

            lines[0] = str(i)+".k\n"

            lines # ["This is the line that's replaced.\n", 'This is the second line.\n']

        with open("gmm01f.in", "w") as f:
            f.writelines(lines)


        outfile = str(i)+".out"
        #check this $... 
        subprocess.run("ifort -o ${outfile} gmm01s.f90")

        with open("run.q") as f:
            lines = f.readlines()

            lines # ['This is the first line.\n', 'This is the second line.\n']

            lines[-1] = "./"+str(outfile)+"\n"

            lines # ["This is the line that's replaced.\n", 'This is the second line.\n']

        with open("run.q", "w") as f:
            f.writelines(lines)
        subprocess.run("sbatch run.q")

                                   
        #capture abs value and append to array maybe this test ^^^
        Agg = pd.read_csv('GMM/gmm01s_00260.out', header=None,skiprows=3, delimiter="\s", engine='python')
        Agg = Agg.apply(pd.to_numeric, errors='coerce')
        Cabs=Agg.at[1,4]
        
        j.append(Cabs)