# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:37:40 2021

@author: chens
"""

import subprocess
import os
import time
import shutil
import configparser
import tqdm

configFile = "config.ini"
config = configparser.ConfigParser()
config.read(configFile)

# XTB=["wsl", "-d", "ubuntu", "-e", "/opt/xtb-6.4.0/bin/xtb"]
XTB=["xtb"]

iF = open(config['IODir']['Template'], "r")
freezeAtoms = list(map(str, range(1, int(iF.readline()) + 1)))
iF.close()
#freezeStrs = config['Properties']['freeze'].split(',')
#freezeAtoms = []
#for freezeStr in freezeStrs:
#    for atomNum in freezeStr.split('-'):
#        if atomNum not in freezeAtoms:
#            freezeAtoms.append(atomNum)
#freezeAtoms.sort()


os.environ['OMP_STACKSIZE'] = '1G'
os.chdir(config["IODir"]["Output"])
if os.path.exists("xtbtmp"):
    shutil.rmtree("xtbtmp")
os.makedirs("xtbtmp")
if os.path.exists("xtbopt"):
    print("Directory xtbopt found. Check if you still need them from the previous runs.")
    if(input("To delete results from previous runs (y/n)? ") != 'y'):
        exit()
    shutil.rmtree("xtbopt")
os.mkdir("xtbopt")
files = []
for file in os.listdir():
    if ".mol" in file:
        files.append(file)
os.chdir("xtbtmp")
oF = open("freeze.inp", "w")
oF.write("$fix\n  atoms: "+",".join(freezeAtoms)+"\n$end")
oF.close()
devNull = open(os.devnull, 'w')
totalItems = len(files)
start = time.time()
pbar = tqdm.tqdm(total = len(files))
for file in files:
    shutil.copyfile("../" + file, "temp.mol")
    subprocess.run(XTB + ["temp.mol", "--opt", "--gfnff", "--input", "freeze.inp"], stdout=devNull, stderr=devNull)
    os.remove("temp.mol")
    if "xtbopt.mol" not in os.listdir():
        pbar.update()
        break
    os.rename("xtbopt.mol", "temp.mol")
    subprocess.run(XTB + ["temp.mol", "--opt", "--input", "freeze.inp"], stdout=devNull, stderr=devNull)
    os.remove("temp.mol")
    if "xtbopt.mol" not in os.listdir():
        pbar.update()
        break
    shutil.move("xtbopt.mol", "../xtbopt/" + file)
    pbar.update()
pbar.close()
print("Finished in", round((time.time() - start) / 3600, 2), "h")
os.chdir("..")
if os.path.exists("xtbtmp"):
    shutil.rmtree("xtbtmp")