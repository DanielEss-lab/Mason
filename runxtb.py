# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:37:40 2021

@author: chens
"""

import subprocess
import os
import shutil
import configparser

configFile = "Pt_Methane.ini"
config = configparser.ConfigParser()
config.read(configFile)

# XTB=["wsl", "-d", "ubuntu", "-e", "/opt/xtb-6.4.0/bin/xtb"]
XTB=["xtb"]

freezeStrs = config['Properties']['freeze'].split(',')
freezeAtoms = []
for freezeStr in freezeStrs:
    for atomNum in freezeStr.split('-'):
        if atomNum not in freezeAtoms:
            freezeAtoms.append(atomNum)
freezeAtoms.sort()


os.environ['OMP_STACKSIZE'] = '1G'
os.chdir(config["IODir"]["Output"])
if os.path.exists("xtbtmp"):
    shutil.rmtree("xtbtmp")
os.makedirs("xtbtmp")
if os.path.exists("xtbopt"):
    print("Directory xtbopt found. Check if you still need them from the previous runs.")
    print("Make sure you don't need it. Please delete the directory by hand.")
    exit()
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
i = 1
for file in files:
    shutil.copyfile("../" + file, "temp.mol")
    subprocess.run(XTB + ["temp.mol", "--opt", "--gfnff", "--input", "freeze.inp"], stdout=devNull)
    os.remove("temp.mol")
    os.rename("xtbopt.mol", "temp.mol")
    subprocess.run(XTB + ["temp.mol", "--opt", "--input", "freeze.inp"], stdout=devNull )
    os.remove("temp.mol")
    shutil.move("xtbopt.mol", "../xtbopt/" + file)
    i += 1
    if i % int(totalItems / 100 + 1) == 0:
        print(str(i / totalItems * 100)+"% finished")
