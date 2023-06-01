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
from openbabel import openbabel as ob


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
conv = ob.OBConversion()            
conv.SetInFormat("mol")
conv.SetOutFormat("mol")
oF = open("failure.txt", "w")
oF.write("filename\treason\n")
for file in files:
    shutil.copyfile("../" + file, "temp.mol")
    subprocess.run(XTB + ["temp.mol", "--opt", "--gfnff", "--input", "freeze.inp"], stdout=devNull, stderr=devNull)
    os.remove("temp.mol")
    if "xtbopt.mol" not in os.listdir():
        pbar.update()
        oF.write(file + "\tgfnff\n")
        continue
    os.rename("xtbopt.mol", "temp.mol")
    if config['Miscs']['dummyAtoms'].split(',')[0].isdigit():
        dummyAtoms = config['Miscs']['dummyAtoms'].split(',')
        mol = ob.OBMol()
        conv.ReadFile(mol, "temp.mol")
        conv.SetInStream(None)
        atoms = []
        for atom in dummyAtoms:
            atoms.append(mol.GetAtom(int(atom)))
        for atom in atoms:
            mol.DeleteAtom(atom)
        conv.WriteFile(mol, "temp.mol")
        conv.SetOutStream(None)
        mol.Clear()
    subprocess.run(XTB + ["temp.mol", "--opt", "--input", "freeze.inp"], stdout=devNull, stderr=devNull)
    os.remove("temp.mol")
    if "xtbopt.mol" not in os.listdir():
        pbar.update()
        oF.write(file + "\tgfn2\n")
        continue
    shutil.move("xtbopt.mol", "../xtbopt/" + file)
    pbar.update()
pbar.close()
oF.close()
print("Finished in", round((time.time() - start) / 3600, 2), "h")
os.chdir("..")
if os.path.exists("xtbtmp"):
    shutil.rmtree("xtbtmp")
'''
if config['Miscs']['dummyAtoms'].split(',')[0].isdigit():
    dummyAtoms = config['Miscs']['dummyAtoms'].split(',')
    if os.path.exists("xtbopt_withoutdummy"):
        print("Directory xtbopt_withoutdummy found. Check if you still need them from the previous runs.")
        if(input("To delete results from previous runs (y/n)? ") != 'y'):
            exit()
        shutil.rmtree("xtbopt_withoutdummy")
    os.mkdir("xtbopt_withoutdummy")
    conv = ob.OBConversion()            
    conv.SetInFormat("mol")
    conv.SetOutFormat("mol")
    mol = ob.OBMol()
    files = []
    for file in os.listdir("xtbopt"):
        if ".mol" in file:
            files.append(file)
    for file in files:
        conv.ReadFile(mol, "xtbopt/" + file)
        atoms = []
        for atom in dummyAtoms:
            atoms.append(mol.GetAtom(int(atom)))
        for atom in atoms:
            mol.DeleteAtom(atom)
        conv.WriteFile(mol, "xtbopt_withoutdummy/" + file)
    mol.Clear()
'''    
