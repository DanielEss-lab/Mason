# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:40:40 2021

@author: chens
"""

import Mason
from openbabel import openbabel as ob
import time
import multiprocessing
import tqdm

#set the file relative to main.py
configFile = "Pt_Methane.ini"

##############################################################################
#                  Initialization                                            #
##############################################################################
# config = configparser.ConfigParser()
# config.read(configFile)
# template, ligandLib = Smarties.GetTemplateAndLigandLib(config)
# #generate the combinations of ligands from ligand library
# ligandCombs = {"": []}
# ligandCombs = Smarties.LigandIter(ligandCombs, ligandLib,\
#                                   int(config['common']['coordNum']))

##############################################################################
#                  Initialization                                            #
##############################################################################


##################for debug#####################
# print(".................")
# for name in ligandCombs:
#     print(name)
#     for ligand in ligandCombs[name]:
#         print(ligand.GetFirstAtom().GetType())
################################################



##############################################################################
#                  MultiProcessing      Begin                                #
##############################################################################
template, ligandCombs, config = Mason.Init(configFile)

def WriteNewTSS(name):
    conv = ob.OBConversion()            
    conv.SetInFormat("mol")
    mol = Mason.AddLigands2Mol(template, ligandCombs[name], config)
    conv.WriteFile(mol, config['IODir']['Output'] + \
                    template.GetAtom(int(config['Properties']['coreNum'])).GetType()\
                        + name + ".mol")
    mol.Clear()

    
if __name__=='__main__':
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    start = time.time()
    pbar = tqdm.tqdm(total = len(ligandCombs))
    update = lambda * args: pbar.update()
    # pool.imap_unordered(WriteNewTSS, ligandCombs)
    for name in ligandCombs:
        results = pool.apply_async(WriteNewTSS, (name, ), callback = update)
    pool.close()
    pool.join()
    pbar.close()
    results.get()
    print("Finished in", round((time.time() - start) / 60, 2), "min")

##############################################################################
#                  MultiProcessing      End                                  #
########################################################################### ###  
xtb = input("Do you wanna run xtb (y/n)?")
if xtb == 'y':
    import runxtb
if xtb == 'no':
    exit()



