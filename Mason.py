# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:18:29 2021

@author: Shusen Chen
"""

from openbabel import openbabel as ob
import random
import os
import configparser


def GetTemplateAndLigandLib(config): 
    #read files and alloc memory to store template and ligands
    # initialize template
    conv = ob.OBConversion()            
    mol = ob.OBMol()
    conv.SetInFormat(config['IODir']['Template'].split('.')[-1])
    conv.ReadFile(mol, config['IODir']['Template'])
    core = mol.GetAtom(int(config['Properties']['coreNum']))
    core.SetFormalCharge(int(config['Properties']['charge']))
    core.SetSpinMultiplicity(int(config['Properties']['multiplicity']))

    # initialize ligands
    ligandLib = {}
    for filename in os.listdir(config['IODir']['LigandLibDir']):
        if '.mol' in filename:
            ligand = ob.OBMol()
            conv.SetInFormat("mol")
            conv.ReadFile(ligand, config['IODir']['LigandLibDir'] + filename)
            ligandLib[filename.split(".")[0]] = ligand
            
    return mol, ligandLib

def Init(configFile):
    config = configparser.ConfigParser()
    config.read(configFile)
    template, ligandLib = GetTemplateAndLigandLib(config)
    if not os.path.exists(config['IODir']['Output']):
        os.makedirs(config['IODir']['Output'])    
    #generate the combinations of ligands from ligand library
    originLigandCombs = {"": []}
    originLigandCombs = LigandIter(originLigandCombs, ligandLib,\
                                  int(config['common']['coordNum']))
    ligandCombs = {}
    ligandCombsNameSplit = []
    for name in originLigandCombs:
        temp = name.split("_")[1:]
        temp.sort()
        if temp not in ligandCombsNameSplit:
            ligandCombsNameSplit.append(temp)
            ligandCombs[name] = originLigandCombs[name]
    return template, ligandCombs, config

def LigandIter(ligandCombs, ligandLib, coordNum):

    ligandMaxCoordNum = 0
    ligandMinCoordNum = 10
    for name in ligandLib:
        ligandCoordNum = len(ligandLib[name].GetTitle().split(","))
        if ligandCoordNum < ligandMinCoordNum:
            ligandMinCoordNum = ligandCoordNum
        if ligandCoordNum > ligandMaxCoordNum:
            ligandMaxCoordNum = ligandCoordNum
    # print(ligandMinCoordNum, ligandMaxCoordNum)
    if coordNum < ligandMinCoordNum:
        # print("type1")
        return ligandCombs
    
    newLigandCombs = {}
    for i in range(ligandMinCoordNum, min(coordNum, ligandMaxCoordNum) + 1): 
        # print("iter1 ligandCoord =", i)
        tempLigandCombs = {}
        for name in ligandLib:
            # print("iter2", name)
            if len(ligandLib[name].GetTitle().split(",")) == i:
                for item in ligandCombs:
                    # print("******************")
                    # print("coordNum = ", coordNum)
                    # print("item = ", item)
                    # print("name = ", name)
                    newitem = item + "_" + name
                    tempLigandCombs[newitem] = []
                    # for a in ligandCombs:
                    #     print("ligandCombs", a)
                    tempLigandCombs[newitem] = list(ligandCombs[item])
                    tempLigandCombs[newitem].append(ligandLib[name])
                    # for a in ligandCombs:
                    #     print("ligandCombs", a)
                    # for a in tempLigandCombs:
                    #     print("tempLigandCombs", a)
                    # print("*************")
                # print(len(newitem))
        # coordNum -= i
        newLigandCombs.update(LigandIter(tempLigandCombs, ligandLib, coordNum - i))
    # print(len(newligands2added))
    # print("Type 2")
    return newLigandCombs



def AddLigand2Mol(mol, ligand, constraints, config):
    #add ligand(OBmol) to mol(OBMol) with constraints by freezing all atoms in mol
    #define which atoms to be connected
    molNConnect = int(config['Properties']['coreNum'])
    ligandNConnects = [int(i) for i in ligand.GetTitle().split(',')]
    # print(ligandNConnects)
    #params which affect the running time, the less the faster
    randomRepeatTimes = int(config['forcefield']['randomRepeatTimes'])
    randomPlacementScale = float(config['forcefield']['randomPlacementScale'])
    stepsPerAtom = int(randomPlacementScale)
    #initial
    vector = ob.vector3()
    minE = float('inf')
    tempLigand = ob.OBMol(ligand)
    minmol = ob.OBMol()
    forcefield = ob.OBForceField.FindForceField(config['forcefield']['ff4AddLigand'])    
    #randomly translating ligand coords to avoid crush
    for i in range(randomRepeatTimes):
        tempmol = ob.OBMol(mol)
        vector.Set(random.random() - .5, random.random() - .5, random.random() - .5)
        vector *= (2 * randomPlacementScale)
        tempLigand.Translate(vector)
        tempmol += tempLigand
        for ligandNConnect in ligandNConnects:
            tempmol.AddBond(molNConnect, mol.NumAtoms() + ligandNConnect, 1)
        #set forcefield with constraints
        forcefield.SetConstraints(constraints)
        forcefield.Setup(tempmol, constraints)
        # print(tempmol.NumAtoms())
        forcefield.ConjugateGradients(stepsPerAtom * tempmol.NumAtoms())
        forcefield.GetCoordinates(tempmol)
        if forcefield.Energy() < minE:
            minE = forcefield.Energy()
            minmol.Clear()
            minmol = tempmol
        else:
            tempmol.Clear()
    tempLigand.Clear()
    return minmol


def AddLigands2Mol(mol, ligands, config):
    # add several ligands(list of OBmols) to mol(OBmol)
    # read config to setup forcefield and steps
    stepsPerAtom = int(config['forcefield']['stepsPerAtom'])
    # set constrains for the atoms from mol
    constraints = ob.OBFFConstraints()
    for i in range(mol.NumAtoms()):
        constraints.AddAtomConstraint(i + 1)
    if config['forcefield']['ff4WholeMol'] != 'NO':
        # print("preforming whole molecule optimazation...")
        forcefield = ob.OBForceField.FindForceField(config['forcefield']['ff4WholeMol'])
        forcefield.SetConstraints(constraints)
    #add ligand one by one with crude optimization
    minE = float('inf')
    minmol = ob.OBMol()
    for i in range(int(config['forcefield']['randomRepeatTimesWholeMol'])):
        tempmol = ob.OBMol(mol)
        # print(i)
        for ligand in ligands:
            tempmol = AddLigand2Mol(tempmol, ligand, constraints, config)
            ## constraints.AddAtomConstraint(mol.NumAtoms() - ligand.NumAtoms() + 1)
        # Optimizing with more expensive method here
        forcefield.Setup(tempmol, constraints)
        # print((mol.NumAtoms() + ligand.NumAtoms()))
        # print(forcefield.Energy())
        forcefield.SteepestDescent(stepsPerAtom * tempmol.NumAtoms())
        # print(forcefield.Energy())
        forcefield.GetCoordinates(tempmol)
        forcefield.Setup(tempmol, constraints)
        forcefield.WeightedRotorSearch(tempmol.NumAtoms(), stepsPerAtom * tempmol.NumAtoms())
        forcefield.GetCoordinates(tempmol)
        forcefield.Setup(tempmol, constraints)
        # print(forcefield.Energy())
        forcefield.ConjugateGradients(stepsPerAtom * tempmol.NumAtoms())
        # print(forcefield.Energy())
        forcefield.GetCoordinates(tempmol)
        if forcefield.Energy() < minE:
            minE = forcefield.Energy()
            minmol.Clear()
            minmol = tempmol
        else:
            tempmol.Clear()        
    return minmol

# def WriteNewTSS(name, template, ligandCombs, config):
#     conv = ob.OBConversion()            
#     conv.SetInFormat("mol")
#     mol = AddLigands2Mol(template, ligandCombs[name], config)
#     conv.WriteFile(mol, config['IODir']['Output'] + \
#                    template.GetAtom(int(config['Properties']['coreNum'])).GetType()\
#                        + name + ".mol")
#     mol.Clear()


    
