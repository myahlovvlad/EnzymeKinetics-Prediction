#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 09:22:29 2023

@author: user
"""

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import pyscf
from pyscf import gto, scf, grad, hessian
import numpy as np

# Константы 
kappa = 1.380649e-23
R = 8.314462618  
T = 298.15

# SMILES
r_smiles = 'CCO'  
p_smiles = 'CC(=O)O'

# Конвертация SMILES -> объекты mol
r_mol = Chem.MolFromSmiles(r_smiles)
p_mol = Chem.MolFromSmiles(p_smiles)

# Генерация координат
AllChem.EmbedMolecule(r_mol)
AllChem.EmbedMolecule(p_mol)

# Конвертация в объекты mole
def mol_to_mole(mol):
  atoms = [a.GetSymbol() for a in mol.GetAtoms()]
  coords = mol.GetConformer().GetPositions()
  
  return gto.M(atom="\n".join([f"{atoms[i]} {x} {y} {z}"  
                               for i, (x, y, z) in enumerate(coords)]), basis='sto-3g')
                               
r_mole = mol_to_mole(r_mol)
p_mole = mol_to_mole(p_mol)                           

# Оптимизация геометрии 
mol_react = r_mole
mol_prod = p_mole

mf_react = scf.RHF(mol_react)
mf_react.kernel()

mf_prod = scf.RHF(mol_prod)
mf_prod.kernel()

# Поиск переходного состояния
mol_ts = pyscf.lib.param.linear_interpolate(mol_react, mol_prod)

conv, e_neb, ts = pyscf.neb.neb(mf_react, mol_react, mol_ts, mf_prod, mol_prod)

# Расчет частот
hessian_ts = hessian.RHF(ts).run()
frequencies = np.linalg.eigh(hessian_ts)[0]

# Расчет константы скорости 
DeltaH = frequencies[0]
DeltaS = np.sum(frequencies[1:])
DeltaG = DeltaH - T*DeltaS

rate = kappa * np.exp(-DeltaG/(R*T))

print("Rate constant is: "+str(rate)+" M^-1 s^-1")