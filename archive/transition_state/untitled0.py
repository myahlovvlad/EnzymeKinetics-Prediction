#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:46:04 2023

@author: user
"""

# Импорты
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

print(f"SMILES реагента: {r_smiles}")
print(f"SMILES продукта: {p_smiles}")

# Конвертация SMILES -> объекты mol
print("Конвертируем SMILES в объекты mol...")
r_mol = Chem.MolFromSmiles(r_smiles)
print("Объект mol реагента:", r_mol)

p_mol = Chem.MolFromSmiles(p_smiles)
print("Объект mol продукта:", p_mol)

# Добавление атомов H
print("Добавляем атомы H...")
r_mol = Chem.AddHs(r_mol)
print("Реагент с H:", r_mol)

p_mol = Chem.AddHs(p_mol)
print("Продукт с H:", p_mol)

# Генерация координат 
print("Генерируем координаты...")
AllChem.EmbedMolecule(r_mol)
print("Координаты реагента сгенерированы")

AllChem.EmbedMolecule(p_mol)
print("Координаты продукта сгенерированы")

# Получение координат
r_conf = r_mol.GetConformer()
r_coords = r_conf.GetPositions()
print("Координаты реагента:", r_coords)

p_conf = p_mol.GetConformer()
p_coords = p_conf.GetPositions()
print("Координаты продукта:", p_coords)


r_coords = np.pad(r_coords, ((0, 0), (0, 1)), mode='constant', constant_values=0)

coords = r_coords + p_coords


# Конвертация в mole объекты
print("Конвертируем в объекты mole...")
def mol_to_mole(mol):
  atoms = [a.GetSymbol() for a in mol.GetAtoms()]  
  coords = mol.GetConformer().GetPositions()

  return gto.M(atom="\n".join([f"{atoms[i]} {x} {y} {z}"  
                               for i, (x, y, z) in enumerate(coords)]), basis='sto-3g')

r_mole = mol_to_mole(r_mol)
p_mole = mol_to_mole(p_mol)
print("Объекты mole созданы")
                           
# Оптимизация геометрии
print("Оптимизация геометрии...")
mol_react = r_mole  
mol_prod = p_mole

mf_react = scf.RHF(mol_react)
mf_react.kernel()

mf_prod = scf.RHF(mol_prod)
mf_prod.kernel()

# Создание объекта переходного состояния
print("Создание объекта TS...")
atoms = [a.GetSymbol() for a in r_mol.GetAtoms()] + [a.GetSymbol() for a in p_mol.GetAtoms()]
coords = r_coords + p_coords
  
mol_ts = gto.M(atom="\n".join(["{} {} {} {}".format(atoms[i], coords[i][0], coords[i][1], coords[i][2]) for i in range(len(atoms))]))
print("Объект TS:", mol_ts)

# Оптимизация TS через NEB
print("Оптимизация TS...")  
neb = pyscf.neb.NEB(mf_react, mol_react, mol_ts, mf_prod, mol_prod)
neb.kernel()