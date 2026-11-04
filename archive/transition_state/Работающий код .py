#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 18:52:27 2023

@author: user
"""
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, MolFromSmiles  
import pyscf
from pyscf import gto, scf, grad, hessian
import numpy as np
r_smiles = 'CCO' # этиловый спирт
p_smiles = 'CC(=O)O' # уксусная кислота

print(f"SMILES реагента: {r_smiles}")
print(f"SMILES продукта: {p_smiles}")

# Конвертация в объекты mol

print("Конвертация SMILES в объекты mol...")

r_mol = Chem.MolFromSmiles(r_smiles)
print("Объект mol реагента создан")

p_mol = Chem.MolFromSmiles(p_smiles) 
print("Объект mol продукта создан")

# Добавление атомов H
print("Добавление атомов H...")

r_mol = Chem.AddHs(r_mol)
print("Атомы H добавлены для реагента")

p_mol = Chem.AddHs(p_mol)
print("Атомы H добавлены для продукта")

# Генерация координат
print("Генерация координат...")

AllChem.EmbedMolecule(r_mol)
print("Координаты сгенерированы для реагента")

AllChem.EmbedMolecule(p_mol)
print("Координаты сгенерированы для продукта")

# Получение координат
conf = r_mol.GetConformer()
r_coords = conf.GetPositions()

conf = p_mol.GetConformer()
p_coords = conf.GetPositions()

print("Координаты получены")

# Вывод объектов mol с координатами
print(r_mol)
print(p_mol)

# Конвертация в объекты mole
def mol_to_mole(mol):
  atoms = [a.GetSymbol() for a in mol.GetAtoms()]
  coords = mol.GetConformer().GetPositions()
  
  return gto.M(atom="\n".join([f"{atoms[i]} {x} {y} {z}" 
                               for i, (x, y, z) in enumerate(coords)]), basis='sto-3g')

r_mole = mol_to_mole(r_mol)
p_mole = mol_to_mole(p_mol)

print(r_mole)
print(p_mole)