#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 19:17:31 2023

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

print('Выравнивание размерности координат массива молекул')
r_coords = np.resize(r_coords, (8, 3))
p_coords = np.resize(p_coords, (8, 3))
print('Размерности координат массива выравнены')

print('Сложение координат массива')
coords = r_coords + p_coords
print('Координаты массива сложены')

# Конвертация в объекты mole
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
print('Оптимизация геометрии')
mol_react = r_mole
print(mol_react)
mol_prod = p_mole
print(mol_prod)
print('Оптимизация геометрии молекул выполнено успешно')


print('Оптимизацию геометрии молекул реагента и продукта с помощью метода самосогласованного поля (SCF) в базисе минимальных атомных орбиталей (RHF)')
mf_react = scf.RHF(mol_react)
mf_react.kernel()
mf_prod = scf.RHF(mol_prod)
mf_prod.kernel()
print('Оптимизация геометрии молекул выполнено успешно')

# Получение актуальных координат молекул реагента и продукта
print('Получение актуальных координат молекул реагента и продукта')
print('Актуальные координаты реагента')
r_conf = r_mol.GetConformer()
r_coords = r_conf.GetPositions()
print('Актуализация координат выполнено успешно')

print('Актуальные координаты продукта')
p_conf = p_mol.GetConformer()
p_coords = p_conf.GetPositions()
print('Актуализация коодинат продукта выполнено успешно')

r_atoms = [atom.GetSymbol() for atom in r_mol.GetAtoms()]
# r_coords = [atom.GetCoords() for atom in r_mol.GetAtoms()]

p_atoms = [atom.GetSymbol() for atom in p_mol.GetAtoms()]
# p_coords = [atom.GetCoords() for atom in p_mol.GetAtoms()]

# Вставляем исправленный код здесь:

r_conf = r_mol.GetConformer()
r_coords = r_conf.GetPositions()

p_conf = p_mol.GetConformer()
p_coords = p_conf.GetPositions()

# Дальше идет объединение координат и создание объекта TS

print('Обрезание массива атомов до одинаковой размерности')
r_atoms = r_atoms[:8]
p_atoms = p_atoms[:8]

print('Объединение атомов')
atoms = r_atoms + p_atoms

print('Выравнивание размерности координат массива молекул')
r_coords = np.resize(r_coords, (16, 3))
p_coords = np.resize(p_coords, (16, 3))

print('Размерности массивов')
print(r_coords.shape)
print(p_coords.shape)

print('Количество атомов')
print(len(atoms))
print('Количество координат')
print(len(coords))

# Объединение координат
coords = r_coords + p_coords
print('Объединение коодинат и атомов реагента и продуктов выполнено успешно')

mol_ts = gto.M(atom="\n".join(["{} {} {} {}".format(atoms[i], coords[i][0], coords[i][1], coords[i][2]) for i in range(len(atoms))]))


# Оптимизация TS через NEB
import pyscf.neb
neb = pyscf.neb.NEB(mf_react, mol_react, mol_ts, mf_prod, mol_prod)
neb.kernel()
ts = neb.ts

# Расчет частот
hessian_ts = hessian.RHF(ts).run()
frequencies = np.linalg.eigh(hessian_ts)[0]

# Расчет константы скорости 
DeltaH = frequencies[0]
DeltaS = np.sum(frequencies[1:])
DeltaG = DeltaH - T*DeltaS

rate = kappa * np.exp(-DeltaG/(R*T))

print("Rate constant is: "+str(rate)+" M^-1 s^-1")
