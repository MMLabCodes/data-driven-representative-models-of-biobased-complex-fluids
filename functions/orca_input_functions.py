# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 14:59:36 2023

@author: danie
"""
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdMolTransforms import CanonicalizeConformer
from rdkit.Chem.rdDistGeom import EmbedMolecule
import ase as ase
from ase.io import read,write
from ase.io.xyz import write_xyz
'''
This file contains functions to prepare and generate orca input files.

Functions contained in this file are (indented functions are used by the one above):
    
    clean_names(list_of_molecule_names)
    
    make_dft_inp(name, smiles, functional, basis_set, input_dir)
        smiles2xyz(smiles, xyz_filepath)
        get_homo_lumo_from_xyz(xyz_filepath)

Further descriptions of each and their usagae can be found after each function.
'''
def clean_names(list_of_molecule_names):
    for i in range(len(list_of_molecule_names)):
        name  = list_of_molecule_names[i]
        name = name.replace('(', '')
        name = name.replace(')', '')
        name = name.replace(',', '')
        name = name.replace(' ', '')
        list_of_molecule_names[i] = name
    return(list_of_molecule_names)
'''
USAGE: clean_names(list_of_molecule_names)

Returns: The same list where any characters uninterperable by ORCA/filepathing are omitted

Note: This removes commas, brackets and blank spaces from the molecule names.
'''
def make_dft_inp(name, smiles, functional, basis_set, input_dir):
    name = name
    if not os.path.isdir(input_dir):
        os.mkdir(input_dir)
    molecule_dir = input_dir + "/" + name
    if not os.path.isdir(molecule_dir):
        os.mkdir(molecule_dir)
    final_xyz_path = input_dir + "/xyz"
    if not os.path.isdir(final_xyz_path):
        os.mkdir(final_xyz_path)
    xyz_name = name + ".xyz"
    inp_name = name + ".inp"
    xyz_filepath = final_xyz_path + "/" + xyz_name
    input_filepath = molecule_dir + "/" + inp_name   
    smiles2xyz(smiles, xyz_filepath)
   
    homo_num = (get_homo_lumo_from_xyz(xyz_filepath))[0]
    lumo_num = (get_homo_lumo_from_xyz(xyz_filepath))[1]
    xyz_line = '*xyzfile 0 1 ' + name + ".xyz"
    homo_cube_name = name + '.homo.cube", ' + str(homo_num)
    lumo_cube_name = name + '.lumo.cube", ' + str(lumo_num)
    dens_name = name
    header = '! ' + functional + ' ' + basis_set + ' keepdens opt CPCM(water)'
    lines = [header, '%scf', ' MaxIter 1000', 'end', '%output', ' Print[ P_Hirshfeld] 1', 'end', '%elprop', ' Polar 1', 'end', '%plots', ' dim1 100', ' dim2 100', ' dim3 100', ' Format Gaussian_Cube', ' ElDens("%s.dens.cube");' % dens_name, ' MO("%s, 0);' % homo_cube_name, ' MO("%s, 0);' % lumo_cube_name, 'end', '%pal', ' nprocs 10', 'end']
    input_filepath = molecule_dir + "/" + name + ".inp"
    lines.append(xyz_line)
    z = open(input_filepath, "w")
    for line in lines:
        z.write(line)
        z.write('\n')
    z.close()
    return()
'''
USAGE: make_dft_inp(name, smiles, functional, basis_set, input_dir)

Returns:
    Nothing explicit is returned but and input file and xyz file are generated for the smiles in the specified directory
'''
def smiles2xyz(smiles, xyz_filepath):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol, nFlipsPerSample = 10, nSample = 500, bondLength = 3.0)
    EmbedMolecule(mol, useRandomCoords = True, forceTol = 0.000075)
    for conf in mol.GetConformers():
        rdMolTransforms.CanonicalizeConformer(conf)
    xyz = Chem.MolToXYZFile(mol, xyz_filepath)
    return(xyz)
'''
USAGE: smiles2xyz(smiles, xyz_filepath)
    
Returns:
    xyz of specified smiles string and generates an xyz file
    
Note:
    This function is used by the 'make_dft_inp' function and xyz filepaths are specified automatically in that function
    If used alone, any xyz filepath can be used
'''
def get_homo_lumo_from_xyz(xyz_filepath):
    atomic_numbers = {'C': 6, 'H': 1, 'O': 8, 'N':14}
    atoms = []
    with open(xyz_filepath, 'r') as xyz_file:
        for i, line in enumerate(xyz_file):
            if i == 0:
                continue
            if i == 1:
                continue
            else:
                atom = (line.split(" ")[0]).strip()
                atoms.append(atom)
    atomic_num = sum(atomic_numbers[atom] for atom in atoms)
    homo_num = int(atomic_num/2 - 1)
    lumo_num = homo_num + 1
    return(homo_num, lumo_num)
'''
USAGE: get_homo_lumo_from_xyz(xyz_filepath)

Returns:
    homo/lumo numbers of the molecule xyz file belongs to   
    
Note:
    More atomic numbers may need to be added to the dictionary as seen below, the current dictionary is an example for pine bark only.
        EXPANDED dictionaty example: atomic_numbers = {'C': 6, 'H': 1, 'O': 8, 'S':16, 'F':9, 'N':14, 'Cl':17}
'''