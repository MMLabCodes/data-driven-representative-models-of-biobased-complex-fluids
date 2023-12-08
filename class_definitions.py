# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 15:08:58 2023

@author: danie
"""
import os
from functions.chem_functions import vol_from_smiles, give_group_label
'''
This file contains class definitions of anything in this code alongside any functions ustilised by them.

Classes contained in this file are:
    
    model_dirs
    
Functions contained in this file are (indented functions are used by the one above):
    
    

Further descriptions of each and their usagae can be found after each function.
'''
class model_dirs:
    def __init__(self, base_dir, model_name):
        self.base_dir = base_dir
        self.dft_dir = os.path.join(base_dir, (model_name  + "_DFT"))
        self.dft_fukui_dir = os.path.join(base_dir, (model_name  + "_DFT"), "fukui")
        self.dft_results_dir = os.path.join(base_dir, (model_name  + "_DFT"), "results")
        self.dft_xyz_dir = os.path.join(base_dir, (model_name  + "_DFT"), "xyz")
        self.model_dir = os.path.join(base_dir, (model_name + "_models"))
'''
This class makes filepath management easier in the program and is used in the next function
'''     
def make_model_dirs(main_dir, model_name):
    model_dir = main_dir + "/" + model_name
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    subdirs = [(model_name  + "_DFT"), ((model_name  + "_DFT") + "/fukui"), ((model_name  + "_DFT") + "/results"), ((model_name  + "_DFT") + "/xyz"), (model_name + "_models")]
    for subdir in subdirs:
        subdir_to_make = os.path.join(model_dir, subdir)
        if not os.path.exists(subdir_to_make):
            os.makedirs(subdir_to_make)
    dirs = model_dirs(model_dir, model_name)
    return(dirs)
'''
USAGE: model_dirs = make_model_dirs(main_dir, model_name)
    main_dir = the main directory where you want different models to be created
    model_name = name of the model
    
Returns:
    A variable called "model_dirs" (in this instance)
    Different paths can be called as so:
        base_dir = model_dirs.base_dir
        dft_dir - model_dirs.dft_dir
        dft_fukui_dir = model_dirs.dft_fukui_dir
        dft_results_dir = model_dirs.dft_results_dir
        dft_xyz_dir = model_dirs.dft_xyz_dir
        
Note: "model_dirs" is an example, this variable can have any name
'''
class orca_molecule:
    def __init__(self, name, smiles, mw, peak_area, total_energy, homo_lumo_gap, chemical_hardness, dipole_moment, polarizability, volume, group): # This is where the properties of the class are specified
       self.name = name
       self.smiles = smiles
       self.mw = mw
       self.peak_area = peak_area
       self.total_energy = total_energy
       self.homo_lumo_gap = homo_lumo_gap
       self.chemical_hardness = chemical_hardness
       self.dipole_moment = dipole_moment
       self.polarizability = polarizability
       self.volume = volume
       self.group = group
'''
This class sets up a molecule object with the specified attributes
'''
def csv_to_orca_class(csv_file):
    import csv
    with open(csv_file, 'r') as file:
      reader = csv.reader(file)
      information = []
      for item in reader:
          information.append(item)
    # self, name, smiles, mw, core, r_groups, tracked_cuts, peak_area, total_energy, homo_lumo_gap, chemical_hardness, dipole_moment):
      molecules = []
      for i in range(len(information)):
          if i == 0:
              continue # skip header row
          else:
              #Molecule,mw,area,smiles,formula,Total electronic energy (eV),HOMO (eV),LUMO(ev),Chemical hardness,Dipole moment,Polarizability, volume, group_lable
              molecule_name = information[i][0].strip() # Strip removes white space from the string
              molecular_weight = round(float(information[i][1]), 2)
              peak_area = information[i][2]
              smiles = information[i][3]
              #formula = information[i][4] # Hashed out as not important for model gen
              tot_energy = information[i][5]
              homo = information[i][6]
              lumo = information[i][7]              
              homo_lumo_gap = abs(float(information[i][7]) - float(information[i][6]))
              chemical_hardness = information[i][8]
              dipole_moment = information[i][9]
              polarizability = information[i][10]
              volume = vol_from_smiles(smiles)             
              m1 = orca_molecule(molecule_name, smiles,  molecular_weight, peak_area, tot_energy, homo_lumo_gap, chemical_hardness, dipole_moment, polarizability, volume, give_group_label(smiles))
              molecules.append(m1)                               
    return(molecules)
'''
USAGE: csv_to_orca_class(csv_file)
    csv_file of DFT results
    
Returns:
    list of molecule objects with the attributes defined in "class Orca_molecule"
'''

