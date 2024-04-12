# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 15:08:58 2023

@author: danie
"""
import os
from rdkit import Chem
from functions.base_functions import *
"""
This file contains class definitions and functions utilized by them.

Classes contained in this file are:

    - model_dirs: Filepath management class for directories related to a model.
        Note: This is a simplified version and can be expanded to include directories for DFT data and other files.

Functions contained in this file are:

    - csv_to_orca_class: Converts information from a CSV file to instances of the 'orca_molecule' class.

Further descriptions of each and their usage can be found after each function.
"""
class model_dirs:
    """
    Filepath management class for directories related to a model.

    Attributes:
        base_dir (str): Base directory path. (master path of the project)
        model_name (str): Name of the model.
        model_dir (str): Path to the model directory. (this is where model output files will be placed after generation)
        results_csv (str): Path to the results CSV file. (this file is not currently generated but will contain information from the)
            Note: All required data is found in the output files.
    """
    def __init__(self, base_dir, model_name):
        """
        Initializes model directory path manager.

        Args:
            base_dir (str): Base directory path. (master path of the project)
            model_name (str): Name of the model. 
        """
        self.base_dir = base_dir
        self.model_dir = os.path.join(self.base_dir, model_name)
        #self.results_csv = os.path.join(self.base_dir, (model_name + "_model_results.csv"))
       
    def gen_model_dirs(self):
        """
        Generate necessary directories for the model if they do not exist.
        """
        if not os.path.exists(self.model_dir):
            os.makedirs(self.model_dir)
        #if not os.path.exists(self.results_csv):
            #with open(self.results_csv, 'w'):
                #pass

class orca_molecule:
    """
    Represents a molecule with additional properties calculated using ORCA.

    Attributes:
        name (str): Name of the molecule.
        smiles (str): SMILES representation of the molecule.
        mw (float): Molecular weight of the molecule.
        peak_area (float): Peak area of the molecule.
        total_energy (float): Total energy of the molecule.
        homo_lumo_gap (float): HOMO-LUMO gap of the molecule.
        chemical_hardness (float): Chemical hardness of the molecule.
        dipole_moment (float): Dipole moment of the molecule.
        polarizability (float): Polarizability of the molecule.
        volume (float): Volume of the molecule.
    """
    def __init__(self, name, smiles, mw, peak_area, total_energy, chemical_hardness, dipole_moment, polarizability, volume): # This is where the properties of the class are specified   
       """
        Initializes molecule properties calculated using ORCA.

        Args:
            name (str): Name of the molecule.
            smiles (str): SMILES representation of the molecule.
            mw (float): Molecular weight of the molecule.
            peak_area (float): Peak area of the molecule.
            total_energy (float): Total energy of the molecule.
            homo_lumo_gap (float): HOMO-LUMO gap of the molecule.
            chemical_hardness (float): Chemical hardness of the molecule.
            dipole_moment (float): Dipole moment of the molecule.
            polarizability (float): Polarizability of the molecule.
            volume (float): Volume of the molecule.
       """
       self.name = name
       self.smiles = smiles
       self.mw = mw
       self.peak_area = peak_area
       self.total_energy = total_energy
       self.chemical_hardness = chemical_hardness
       self.dipole_moment = dipole_moment
       self.polarizability = polarizability
       self.volume = volume

def csv_to_orca_class(csv_file):
    """
    Converts information from a CSV file containing calculated DFT properties to instances of the 'orca_molecule' class.

    Args:
        csv_file (str): Path to the CSV file containing molecule information and calculated DFT properties.

    Returns:
        list: List of 'orca_molecule' instances.
    """
    import csv
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)  # Read the header row

        # Map column names to their indices
        column_indices = {column_name: index for index, column_name in enumerate(header)}
        
        # Initiate a list for objects to be appended to
        molecules = []
        
        # Extract relevant info for each line in the csv file
        for row in reader:
            molecule_name = row[column_indices['Molecule']].strip()
            molecular_weight = round(float(row[column_indices['mw']]), 2)
            peak_area = row[column_indices['area']]
            smiles = row[column_indices['smiles']]
            tot_energy = row[column_indices['Total electronic energy (eV)']]
            homo = row[column_indices['HOMO (eV)']]
            lumo = row[column_indices['LUMO(ev)']]
            chemical_hardness = row[column_indices['Chemical hardness']]
            dipole_moment = row[column_indices['Dipole moment']]
            polarizability = row[column_indices['Polarizability']]
            
            # Calculate volume for each molecule - not used in model generation, but is useful for generated input files for molecular dynamics
            mol = Chem.MolFromSmiles(smiles)
            volume = vol_from_smiles(smiles)
            
            # Create molecule objects frome using extracted data from the csv
            m1 = orca_molecule(molecule_name, smiles, molecular_weight, peak_area, tot_energy,
                                chemical_hardness, dipole_moment, polarizability, volume)
            molecules.append(m1)

    return molecules
