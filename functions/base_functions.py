# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 10:46:41 2024

@author: danie
"""
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
"""
This file contains functions for the generation of data driven representative models.
These are general functions that can be used with RDkit objects and smilesstrings,
    they have been developed for use with the functions that gnerate data driven representative models.
    
The functions contained in this file are:
    vol_from_smiles
    vol_from_mol
    has_heteroatoms
    has_rings
    get_group_area
    calculate_heteroatom_percentages
    find_minimum
    min_mols_4_simulation
"""
def vol_from_smiles(smiles):
    """
    Calculate the volume of a molecule from its SMILES representation.
    
    Args:
        - smiles (str): The SMILES representation of the molecule.
    
    Returns:
        - volume (float): The volume of the molecule.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    volume = AllChem.ComputeMolVolume(mol) 
    return volume

def vol_from_mol(mol):
    """
    Calculate the volume of a molecule from its RDKit Mol object.
    
    Args:
        - mol (RDKit Mol): The RDKit Mol object representing the molecule.
    
    Returns:
        - volume (float): The volume of the molecule.
    """
    AllChem.EmbedMolecule(mol)
    volume = AllChem.ComputeMolVolume(mol)   
    return volume

def has_heteroatoms(mol):
    """
    Check if a molecule contains heteroatoms and return their names.
    
    Args:
        - mol (RDKit Mol): The RDKit Mol object representing the molecule.
    
    Returns:
        - heteroatoms_in_mol (list): List of heteroatoms present in the molecule.
        
    Note:
        - The list of heteroatoms can be expanded with more SMARTS strings and names of heteroatoms, these are merely examples.
    """
    heteroatom_smarts = ["[#7]", "[#8]", "[#16]"]
    heteroatom_names = ["nitrogen", "oxygen", "sulfur"]
    heteroatoms_in_mol = []
    for i in range(len(heteroatom_smarts)):
                   pattern_smarts = Chem.MolFromSmarts(heteroatom_smarts[i])
                   if mol.HasSubstructMatch(pattern_smarts) == True:
                       heteroatoms_in_mol.append(heteroatom_names[i])
    if len(heteroatoms_in_mol) == 0:
        heteroatoms_in_mol.append("No heteroatoms")
    return(heteroatoms_in_mol)

def has_rings(mol):
    """
    Check if a molecule contains rings and return their names.
    
    Args:
        - mol (RDKit Mol): The RDKit Mol object representing a molecule.
    
    Returns:
        - ring_groups_in_mol (list): List of ring names present in the molecule.
    """
    ring_functionals = ["[r5]", "[r6]", "[r7]", "[r8]"]
    ring_names = ["5-membered ring", "6-membered ring", "7-membered ring", "8-membered ring"]
    ring_groups = []
    ring_groups.append(ring_names)
    ring_groups.append(ring_functionals)
    ring_groups_in_mol = []
    for i in range(len(ring_groups[0])):
        pattern_smarts = Chem.MolFromSmarts(ring_groups[1][i])
        if mol.HasSubstructMatch(pattern_smarts) == True:
            ring_groups_in_mol.append(ring_groups[0][i])
    if len(ring_groups) != 0:
        return(ring_groups_in_mol)
    if len(ring_groups) == 0:
        no_list = ["N"]
        return(no_list)
    
def get_group_area(list_of_molecule_class_objects):
    """
    Calculate the total peak area of a list of molecule objects.
    
    Args:
        - list_of_molecule_class_objects (list): List of molecule objects generated with "csv_to_orca_class(csv_file)".
    
    Returns:
        - tot_peak_area (float): Total peak area of the molecules.
    """
    tot_peak_area = 0
    for m in list_of_molecule_class_objects:
        tot_peak_area = tot_peak_area + float(m.peak_area)
    return(tot_peak_area)

def calculate_heteroatom_percentages(model):
    """
    Calculate the percentage of heteroatoms in a model of molecules.
    
    Args:
        - model (list): List of molecule objects generated with "csv_to_orca_class(csv_file)" and processed with one of the model generation functions found in "model_functions.py"
    
    Returns:
        - atom_weights (dict): Dictionary containing the percentage of each heteroatom in the model.
    """
    mixture = get_mixture(model) 
    total_peak_area = get_group_area(model)
    
    atomic_weights = {'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'P': 30.97, 'S': 32.07, 'F': 18.99, 'Cl': 35.45, 'Br': 79.90, 'I': 126.90}
    heteroatoms = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']  # List of heteroatoms
      
    atom_weights = {}  
    for atom in heteroatoms:
        content = 0.0
        for compound, weight in mixture.items():
            count = compound.count(atom)
            content += (weight/total_peak_area) * count * atomic_weights[atom]
        atom_weights[atom] = content
     # Now need to get H content
    tot_h_wt = 0.0
    for compound, weight in mixture.items():
        hydrogen_weight = 0.0
        mol = Chem.MolFromSmiles(compound)
        mol_with_hydrogens = Chem.AddHs(mol)
        count = 0
        for atom in mol_with_hydrogens.GetAtoms():
                if atom.GetSymbol() == 'H':
                    count += 1
        content = (weight/total_peak_area) * count * 1.0
        tot_h_wt += content
    atom_weights['H'] = tot_h_wt
   
    heteroatoms.append('H') 
    total_weight = sum(atom_weights.values())
    for atom in heteroatoms:
        if total_weight == 0:
            atom_weights[atom] = 0
        else:
            atom_weights[atom] = (atom_weights[atom]/total_weight)*100
   
    return atom_weights

def get_mixture(model):   
    """
    Get a dictionary representing the mixture of molecules in a model.
    
    Args:
        - model (list): List of molecule objects generated with "csv_to_orca_class(csv_file)".
    
    Returns:
        - mixture (dict): Dictionary representing the mixture of molecules.
    """
    mixture = {obj.smiles: float(obj.peak_area) for obj in model}
    return(mixture)

def calculate_heteroatom_percentages_sing_mol(molecule):
    """
    Calculate the percentage of heteroatoms in a single molecule.
    
    Args:
        - molecule (Molecule): Single molecule object generated with "csv_to_orca_class(csv_file)".
    
    Returns:
        - atom_weights (dict): Dictionary containing the percentage of each heteroatom in the molecule.
    """
    smiles = molecule.smiles
    weight = molecule.mw
    atomic_weights = {'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'P': 30.97, 'S': 32.07, 'F': 18.99, 'Cl': 35.45, 'Br': 79.90, 'I': 126.90}
    heteroatoms = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']  # List of heteroatoms
    
    atom_weights = {}  
    for atom in heteroatoms:
        content = 0.0    
        count = smiles.count(atom)
        content += (count * atomic_weights[atom])/weight
        atom_weights[atom] = content
        
    tot_h_wt = 0.0
    mol = Chem.MolFromSmiles(smiles)
    mol_with_hydrogens = Chem.AddHs(mol)
    count = 0
    for atom in mol_with_hydrogens.GetAtoms():
            if atom.GetSymbol() == 'H':
                count += 1
    content = (count * 1.0)/weight
    tot_h_wt += content
    atom_weights['H'] = tot_h_wt
    
    for key in atom_weights:
        if isinstance(atom_weights[key], (int, float)):
            atom_weights[key] *= 100
    
    return atom_weights

def find_minimum(list_of_molecule_class_objects, attribute):
    """
    Find the minimum value of a given attribute among a list of molecules.
    
    Args:
        - list_of_molecule_class_objects (list): List of molecule objects generated with "csv_to_orca_class(csv_file)".
        - attribute (str): Attribute of the molecule objects to compare.
    
    Returns:
        - min_value (float): Minimum value of the attribute among the molecules.
    """
    min_value = float('inf')  # Initialize with a large value, should be large enough 0_o
    
    for obj in list_of_molecule_class_objects:
        attr_value = getattr(obj, attribute)
        try:
            attr_value = float(attr_value)
            min_value = min(min_value, attr_value)
        except ValueError:
            print(f"Skipping object: Invalid value for attribute '{attribute}'")
    return min_value

def min_mols_4_simulation(model, model_type):
    """
    Estimate the number of molecules required for a simulation based on peak area and model type.
    
    Args:
        - model (list): List of molecule objects generated with "csv_to_orca_class(csv_file)".
        - model_type (str): Type of the model.
    
    Returns:
        - num_compounds (int): Estimated number of molecules required for a simulation.
    """
    if model_type == "all_model" or model_type == "five_percent_model":
        min_peak = find_minimum(model, 'peak_area')
        num_compounds = 0
        for i in range(len(model)):
            num_compounds += (float(model[i].peak_area) / float(min_peak)) * (len(model) + 1)
    elif model_type == "group_model" or model_type == "score_model":
        peaks = model[2]
        min_peak = min(peaks)
        num_compounds = 0
        for i in range(len(model[0])):
            num_compounds += (float(peaks[i]) / float(min_peak)) * (len(model[0]) + 1)
    elif model_type == "propr_rep":  # Assuming "propr_rep" is a string
        model = model[0]
        min_peak = find_minimum(model, 'peak_area')
        num_compounds = 0
        for i in range(len(model)):
            num_compounds += (float(model[i].peak_area) / float(min_peak)) * (len(model) + 1)

    return int(num_compounds)