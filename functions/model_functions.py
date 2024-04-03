# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 14:15:20 2023

@author: danie
"""
from rdkit import Chem
import math
import pandas as pd
from functions.base_functions import *
#from functions.chem_functions import *
'''
This file contains functions for
    the generation of data driven representative models

Functions contained in this file are (indented functions are used by the one above):
    
    all_model(molecules)
    FT_model(molecules, threshold)
    PT_model(molecules)
        calculate_normalised_total(all_info)
    AG_model(molecules)
    SG_model(group_model)
    rank_models(df)
        generate_model_df(molecules, models, model_types)

Further descriptions of each and their usagae can be found after each function.
'''
def get_heteroatom_content(smiles, class_attribute):
    """
    Calculates the weighted heteroatom content based on the specified class attribute.

    Args:
        smiles (str): SMILES representation of the molecule.
        class_attribute (str): The class attribute specifying the heteroatom content.

    Returns:
        float: Weighted heteroatom content.
    """
    het_weights = {"Heteroatom_content_O": 16, "Heteroatom_content_N": 14, "Heteroatom_content_S": 32}

    het_wt = sum(smiles.count(atom) * weight for atom, weight in {"O": 16, "N": 14, "S": 32}.items() if atom.lower() in smiles.lower())

    return het_wt * het_weights.get(class_attribute, 0)

def get_weighted_average(molecules, class_attribute, model_type):
    """
    Calculates the weighted average of a given attribute for a list of molecules.

    Args:
        molecules (list): List of molecules.
        class_attribute (str): The attribute for which the weighted average is calculated.
        model_type (str): The type of model.

    Returns:
        float: Weighted average of the specified attribute.
    """
    if isinstance(molecules[0], list):  # Check if molecules is a list of lists
        num_mols = len(molecules[0]) if model_type in {"group_model", "score_model", "propr_rep"} else len(molecules)
        peak_area = [get_group_area(molecules[1][i]) if model_type in {"group_model", "score_model"} else float(molecules[1][i]) for i in range(num_mols)]
        
        if "Heteroatom_content" in class_attribute:
            weights = [get_heteroatom_content(molecules[0][i].smiles) if model_type in {"group_model", "score_model"} else get_heteroatom_content(molecules[0][i].smiles) for i in range(num_mols)]
        else:
            weights = [float(getattr(molecules[0][i], class_attribute)) if model_type in {"group_model", "score_model"} else float(getattr(molecules[0][i], class_attribute)) for i in range(num_mols)]
    else:  # molecules is a list of individual objects
        num_mols = len(molecules)
        peak_area = [float(m.peak_area) for m in molecules] if model_type != "propr_rep" else [float(m) for m in molecules]
        
        if "Heteroatom_content" in class_attribute:
            weights = [get_heteroatom_content(m.smiles) for m in molecules]
        else:
            weights = [float(getattr(m, class_attribute)) for m in molecules]
    
    sum_peak_area = sum(peak_area)
    normalized_peak_area = [area / sum_peak_area for area in peak_area]
    weighted_average = sum(weight * area for weight, area in zip(weights, normalized_peak_area))
    
    return weighted_average


def gen_all_model(molecules):
    return(molecules)
'''
USAGE: all_model(molecules)
    molecules = list of orca class molecules

Returns:
    The same list of molecules (this is the benchmark "all-molecule model")
'''
def gen_FT_model(molecules, threshold):
    five_percent_model = []
    for molecule in molecules:
        if float(molecule.peak_area) >= threshold:
            five_percent_model.append(molecule)
    return(five_percent_model)
'''
USAGE: FT_model(molecules, threshold)
    molecules = list of orca class molecules
    threshold = integer to define the desired selection threshold
'''
def gen_PT_model(molecules):
       all_info = []
       peak_areas = [] # Just used for inital normalising calculation
       #print(molecules)
       for molecule in molecules:
           peak_areas.append(float(molecule.peak_area))
       for molecule in molecules:
           info = [] # mol object, Smiles, peak area, normalised
           info.append(molecule)
           info.append(molecule.smiles)
           info.append(molecule.peak_area)
           info.append(float(molecule.peak_area)/min(peak_areas))
           all_info.append(info)
       max_mols = len(molecules)
       
       def calculate_normalised_total(all_info):
           normalised_values = []
           for i in range(len(all_info)):
               normalised_values.append(all_info[i][3])
           return(sum(normalised_values))
       '''
       USAGE: calculate_normalised_total(all_info)
           all_info = [[molecule, molecule.smiles, molecule.peak_area, normalised_area], [as previous for next molecule and so forth...]]
           
       Returns:
           The sum of all normalised values
       '''
       for i in range(len(molecules)):
           normalised_total = calculate_normalised_total(all_info)
           #print(normalised_total)
           #print("Normalised total = ", normalised_total)
           if normalised_total > max_mols:
               min_peak_area = peak_areas.index(min(peak_areas))
               peak_areas.pop(min_peak_area)
               all_info.pop(min_peak_area)
               #print("Peak areas", peak_areas)
               #print("min peak area is", min(peak_areas))
           for j in range(len(all_info)):
               all_info[j][3] = float(all_info[j][2])/min(peak_areas)
       mols_to_return = [[], []] # mol objects, normalised values
       for i in range(len(all_info)):
           mols_to_return[0].append(all_info[i][0])
           mols_to_return[1].append(all_info[i][3])
       molecule_ratios = mols_to_return[1]
       summ = sum(molecule_ratios)
       for i in range(len(molecule_ratios)):
           molecule_ratios[i] = molecule_ratios[i]/summ
       return(mols_to_return)
'''
USAGE: PT_model(molecules)
    molecules = list of orca class molecules
    
Returns:
    A 2d array in this format:
        [[molecule_1, molecule_2, ...], [proportion_1, proportion_2, ...]]
        Where molecule_x is a molecule object and proportion_x is its relative proportion in the model
'''   
def group_molecules(molecule_class_list):
    import csv
    smiles_list = []
    for thing in molecule_class_list:
        a = thing.smiles
        smiles_list.append(thing.smiles)
    mols_with_n_and_o = [[], [], [], []]       
    mols_with_n = [[], [], [], []] # 5 RING, 6 RING, OTHER RING(MAYBE COMBO CORES), NO RING
    mols_with_o = [[], [], [], []]
    no_hetero = [[], [], [], []]       
    #print("SMILES ARE", smiles_list)
    N = "nitrogen"
    O = "oxygen"
    S = "sulfur"

    five = "5-membered ring"
    six = "6-membered ring"
     
    for molecule in molecule_class_list:
        mol = Chem.MolFromSmiles(molecule.smiles)
        
   # for i in range(len(smiles_list)):
    #    mol = Chem.MolFromSmiles(smiles_list[i])
        a = has_heteroatoms(mol)
        #print("a is", a)
        if N in a:
            if O in a:
                b = has_rings(mol)
               #print(b)
                if five in b and six not in b:
                    mols_with_n_and_o[0].append(molecule) # Pulls out 5 membered rings
                elif six in b and five not in b:
                    mols_with_n_and_o[1].append(molecule) # Pulls out 6 membered rings
                elif len(b) != 0:
                    mols_with_n_and_o[2].append(molecule) # Pulls out any rings that aren't just 5 or 6 (ie. indoles which have both five and six rings)
                else:
                    mols_with_n_and_o[3].append(molecule)
            elif O not in a:
                b = has_rings(mol)
               #print(b)
                if five in b and six not in b:
                    mols_with_n[0].append(molecule) # Pulls out 5 membered rings
                elif six in b and five not in b:
                    mols_with_n[1].append(molecule) # Pulls out 6 membered rings
                elif len(b) != 0:
                    mols_with_n[2].append(molecule) # Pulls out any rings that aren't just 5 or 6 (ie. indoles which have both five and six rings)
                else:
                    mols_with_n[3].append(molecule)           
            else:
                pass
        elif O in a:
            b = has_rings(mol)
           #print(b)
            if five in b and six not in b:
                mols_with_o[0].append(molecule) # Pulls out 5 membered rings
            elif six in b and five not in b:
                mols_with_o[1].append(molecule) # Pulls out 6 membered rings
            elif len(b) != 0:
                mols_with_o[2].append(molecule) # Pulls out any rings that aren't just 5 or 6 (ie. indoles which have both five and six rings)
            else:
                mols_with_o[3].append(molecule)
        else:
            b = has_rings(mol)
           #print(b)
            if five in b and six not in b:
                no_hetero[0].append(molecule) # Pulls out 5 membered rings
            elif six in b and five not in b:
                no_hetero[1].append(molecule) # Pulls out 6 membered rings
            elif len(b) != 0:
                no_hetero[2].append(molecule) # Pulls out any rings that aren't just 5 or 6 (ie. indoles which have both five and six rings)
            else:
                no_hetero[3].append(molecule)
    grouped_mols = []
    grouped_mols.append(mols_with_n_and_o)
    grouped_mols.append(mols_with_n)
    grouped_mols.append(mols_with_o)
    grouped_mols.append(no_hetero)
    return(grouped_mols)

def gen_AG_model(molecules):
    grouped_molecules = group_molecules(molecules)
    group_model = []
    all_mols_in_group = []
    groups_in_model = 0
    for j in range(len(grouped_molecules)): # for i in range number of groups
        for k in range(len(grouped_molecules[j])): # for i in range number of subgroups
            group = grouped_molecules[j][k]
            if len(group) == 0:
                pass
            else:
                all_mols_in_group.append(group)
                groups_in_model += 1
                molecule_to_model = group[0]
                for i in range(len(group)): # for i in range mols in group
                    if i == 0: # No mols = skip group
                        pass
                    else:
                        if float(group[i].peak_area) > float(molecule_to_model.peak_area):
                            molecule_to_model = group[i]
                        else:
                            pass
                        #print("MOL to model", molecule_to_model.name)
                group_model.append(molecule_to_model)
    group_peaks = [] # this needs to be the group area
    mw = []
    # This loop just get group peak areas
    for j in range(len(grouped_molecules)):
        for k in range(len(grouped_molecules[j])):
            group = grouped_molecules[j][k]
            if len(group) == 0:
                pass
            else:
                group_area = get_group_area(group)
                group_peaks.append(group_area)
    summ = sum(group_peaks)
    for i in range(len(group_peaks)):
        group_peaks[i] = group_peaks[i]/summ
    return(group_model, all_mols_in_group, group_peaks)
'''
USAGE: AG_model(molecules)
    molecules = list of orca class molecules
    
Returns:
    A 2d array in this format:
        [[molecule_1, molecule_2, ...], [[group_1], [group_2], ...], [proporion_1, proportion_2, ...]]
        Where molecule_x is a molecule object, group_x is a molecular subclass and proportion_x is its relative proportion of that molecular subclass the model
'''  
def gen_SG_model(group_model): # Takes group_model as an input as it is a avariant of that model
    selected_molecules = []
    all_groups = []
    group_peaks = []
    for i in range(len(group_model[1])): # for i in range number of groups
        scored_molecules= [[], []]
        model = group_model[1][i]
        all_groups.append(model)
        group_peak = get_group_area(model) # gets the total area of the group
        group_peaks.append(group_peak)
        # This part bypasses the scoring function if there is only one compound in the group
        # There is no need to score and compare values of one compound
        if len(model) == 1:
            scored_molecules[0].append(model[0])
            scored_molecules[1].append(1.0)
            selected_molecules.append(model[0])
            continue
        # Gets group average attributes
        oxygen_content = calculate_heteroatom_percentages(model)['O']
        wa_mw = get_weighted_average(model, "mw", gen_all_model) # wq = weighted average
        wa_tot_en = get_weighted_average(model, "total_energy", gen_all_model)
        wa_polar = get_weighted_average(model, "polarizability", gen_all_model)
        wa_dipole = get_weighted_average(model, "dipole_moment", gen_all_model)
        wa_chem_hard = get_weighted_average(model, "chemical_hardness", gen_all_model)
        mws = [float(mol.mw) for mol in model]
        min_mw = min(mws)
        max_mw = max(mws)
        #print("MIN", min_mw, "MAX", max_mw)
        normalized_mws = [(mw - min_mw) / (max_mw - min_mw) for mw in mws]
        avg_norm_mw = (wa_mw - min_mw) / (max_mw - min_mw)
        tot_ens = [float(mol.total_energy) for mol in model]
        min_tot_en = min(tot_ens)
        max_tot_en = max(tot_ens)
        normalized_tot_ens = [(tot_en - min_tot_en) / (max_tot_en - min_tot_en) for tot_en in tot_ens]
        avg_norm_tot_en = (wa_tot_en - min_tot_en) / (max_tot_en - min_tot_en)
        polars = [float(mol.polarizability) for mol in model]
        min_polar = min(polars)
        max_polar = max(polars)
        normalized_polars = [(polar - min_polar) / (max_polar - min_polar) for polar in polars]
        avg_norm_polar = (wa_polar - min_polar) / (max_polar - min_polar)
        dipoles = [float(mol.dipole_moment) for mol in model]
        min_dipole = min(dipoles)
        max_dipole = max(dipoles)
        normalized_dipoles = [(dipole - min_dipole) / (max_dipole - min_dipole) for dipole in dipoles]
        avg_norm_dipole = (wa_dipole - min_dipole) / (max_dipole - min_dipole)
        chem_hards = [float(mol.chemical_hardness) for mol in model]
        min_chem_hard = min(chem_hards)
        max_chem_hard = max(chem_hards)
        normalized_chem_hards = [(chem_hard - min_chem_hard) / (max_chem_hard - min_chem_hard) for chem_hard in chem_hards]
        avg_norm_chem_hard = (wa_chem_hard - min_chem_hard) / (max_chem_hard - min_chem_hard)
        oxys = [calculate_heteroatom_percentages_sing_mol(mol)['O'] for mol in model]
        min_oxy = min(oxys)
        max_oxy = max(oxys)
        if min_oxy == 0:
            if max_oxy ==0:
                key="no_oxy"
        else:
            key="oxy_present"
            normalized_oxys = [(oxy - min_oxy ) / (max_oxy  - min_oxy ) for oxy  in oxys]
            avg_norm_oxy = (oxygen_content - min_oxy) / (max_oxy - min_oxy)
        for j in range(len(model)):
            if key == "oxy_present":
                score_value = 1/(1+(math.sqrt((normalized_mws[j]-float(avg_norm_mw))**2 + (normalized_tot_ens[j]-float(avg_norm_tot_en))**2 +(normalized_polars[j]-float(avg_norm_polar))**2 + (normalized_dipoles[j]-float(avg_norm_dipole))**2 +(normalized_chem_hards[j]-float(avg_norm_chem_hard))**2 +(normalized_oxys[j]-float(avg_norm_oxy))**2)))
            else:
                score_value = 1/(1+(math.sqrt((normalized_mws[j]-float(avg_norm_mw))**2 + (normalized_tot_ens[j]-float(avg_norm_tot_en))**2 +(normalized_polars[j]-float(avg_norm_polar))**2 + (normalized_dipoles[j]-float(avg_norm_dipole))**2 +(normalized_chem_hards[j]-float(avg_norm_chem_hard))**2)))
            scored_molecules[0].append(model[j])
            scored_molecules[1].append(score_value)
        best_score = max(scored_molecules[1])
        best_score_index = scored_molecules[1].index(best_score)
        selected_molecules.append(scored_molecules[0][best_score_index])
    #print(scored_molecules[0][0].name, scored_molecules[1])
    return(selected_molecules, all_groups, group_peaks)
'''
USAGE: SG_model(group_model)
    group_model = AG_model
    AG_model can be obtained this way:
        AG_model(molecules)
        where, molecules = list of orca class molecules
    
Returns:
    A 2d array in this format:
        [[molecule_1, molecule_2, ...], [[group_1], [group_2], ...], [proporion_1, proportion_2, ...]]
        Where molecule_x is a molecule object, group_x is a molecular subclass and proportion_x is its relative proportion of that molecular subclass the model
'''
def generate_model_df(molecules, models, model_types): # molecules is all set, models are subsets of molecules
    mws = []
    chem_hards = []
    polars = []
    dipoles = []
    tot_ens = []
    oxygen_content = []
    for i in range(len(models)):
        mws.append(get_weighted_average(models[i], "mw", model_types[i]))
        chem_hards.append(get_weighted_average(models[i], "chemical_hardness", model_types[i]))
        polars.append(get_weighted_average(models[i], "polarizability", model_types[i]))
        dipoles.append(get_weighted_average(models[i], "dipole_moment", model_types[i]))
        tot_ens.append(get_weighted_average(models[i], "total_energy", model_types[i]))
        if model_types[i] == "group_model":
            oxy = calculate_heteroatom_percentages(models[i][0])["O"]
        elif model_types[i] == "score_model":
            oxy = calculate_heteroatom_percentages(models[i][0])["O"]
        elif model_types[i] == "propr_rep":
            oxy = calculate_heteroatom_percentages(models[i][0])["O"]
        else:
            oxy = calculate_heteroatom_percentages(models[i])["O"]
        oxygen_content.append(oxy)
    # Create a DataFrame with model properties
    df = pd.DataFrame({
        "Model Type": model_types,
        "Mw": mws,
        "Chem_hard": chem_hards,
        "polarizability": polars,
        "Dipole": dipoles,
        "total_energy": tot_ens,
        "oxygen_content": oxygen_content
    })
    return df
'''
USAGE: generate_model_df(molecules, models, model_types)
    molecules = all_molecules from the DFT csv
    models = a list of generated models
        i.e. models = [all_model, FT_model, PT_model, AG_model, SG_model]
    model_types = a list of model types in string format
        i.e. model_types = ["all_model", "five_percent_model", "propr_rep", "group_model", "score_model"]

Returns:
    A dataframe with the average properties from each model passed to the function
'''
def rank_models(df):
    properties = ["Mw", "Chem_hard", "polarizability", "Dipole", "total_energy", "oxygen_content"]
    # Create a new dataframe to store the calculated differences and ranks
    differences_df = pd.DataFrame()
    ranks_df = pd.DataFrame()
    # Copy the "Model Type" column from the input dataframe to the output dataframes
    differences_df["Model Type"] = df["Model Type"]
    ranks_df["Model Type"] = df["Model Type"]
    # Calculate absolute differences for each property against the "all_model" row
    for prop in properties:
        differences_df[prop + "_Difference"] = df[prop].sub(df.loc[df["Model Type"] == "all_model", prop].iloc[0]).abs()
        # Calculate ranks for each descriptor against the "all_model" benchmark
        ranks_df[prop + "_Rank"] = differences_df[prop + "_Difference"].rank()
    # Calculate cumulative difference for each model based on the property differences
    differences_df["Cumulative_Difference"] = differences_df[[prop + "_Difference" for prop in properties]].sum(axis=1)
    # Calculate final score for each model
    min_diff = differences_df["Cumulative_Difference"].min()
    differences_df["Final_Score"] = differences_df["Cumulative_Difference"].apply(lambda diff: (min_diff - diff) + 1)
    # Calculate the sum of ranks for each model
    ranks_df["Final_Rank"] = ranks_df[[prop + "_Rank" for prop in properties]].sum(axis=1)
    return ranks_df
'''
USAGE: rank_models(df)
    df = dataframe constructed from "generate_model_df(molecules, models, model_types)"

Returns:
    A ranked dataframe on how well the models perform against the all-molecules model
'''
def model_output_block(molecules, model, model_type, ranked_data):
    supported_models = ["five_percent_model", "all_model", "group_model", "propr_rep", "score_model"]
    # Set up a dictionary for the benchmark so I can score each model versus the benchmark
    benchmark = gen_all_model(molecules)
    benchmark_values = {"Mw":0, "Chem_hard":0, "polarizability":0, "Dipole":0, "total_energy":0, "oxygen_content":0}
    benchmark_values["Mw"] = get_weighted_average(benchmark, "mw", "all_model")
    benchmark_values["Chem_hard"] = get_weighted_average(benchmark, "chemical_hardness", "all_model")
    benchmark_values["polarizability"] = get_weighted_average(benchmark, "polarizability", "all_model")
    benchmark_values["Dipole"] = get_weighted_average(benchmark, "dipole_moment", "all_model")
    benchmark_values["total_energy"] = get_weighted_average(benchmark, "total_energy", "all_model")
    benchmark_values["oxygen_content"] = (calculate_heteroatom_percentages(benchmark))["O"]
    # Creation of a dictionary for the model_values
    model_values = {"Mw":0, "Chem_hard":0, "polarizability":0, "Dipole":0, "total_energy":0, "oxygen_content":0}
    if model_type == "group_model":
        num_mols = len(model[0])
    elif model_type == "score_model":
        num_mols = len(model[0])
    elif model_type == "propr_rep":
        num_mols = len(model[0])
    else:
        num_mols = len(model)
    spacer = ""
    molecule_header = "The molecules in this model are:"
    smiles_header = "The smiles for these molecules are:"
    model_lines = []
    recorded_peak_area = 0
    # This loop gets total peak area of the entire dataset passed to it
    for i in range(len(molecules)):
        recorded_peak_area = recorded_peak_area + float(molecules[i].peak_area)
    if model_type not in supported_models:
        print("This function currently does not support: ", model_type)
    else:
        model_line = "Model: " + model_type
        # This gets the peak area of the compounds in the model (or subset)
        model_area = 0
        for i in range(num_mols):
            if model_type == "group_model":
                model_area = model_area + float(model[0][i].peak_area)
            elif model_type == "score_model":
                model_area = model_area + float(model[0][i].peak_area)
            elif model_type == "propr_rep":
                model_area = model_area + float(model[0][i].peak_area)
            else:
                model_area = model_area + float(model[i].peak_area)
        # Representation of model in a string to write to file
        representation_line = "This model represents " + str((model_area/recorded_peak_area)*100) + "% of the charecterised biocrude"
        # Other properties of model to write to line
        weighted_line = "Weighted values for molecular weight, chemical hardness, polarizability, chemical hardness and dipole moment are listed below:"
        mw = "Molecular weight  = " + str(get_weighted_average(model, "mw", model_type))
        chem_hard = "Chemical hardness = " + str(get_weighted_average(model, "chemical_hardness", model_type))
        polar = "Polarizability    = " + str(get_weighted_average(model, "polarizability", model_type))
        dipole = "Dipole moment     = " + str(get_weighted_average(model, "dipole_moment", model_type))
        tot_energy = "Total energy      = " + str(get_weighted_average(model, "total_energy", model_type))
        if model_type == "group_model":
            model_values["oxygen_content"] = (calculate_heteroatom_percentages(model[0]))["O"]
        elif model_type == "score_model":
            model_values["oxygen_content"] = (calculate_heteroatom_percentages(model[0]))["O"]
        elif model_type == "propr_rep":
            model_values["oxygen_content"] = (calculate_heteroatom_percentages(model[0]))["O"]
        else:
            model_values["oxygen_content"] = (calculate_heteroatom_percentages(model))["O"]
        model_values["Mw"] = get_weighted_average(model, "mw", model_type)
        model_values["Chem_hard"] = get_weighted_average(model, "chemical_hardness", model_type)
        model_values["polarizability"] = get_weighted_average(model, "polarizability", model_type)
        model_values["Dipole"] = get_weighted_average(model, "dipole_moment", model_type)
        model_values["total_energy"] = get_weighted_average(model, "total_energy", model_type)
        # Functions to calculate how good a model is versus the benchmark
        
        model_score = ranked_data.loc[ranked_data["Model Type"] == model_type, "Final_Rank"].values[0]
        score_line = "This model has a score of: " + str(model_score/6) + " versus the benchmark"
        model_lines.append(model_line)
        model_lines.append(representation_line)
        model_lines.append(score_line)
        model_lines.append(spacer)
        model_lines.append(weighted_line)
        model_lines.append(mw)
        model_lines.append(chem_hard)
        model_lines.append(polar)
        model_lines.append(dipole)
        model_lines.append(tot_energy)
        # Now we write heteroatoms content to file
        model_lines.append(spacer)
        model_lines.append("The heteroatom content of this model is:")
        if model_type == "group_model":
            heteroatom_content = calculate_heteroatom_percentages(model[0])
        elif model_type == "score_model":
            heteroatom_content = calculate_heteroatom_percentages(model[0])
        elif model_type == "propr_rep":
            heteroatom_content = calculate_heteroatom_percentages(model[0])
        else:
            heteroatom_content = calculate_heteroatom_percentages(model)
        for atom, weight in heteroatom_content.items():
            atom_content = f"{atom}: {weight}"
            model_lines.append(atom_content)
        num_atoms_req = min_atoms_4_simulation(model, model_type)
        model_lines.append(spacer)
        atoms_line = "This model will require the following minimum number of atoms for unbiased simulation:"
        model_lines.append(atoms_line)
        model_lines.append(str(num_atoms_req))
        model_lines.append(spacer)
        atoms_line = "This model will require the following minimum number of molecules for unbiased simulation:"
        model_lines.append(atoms_line)
        num_mols_req = min_mols_4_simulation(model, model_type)
        model_lines.append(str(num_mols_req))
        model_lines.append(spacer)
        model_lines.append(molecule_header)
        for i in range(num_mols):
            if model_type == "group_model":
                group_area = get_group_area((model[1][i]))
                line = model[0][i].name + " and this group accounts for: " + str(group_area) + " % of the biocrude"
                model_lines.append(line)
            elif model_type == "score_model":
                group_area = get_group_area((model[1][i]))
                line = model[0][i].name + " and this group accounts for: " + str(group_area) + " % of the biocrude"
                model_lines.append(line)
            elif model_type == "propr_rep":
                model_lines.append(model[0][i].name)
            else:
                model_lines.append(model[i].name)
        model_lines.append(spacer)
        model_lines.append(smiles_header)
        for i in range(num_mols):
            if model_type == "group_model":
                model_lines.append(model[0][i].smiles)
            elif model_type == "score_model":
                model_lines.append(model[0][i].smiles)
            elif model_type == "propr_rep":
                model_lines.append(model[0][i].smiles)
            else:
                model_lines.append(model[i].smiles)
        return(model_lines)
'''
USAGE: model_output_block(molecules, model, model_type, ranked_data)
    molecules = list of molecules objects from initial csv
    model = model (list of molecules objects)
    model_type = string
    ranked_data = ranked dataframe coming from "rank_models(df)"
'''
def write_output(filepath, lines):
    f = open(filepath, "w")
    for line in lines:
        f.write(line)
        f.write('\n')
    f.close()
    return()
'''
USAGE: write_output(filepath, lines)
    filepath = where you want to write a file and filename
    lines = list of things to write to file from "model_output_block(molecules, model, model_type, ranked_data)"
    Writes a list of strings to a file
    
Returns:
    information written to an output file
'''
def min_atoms_4_simulation(model, model_type):
    if model_type == "all_model" or model_type == "five_percent_model":
        min_peak = find_minimum(model, 'peak_area')   
        num_atoms = 0
        for i in range(len(model)):           
            mol = Chem.MolFromSmiles(model[i].smiles)
            mol = Chem.AddHs(mol)          
            num_atoms = num_atoms + (mol.GetNumAtoms()*(float(model[i].peak_area) / float(min_peak))*(len(model) + 1))
            
    elif model_type == "group_model" or model_type == "score_model":
        peaks = model[2]
        peaks = [get_group_area(thing) for thing in model[1]]
        
        compounds = model[0]
        min_peak = min(peaks)

        num_atoms = 0
        for i in range(len(compounds)):
            mol = Chem.MolFromSmiles(compounds[i].smiles)
            mol = Chem.AddHs(mol)
            num_atoms = num_atoms + (mol.GetNumAtoms()*(float(peaks[i]) / float(min_peak))*(len(compounds) + 1))
            
    elif model_type == "propr_rep":  # Assuming "propr_rep" is a string
        model = model[0]
        min_peak = find_minimum(model, 'peak_area')
        num_atoms = 0
        for i in range(len(model)):
            mol = Chem.MolFromSmiles(model[i].smiles)
            mol = Chem.AddHs(mol)
            num_atoms = num_atoms + (mol.GetNumAtoms()*(float(model[i].peak_area) / float(min_peak))*(len(model) + 1))

    return int(num_atoms)
