# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:18:48 2024

@author: danie
"""
import os 
import csv
import ast
import time
import re
import subprocess
import shutil
import pandas as pd
import numpy as np

import ipywidgets as widgets
from IPython.display import display

'''
SnippetSim Directories Module

This module defines the SnippetSimManage class, which is used for initializing
and organizing directories for polymer simulation setup. It creates specific
directories such as 'python_scripts', 'pdb_files', 'molecules', 'systems', and manages
files like 'residue_codes.csv'.

Attributes:
    main_dir (str): The main directory for snippet simulation setup.

Example:
    >>> from PolymerSimulatorDirs import PolymerSimulatorDirs
    >>> polymer_dirs = PolymerSimulatorDirs('/path/to/main/dir/')
    >>> print(polymer_dirs.pdb_file_dir)
    '/path/to/main/dir/pdb_files'

Note:
    The main_dir must be provided as a valid path and must already exist.
'''
class SnippetSimManage:
    packmol_path = "/home/dan/packmol-20.14.4-docs1/packmol-20.14.4-docs1/packmol"
    # Class for intialsing and building files and folders for polymer simulation set up and the like
    def __init__(self, main_dir): 
        """
        Initialize SnippetSimManage object.

        Args:
            main_dir (str): The main directory for polymer simulation setup.
                            Example: '/path/to/main/dir/'

        Raises:
            FileNotFoundError: If main_dir does not exist.
        """
        if not os.path.exists(main_dir):
            raise FileNotFoundError(f"The specified main directory '{main_dir}' does not exist.")       
        self.main_dir = main_dir 
        
        # Get path to python script directory and check if exists
        self.python_script_dir = os.path.join(main_dir, 'python_scripts')        
        if not os.path.exists(self.python_script_dir):
            os.makedirs(self.python_script_dir)
         
        # Get path to pdb file directory and check if exists
        self.pdb_file_dir = os.path.join(main_dir, 'pdb_files')      
        if not os.path.exists(self.pdb_file_dir):
            os.makedirs(self.pdb_file_dir)     
         
        # Get path to csv to pdb file directory and check if exists
        self.csv_to_pdb_dir = os.path.join(main_dir, 'csvs_to_pdb')
        if not os.path.exists(self.csv_to_pdb_dir):
            os.makedirs(self.csv_to_pdb_dir)
         
        # Get path to residue code csv and check if exists
        self.residue_code_csv = os.path.join(self.pdb_file_dir, 'residue_codes.csv')       
        if not os.path.exists(self.residue_code_csv):
            with open(self.residue_code_csv, 'w') as file:
                pass # Don't want to do anything with it, just generate it
        else:
            pass # Pass as the file already exists
        
        # Get path to molecules directory and check if exists - this will contain the parameters and pdb for individual molecules
        self.molecules_dir = os.path.join(self.pdb_file_dir, 'molecules')     
        if not os.path.exists(self.molecules_dir):
            os.makedirs(self.molecules_dir)

        
        # Get path to systems directory and check if exists - this will contain the parameters and pdb for systems for md simulations
        self.systems_dir = os.path.join(self.pdb_file_dir, 'systems')   
        if not os.path.exists(self.systems_dir):
            os.makedirs(self.systems_dir)

    def mol2_files_avail(self):
        mol2_avail = []
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".mol2"):
                    # Construct the full path to the .pdb file
                    mol2_filepath = os.path.join(root, file)
                    mol2_avail.append(mol2_filepath)
                    # Extract molecule name
                    mol2_file = mol2_filepath.split("/")[-1]
                    #print(mol2_file)
        return(mol2_avail)

    def load_mol2_filepath(self, molecule_name=None):
        if molecule_name == None:
            print("Please provide the name of the system you are retrieving files as follows: 'mol2_file = directories.load_mol2_filepath('ethane')")
            print("Change ethane for the name of the desired system")
            print("NOTE: The mol2 file of the requested molecule must be generated with tleap prior to use of this function")
            return(None)
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".mol2"):
                    if (molecule_name + ".mol2") == file:
                        #print(file)
                        mol2_file_path = os.path.join(root, file)
                        return mol2_file_path 

    def load_pckml_filepath(self, system_name=None):
        if system_name == None:
            print("Please provide the name of the system you are retrieving files as follows: 'pckml_file = directories.load_pckml_filepath('pb_ph_41')")
            print("NOTE: Packmol input files must be written manually.")
            return(None)
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".pckml"):
                    if (system_name + ".pckml") == file:
                        pckml_file_path = os.path.join(root, file)
                        return pckml_file_path 
                        
    def pckml_files_avail(self):
        pckml_avail = []
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".pckml"):
                    # Construct the full path to the .pdb file
                    pckml_filepath = os.path.join(root, file)
                    pckml_avail.append(pckml_filepath)
                    # Extract molecule name
                    pckml_file = pckml_filepath.split("/")[-1]
                    #print(pckml_file)
        return(pckml_avail)

    def pdb_files_avail(self):
        pdb_avail = []
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".pdb"):
                    # Construct the full path to the .pdb file
                    pdb_filepath = os.path.join(root, file)
                    pdb_avail.append(pdb_filepath)
                    # Extract molecule name
                    pdb_file = pdb_filepath.split("/")[-1]
                    #print(pdb_file)
        return(pdb_avail)

    def load_pdb_filepath(self, molecule_name=None):
        if molecule_name == None:
            print("Please provide the name of the system you are retrieving files as follows: 'pdb_file = directories.load_pdb_filepath('ethane')")
            print("Change ethane for the name of the desired system")
            print("NOTE: If requesting a system for molecular dynamics - PDB files of a system must be generated using tleap prior to this step")
            return(None)

         
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                #if file.endswith(".pdb"):
                if file == (molecule_name + ".pdb"):
                    #if molecule_name in file:
                    pdb_file_path = os.path.join(root, file)
                    return pdb_file_path

        
                    
    def ac_files_avail(self):
        ac_avail = []
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.pdb_file_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .ac extension
                if file.endswith(".ac"):
                    # Construct the full path to the .ac file
                    ac_filepath = os.path.join(root, file)
                    ac_avail.append(ac_filepath)
                    # Extract molecule name
                    ac_file = ac_filepath.split("/")[-1]
                    #print(ac_file)
        return(ac_avail)
                
            
    def amber_systems_avail(self):
        a = False
        amber_system_avail =[]
        for root, dirs, files in os.walk(self.systems_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".prmtop"):
                    a = True
                    # Construct the full path to the .pdb file
                    prmtop_filepath = os.path.join(root, file)
                    # Extract molecule name
                    prmtop_file = prmtop_filepath.split("/")[-1]
                    amber_system_avail.append(prmtop_file)
                    #print(pdb_file)
                if file.endswith(".rst7"):
                    a = True
                    # Construct the full path to the .pdb file
                    rst7_filepath = os.path.join(root, file)
                    amber_system_avail.append(rst7_filepath)
                    # Extract molecule name
                    rst7_file = rst7_filepath.split("/")[-1]
                    #print(pdb_file)
        if a == True:
            print("")
            print("Remember you need both .prmtop and .rst7 files to run a simulation")    
            return(amber_system_avail)
        if a == False:
            print("No parametrized molecules.")
            return(None)

    def load_amber_filepaths(self, system_name=None):
        prmtop_file_path = None
        coord_file_path = None
        if system_name == None:
            print("Please provide the name of the system you are retrieving files as follows: 'topology_file, coordinate_file = directories.retrieve_top_crds('ethane')")
            print("Change ethane for the name of the desired system")
            print("NOTE: Amber files must be generated using tleap prior to this step")
            return(None)
        for root, dirs, files in os.walk(self.systems_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                #if file.endswith(".prmtop") and system_name in file: 
                if file == (system_name + ".prmtop"):
                    # Construct the full path to the .pdb file
                    prmtop_file_path = os.path.join(root, file)
                #if file.endswith(".rst7") and system_name in file:
                if file == (system_name + ".rst7") or file == (system_name + ".inpcrd"):
                    # Construct the full path to the .pdb file
                    coord_file_path = os.path.join(root, file)
        if (prmtop_file_path is not None) and (coord_file_path is not None):
            return(prmtop_file_path, coord_file_path)
        else:
            print("Files not found. Check name of molecule/system and if files have been generated.")
            return(None)

    def retrieve_files_for_MDanalysis(self, system_name=None):
        prmtop_file_path = None
        if system_name == None:
            print("Please provide the name of the system you are retrieving files for.")
            return(None)
        for root, dirs, files in os.walk(self.systems_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            for file in files:
                if file == (system_name + ".prmtop"):
                    prmtop_file_path = os.path.join(root, file)
                    
        folders = os.listdir(self.systems_dir)
        anneal_files, equili_files, prod_files = [], [], []
        for folder in folders:
            if folder == system_name:
                print(folder)
                folder_path = os.path.join(self.systems_dir, folder)
                folder_contents = os.listdir(folder_path)
                for item in folder_contents:
                        item_path = os.path.join(folder_path, item)
                        if os.path.isdir(item_path):
                            output_contents = os.listdir(item_path)
                            for file in output_contents:
                                print(item_path)
                                if ".dcd" in file:
                                    if "anneal" in file:
                                        anneal_files.append(os.path.join(item_path, file))
                                    if "atm" in file:
                                        equili_files.append(os.path.join(item_path, file))
                                    if "prod" in file:
                                        prod_files.append(os.path.join(item_path, file))
                                if ".pdb" in file:
                                    if "anneal" in file:
                                        anneal_files.append(os.path.join(item_path, file))
                                    if "atm" in file:
                                        equili_files.append(os.path.join(item_path, file))
                                    if "prod" in file:
                                        prod_files.append(os.path.join(item_path, file))
            
        
        return(prmtop_file_path, anneal_files, equili_files, prod_files)
        
    
    def bash_submission(self):
        pass
        # Need an if instance for whether its amber or ani
    
    @staticmethod  # Dont take instance or class - can run this without an instance   
    def unpack_csv(csv_with_mol_info):
        names = []
        smiles = []

        with open(csv_with_mol_info, 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                if len(row) >= 2:
                    names.append(row[0])
                    smiles.append(row[1])

        return names, smiles
        
    def retrieve_polymeric_rescodes(self, molecule_name):
        with open(self.residue_code_csv, 'r') as file:
            head_code, mainchain_code, tail_code = None, None, None
            for line in file:
                parts = (line.strip().split('\t'))[0]
                if len(parts.split(",")) >= 3 and parts.split(",")[0] == ("head_" + molecule_name):
                    head_code = parts.split(",")[2]
                if len(parts.split(",")) >= 3 and parts.split(",")[0] == ("mainchain_" + molecule_name):
                    mainchain_code = parts.split(",")[2]
                if len(parts.split(",")) >= 3 and parts.split(",")[0] == ("tail_" + molecule_name):
                    tail_code = parts.split(",")[2]    
        return(head_code, mainchain_code, tail_code)


    def retrieve_rescode(self, molecule_name):
        with open(self.residue_code_csv, 'r') as file:
            rescode = None
            for line in file:
                parts = (line.strip().split('\t'))[0]
                if len(parts.split(",")) >= 3 and parts.split(",")[0] == (molecule_name):
                    rescode = parts.split(",")[2]
        return(rescode)


        
    # Need to add stuff for md results and the like
    # Need to add a way to import the bash submission scripts - will get to this at some point

        if not os.path.exists(main_dir):
            raise FileNotFoundError(f"The specified main directory '{main_dir}' does not exist.")       
        self.main_dir = main_dir 
        
        # Get path to python script directory and check if exists
        self.python_script_dir = os.path.join(main_dir, 'python_scripts')        
        if not os.path.exists(self.python_script_dir):
            os.makedirs(self.python_script_dir)
         
        # Get path to pdb file directory and check if exists
        self.pdb_file_dir = os.path.join(main_dir, 'pdb_files')      
        if not os.path.exists(self.pdb_file_dir):
            os.makedirs(self.pdb_file_dir)     
         
        # Get path to csv to pdb file directory and check if exists
        self.csv_to_pdb_dir = os.path.join(main_dir, 'csvs_to_pdb')
        if not os.path.exists(self.csv_to_pdb_dir):
            os.makedirs(self.csv_to_pdb_dir)
         
        # Get path to residue code csv and check if exists
        self.residue_code_csv = os.path.join(self.pdb_file_dir, 'residue_codes.csv')       
        if not os.path.exists(self.residue_code_csv):
            with open(self.residue_code_csv, 'w') as file:
                pass # Don't want to do anything with it, just generate it
        else:
            pass # Pass as the file already exists
        
        # Get path to molecules directory and check if exists - this will contain the parameters and pdb for individual molecules
        self.molecules_dir = os.path.join(self.pdb_file_dir, 'molecules')     
        if not os.path.exists(self.molecules_dir):
            os.makedirs(self.molecules_dir)

        
        # Get path to systems directory and check if exists - this will contain the parameters and pdb for systems for md simulations
        self.systems_dir = os.path.join(self.pdb_file_dir, 'systems')   
        if not os.path.exists(self.systems_dir):
            os.makedirs(self.systems_dir)

    def simulations_avail(self, system_name):
        simulation_dir = os.path.join(self.systems_dir, system_name)
        avail_sims = []
        for item in os.listdir(simulation_dir):
            item_path = os.path.join(simulation_dir, item)
            if os.path.isdir(item_path):
                avail_sims.append(item_path)      
        if avail_sims == []:
            print("No simulations found for the system.")
            print("")
            print("Please ensure simulation files are avaialable and system name is correct.")
            return(None)
        else:
            print("Output contains paths to simulation directories.")
            return(avail_sims)

class BioOilDirs(SnippetSimManage):

    def __init__(self, main_dir, *args, **kwargs):
        """
        Initialize PolymerSimulatorDirs object.

        Args:
            main_dir (str): The main directory for polymer simulation setup.
                            Example: '/path/to/main/dir/'

        Raises:
            FileNotFoundError: If main_dir does not exist.
        """
        # Call the parent class's __init__ method with all arguments
        super().__init__(main_dir, *args, **kwargs)
        
        if not os.path.exists(main_dir):
            raise FileNotFoundError

        self.bio_oil_dir = os.path.join(main_dir, 'bio_oil')
        if not os.path.exists(self.bio_oil_dir):
            os.makedirs(self.bio_oil_dir)

        self.bio_oil_systems_dir = os.path.join(self.bio_oil_dir, 'systems')
        if not os.path.exists(self.bio_oil_systems_dir):
            os.makedirs(self.bio_oil_systems_dir)

        self.bio_oil_GC_data = os.path.join(self.bio_oil_dir, 'GC_data')
        if not os.path.exists(self.bio_oil_GC_data):
            os.makedirs(self.bio_oil_GC_data)

        self.bio_oil_models_dir = os.path.join(self.bio_oil_dir, 'models')
        if not os.path.exists(self.bio_oil_models_dir):
            os.makedirs(self.bio_oil_models_dir)

        self.packmol_inputs = os.path.join(self.bio_oil_dir, 'packmol_inputs')
        if not os.path.exists(self.packmol_inputs):
            os.makedirs(self.packmol_inputs)

        self.packmol_systems = os.path.join(self.bio_oil_dir, 'packmol_systems')
        if not os.path.exists(self.packmol_systems):
            os.makedirs(self.packmol_systems)

    def GC_data_avail(self):
        import pandas as pd
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.bio_oil_GC_data):
            # Exclude 'depreceated' directories
            dirs[:] = [d for d in dirs if d != 'depreceated']
        
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .csv extension
                if file.endswith(".csv"):
                    # Construct the full path to the .csv file
                    csv_file_path = os.path.join(root, file)
                    # Extract filename
                    csv_file = os.path.basename(csv_file_path)

                    try:
                        # Read the CSV file
                        df = pd.read_csv(csv_file_path)

                        # Drop any columns without headers (unnamed columns)
                        df = df.loc[:, ~df.columns.str.startswith('Unnamed')]

                        # Check for null values in the remaining columns
                        if df.isnull().any().any():  # Check for any NaN in valid columns
                            print(f"{csv_file} - Incomplete")
                        else:
                            print(csv_file)
                    except Exception as e:
                        print(f"Error reading {csv_file}: {e}")

    def bio_oil_pckml_inputs_avail(self):
        avail_files = []
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.bio_oil_systems_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".pdckml") or file.endswith(".inp"):
                    # Construct the full path to the .pdb file
                    pdb_file_path = os.path.join(root, file)
                    # Extract molecule name
                    pdb_file = pdb_file_path.split("/")[-1]
                    avail_files.append(pdb_file_path)
                    print(pdb_file)
        return(avail_files)

    def bio_oil_pckml_systems_avail(self):
        avail_files = []
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.packmol_systems):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".pdb"):
                    # Construct the full path to the .pdb file
                    pdb_file_path = os.path.join(root, file)
                    # Extract molecule name
                    pdb_file = pdb_file_path.split("/")[-1]
                    avail_files.append(pdb_file_path)
                    print(pdb_file)
        return(avail_files)

    def bio_oil_amber_systems_avail(self):
        a = False
        amber_system_avail =[]
        for root, dirs, files in os.walk(self.bio_oil_systems_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".prmtop"):
                    a = True
                    # Construct the full path to the .pdb file
                    prmtop_filepath = os.path.join(root, file)
                    # Extract molecule name
                    prmtop_file = prmtop_filepath.split("/")[-1]
                    amber_system_avail.append(prmtop_file)
                    #print(pdb_file)
                if file.endswith(".rst7"):
                    a = True
                    # Construct the full path to the .pdb file
                    rst7_filepath = os.path.join(root, file)
                    amber_system_avail.append(rst7_filepath)
                    # Extract molecule name
                    rst7_file = rst7_filepath.split("/")[-1]
                    #print(pdb_file)
        if a == True:
            print("")
            print("Remember you need both .prmtop and .rst7 files to run a simulation")    
            return(amber_system_avail)
        if a == False:
            print("No parametrized molecules.")
            return(None)

    def load_bio_oil_amber_filepaths(self, system_name=None):
        prmtop_file_path = None
        coord_file_path = None
        if system_name == None:
            print("Please provide the name of the system you are retrieving files as follows: 'topology_file, coordinate_file = directories.retrieve_top_crds('ethane')")
            print("Change ethane for the name of the desired system")
            print("NOTE: Amber files must be generated using tleap prior to this step")
            return(None)
        for root, dirs, files in os.walk(self.bio_oil_systems_dir):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                #if file.endswith(".prmtop") and system_name in file: 
                if file == (system_name + ".prmtop"):
                    # Construct the full path to the .pdb file
                    prmtop_file_path = os.path.join(root, file)
                #if file.endswith(".rst7") and system_name in file:
                if file == (system_name + ".rst7"):
                    # Construct the full path to the .pdb file
                    coord_file_path = os.path.join(root, file)
        if (prmtop_file_path is not None) and (coord_file_path is not None):
            return(prmtop_file_path, coord_file_path)
        else:
            print("Files not found. Check name of molecule/system and if files have been generated.")
            return(None)


# This class should be inititated when you start handing specific data for a model
# i.e. there will be an instance of this specific to each model
class complex_model_dirs(BioOilDirs):
    def __init__(self, main_dir, model_name, *args, **kwargs):
        super().__init__(main_dir, *args, **kwargs)
        
        self.complex_model_dir = os.path.join(self.bio_oil_dir, model_name)
        if not os.path.exists(self.complex_model_dir):
            os.makedirs(self.complex_model_dir)

        self.dft_input_dir = os.path.join(self.complex_model_dir, "DFT_inputs")
        if not os.path.exists(self.dft_input_dir):
            os.makedirs(self.dft_input_dir)

        self.packmol_inputs = os.path.join(self.complex_model_dir, "packmol_inputs")
        if not os.path.exists(self.packmol_inputs):
            os.makedirs(self.packmol_inputs)

        self.packmol_systems = os.path.join(self.complex_model_dir, "packmol_systems")
        if not os.path.exists(self.packmol_systems):
            os.makedirs(self.packmol_systems)

        self.output_files = os.path.join(self.complex_model_dir, "output_files")
        if not os.path.exists(self.output_files):
            os.makedirs(self.output_files)

    def packmol_systems_avail(self):
        avail_files = []
        # Walk through the directory tree recursively
        for root, dirs, files in os.walk(self.packmol_systems):
            dirs[:] = [d for d in dirs if d != 'depreceated']
            # Check each file in the current directory
            for file in files:
                # Check if the file has a .pdb extension
                if file.endswith(".pdb"):
                    # Construct the full path to the .pdb file
                    pdb_file_path = os.path.join(root, file)
                    # Extract molecule name
                    pdb_file = pdb_file_path.split("/")[-1]
                    avail_files.append(pdb_file_path)
                    #print(pdb_file)
        return(avail_files)


