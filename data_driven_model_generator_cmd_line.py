# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 11:18:20 2023

@author: danie
"""
import sys
from functions.class_definitions import *
from functions.model_functions import *
"""
This script generates data-driven representative models from a final CSV containing DFT data.
"""
def main():
    """
    Specify filepaths for the csv_file and master_path in the command line when executing this script.
    
    csv_file = csv_file containing DFT data in the format given below:
        
        Molecule, mw, area, smiles, formula, Total electronic energy (eV), HOMO (eV), LUMO(ev), Chemical hardness, Dipole moment
        
        Note: many of the information in this csv_file comes from DFT, such as: Total electronic energy (eV), HOMO (eV), LUMO(ev), Dipole moment
        
    master_path = main directory where model information will be generated
    """
    csv_file = sys.argv[1]
    master_path = sys.argv[2]
    """
    Extract model name from the csv_file 
    
    Note: Our example datafiles contain the 'model name' and 'DFT_results'
        We remove this 'DFT_results' string to get the basename of the csv we are modelling; this may not be necessary in some cases. 
    """
    model_name = os.path.basename(csv_file).split(".")[0]
    if "_DFT_results" in model_name:
        model_name = model_name.replace("_DFT_results", "")
    """
    Generate directories for the models by calling the class 'model_dirs' class
    Return the output directory where we can send our outputs
    
    Note: The class "model_dirs" can be editted to aid with other filepath management during your own DFT calculations
    """
    directory_manager = model_dirs(master_path, model_name)
    directory_manager.gen_model_dirs()
    output_directory = directory_manager.model_dir
    """
    Extract the data into 'orca_molecule' objects to allow for model generation
    """
    molecules = csv_to_orca_class(csv_file)
    """
    Generate all required models:
        1) Benchmark model (containing all molecules)
        2) FT model with 5% selection threshold
            Note: The selection threshold can altered
        3) PT model
        4) AG model
        5) SG model
            Note: The SG model is generated from the AG model
    """
    all_molecules_model = gen_all_model(molecules)
    FT_model = gen_FT_model(molecules, 5)
    PT_model = gen_PT_model(molecules)
    AG_model = gen_AG_model(molecules)
    SG_model = gen_SG_model(AG_model)
    """
    Place models and model types in a list so they can be ranked
    
    Note: A new ranking system will be implemented in the future using methods such as AHP analysis.
        This current ranking system provides a general performance metric, but is only based on rankings against the benchmark model.
   
    Generate a dataframe for model properties and rank them
    """
    models = [all_molecules_model, FT_model, PT_model, AG_model, SG_model]
    model_types = ["all_model", "five_percent_model", "propr_rep", "group_model", "score_model"]
    score_df = generate_model_df(molecules, models, model_types) 
    ranks = rank_models(score_df)
    """
    Write output files for each model which is structured as below (AG model output as an example):
        
        Model: group_model
        This model represents 18.397093713247152% of the charecterised biocrude
        This model has a score of: 4.5 versus the benchmark

        Weighted values for molecular weight, chemical hardness, polarizability, chemical hardness and dipole moment are listed below:
        Molecular weight  = 210.42959924776352
        Chemical hardness = 4.232922624242967
        Polarizability    = 165.15302923925347
        Dipole moment     = 1.5091704133772004
        Total energy      = -17157.37626977363

        The heteroatom content of this model is:
        B: 0.0
        C: 79.66551336601214
        N: 0.0
        O: 8.36950139125679
        P: 0.0
        S: 0.0
        F: 0.0
        Cl: 0.0
        Br: 0.0
        I: 0.0
        H: 11.964985242731068

        This model will require the following minimum number of atoms for unbiased simulation:
        11816

        This model will require the following minimum number of molecules for unbiased simulation:
        299

        The molecules in this model are:
        Furfural and this group accounts for: 1.8847800000000001 % of the biocrude
        2-methylPhenol and this group accounts for: 53.52037999999999 % of the biocrude
        Stigmast-5-ene3-methoxy-3 and this group accounts for: 1.87152 % of the biocrude
        Tetracosanoicacid and this group accounts for: 22.116580000000003 % of the biocrude
        1-Tetracosene and this group accounts for: 14.11988 % of the biocrude

        The smiles for these molecules are:
        C1=COC(=C1)C=O
        CC1=CC=CC=C1O
        CCC(CCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)OC)C)C)C(C)C
        CCCCCCCCCCCCCCCCCCCCCCCC(=O)O
        CCCCCCCCCCCCCCCCCCCCCCC=C
        
    The output files are sent to the filepath specified in 'output_directory'
    """
    model_names = ["all_model", "FT_model", "PT_model", "AG_model", "SG_model"]
    for i in range(len(models)):
        output_filename = output_directory + "/" + model_names[i] + ".out"
        write_output(output_filename, model_output_block(molecules, models[i], model_types[i], ranks))
        
if __name__ == "__main__":
    main()"""
    This script generates data-driven representative models from a final CSV containing DFT data.
    """
    def main():
        """
        Specify filepaths for the csv_file and master_path
        
        csv_file = csv_file containing DFT data in the format given below:
            
            Molecule, mw, area, smiles, formula, Total electronic energy (eV), HOMO (eV), LUMO(ev), Chemical hardness, Dipole moment
            
            Note: many of the information in this csv_file comes from DFT, such as: Total electronic energy (eV), HOMO (eV), LUMO(ev), Dipole moment
            
        master_path = main directory where model information will be generated
        """
        csv_file = "C:/Users/danie/Swansea University/Francisco Martin-Martinez - York_Dan/05_Others/97_bio_oil_test/data-driven-representative-models-of-biobased-complex-fluids/pb_cp_DFT_results.csv"
        master_path = "C:/Users/danie/Swansea University/Francisco Martin-Martinez - York_Dan/05_Others/97_bio_oil_test/data-driven-representative-models-of-biobased-complex-fluids"
        """
        Extract model name from the csv_file 
        
        Note: Our example datafiles contain the 'model name' and 'DFT_results'
            We remove this 'DFT_results' string to get the basename of the csv we are modelling; this may not be necessary in some cases. 
        """
        model_name = os.path.basename(csv_file).split(".")[0]
        if "_DFT_results" in model_name:
            model_name = model_name.replace("_DFT_results", "")
        """
        Generate directories for the models by calling the class 'model_dirs' class
        Return the output directory where we can send our outputs
        
        Note: The class "model_dirs" can be editted to aid with other filepath management during your own DFT calculations
        """
        directory_manager = model_dirs(master_path, model_name)
        directory_manager.gen_model_dirs()
        output_directory = directory_manager.model_dir
        """
        Extract the data into 'orca_molecule' objects to allow for model generation
        """
        molecules = csv_to_orca_class(csv_file)
        """
        Generate all required models:
            1) Benchmark model (containing all molecules)
            2) FT model with 5% selection threshold
                Note: The selection threshold can altered
            3) PT model
            4) AG model
            5) SG model
                Note: The SG model is generated from the AG model
        """
        all_molecules_model = gen_all_model(molecules)
        FT_model = gen_FT_model(molecules, 5)
        PT_model = gen_PT_model(molecules)
        AG_model = gen_AG_model(molecules)
        SG_model = gen_SG_model(AG_model)
        """
        Place models and model types in a list so they can be ranked
        
        Note: A new ranking system will be implemented in the future using methods such as AHP analysis.
            This current ranking system provides a general performance metric, but is only based on rankings against the benchmark model.
       
        Generate a dataframe for model properties and rank them
        """
        models = [all_molecules_model, FT_model, PT_model, AG_model, SG_model]
        model_types = ["all_model", "five_percent_model", "propr_rep", "group_model", "score_model"]
        score_df = generate_model_df(molecules, models, model_types) 
        ranks = rank_models(score_df)
        """
        Write output files for each model which is structured as below (AG model output as an example):
            
            Model: group_model
            This model represents 18.397093713247152% of the charecterised biocrude
            This model has a score of: 4.5 versus the benchmark

            Weighted values for molecular weight, chemical hardness, polarizability, chemical hardness and dipole moment are listed below:
            Molecular weight  = 210.42959924776352
            Chemical hardness = 4.232922624242967
            Polarizability    = 165.15302923925347
            Dipole moment     = 1.5091704133772004
            Total energy      = -17157.37626977363

            The heteroatom content of this model is:
            B: 0.0
            C: 79.66551336601214
            N: 0.0
            O: 8.36950139125679
            P: 0.0
            S: 0.0
            F: 0.0
            Cl: 0.0
            Br: 0.0
            I: 0.0
            H: 11.964985242731068

            This model will require the following minimum number of atoms for unbiased simulation:
            11816

            This model will require the following minimum number of molecules for unbiased simulation:
            299

            The molecules in this model are:
            Furfural and this group accounts for: 1.8847800000000001 % of the biocrude
            2-methylPhenol and this group accounts for: 53.52037999999999 % of the biocrude
            Stigmast-5-ene3-methoxy-3 and this group accounts for: 1.87152 % of the biocrude
            Tetracosanoicacid and this group accounts for: 22.116580000000003 % of the biocrude
            1-Tetracosene and this group accounts for: 14.11988 % of the biocrude

            The smiles for these molecules are:
            C1=COC(=C1)C=O
            CC1=CC=CC=C1O
            CCC(CCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)OC)C)C)C(C)C
            CCCCCCCCCCCCCCCCCCCCCCCC(=O)O
            CCCCCCCCCCCCCCCCCCCCCCC=C
            
        The output files are sent to the filepath specified in 'output_directory'
        """
        model_names = ["all_model", "FT_model", "PT_model", "AG_model", "SG_model"]
        for i in range(len(models)):
            output_filename = output_directory + "/" + model_names[i] + ".out"
            write_output(output_filename, model_output_block(molecules, models[i], model_types[i], ranks))
            
    if __name__ == "__main__":
        main()