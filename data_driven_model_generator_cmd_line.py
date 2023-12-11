# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 11:18:20 2023

@author: danie
"""
import sys
from functions.class_definitions import *
from functions.model_functions import *
'''
This generates data driven representative models from the final csv containing DFT data
'''
'''
Part one:
    
    1) Obtain csv_file and master_path paths from command line arguments
    
    2) Extract model name from csv file
    
    3) Unpack csv_file into molecule objects using csv_to_orca_class
'''
csv_file = sys.argv[1]
master_path = sys.argv[2]

model_name = (os.path.basename(csv_file)).split(".")[0]
model_name = model_name.replace("_DFT_results", "")
model_dirs = make_model_dirs(master_path, model_name)
output_directory = model_dirs.model_dir

molecules = csv_to_orca_class(csv_file)
'''
Part two:
    
    1) Generate benchmark model
    2) Generate FT model with 5% selection threshold
        Note: Changing the number will alter the selection threshold
    3) Generate PT model
    4) Generate AG model
    5) Generate SG model
        Note: The SG model can only be generated from the AG model
'''
all_molecules_model = all_model(molecules)
FT_model = FT_model(molecules, 5)
PT_model = PT_model(molecules)
AG_model = AG_model(molecules)
SG_model = SG_model(AG_model)
'''
Part three: 
    
    1) Put models in a list
    2) Enter model types:
        Note: This version uses old names, model_types = ["all_model", "five_percent_model", "propr_rep", "group_model", "score_model"]
            These function the same 
    3) Rank the models
    
'''
models = [all_molecules_model, FT_model, PT_model, AG_model, SG_model]
model_types = ["all_model", "five_percent_model", "propr_rep", "group_model", "score_model"]
score_df = generate_model_df(molecules, models, model_types) # get data as pandas datafrmae
ranks = rank_models(score_df) # creates a ranked dataframe
'''
Part four:
    
    1) Write model output files
        1.1) 
'''
model_names = ["all_model", "FT_model", "PT_model", "AG_model", "SG_model"]
for i in range(len(models)):
        output_filename = output_directory + "/" + model_names[i] + ".out"
        write_output(output_filename, model_output_block(molecules, models[i], model_types[i], ranks))
   

