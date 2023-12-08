# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 15:08:58 2023

@author: danie
"""
import os
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
'''
This class makes filepath management easier in the program and is used in the next function
'''     
def make_model_dirs(main_dir, model_name):
    model_dir = main_dir + "/" + model_name
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    subdirs = [(model_name  + "_DFT"), ((model_name  + "_DFT") + "/fukui"), ((model_name  + "_DFT") + "/results"), ((model_name  + "_DFT") + "/xyz")]
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

