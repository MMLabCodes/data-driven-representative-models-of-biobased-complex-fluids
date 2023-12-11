# Data driven representative models of biobased complex fluids

An open-source code that generates a range of data driven models to represent complex biobased fluids utilising data that 
quantifies molecular components within a system.

# Info

Authors: Daniel York, Isaac Vidal-Daza, Francisco J. Martin-Martinez 

Swansea University, 2022 - 2023

Email correspondence: 983045@swansea.ac.uk

# Publications

York D., Vidal-Daza I., Segura C., _et al._ Data-driven representative models to accelerate scaled-up atomistic simulations of bitumen and biobased complex fluids. _Submitted to Digital discovery._

# Requirements

An anaconda python environment is recommended 


# Usage

The script to generate data driven representative models of complex organic fluids is called;
	"data_driven_model_generator.py"

There are 2 versions:

	1) data_driven_model_generator.py
		This is the interactive script for use in python IDE programmes such as spyder or vscode.

	2) data_driven_model_generator_cmd_line.py
		This is the command line version of the script

Before running these scripts;
	1) The functions folder from this associated github must be in the same location as the above script/scripts.
	2) A csv file containing the following information must be available in the format below;
		Note: 'area' is the peak area coming from GCMS data and could be substituted for any quantifiable data 	.
		Note: A lot of information in the csv is obtained from DFT calculations and should be collated prior to executing these scripts.		

		Molecule, mw, area, smiles, formula, Total electronic energy (eV), HOMO (eV), LUMO(ev), Chemical hardness, Dipole moment, Polarizability

Running these scripts;
	
	Interactive version:
		1) Load script into python IDE and change the following 2 variables
			
			csv_file = "path/to/csv/file" # This should be the path to the csv file in the format specified above.
			master_path = "path/to/master_path" # This will be where results files for the different models are generated. RECOMMENDED: specify this to be the directory of where the script is being run from, a new directory will be generated for each complex organic mixture.

		
		2) Run script in python IDE

	Command line version:
		1) Move to directory where this script is located
		
		2) Execute as follows;
			python3 data_driven_model_generator_cmd_line.py <csv_file> <master_path>
			
			The variables should be defined as:
				<csv_file> = Path to csv file in the format specified above.
				<master_path> = This will be where results files for the different models are generated. RECOMMENDED: specify this to be the directory of where the script is being run from, a new directory will be generated for each complex organic mixture.
