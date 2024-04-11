# Data driven representative models of biobased complex fluids :star:

An open-source code that generates a range of data driven models to represent complex biobased fluids from data that quantifies molecular components within a system.

## Info :information_source:

Authors: Daniel York, Isaac Vidal-Daza, Francisco J. Martin-Martinez 

Swansea University, 2022 - 2023

Email correspondence: [f.j.martin-martinez@swansea.ac.uk](mailto:f.j.martin-martinez@swansea.ac.uk)

## Publications :books:

York D., Vidal-Daza I., Segura C., _et al._ Data-driven representative models to accelerate scaled-up atomistic simulations of bitumen and biobased complex fluids. _Submitted to Digital discovery._

# Installation
## Installation of a virtual environment

We first create a virtual environment to avoid dealing with other packages installed in the system:

Test folder:
```shell
'$ mkdir Test && cd Test'
'$ git clone https://eur03.safelinks.protection.outlook.com/?url=https%3A%2F%2Fgithub.com%2FMMLabCodes%2Fdata-driven-representative-models-of-biobased-complex-fluids.git%2560&data=05%7C02%7CF.J.Martin-Martinez%40Swansea.ac.uk%7C919096b1d4cb426173af08dc4834d5fb%7Cbbcab52e9fbe43d6a2f39f66c43df268%7C0%7C0%7C638464639261410676%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C0%7C%7C%7C&sdata=OvOrJ8pG0DNPvXNc%2BMaM7NVOk7PybEykdl%2FQoQ2EdDU%3D&reserved=0'
```
We create a virtual environment:
```shell
'$ python -m venv biobased_venv'
```
To activate this environment, we use:
```shell
'$ source ./biobased_venv/bin/activate'
```
To deactivate the environment:
```shell
'$ deactivate'
```
## Install requirements :wrench:

An anaconda python environment is recommended and the following packeages are required (python 3.7 or later):

```shell
$ pip install -r requirements.txt
```
This will install pandas and rdkit.

![Python version](https://img.shields.io/badge/python-3.7+-blue)

## Usage :computer:

The script to generate data driven representative models of complex organic fluids is called `data_driven_model_generator.py`.
There are 2 versions:

1. `data_driven_model_generator.py`: This is the interactive script for use in python IDE programmes such as spyder or vscode.
2. `data_driven_model_generator_cmd_line.py`: This is the command line version of the script

Before running these scripts:

1. The functions folder from this associated [github](https://github.com/dyork1/data_driven_model_generator) must be in the same location as the above script/scripts.
2. A csv file containing the following information must be available in the format below:

| Molecule | mw | area | smiles | formula | Total electronic energy (eV) | HOMO (eV) | LUMO(ev) | Chemical hardness | Dipole moment | Polarizability |
|----------|----|------|--------|---------|------------------------------|-----------|----------|-------------------|---------------|---------------|

> Note 1: 'area' is the peak area coming from GCMS data and could be substituted for any quantifiable data.
>
> Note 2: A lot of information in the csv is obtained from DFT calculations and should be collated prior to executing these scripts.

### Running the scripts
- Interactive version:
	1. Load script into python IDE and change the following 2 variables:
	    - csv_file = "path/to/csv/file" # This should be the path to the csv file in the format specified above.
	    - master_path = "path/to/master_path" # This will be where results files for the different models are generated. RECOMMENDED: specify this to be the directory of where the script is being run from, a new directory will be generated for each complex organic mixture.
	2. Run script in python IDE
- Command line version:
	1. Move to directory where this script is located
	2. Execute as follows:
	   ```shell
  		$ python3 data_driven_model_generator_cmd_line.py <csv_file> <master_path>
    	```
	   The variables should be defined as:
	    - <csv_file> = Path to csv file in the format specified above.
	    - <master_path> = This will be where results files for the different models are generated. RECOMMENDED: specify this to be the directory of where the script is being run from, a new directory will be generated for each complex organic mixture.
### Output data
- The output files for each model will be located in directory generated when the script runs.
- An output file for each model (i.e. all molecules, FT, PT, AG and SG) will be available.
- A "_results.csv" file is generated to compile all data coming from the models.
    - Currently this is empty but all results can be found from "...out" files.
    