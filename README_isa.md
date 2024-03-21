# Data driven representative models of biobased complex fluids :star:

## Description

An open-source code that generates a range of data driven models to represent complex biobased fluids from data that quantifies molecular components within a system.

**Authors**: Daniel York, Isaac Vidal-Daza, Francisco J. Martin-Martinez *Swansea University*, 2022-2023

## Publications :books:

York D., Vidal-Daza I., Segura C., _et al._ Data-driven representative models to accelerate scaled-up atomistic simulations of bitumen and biobased complex fluids. _Submitted to Digital discovery._

## Installation
Instructions for setting up your project, including required software and dependencies.

### Installation of virtual environment

We first create a virtual environment to avoid dealing with other packages installed in the system:

```shell
$ mkdir Test && cd Test
$ git clone https://github.com/MMLabCodes/data-driven-representative-models-of-biobased-complex-fluids.git
```

A virtual environment is created: 

```shell
$ python -m venv biobased_venv
``` 

To activate this environment, use:

```shell
$ source ./biobased_venv/bin/activate
``` 

To deactivate the environment: 

```shell
$ deactivate
```

### Installation of dependencies and packages

```shell
(biobased_venv) $ python -m pip install -r requirements.txt 
```

## Usage

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

### Running the code

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

## Data

(Zenodo link to data)

## Results
Summarize the results of your data analysis or model. This could include charts or other visuals.

## Contributing
If this is a collaborative project, explain how others can contribute. Include any rules or guidelines you want collaborators to follow.

## License
Information about the license (if any).

## Contact

Email correspondence: [f.j.martin-martinez@swansea.ac.uk](mailto:f.j.martin-martinez@swansea.ac.uk)