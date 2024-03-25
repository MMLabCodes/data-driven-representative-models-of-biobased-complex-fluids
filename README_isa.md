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

Creating the Test directory and cloning the repository:

```shell
$ mkdir Test && cd Test
$ git clone https://github.com/MMLabCodes/data-driven-representative-models-of-biobased-complex-fluids.git
```

### Installation of dependencies and packages

```shell
(biobased_venv) $ cd data-driven-representative-models-of-biobased-complex-fluids
(biobased_venv) $ python -m pip install -r requirements.txt 
```

## Usage

The script `data_driven_model_generator.py` generates data-driven representative models of complex organic fluids. There are two versions of the script:

1. `data_driven_model_generator.py`: Interactive version for use in Python IDEs like Spyder or VSCode.
2. `data_driven_model_generator_cmd_line.py`: Command line version of the script.

### Preparation Before Running the Scripts

Before running these scripts, ensure the following:

1. The functions folder from this associated [GitHub repository](https://github.com/dyork1/data_driven_model_generator) must be in the same location as the scripts.

2. A CSV file with the following information must be available:
   
   - **Molecule**: The name of the molecule.
   - **mw**: Molecular weight of the molecule, typically measured in atomic mass units (amu).
   - **area**: The peak area from GCMS data. This could be substituted for any quantifiable data.
   - **smiles**: Simplified molecular-input line-entry system (SMILES), a specification in the form of a line notation for describing the structure of chemical species.
   - **formula**: The chemical formula of the molecule.
   - **Total electronic energy (eV)**: The total energy of the electrons in the molecule, typically measured in electron volts (eV).
   - **HOMO (eV)**: Energy of the highest occupied molecular orbital, typically measured in electron volts (eV).
   - **LUMO (eV)**: Energy of the lowest unoccupied molecular orbital, typically measured in electron volts (eV).
   - **Chemical hardness**: A measure of the stability of a chemical species. It is defined as the average energy gap between the highest occupied molecular orbital (HOMO) and the lowest unoccupied molecular orbital (LUMO).
   - **Dipole moment**: A measure of the polarity of a molecule. It is the product of the charge and the distance between the charges.
   - **Polarizability**: A measure of how much the electron cloud around an atom or molecule can be distorted by an external electric field.
  
> *Much of the information in the CSV is obtained from DFT calculations and should be collated prior to running these scripts.*

The table header must be:

| Molecule | mw | area | smiles | formula | Total electronic energy (eV) | HOMO (eV) | LUMO(ev) | Chemical hardness | Dipole moment | Polarizability |
|----------|----|------|--------|---------|------------------------------|-----------|----------|-------------------|---------------|---------------|


### Running the Code

- Interactive Version:
	1. Load the script into your Python IDE and change the following two variables:
	    - csv_file = "path/to/csv/file" 
	    - master_path = "path/to/master_path" 
	2. Run the script in your Python IDE.

- Command Line Version:
	1. Navigate to the directory where this script is located.
	2. Execute the script as follows:
	   ```shell
  		$ python3 data_driven_model_generator_cmd_line.py <csv_file> <master_path>
    	```

The variables should be defined as:

- <csv_file>: Path to the CSV file in the format specified above.
- <master_path>: This will be where result files for the different models are generated. RECOMMENDED: specify this to be the directory from where the script is being run, a new directory will be generated for each complex organic mixture.

---

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