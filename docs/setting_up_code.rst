⚙️ Project Setup Guide
=======================

This guide goes through some prerequisites before any data driven representative models can be generated:

- Importing packages
- Intialising a filepath manager
- Data structure (.csv file)
- Extracting data from csv file


⬇️ Importing python modules
==========================

📁 **modules.sw_directories**
-----------------------------

This module functions as a **file path manager and organizer** 🗂️.

It ensures that:
- 📌 All saved files are directed to appropriate folders based on the assigned **model name**.
- 🗃️ Project files remain clean, well-structured, and easy to navigate.

By automating directory creation and file sorting, it helps maintain an organized workspace throughout your modeling workflow.

.. code-block:: python

   # Import filepath manager
   from modules.sw_directories import *

⚛️ **modules.sw_orca**
----------------------

This module contains one core class that:
- 🔄 Converts raw **DFT data** (from your input `.csv`)  
- 🧱 Transforms each molecule’s data into a structured **Python object**  

This greatly simplifies handling, accessing, and integrating quantum chemistry results throughout the modeling process.

.. code-block:: python

   # Import filepath manager
   from modules.sw_orca import *

🧪 **modules.sw_complex_fluid_models**
--------------------------------------

This is the **core module** of the project 🔧.

It processes groups of molecules from your raw **DFT dataset** and generates a range of **complex fluid models** 🧬.  
The primary output is an instance of the `complex_fluid_model` class, which includes several attributes — the most important being the **group of selected molecules** (more on this later).

🧠 Model Types Generated:

- 🧱 **FT_model** — Fixed Threshold Model  
- ⚖️ **PT_model** — Proportional Threshold Model  
- 🌐 **AG_model** — Abundancy Grouping Model  
- 🎯 **SG_model** — Scored Grouping Model  

.. code-block:: python

   # Import filepath manager
   from modules.sw_complex_fluid_models import *

📁 Setting up filepath manager
==============================

Now all of the modules are imported first step is to set up the **file system** structure using the command below.  
This will automatically generate the necessary folders to organize your raw and processed data.

📌 *Note:* You must always run the initialization command — but folders will only be created the **first time** the `manager` object is instantiated.

🗃️ Once initialized, a target folder will be available for you to **copy your raw data files into** (details and examples below).

.. code-block:: python

    # Import os 
    import os as os

    manager = BioOilDirs(os.getcwd())

**BioOilDirs** is a class that contains different filepaths as variables and creates a filetree that looks like this:

.. image:: _static/file_tree.png
   :alt: Directory structure for project
   :align: center
   :width: 600px

🎯 **Key Areas Highlighted**

- 🔴 **Red**: The main project directory — models are stored here.  
  ⤷ Place all **raw data** files in the `GC_data/` folder.

- 🟡 **Yellow**: The `modules/` directory — contains all core source code.

- 🟢 **Green**: The `molecules/` directory — holds individual molecule data (not critical at this stage).

.. note::

   💡 The folders and class names (e.g., ``bio_oil``) are based on the published research example.  
   You can adapt them to suit your own project or dataset structure.


📦 Uploading raw data
=====================

Raw data should be placed inside the ``GC_data/`` folder (as shown above) and must include:

1. 🧪 **Quantification data** — e.g., GC-MS or LC-MS output  
2. 🧬 **SMILES strings** — used for molecular analysis with RDKit  
3. ⚛️ **DFT data** — optional but required for post-analysis and model scoring  
   *(Note: FT, PT, and AG models do not require DFT data)*

📁 Example Raw Data Structure
-----------------------------

Our provided datasets (already inside ``GC_data/``) are an ideal reference:

.. image:: _static/raw_data_example.png
   :alt: Example structure of raw data folder
   :align: center
   :width: 600px

.. warning::

   ⚠️ **Important:**  
   The raw data file **must match the format** (column number and order) of the provided example.

🧬 Unpacking the raw data file
==============================

Here the raw data file is unpacked into a list of python objects (*i.e. each molecule is a python object*) that are passed to a series of functions
to generate the models.

For this example, we will use the ``pb_cp.csv`` file — containing **GC-MS** and **DFT data** for each molecule in the pine-bark derived bio-oil fraction.

🗂️ Obtain raw data filepath
---------------------------

The first task is to obtain the **filepath** of the folder containing the raw data.  
We will use the **manager module** for this, which has an attribute pointing to the **GC_data/** folder as shown below.

.. code-block:: python

    # data file name
    data_filename = "pb_cp.csv"

    # data folder path
    data_folder = manager.bio_oil_GC_data

    # create filepath to data file
    data_filepath = os.path.join(data_folder, data_filename) 

    # show data filepath
    print(data_filepath)

🔄 Convert Raw Data to Molecule Objects
---------------------------------------

Now that we have the full **data file path**, we can pass this to the function ``csv_to_orca_class``, which will convert the raw data into a list of **molecule objects**.  
These objects can be easily manipulated and analyzed using Python.

Pass the ``data_filepath`` to the function as shown below:

.. code-block:: python

   from modules.molecule_parser import csv_to_orca_class

   # Convert the CSV file into molecule objects
   molecules = csv_to_orca_class(data_filepath)

🧬 What Information is Contained in a Molecule Object?
------------------------------------------------------

Now that we've created the list of molecule objects using the ``csv_to_orca_class`` function, let's explore what information each molecule object contains.

In Python, an **object** is a collection of data (attributes) and methods (functions) that act on the data. In our case, each **molecule** is represented as an object, and these objects are stored in a **list**.

Each molecule object holds a variety of associated information, such as:

- **Molecular properties** (e.g., molecular weight, structure)
- **GC-MS data** (e.g., retention time, intensity)
- **DFT data** (e.g., energy, optimization results)

🧪 Inspecting a Single Molecule Object
--------------------------------------

To understand how the molecule object is structured, let’s inspect the attributes of one molecule.

.. code-block:: python

    # define a single molecule from the list of molecules
    molecule = molecules[0]

    # show all attributes
    molecule.__dict__

🔑 Deepdive: Molecule Object Attributes
---------------------------------------

Each molecule object contains a variety of attributes that describe its physical, structural, and electronic properties:

- **``name``**: The name of the molecule  
- **``smiles``**: The SMILES notation of the molecule  
- **``mw``**: **Molecular Weight** *(g/mol)*  
- **``peak_area``**: The proportion of the molecule in the characterization (e.g., peak area in GC-MS)  
- ***``homo_lumo_gap``***: **HOMO-LUMO gap** *(energy difference between highest occupied and lowest unoccupied molecular orbitals)*  
- ***``chemical_hardness``***: **Chemical hardness** *(related to molecular stability)*  
- ***``dipole_moment``***: The **dipole moment** *(measure of polarity)*  
- ***``polarizability``***: The **polarizability** *(how easily the electron cloud is distorted)*  
- **``volume``**: The **molecular volume**

Some of the attributes listed above are derived from **DFT data** and are marked *in italics* below:

Other attributes such as **``mw``** and **``peak_area``** typically come directly from the raw characterization data (e.g., GC-MS).

🔍 Accessing Attributes
-----------------------

You can access these attributes directly from the molecule object. For example, to print the molecular weight of the first molecule in the list:

.. code-block:: python

    print(molecules[0].mw)

The methodology is the same for showing other attributes and print statements can also be used here:

.. code-block:: python

    # extract data for the initial molecule
    print(f"The molecule is called {molecule.name} and has a molecular weight of {molecule.mw} g/mol.")
    print(f"It has an estimated volume of {molecule.volume:.3f} A^3 and a peak area in the charecterization of {float(molecule.peak_area):.3f} %.")



