⚛️ Complex Fluid Molecular Models: Generation
===============================================

Overview
========

This guide demonstrates how to generate and inspect **complex fluid models** from data (that is supplemented with DFT calculations) that quantifies the
porportion of each molecule in a mixture.

This section of the documentation covers working with:

- A list of molecule objects with peak area and DFT attributes.
- Five model types:
  
  - 🧩 **FT** — Fixed Threshold  
  - ⚖️ **PT** — Proportional Threshold  
  - 🧬 **AG** — Abundancy Grouping  
  - 🏆 **SG** — Scored Grouping  
  - 🔍 **All** — All Molecules

Each of these models work by taking a list of molecule objects and applying predetermined algorithms to create a new list of molecules.

Prerequisites
-------------

You must have set up the code and extracted the raw data from the csv file into a list of molecule objects as shown in the documentation.

Common Model Structure
======================

All model instances are of the same class: ``complex_fluid_model``. This consistency makes it easy to analyze and compare across selection methods.

**Key Attributes:**

===========================  ==============================================================
Attribute                    Description
===========================  ==============================================================
``model_name``               Name of the model (e.g., "pb_cp_FT_model")
``molecules``                List of selected molecule objects
``group_molecules``          Groups of molecules (if any), or ``None`` if no grouping
``molecule_ratios``          Relative contribution of each molecule (peak-area-based) - this is adjusted so representative models are '100%'
``model_type``               Model identifier ("FT", "PT", "AG", "SG", "All")
``wa_mw``                    Weighted average molecular weight (g/mol)
``wa_chemical_hardness``     Weighted average chemical hardness (eV)
``wa_dipole_moment``         Weighted average dipole moment (D)
``wa_polarizability``        Weighted average polarizability (a.u.)
``wa_total_energy``          Weighted average total DFT energy (eV)
``wa_oxygen_content``        Weighted average oxygen content (%)
``wa_nitrogen_content``      Weighted average nitrogen content (%)
``wa_sulfur_content``        Weighted average sulfur content (%)
``min_mols_for_sim``         Minimum number of molecules needed for a simulation
``min_atoms_for_sim``        Minimum number of atoms for unbiased simulation
``min_vol_for_sim``          Minimum volume required
===========================  ==============================================================

Inspecting Model Output
=======================

As detailed above, all outputs from generating models are instances of the **complex_fluid_model** class attached to the assigned variable.

This means they can be inspected in the same manner. Using **__dict__** will show all model attributes.

.. code-block:: python

    model.__dict__

This means attributes of each model can also be printed.

.. code-block:: python

    # Example: Inspecting the FT model's attributes
    print(f"Model Name: {model.model_name}")
    print(f"Number of Selected Molecules: {model.min_mols_for_sim}")

    # Accessing the list of selected molecules
    print(f"Selected Molecules: {model.molecules}")
  
Examples for each generated model will be shown in the associated jupyter notebook.

Model Generation
================

🚀 Output: Model-Specific Molecule Lists
----------------------------------------

After running the functions, you will receive a **new list** of molecule objects for each model. These objects will carry attributes specific to that model type, such as:

- **FT Model**: Molecules sorted by a fixed threshold
- **PT Model**: Molecules sorted proportionally based on a threshold
- **AG Model**: Molecules grouped by their abundance in the mixture
- **SG Model**: Molecules grouped and scored according to certain criteria

Each model allows for different levels of detail depending on how the molecules are categorized and 

2 Arguments are passed to every model, and these inputs are consistent among each model:
 
1. **`model_name`**: The name should match that of your CSV file (without the .csv extension). The resulting model files will be named as `model_name_FT`, `model_name_PT`, etc., based on the model type.
2. **`orca_molecules`**: The list of molecule objects containing the associated data. These are called **"orca_molecules"** because the DFT data originates from **ORCA** calculations.

The FT model is the only model with a third argument - **selection_threshold** which will be explained below.

🧩 FT Model (Fixed Threshold)
-----------------------------

Selects molecules whose **peak area** is above a user-defined fixed threshold.
This allows us to filter out molecules based on their relative abundance in the mixture, making the model more focused on abundant components.

**selection threshold** is an argument unique to this model and is simply an integer passed to the function
The selection process for the FT model is governed by the following equation, as detailed in the publication:

$$ a_i > X $$

**Code Example:**

.. code-block:: python

   pb_cp_FT_model = complex_fluid_models.fixed_threshold_model(
       model_name="pb_cp",
       orca_molecules=molecules,
       selection_threshold=5
   )

**Inspecting FT Model:**

.. code-block:: python

   print(f"Model Name: {pb_cp_FT_model.model_name}")
   print(f"Model Type: {pb_cp_FT_model.model_type}")
   print(f"Selected Molecules: {len(pb_cp_FT_model.molecules)}")
   print(f"Weighted MW: {pb_cp_FT_model.wa_mw:.3f} g/mol")

⚖️ PT Model (Proportional Threshold)
------------------------------------

The inputs to the functions will not be covered from herein as they are the same between all types of model.

For a molecule to be selected is proportion in the mixture must exceed the selection threshold which is governed by the following equation:

.. image:: _static/pt_criteria.png
   :alt: Directory structure for project
   :align: center
   :width: 600px

**Code Example:**

.. code-block:: python

   pb_cp_PT_model = complex_fluid_models.proportional_threshold_model(
       model_name="pb_cp",
       orca_molecules=molecules
   )

The resulting model can be inspected in the ways detailed previously.

🔄 AG Model (Abundancy Grouping Model)
--------------------------------------

The **AG model** is generated by grouping molecules based on their **structural similarities** and **heteroatom content**. After grouping, the most **abundant molecule** in each group is selected. This model is useful for identifying dominant species within certain structural or functional classes.

The selection Criteria for AG Model is detailed in the flowchart below

The **AG model** follows these steps:
1. **Grouping**: Molecules are grouped based on similarities in **structure** and **heteroatom content** (i.e., the type and number of non-carbon atoms present in the molecule).
2. **Most Abundant Selection**: For each group, the molecule with the **highest peak area** (i.e., the most abundant molecule) is selected for inclusion in the model.

This approach ensures that the most **representative molecules** from each group are included in the final model, allowing for a more balanced representation of the dataset and visual representaiton of how
the groups are sorted is shown below.

.. image:: _static/ag_sg_criteria.png
   :alt: Directory structure for project
   :align: center
   :width: 600px
**Code Example:**

.. code-block:: python

   pb_cp_PT_model = complex_fluid_models.abundancy_grouped_model(
       model_name="pb_cp",
       orca_molecules=molecules
   )

The resulting model can be inspected in the ways detailed previously.

🏆 SG Model (Scored Grouping)
-----------------------------

Similar to AG, but selects molecule in each group with the **highest score**, based on this formula:

.. math::

   \text{Score} = \sum \left(\frac{X_i}{X_{\text{group average}}}\right)

- ``group_molecules`` populated with structure-based groups.

**Code Example:**

.. code-block:: python

   pb_cp_SG_model = complex_fluid_models.scored_grouped_model(
       model_name="pb_cp",
       orca_molecules=molecules
   )

The resulting model can be inspected in the ways detailed previously.

🔍 All Model (Benchmark)
------------------------

Includes **every molecule** in the dataset — no filtering or grouping.
This means that, regardless of the algorithm applied, the **ALL model** will always include every molecule, 
ensuring that you can compare the results from the other models with a comprehensive, baseline dataset.

**Code Example:**

.. code-block:: python

   pb_cp_ALL_model = complex_fluid_models.all_model(
       model_name="pb_cp",
       orca_molecules=molecules
   )

The resulting model can be inspected in the ways detailed previously.

Detailed model output inspection
================================

Use the following helper function to print key attributes of any model:

.. code-block:: python

   def print_model_info(model):
       print(f"Information about the {model.model_name}")
       print(f"Model type: {model.model_type}")
       print(f"Minimum atoms for simulation: {model.min_atoms_for_sim}")
       print(f"Minimum molecules for simulation: {model.min_mols_for_sim}\n")

       print("Weighted averages of molecular properties:")
       print(f"  MW: {model.wa_mw:.3f} g/mol")
       print(f"  Chemical hardness: {model.wa_chemical_hardness:.3f} eV")
       print(f"  Dipole moment: {model.wa_dipole_moment:.3f} D")
       print(f"  Polarizability: {model.wa_polarizability:.3f} a.u.")
       print(f"  Total energy: {model.wa_total_energy:.3f} eV")
       print(f"  Oxygen content: {model.wa_oxygen_content:.3f} %")
       print(f"  Nitrogen content: {model.wa_nitrogen_content:.3f} %")
       print(f"  Sulfur content: {model.wa_sulfur_content:.3f} %\n")

       if model.group_molecules is not None:
           print(f"Number of molecule groups: {len(model.group_molecules)}")
       else:
           print("No molecule grouping applied for this model.")

**Example Usage:**

.. code-block:: python

   print_model_info(pb_cp_AG_model)
   print_model_info(pb_cp_SG_model)
