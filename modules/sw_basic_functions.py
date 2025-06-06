from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import math
import os

def vol_from_smiles(smiles):
    """
    Calculates the molecular volume from a SMILES string using RDKit.

    Args:
        smiles (str): A valid SMILES string.

    Returns:
        float: Estimated molecular volume in cubic angstroms.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    volume = AllChem.ComputeMolVolume(mol) 
    return volume

def vol_from_mol(mol):
    """
    Calculates molecular volume from an RDKit molecule object.

    Args:
        mol (Chem.Mol): RDKit molecule object.

    Returns:
        float: Estimated molecular volume.
    """
    AllChem.EmbedMolecule(mol)
    volume = AllChem.ComputeMolVolume(mol)   
    return volume

def estimated_volume(pdb_file_path):
    """
    Estimates molecular volume from a PDB file using van der Waals radii.

    Args:
        pdb_file_path (str): Path to the PDB file.

    Returns:
        float: Estimated molecular volume (Å³), or None if loading fails.
    """
    # Van der Waals radii (in angstroms) for common elements
    vdw_radii = {
        'H': 1.20,  'C': 1.70,  'N': 1.55,  'O': 1.52,
        'P': 1.80,  'S': 1.80,  'F': 1.47,  'Cl': 1.75,
        'Br': 1.85, 'I': 1.98
    }

    # Load the PDB file without sanitizing
    mol = Chem.MolFromPDBFile(pdb_file_path, removeHs=False, sanitize=False)
    
    if not mol:
        print("Failed to load molecule from PDB.")
        return None

    # Try sanitizing, but continue if it fails
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        print(f"Warning: Molecule sanitization failed due to valence issues. Proceeding anyway. ({e})")

    # Calculate total volume
    total_volume = 0.0

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in vdw_radii:
            radius = vdw_radii[symbol]
            atom_volume = (4/3) * math.pi * (radius ** 3)
            total_volume += atom_volume
        else:
            print(f"Warning: Van der Waals radius not found for atom type '{symbol}'.")

    return total_volume / 2

# Function to calculate volume from a PDB file
def vol_from_pdb(pdb_file):
    """
    Calculates the volume of a molecule from a PDB file using RDKit's 3D modeling.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        float: Calculated molecular volume.
    """
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)  # Load the molecule from PDB file
    if mol is None:
        raise ValueError(f"Could not load molecule from PDB file: {pdb_file}")
    
    mol = Chem.AddHs(mol)  # Ensure hydrogens are added for volume calculation
    AllChem.EmbedMolecule(mol)  # Generate 3D coordinates for the molecule
    volume = AllChem.ComputeMolVolume(mol)  # Calculate the molecular volume
    return volume

def has_heteroatoms(mol):
    """
    Checks for presence of nitrogen, oxygen, or sulfur in a molecule.

    Args:
        mol (Chem.Mol): RDKit molecule object.

    Returns:
        list: Names of heteroatoms found, or 'No heteroatoms'.
    """
    heteroatom_smarts = ["[#7]", "[#8]", "[#16]"]
    heteroatom_names = ["nitrogen", "oxygen", "sulfur"]
    heteroatoms_in_mol = []
    for i in range(len(heteroatom_smarts)):
                   pattern_smarts = Chem.MolFromSmarts(heteroatom_smarts[i])
                   if mol.HasSubstructMatch(pattern_smarts) == True:
                       heteroatoms_in_mol.append(heteroatom_names[i])
    if len(heteroatoms_in_mol) == 0:
        heteroatoms_in_mol.append("No heteroatoms")
    return(heteroatoms_in_mol)

def has_rings(mol):
    """
    Detects common ring sizes in a molecule (5- to 8-membered).

    Args:
        mol (Chem.Mol): RDKit molecule object.

    Returns:
        list: Names of detected ring types, or ['N'] if none.
    """
    ring_functionals = ["[r5]", "[r6]", "[r7]", "[r8]"]
    ring_names = ["5-membered ring", "6-membered ring", "7-membered ring", "8-membered ring"]
    ring_groups = []
    ring_groups.append(ring_names)
    ring_groups.append(ring_functionals)
    ring_groups_in_mol = []
    for i in range(len(ring_groups[0])):
        pattern_smarts = Chem.MolFromSmarts(ring_groups[1][i])
        #print(type(mol))
        if mol.HasSubstructMatch(pattern_smarts) == True:
            ring_groups_in_mol.append(ring_groups[0][i])
    if len(ring_groups) != 0:
        return(ring_groups_in_mol)
    if len(ring_groups) == 0:
        no_list = ["N"]
        return(no_list)

def write_output(filepath, lines):
    """
    Writes a list of strings to a file, one per line.

    Args:
        filepath (str): Output file path.
        lines (list): List of string lines to write.
    """
    f = open(filepath, "w")
    for line in lines:
        f.write(line)
        f.write('\n')
    f.close()   
    return()

def get_homo_lumo_from_xyz(xyz_filepath):
    """
    Estimates HOMO and LUMO indices from an XYZ file based on total electron count.

    Args:
        xyz_filepath (str): Path to the .xyz file.

    Returns:
        tuple: (HOMO index, LUMO index)
    """
    atomic_numbers = {'C': 6, 'H': 1, 'O': 8, 'S': 16, 'F': 9, 'N': 7, 'Cl': 17}
    atoms = []

    with open(xyz_filepath, 'r') as xyz_file:
        for i, line in enumerate(xyz_file):
            # Skip the first two lines (header lines)
            if i < 2:
                continue
            
            # Split the line and extract the atomic symbol
            parts = line.split()
            if parts:  # Ensure the line isn't empty
                atom = parts[0].strip()
                # Verify that the atom is in the atomic_numbers dictionary
                if atom in atomic_numbers:
                    atoms.append(atom)
                else:
                    raise KeyError(f"Unknown atomic symbol '{atom}' in line: {line.strip()}")
    
    # Calculate the total atomic number and derive HOMO and LUMO indices
    atomic_num = sum(atomic_numbers[atom] for atom in atoms)
    homo_num = int(atomic_num / 2 - 1)
    lumo_num = homo_num + 1
    return homo_num, lumo_num

def SmilesToPDB(smiles_string, output_file):
    """
    Converts a SMILES string into a 3D PDB structure using Pybel.

    Args:
        smiles_string (str): Input SMILES string.
        output_file (str): Destination PDB file path.

    Returns:
        None
    """
    # Create a molecule from the SMILES string
    molecule = pybel.readstring("smi", smiles_string)

    # Generate the 3D structure
    molecule.make3D()

    # Convert the molecule to PDB format
    pdb_string = molecule.write("pdb")

    # Write the PDB string to a file
    with open(output_file, "w") as f:
        f.write(pdb_string)

def pdb_to_mol(pdb_filename):
    """
    Converts a PDB file to an RDKit Mol object.

    Args:
        pdb_filename (str): Path to PDB file.

    Returns:
        Chem.Mol: RDKit molecule, or None if failed.
    """
    with open(pdb_filename, 'r') as pdb_file:
        pdb_content = pdb_file.read()

    mol = Chem.MolFromPDBBlock(pdb_content)
    if mol is not None:
        AllChem.AssignAtomChiralTagsFromStructure(mol)

    return mol

def clean_mol_name(molecule_name):
    """
    Sanitizes molecule names by removing problematic characters.

    Args:
        molecule_name (str): Input name.

    Returns:
        str: Cleaned name.
    """
    chars_to_replace = ["(", ")", " ", ","] 
    for char in chars_to_replace:
        molecule_name = molecule_name.replace(char, '')
    return(molecule_name)

def remove_conect_master_lines(file_path):
    """
    Removes 'CONECT' and 'MASTER' lines from a PDB file.

    Args:
        file_path (str): Path to the PDB file.

    Returns:
        None
    """
    try:
        # Read the content of the file
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Filter out lines containing 'CONECT' or 'MASTER'
        filtered_lines = [line for line in lines if not ('CONECT' in line or 'MASTER' in line)]

        # Overwrite the file with the filtered content
        with open(file_path, 'w') as file:
            file.writelines(filtered_lines)
    except Exception as e:
        print(f"An error occurred when removing lines from PDB file: {e}")

def volume_model(T, a, b, c):
    """
    Predicts material expansion volume using a quadratic model.

    Args:
        T (float): Temperature.
        a (float): Constant term.
        b (float): Linear coefficient.
        c (float): Quadratic coefficient.

    Returns:
        float: Estimated volume.
    """
    # used for predicting expansion of a material
    return a + b*T + c*T**2

# Define atomic weights for common elements in PDB files
ATOMIC_WEIGHTS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "Na": 22.99, "K": 39.10, "Ca": 40.08, "Mg": 24.305, "Fe": 55.845, "F":19
}

# Extracts element symbol from a PDB line intelligently
def get_element_from_pdb_line(line):
    """
    Extracts atomic element from a line in a PDB file using various heuristics.

    Args:
        line (str): Line from a PDB file.

    Returns:
        str or None: Element symbol if found, else None.
    """
    element = line[76:78].strip()  # Try the standard element column first

    if not element or element not in ATOMIC_WEIGHTS:
        atom_name = line[12:16].strip()  # Extract atomic name (cols 13-16)
        
        # Infer element from atomic name (handling multi-character elements like 'CL', 'CA')
        possible_element = atom_name[:2].capitalize() if atom_name[:2].capitalize() in ATOMIC_WEIGHTS else atom_name[0]
        
        if possible_element in ATOMIC_WEIGHTS:
            element = possible_element

    return element if element in ATOMIC_WEIGHTS else None  # Ensure it's valid

def count_elements_in_pdb(pdb_filepath):
    """
    Counts occurrences of each element in a PDB file.

    Args:
        pdb_filepath (str): Path to the PDB file.

    Returns:
        dict: Element counts.
    """
    element_counts = {}

    with open(pdb_filepath, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith(("ATOM", "HETATM")):  # Look for atomic records
                element = get_element_from_pdb_line(line)
                if element:
                    element_counts[element] = element_counts.get(element, 0) + 1

    return element_counts

def calculate_molecular_weight(pdb_filepath):
    """
    Calculates molecular weight from a PDB file based on known atomic weights.

    Args:
        pdb_filepath (str): Path to the PDB file.

    Returns:
        tuple: (Total molecular weight, dict of element counts)
    """
    if not os.path.exists(pdb_filepath):
        raise FileNotFoundError(f"File '{pdb_filepath}' not found!")

    element_counts = count_elements_in_pdb(pdb_filepath)
    
    if not element_counts:
        print("Warning: No valid elements found in the PDB file!")
    
    molecular_weight = sum(ATOMIC_WEIGHTS[el] * count for el, count in element_counts.items())

    return molecular_weight, element_counts