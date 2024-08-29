from pymatgen.core import Structure
import sisl

from ase import Atoms
import numpy as np

def pymatgen_to_ase(structure, ghost_tags=None):
    """
    Converts a Pymatgen Structure object to an ASE Atoms object.
    
    Parameters:
    structure (pymatgen.core.Structure): A Pymatgen Structure object.
    ghost_tags (list of bool, optional): List indicating which atoms are ghost atoms. 
                                         If None, the function will check for 'ghost_tags' in structure.site_properties.
    
    Returns:
    ase.Atoms: An ASE Atoms object with adjusted atomic numbers for ghost atoms.
    """
    # Add ghost_tags as a site property if provided
    if ghost_tags is not None:
        structure.add_site_property("ghost_tags", ghost_tags)        

    # Extract lattice and atomic positions from Pymatgen structure
    lattice_matrix = structure.lattice.matrix  # 3x3 array representing the lattice vectors
    atomic_positions = structure.frac_coords   # Fractional atomic coordinates
    species = structure.species                # Atomic species
    
    # Get atomic numbers
    atomic_numbers = [site.specie.Z for site in structure]
    
    # Check if 'ghost_tags' is available in site properties 
    if "ghost_tags" in structure.site_properties:
        ghost_tags = structure.site_properties["ghost_tags"]
        # Adjust atomic numbers based on ghost tags
        atomic_numbers = [-num if ghost else num for num, ghost in zip(atomic_numbers, ghost_tags)]
    else:
        # If ghost_tags are not available, show a message and proceed without ghost atoms
        print("Note: The structure does not have 'ghost_tags' site property. Proceeding without ghost atoms.")
    
    # Convert fractional coordinates to Cartesian coordinates
    cartesian_positions = np.dot(atomic_positions, lattice_matrix)

    # Create ASE Atoms object
    ase_atoms = Atoms(numbers=atomic_numbers, positions=cartesian_positions, cell=lattice_matrix, pbc=True)
    
    return ase_atoms

# Example usage
# Assuming `structure` is a Pymatgen Structure object with or without the 'ghost_tags' property
# ase_atoms = pymatgen_to_ase(structure)
# print(ase_atoms)


def pymatgen_to_sisl(structure,ghost_tags=None):
    """
    Converts a Pymatgen Structure object to a sisl Geometry object.
    
    Parameters:
    structure (pymatgen.core.Structure): A Pymatgen Structure object.

    Returns:
    sisl.Geometry: A sisl Geometry object with adjusted atomic numbers for ghost atoms.
    """
    if ghost_tags is not None:
        structure.add_site_property("ghost_tags", ghost_tags)        

    # Extract lattice and atomic positions from Pymatgen structure
    lattice_matrix = structure.lattice.matrix  # 3x3 array representing the lattice vectors
    atomic_positions = structure.frac_coords   # Fractional atomic coordinates
    species = structure.species                # Atomic species
    
    # Get atomic numbers and adjust based on ghost tags
    atomic_numbers = [site.specie.Z for site in structure]
    
    # Check if 'ghost_tags' is available in site properties 
    if "ghost_tags" in structure.site_properties:
        ghost_tags = structure.site_properties["ghost_tags"]
        atomic_numbers = [-num if ghost else num for num, ghost in zip(atomic_numbers, ghost_tags)]
    else:
        print("Note: The structure does not have 'ghost_tags' site property. Proceeding without ghost atoms.")
        # Optionally, you could also use a warning
        # warnings.warn("The structure does not have 'ghost_tags' site property. Proceeding without ghost atoms.")
        ghost_tags = [False] * len(atomic_numbers)  # No ghost atoms
        #raise ("The structure does not have 'ghost_tags' site property.")
        #pass 

 
    
    # Create a list of sisl Atom objects
    atoms = [sisl.Atom(an) for an in atomic_numbers]

    # Create sisl Geometry object
    geo = sisl.Geometry(atomic_positions, atoms, lattice=lattice_matrix)
    
    return geo

# Example usage
# Assuming `structure` is a Pymatgen Structure object with the 'ghost_tags' property
# geo = pymatgen_to_sisl(structure)
# print(geo)

import yaml

def read_outvars(file_path):
    """
    Reads the OUTVARS.yml file and returns its contents as a dictionary.
    
    Parameters
    ----------
    file_path : str
        The path to the OUTVARS.yml file.
        
    Returns
    -------
    dict
        A dictionary containing the contents of the YAML file.
    """
    try:
        with open(file_path, 'r') as file:
            data = yaml.safe_load(file)
        return data
    except FileNotFoundError:
        print(f"The file {file_path} does not exist.")
    except yaml.YAMLError as e:
        print(f"Error reading the YAML file: {e}")

# Usage example
#file_path = '/mnt/data/OUTVARS.yml'  # Path to your file
#outvars_data = read_outvars(file_path)

# Print the contents of the YAML file
#print(outvars_data)


import sisl
import json

def siesta_fdf_to_json(siesta_fdf_path, json_output_path):
    # Read the FDF file using sisl
    siesta_fdf = sisl.get_sile(siesta_fdf_path)
    
    # Create a dictionary to store the extracted information
    fdf_data = {}

    # Extracting simple key-value pairs
    fdf_data['SCFMustConverge'] = siesta_fdf.get("SCFMustConverge")
    fdf_data['Spin'] = siesta_fdf.get("Spin")
    fdf_data['XC.functional'] = siesta_fdf.get("XC.functional")
    fdf_data['XC.authors'] = siesta_fdf.get("XC.authors")
    fdf_data['MeshCutoff'] = siesta_fdf.get("MeshCutoff")
    fdf_data['PAO.EnergyShift'] = siesta_fdf.get("PAO.EnergyShift")
    fdf_data['NumberOfSpecies'] = siesta_fdf.get("NumberOfSpecies")
    fdf_data['NumberOfAtoms'] = siesta_fdf.get("NumberOfAtoms")
    fdf_data['LatticeConstant'] = siesta_fdf.get("LatticeConstant")
    fdf_data['AtomicCoordinatesFormat'] = siesta_fdf.get("AtomicCoordinatesFormat")
    fdf_data['DM.UseSaveDM'] = siesta_fdf.get("DM.UseSaveDM")
    fdf_data['SaveRho'] = siesta_fdf.get("SaveRho")
    fdf_data['WriteFoces'] = siesta_fdf.get("WriteFoces")
    fdf_data['LongOutput'] = siesta_fdf.get("LongOutput")

    # Extracting blocks
    fdf_data['ChemicalSpeciesLabel'] = siesta_fdf.get("ChemicalSpeciesLabel")
    fdf_data['PAO.BasisSizes'] = siesta_fdf.get("PAO.BasisSizes")
    fdf_data['LatticeVectors'] = siesta_fdf.get("LatticeVectors")
    fdf_data['AtomicCoordinatesAndAtomicSpecies'] = siesta_fdf.get("AtomicCoordinatesAndAtomicSpecies")
    fdf_data['DM.InitSpin'] = siesta_fdf.get("DM.InitSpin")
    fdf_data['kgrid_Monkhorst_Pack'] = siesta_fdf.get("kgrid_Monkhorst_Pack")

    # Convert the dictionary to a JSON formatted string
    json_data = json.dumps(fdf_data, indent=4)

    # Write the JSON string to a file
    with open(json_output_path, 'w') as json_file:
        json_file.write(json_data)

    print(f"JSON file '{json_output_path}' created successfully.")

# Example usage
#siesta_fdf_to_json('siesta.fdf', 'siesta_input.json')
