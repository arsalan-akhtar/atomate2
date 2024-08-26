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
