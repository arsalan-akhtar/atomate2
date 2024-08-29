from ase import Atoms
from ase.atom import Atom
from ase.calculators.siesta import Siesta
from ase.units import Ry
from ase.calculators.siesta.parameters import Species
from dataclasses import dataclass, field
from atomate2.siesta.sets.utils import pymatgen_to_ase 
from atomate2.siesta.sets.utils import siesta_fdf_to_json 
from typing import Any, Dict, List
from pymatgen.core import Structure
#from ase.calculators.siesta.parameters import Specie
from pymatgen.io.core import InputGenerator, InputSet


@dataclass
class SiestaInputGenerator(InputGenerator):
    """ """
    #structure : Structure
    xc : str = field (default='PBE') # Optional xc str
    mesh_cutoff : float = field(default=100.0) # Optional mesh_cutoff float in eV 
    energy_shift : float =field(default=0.01) # Optional energy_shift float in eV
    basis_set : str = field (default='DZ') # Optional basis_set str
    kpts : List[int] = field(default_factory=lambda: [1, 1, 1]) # Optional kpts list
    fdf_arguments : Dict[str, Any] = field(default_factory=dict) # Optional fdf_arguments
    species: List[Species] = field(default_factory=list)  # Optional species list

    basis_set_allowed = ['SZ','SZP','DZ','DZP','TZP']  # Allowed basis sets

    def __post_init__(self):
        """ """
        # Validate the basis set
        if self.basis_set not in self.basis_set_allowed:
            raise ValueError(f"Invalid basis set '{self.basis_set}'. Allowed values are: {self.basis_set_allowed}")

        self.fdf_arguments = {"LongOutput": "True",
                              "WriteForces":"True",}
        #self.ase_atoms = pymatgen_to_ase(structure=self.structure)
        self.siesta_input_generator = Siesta(
              #pseudo_path = "/home/aakhtar/1-Project/Hoshtane/calculations/siesta/2-test-siesta-ase/",
              label='siesta',
              xc=self.xc,
              mesh_cutoff=self.mesh_cutoff * Ry,
              energy_shift=self.energy_shift * Ry,
              basis_set=self.basis_set,
              kpts=self.kpts,
              fdf_arguments=self.fdf_arguments,
              species = self.species,
              )
        # Some Debug Prints 
        print(f"{self.xc=}")
        print(f"{self.mesh_cutoff=} Ry")
        print(f"{self.energy_shift=} Ry")
        print(f"{self.basis_set=}")
        print(f"{self.kpts=}")
        print(f"{self.fdf_arguments=}")
        print(f"{self.species=}")


    #def write(self):
    #    """ """
    #    self.siesta_input_generator.write_input(self.ase_atoms,'density')
    
    def write_siesta_fdf(self,structure: Structure,directory=None):
        """ """
        ase_atoms = pymatgen_to_ase(structure=structure)
        self.siesta_input_generator.write_input(ase_atoms,'density')
        siesta_fdf_to_json("siesta.fdf",json_output_path="siesta_input.json")
    
    #def __repr__(self):
    #    return f"SiestaInputGenerator(structure={self.structure}, ...)"