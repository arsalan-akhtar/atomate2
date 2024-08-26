from ase import Atoms
from ase.atom import Atom
from ase.calculators.siesta import Siesta
from ase.units import Ry
from ase.calculators.siesta.parameters import Species
from dataclasses import dataclass, field
from atomate2.siesta.sets.utils import pymatgen_to_sisl,pymatgen_to_ase


@dataclass
class SiestaInputGenerator:
    """ """
    #user_params: dict[str, Any] = field(default_factory=dict)
    #user_kpoints_settings: dict[str, Any] = field(default_factory=dict)

    siesta_input_generator = Siesta(#species = None,
              #pseudo_path = "/home/aakhtar/1-Project/Hoshtane/calculations/siesta/2-test-siesta-ase/",
              label='siesta',
              xc='PBE',
              mesh_cutoff=200 * Ry,
              energy_shift=0.01 * Ry,
              basis_set='DZ',
              kpts=[10, 10, 10],
              fdf_arguments={'DM.MixingWeight': 0.1,
                             'MaxSCFIterations': 100},
              )
    
    def write(self,structure):
        ase_atoms = pymatgen_to_ase(structure=structure)
        self.siesta_input_generator.write_input(ase_atoms,'density')
