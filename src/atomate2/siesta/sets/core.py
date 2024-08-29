from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field

from atomate2.siesta.sets.base import SiestaInputGenerator
from typing import TYPE_CHECKING
from pymatgen.core import Structure


from ase.calculators.siesta import Siesta
from ase.units import Ry

if TYPE_CHECKING:
    from typing import Any

    from pymatgen.core import Molecule

@dataclass
class RelaxSetGenerator(SiestaInputGenerator):
    """RelaxSetGenerator class for relaxation-specific calculations."""

    relax_cell: bool = True
    max_force: float = 1e-3  # Convergence criteria for relaxation in eV/Å
    method: str = "CG"  # Default method for relaxation (conjugate gradient)
    number_cg_steps = 100
    
    def __init__(self, **kwargs):
        """Initialize RelaxSetGenerator with a structure and optional parameters."""
        super().__init__(**kwargs)

        # Set relaxation-specific FDF arguments
        relaxation_fdf_arguments = {
            "MD.TypeOfRun": self.method,
            "MD.NumCGsteps":self.number_cg_steps,
            "MD.MaxForceTol": f"{self.max_force}  eV/Ang",
            "MD.MaxStressTol": f"{self.max_force if self.relax_cell else 0} GPa" ,
            "MD.VariableCell": "true" if self.relax_cell else "false",
            "MD.UseSaveXV": "false",
        }

        # Combine user-defined FDF arguments with relaxation-specific ones
        self.fdf_arguments = {**self.fdf_arguments, **relaxation_fdf_arguments}


        # Call the parent __post_init__ to initialize the Siesta calculator
        super().__post_init__()



@dataclass
class RelaxSetGeneratorOLD2(SiestaInputGenerator):
    """RelaxSetGenerator class for relaxation-specific calculations."""

    relax_cell: bool = True
    max_force: float = 1e-3  # Convergence criteria for relaxation in eV/Å
    method: str = "CG"  # Default method for relaxation (conjugate gradient)
    number_cg_steps = 100
    
    def __init__(self, structure: Structure, **kwargs):
        """Initialize RelaxSetGenerator with a structure and optional parameters."""
        super().__init__(structure=structure, **kwargs)

        # Set relaxation-specific FDF arguments
        relaxation_fdf_arguments = {
            "MD.TypeOfRun": self.method,
            "MD.NumCGsteps":self.number_cg_steps,
            "MD.MaxForceTol": self.max_force + "eV/Ang",
            "MD.MaxStressTol": self.max_force if self.relax_cell else 0 ,
            "MD.VariableCell": "true" if self.relax_cell else "false",
            "MD.UseSaveXV": "false",
        }

        # Combine user-defined FDF arguments with relaxation-specific ones
        self.fdf_arguments = {**self.fdf_arguments, **relaxation_fdf_arguments}

        # Initialize the Siesta calculator again with the updated FDF arguments
        self.siesta_input_generator = Siesta(
            label='siesta',
            xc=self.xc,
            mesh_cutoff=self.mesh_cutoff * Ry,
            energy_shift=self.energy_shift * Ry,
            basis_set=self.basis_set,
            kpts=self.kpts,
            fdf_arguments=self.fdf_arguments,
            species=self.species,
        )

    #@classmethod
    #def from_structure(cls, structure: Structure, **kwargs):
    #    """Create a RelaxSetGenerator instance from a pymatgen Structure."""
    #    return cls(structure=structure, **kwargs)

    #def generate(self):
    #    """Method to explicitly set up for relaxation if needed."""
    #    print("Generating relaxation input settings...")


@dataclass
class RelaxSetGeneratorOLD(SiestaInputGenerator):
    """RelaxSetGenerator class for relaxation-specific calculations."""

    relax_cell: bool = True
    max_force: float = 1e-3  # Convergence criteria for relaxation in eV/Å
    method: str = "CG"  # Default method for relaxation (conjugate gradient)

    def __post_init__(self):
        """Initialize the Siesta calculator with relaxation-specific settings."""
        # Set relaxation-specific FDF arguments
        relaxation_fdf_arguments = {
            "MD.TypeOfRun": self.method,
            "MD.MaxForceTol": self.max_force,
            "MD.MaxStressTol": self.max_force if self.relax_cell else 0,
            "MD.VariableCell": "true" if self.relax_cell else "false",
            "MD.UseSaveXV": "false",
        }

        # Combine user-defined FDF arguments with relaxation-specific ones
        self.fdf_arguments = {**self.fdf_arguments, **relaxation_fdf_arguments}

        # Call the parent __post_init__ to initialize the Siesta calculator
        super().__post_init__()
        

    @classmethod
    def from_structure(cls, structure: Structure, **kwargs):
        """Create a RelaxSetGenerator instance from a pymatgen Structure."""
        return cls(structure=structure, **kwargs)

    def generate(self):
        """Method to explicitly set up for relaxation if needed."""
        print("Generating relaxation input settings...")
    # def generate(self) -> None:
    #     """Generate settings specific for relaxation."""

    #     # Define additional FDF arguments specific to relaxation
    #     relaxation_fdf_arguments = {
    #         "MD.TypeOfRun": "CG",
    #         "MD.MaxForceTol": self.max_force,
    #         "MD.MaxStressTol": self.max_force if self.relax_cell else 0,
    #         "MD.VariableCell": "true" if self.relax_cell else "false",
    #         "MD.UseSaveXV": "false",
    #     }

    #     # Combine user-defined FDF arguments with relaxation-specific ones
    #     self.fdf_arguments = {**self.fdf_arguments, **relaxation_fdf_arguments}

    #     # Update the Siesta calculator with new FDF arguments
    #     self.siesta_input_generator.fdf_arguments.update(self.fdf_arguments)

    # def generate(self,structure: Structure | Molecule) -> SiestaInputGenerator:
    #     """Generate a SiestaInputGenerator with relaxation settings."""

    #     # Define additional FDF arguments specific to relaxation
    #     relaxation_fdf_arguments = {

    #         "MD.TypeOfRun": "CG",
    #         "MD.MaxForceTol": self.max_force,
    #         "MD.MaxStressTol": self.max_force if self.relax_cell else 0,
    #         "MD.VariableCell": "true" if self.relax_cell else "false",
    #         "MD.UseSaveXV": "false",
    #     }

    #     # Combine user-defined FDF arguments with relaxation-specific ones
    #     combined_fdf_arguments = {**self.fdf_arguments, **relaxation_fdf_arguments}

    #     # Create a new SiestaInputGenerator instance with combined settings
    #     return SiestaInputGenerator(
    #         structure=self.structure,
    #         mesh_cutoff=self.mesh_cutoff,
    #         fdf_arguments=combined_fdf_arguments,
    #         kpts=self.kpts,
    #         basis_set=self.basis_set,
    #         species=self.species
    #     )
    
    # def get_parameter_updates(self, structure: Structure | Molecule, prev_parameters: dict[str, Any]) -> dict:
    #     """Get the parameter updates for the calculation.

    #     Args:
    #         structure (Structure or Molecule): The structure to calculate the bands for
    #     prev_parameters (Dict[str, Any]): The previous parameters

    #     Returns:
    #         dict: The updated for the parameters for the output section of FHI-aims
    #     """
    #     updates = {"relax_geometry": f"{self.method} {self.max_force:e}"}
    #     if isinstance(structure, Structure) and self.relax_cell:
    #         updates["relax_unit_cell"] = "full"
    #     elif isinstance(structure, Structure):
    #         updates["relax_unit_cell"] = "none"

    #     return updates


@dataclass
class SocketIOSetGenerator(SiestaInputGenerator):
    """ """
    pass 

@dataclass
class StaticSetGenerator(SiestaInputGenerator):
    """ """
    pass 

@dataclass
class BandStructureSetGenerator(SiestaInputGenerator):
    """ """
    pass 