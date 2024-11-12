"""Core jobs for running ABINIT calculations."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, ClassVar

from abipy.flowtk.events import (
    AbinitCriticalWarning,
    NscfConvergenceWarning,
    RelaxConvergenceWarning,
    ScfConvergenceWarning,
)
from jobflow import Job, job

from atomate2.abinit.jobs.base import BaseAbinitMaker
from atomate2.abinit.sets.core import (
    LineNonSCFSetGenerator,
    NonSCFSetGenerator,
    NonScfWfqInputGenerator,
    RelaxSetGenerator,
    StaticSetGenerator,
    UniformNonSCFSetGenerator,
)

import os
import subprocess
import pandas as pd
from jobflow import Maker

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pymatgen.core.structure import Structure

    from atomate2.abinit.sets.base import AbinitInputGenerator
    from atomate2.abinit.utils.history import JobHistory

logger = logging.getLogger(__name__)


@dataclass
class StaticMaker(BaseAbinitMaker):
    """Maker to create ABINIT scf jobs.

    Parameters
    ----------
    name : str
        The job name.
    """

    calc_type: str = "scf"
    name: str = "Scf calculation"
    input_set_generator: AbinitInputGenerator = field(
        default_factory=StaticSetGenerator
    )

    CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
        ScfConvergenceWarning,
    )


@dataclass
class LineNonSCFMaker(BaseAbinitMaker):
    """Maker to create a jobs with a non-scf ABINIT calculation along a line.

    Parameters
    ----------
    name : str
        The job name.
    """

    calc_type: str = "nscf_line"
    name: str = "Line non-Scf calculation"
    input_set_generator: AbinitInputGenerator = field(
        default_factory=LineNonSCFSetGenerator
    )

    CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
        NscfConvergenceWarning,
    )


@dataclass
class UniformNonSCFMaker(BaseAbinitMaker):
    """Maker to create a jobs with a non-scf ABINIT calculation along a line.

    Parameters
    ----------
    name : str
        The job name.
    """

    calc_type: str = "nscf_uniform"
    name: str = "Uniform non-Scf calculation"
    input_set_generator: AbinitInputGenerator = field(
        default_factory=UniformNonSCFSetGenerator
    )

    CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
        NscfConvergenceWarning,
    )


@dataclass
class NonSCFMaker(BaseAbinitMaker):
    """Maker to create non SCF calculations."""

    calc_type: str = "nscf"
    name: str = "non-Scf calculation"

    input_set_generator: AbinitInputGenerator = field(
        default_factory=NonSCFSetGenerator
    )

    # Non dataclass variables:
    CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
        NscfConvergenceWarning,
    )

    @job
    def make(
        self,
        structure: Structure | None = None,
        prev_outputs: str | list[str] | None = None,
        restart_from: str | list[str] | None = None,
        history: JobHistory | None = None,
        mode: str = "uniform",
    ) -> Job:
        """Run a non-scf ABINIT job.

        Parameters
        ----------
        structure : .Structure
            A pymatgen structure object.
        mode : str
            Type of band structure calculation. Options are:
            - "line": Full band structure along symmetry lines.
            - "uniform": Uniform mesh band structure.
        """
        self.input_set_generator.mode = mode

        return super().make.original(
            self,
            structure=structure,
            prev_outputs=prev_outputs,
            restart_from=restart_from,
            history=history,
        )


@dataclass
class NonSCFWfqMaker(NonSCFMaker):
    """Maker to create non SCF calculations for the WFQ."""

    calc_type: str = "nscf_wfq"
    name: str = "non-Scf calculation"

    input_set_generator: AbinitInputGenerator = field(
        default_factory=NonScfWfqInputGenerator
    )

    # Non dataclass variables:
    CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
        NscfConvergenceWarning,
    )


@dataclass
class RelaxMaker(BaseAbinitMaker):
    """Maker to create relaxation calculations."""

    calc_type: str = "relax"
    input_set_generator: AbinitInputGenerator = field(default_factory=RelaxSetGenerator)
    name: str = "Relaxation calculation"

    # non-dataclass variables
    CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
        RelaxConvergenceWarning,
    )

    @classmethod
    def ionic_relaxation(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        return cls(
            input_set_generator=RelaxSetGenerator(*args, relax_cell=False, **kwargs),
            name=cls.name + "-ions-only",
        )

    @classmethod
    def full_relaxation(cls, *args, **kwargs) -> Job:
        """Create a full relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        return cls(
            input_set_generator=RelaxSetGenerator(*args, relax_cell=True, **kwargs),
            name=cls.name + "-ions-and-cells",
        )

#AA
#-------------------------------------------------
# @dataclass
# class AARelaxMaker(BaseAbinitMaker):
#     """Maker to create relaxation calculations."""

#     calc_type: str = "relax"
#     input_set_generator: AbinitInputGenerator = field(default_factory=RelaxSetGenerator)
#     name: str = "Relaxation calculation"

#     # non-dataclass variables
#     CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
#         RelaxConvergenceWarning,
#     )

#     @classmethod
#     def ionic_relaxation(cls, *args, **kwargs) -> Job:
#         """Create an ionic relaxation maker."""
#         # TODO: add the possibility to tune the RelaxInputGenerator options
#         #  in this class method.
#         return cls(
#             input_set_generator=RelaxSetGenerator(*args, relax_cell=False, **kwargs),
#             name=cls.name + "-ions-only",
#         )

#     @classmethod
#     def full_relaxation(cls, *args, **kwargs) -> Job:
#         """Create a full relaxation maker."""
#         # TODO: add the possibility to tune the RelaxInputGenerator options
#         #  in this class method.
#         return cls(
#             input_set_generator=RelaxSetGenerator(*args, relax_cell=True, **kwargs),
#             name=cls.name + "-ions-and-cells)",
#         )

#TODO: Move All ML stuff to ml.py
@dataclass
class AAMLRelaxMaker(BaseAbinitMaker):
    """ Maker to create ml relaxation calculations.
        TODO: Fixing/Registering nn_name in abipy
    """


    calc_type: str = "relax"
    input_set_generator: AbinitInputGenerator = field(default_factory=RelaxSetGenerator)
    name: str = "ML-Abinit Relaxation"

    # non-dataclass variables
    CRITICAL_EVENTS: ClassVar[Sequence[AbinitCriticalWarning]] = (
        RelaxConvergenceWarning,
    )

    @classmethod
    def ionic_relaxation(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("ionic-relaxation-3")
        return cls(
            input_set_generator=RelaxSetGenerator(*args,
                                                  relax_cell=False, 
                                                   user_abinit_settings={"prtwf" : 0,
                                                                         "prtden": 0,
                                                                         "prteig": 0,
                                                                         "prtddb": 0,},
                                                   **kwargs),
            name=cls.name + "-3-ions-only",
        )

    @classmethod
    def ionic_relaxation_30_none(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("ionic-relaxation-30-none")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=False,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov":30,
                                                                        "nn_correction":"none",
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-30-none-ions-only",
        )
    
    @classmethod
    def ionic_relaxation_30_delta(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("ionic-relaxation-30-delta")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=False,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov":30,
                                                                        "nn_correction":'"delta"',
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-30-delta-ions-only",
        )

    @classmethod
    def ionic_relaxation_31_none(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("ionic-relaxation-31-none")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=False,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov": 31,
                                                                        "nn_correction": "none",
                                                                        "prtwf":  0,
                                                                        "prtden": 0,
                                                                        "prteig": 0,
                                                                        "prtddb": 0,
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-31-none-ions-only",
        )
    
    @classmethod
    def ionic_relaxation_31_delta(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("ionic-relaxation-31-delta")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=False,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov":31,
                                                                        "nn_correction":'"delta"',
                                                                        "prtwf":  0,
                                                                        "prtden": 0,
                                                                        "prteig": 0,
                                                                        "prtddb": 0, 
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-31-delta-ions-only",
        )

    @classmethod
    def ionic_relaxation_32_none(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("ionic-relaxation-32-none")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=False,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov": 32,
                                                                        "nn_correction": "none",
                                                                        "prtwf":  0,
                                                                        "prtden": 0,
                                                                        "prteig": 0,
                                                                        "prtddb": 0,
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-32-none-ions-only",
        )
    
    @classmethod
    def ionic_relaxation_32_delta(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("ionic-relaxation-32-delta")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=False,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov":32,
                                                                        "nn_correction":'"delta"',
                                                                        "prtwf":  0,
                                                                        "prtden": 0,
                                                                        "prteig": 0,
                                                                        "prtddb": 0, 
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-32-delta-ions-only",
        )


    @classmethod
    def full_relaxation(cls, *args, **kwargs) -> Job:
        """Create a full relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("full_relaxation")
        return cls(
            input_set_generator=RelaxSetGenerator(*args,
                                                  relax_cell=True,
                                                  user_abinit_settings={"prtwf":  0,
                                                                        "prtden": 0,
                                                                        "prteig": 0,
                                                                        "prtddb": 0,
                                                                        },
                                                  **kwargs),
            name=cls.name + "-full",                                                  
        )


    @classmethod
    def full_relaxation_31_none(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("full-relaxation-31-none")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=True,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov": 31,
                                                                        "nn_correction": "none",
                                                                        "prtwf":  0,
                                                                        "prtden": 0,
                                                                        "prteig": 0,
                                                                        "prtddb": 0,
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-31-none-full",
        )
    
    @classmethod
    def full_relaxation_31_delta(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("full-relaxation-31-delta")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=True,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov":31,
                                                                        "nn_correction":'"delta"',
                                                                        "prtwf":  0,
                                                                        "prtden": 0,
                                                                        "prteig": 0,
                                                                        "prtddb": 0, 
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-31-delta-full",
        )

    @classmethod
    def full_relaxation_32_none(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("full-relaxation-32-none")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=True,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov": 32,
                                                                        "nn_correction": "none",
                                                                        "prtwf":  0,
                                                                        "prtden": 0,
                                                                        "prteig": 0,
                                                                        "prtddb": 0,
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-32-none-full",
        )
    
    @classmethod
    def full_relaxation_32_delta(cls, *args, **kwargs) -> Job:
        """Create an ionic relaxation maker."""
        # TODO: add the possibility to tune the RelaxInputGenerator options
        #  in this class method.
        print("full-relaxation-32-delta")
        return cls(
            input_set_generator=RelaxSetGenerator(*args, 
                                                  relax_cell=True,
                                                  user_abinit_settings={"nn_name":'"chgnet"',
                                                                        "ionmov":32,
                                                                        "nn_correction":'"delta"',
                                                                        "prtwf":  0,
                                                                        "prtden": 0,
                                                                        "prteig": 0,
                                                                        "prtddb": 0, 
                                                                        },
                                                  **kwargs),
            #input_set_generator=update_user_abinit_settings(input_set_generator,{"nn_name":"'chgnet'"}),
            name=cls.name + "-32-delta-full",
        )




@dataclass
#class DeepRelaxMaker(BaseAbinitMaker):
class DeepRelaxMaker(Maker):
    """Maker to create DeepRelax relaxation calculations."""

    @classmethod
    @job
    def prediction(cls, structure, restart_from=None)-> Structure:
        """Create a basic ionic relaxation using DeepRelax."""
        # Save structure to a CIF file (no relaxation type)
        cif_filename = cls.save_structure_to_cif(structure)
        #print(f"DEBUG: {cif_filename=}")
        # Create test.csv file
        cls.create_test_csv(cif_filename, structure)
        
        # Run DeepRelax
        #cls.run_deep_relax(cif_filename)
        
        #return f"DeepRelax ionic relaxation completed for {structure.composition.formula}"
            # Run DeepRelax and get the relaxed structure
        relaxed_structure = cls.run_deep_relax(cif_filename)

        # Return the relaxed structure
        return {"structure": relaxed_structure}
        #return relaxed_structure

    @classmethod
    @job
    def ionic_relaxation_32_none(cls, structure, restart_from=None):
        """Create a 32-none ionic relaxation using DeepRelax."""
        cif_filename = cls.save_structure_to_cif(structure)
        cls.run_deep_relax(cif_filename)
        return f"DeepRelax ionic relaxation 32-none completed for {structure.composition.formula}"

    @classmethod
    @job
    def ionic_relaxation_32_delta(cls, structure, restart_from=None):
        """Create a 32-delta ionic relaxation using DeepRelax."""
        cif_filename = cls.save_structure_to_cif(structure)
        cls.run_deep_relax(cif_filename)
        return f"DeepRelax ionic relaxation 32-delta completed for {structure.composition.formula}"

    @staticmethod
    def save_structure_to_cif(structure):
        """Save structure to a CIF file for DeepRelax."""
        # Create the CIF directory (make sure no duplicate 'CIF')
        cif_dir = os.path.join("deep_relax", "CIF")
        os.makedirs(cif_dir, exist_ok=True)
        
        # Ensure that the CIF filename is correct without extra suffix
        cif_filename = os.path.join(cif_dir, f"{structure.composition.reduced_formula}_unrelaxed.cif")
        
        # Save the structure in CIF format
        structure.to(fmt="cif", filename=cif_filename)
        return cif_filename

   
    @staticmethod
    def create_test_csv(cif_filename, structure):
        """Create a test.csv file with basic structure information."""
        #data_root = os.path.dirname(cif_filename)
        data_root = os.path.dirname(os.path.dirname(cif_filename))  # Move up one directory to 'deep_relax'
        csv_path = os.path.join(data_root, "test.csv")
        
        # Example structure of test.csv (without relaxation type)
        test_data = {
            "atoms_id": [f"{structure.composition.reduced_formula}"],
            #"atoms_id": [f"{structure.composition.reduced_formula}_unrelaxed.cif"],
            #"atoms_id": [os.path.basename(cif_filename)],  # Just the filename, no path or double extensions
            }
        test_df = pd.DataFrame(test_data)
        
        # Save the DataFrame to test.csv
        test_df.to_csv(csv_path, index=False)
        #print(f"Created test.csv at {csv_path}")

    @staticmethod
    def run_deep_relax(cif_filename):
        """Run the DeepRelax script using the given CIF file."""
        #print("DEBUG: RUN_DEEP_REALX")
        from pymatgen.io.cif import CifParser

        #data_root = os.path.dirname(cif_filename)
        data_root = os.path.dirname(os.path.dirname(cif_filename))  # One level up to 'deep_relax'
        #print(f"DEBUG:{data_root=}")
        model_path = "/home/akhtar/2-Areas/Miniconda3-apps-abinit/DeepRelax/trained_model/model.pt"

        command = [
            "python",
            "/home/akhtar/2-Areas/Miniconda3-apps-abinit/DeepRelax/predict_relaxed_structure.py",
            "--data_root", data_root,
            "--model_path", model_path,
        ]
        
        #print(f"Running DeepRelax with CIF: {cif_filename}")
        #print(f"{command=}")
        subprocess.run(command, check=True)
        print(f"DeepRelax completed.")

        # Construct the predicted CIF filename by replacing '_unrelaxed' with '_predicted'
        predicted_dir = os.path.join(os.path.dirname(data_root), "predicted_structures")  # Move one level up from 'deep_relax'
        predicted_cif = cif_filename.replace("_unrelaxed", "_predicted")
        predicted_cif = os.path.join(predicted_dir, os.path.basename(predicted_cif))  # Keep the same name pattern


        # Check if the predicted CIF file exists
        if os.path.exists(predicted_cif):
            #print(f"Predicted structure CIF found: {predicted_cif}")

        # Read the predicted CIF file using pymatgen's CifParser
            parser = CifParser(predicted_cif)
            #structure = parser.get_structures()[0]  # Get the first structure from the CIF file
            structure = parser.parse_structures(primitive=False)[0]
            #print(f"Read predicted structure from {predicted_cif} as a pymatgen Structure")
            #print(f"{structure=}")
            print(f"{type(structure)=}")
            # Optionally, return or process the structure further
            return structure
        else:
            raise FileNotFoundError(f"Predicted CIF file not found: {predicted_cif}")

