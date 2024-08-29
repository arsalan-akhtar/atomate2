"""A definition of a MSON document representing an Abinit task."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional, Union

#from abipy.abio.inputs import AbinitInput
#from abipy.flowtk import events
from emmet.core.math import Matrix3D, Vector3D
from emmet.core.structure import StructureMetadata
from pydantic import BaseModel, Field
from pymatgen.core import Structure
from typing_extensions import Self

from atomate2.siesta.files import  load_siesta_input

#from atomate2.abinit.utils.common import LOG_FILE_NAME, MPIABORTFILE

from atomate2.siesta.schemas.calculation import Calculation 
#from atomate2.siesta.schemas.calculation import SiestaObject
from atomate2.siesta.schemas.calculation import TaskState
from atomate2.utils.datetime import datetime_str
from atomate2.utils.path import get_uri, strip_hostname

logger = logging.getLogger(__name__)


class InputDoc(BaseModel):
    """Summary of the inputs for an SIESTA calculation.

    Parameters
    ----------
    structure: Structure
        The final pymatgen Structure of the final system
    """

    #structure: Union[Structure] = Field(None, description="The input structure object")
    #abinit_input: SiestaInput = Field(
    #    None, description="SiestaInput used to perform calculation."
    #)
    xc_functional: str = Field(
        None, description="Exchange-correlation functional used if not the default"
    )
    xc_authors: str = Field( None, description="Exchange-correlation authors used if not the default")
    @classmethod
    def from_siesta_calc_doc(cls, calc_doc: Calculation) -> Self:
        """Create a summary from an SIESTA CalculationDocument.

        Parameters
        ----------
        calc_doc: .Calculation
            An SIESTA calculation document.

        Returns
        -------
        .InputDoc
            The calculation input summary.
        """
        siesta_input = load_siesta_input(calc_doc.dir_name)
        print(f"{siesta_input}")
        return cls(
            #structure=abinit_input.structure,
            #abinit_input=abinit_input,
            xc_functional=siesta_input["XC.functional"],
            xc_authors=siesta_input["XC.authors"],
        )


class OutputDoc(BaseModel):
    """Summary of the outputs for an Abinit calculation.

    Parameters
    ----------
    structure: Structure
        The final pymatgen Structure of the final system
    trajectory: List[Structure]
        The trajectory of output structures
    energy: float
        The final total DFT energy for the last calculation
    energy_per_atom: float
        The final DFT energy per atom for the last calculation
    bandgap: float
        The DFT bandgap for the last calculation
    cbm: float
        CBM for this calculation
    vbm: float
        VBM for this calculation
    forces: List[Vector3D]
        Forces on atoms from the last calculation
    stress: Matrix3D
        Stress on the unit cell from the last calculation
    """

    structure: Union[Structure] = Field(None, description="The output structure object")
    trajectory: Optional[Sequence[Union[Structure]]] = Field(
        None, description="The trajectory of output structures"
    )
    energy: float = Field(
        None, description="The final total DFT energy for the last calculation"
    )
    energy_per_atom: float = Field(
        None, description="The final DFT energy per atom for the last calculation"
    )
    bandgap: Optional[float] = Field(
        None, description="The DFT bandgap for the last calculation"
    )
    cbm: Optional[float] = Field(None, description="CBM for this calculation")
    vbm: Optional[float] = Field(None, description="VBM for this calculation")
    forces: Optional[list[Vector3D]] = Field(
        None, description="Forces on atoms from the last calculation"
    )
    stress: Optional[Matrix3D] = Field(
        None, description="Stress on the unit cell from the last calculation"
    )

    @classmethod
    def from_siesta_calc_doc(cls, calc_doc: Calculation) -> Self:
        """Create a summary from an abinit CalculationDocument.

        Parameters
        ----------
        calc_doc: .Calculation
            An Abinit calculation document.

        Returns
        -------
        .OutputDoc
            The calculation output summary.
        """
        return cls(
            structure=calc_doc.output.structure,
            #energy=calc_doc.output.energy,
            #energy_per_atom=calc_doc.output.energy_per_atom,
            bandgap=calc_doc.output.bandgap,
            cbm=calc_doc.output.cbm,
            vbm=calc_doc.output.vbm,
            forces=calc_doc.output.forces,
            stress=calc_doc.output.stress,
        )


class SiestaTaskDoc(StructureMetadata):
    """Definition of Abinit task document.

    Parameters
    ----------
    dir_name: str
        The directory for this Abinit task
    last_updated: str
        Timestamp for when this task document was last updated
    completed_at: str
        Timestamp for when this task was completed
    input: .InputDoc
        The input to the first calculation
    output: .OutputDoc
        The output of the final calculation
    structure: Structure
        Final output structure from the task
    state: .TaskState
        State of this task
    included_objects: List[.AbinitObject]
        List of Abinit objects included with this task document
    abinit_objects: Dict[.AbinitObject, Any]
        Abinit objects associated with this task
    task_label: str
        A description of the task
    tags: List[str]
        Metadata tags for this task document
    author: str
        Author extracted from transformations
    icsd_id: str
        International crystal structure database id of the structure
    calcs_reversed: List[.Calculation]
        The inputs and outputs for all Abinit runs in this task.
    transformations: Dict[str, Any]
        Information on the structural transformations, parsed from a
        transformations.json file
    custodian: Any
        Information on the custodian settings used to run this
        calculation, parsed from a custodian.json file
    additional_json: Dict[str, Any]
        Additional json loaded from the calculation directory
    """

    dir_name: Optional[str] = Field(
        None, description="The directory for this Abinit task"
    )
    last_updated: Optional[str] = Field(
        default_factory=datetime_str,
        description="Timestamp for when this task document was last updated",
    )
    completed_at: Optional[str] = Field(
        None, description="Timestamp for when this task was completed"
    )
    input: Optional[InputDoc] = Field(
        None, description="The input to the first calculation"
    )
    output: Optional[OutputDoc] = Field(
        None, description="The output of the final calculation"
    )
    structure: Union[Structure] = Field(
        None, description="Final output atoms from the task"
    )
    state: Optional[TaskState] = Field(None, description="State of this task")
    
    task_label: Optional[str] = Field(None, description="A description of the task")
    tags: Optional[list[str]] = Field(
        None, description="Metadata tags for this task document"
    )
    author: Optional[str] = Field(
        None, description="Author extracted from transformations"
    )
    icsd_id: Optional[str] = Field(
        None, description="International crystal structure database id of the structure"
    )
    calcs_reversed: Optional[list[Calculation]] = Field(
        None, description="The inputs and outputs for all Abinit runs in this task."
    )
    transformations: Optional[dict[str, Any]] = Field(
        None,
        description="Information on the structural transformations, parsed from a "
        "transformations.json file",
    )
    custodian: Any = Field(
        None,
        description="Information on the custodian settings used to run this "
        "calculation, parsed from a custodian.json file",
    )
    additional_json: Optional[dict[str, Any]] = Field(
        None, description="Additional json loaded from the calculation directory"
    )

    @classmethod
    def from_directory(
        cls,
        dir_name: Path | str,
        additional_fields: dict[str, Any] = None,
        **siesta_calculation_kwargs,
    ) -> Self:
        """Create a task document from a directory containing SIESTA files.

        Parameters
        ----------
        dir_name: Path or str
            The path to the folder containing the calculation outputs.
        additional_fields: Dict[str, Any]
            Dictionary of additional fields to add to output document.
        **siesta_calculation_kwargs
            Additional parsing options that will be passed to the
            :obj:`.Calculation.from_abinit_files` function.

        Returns
        -------
        .AbinitTaskDoc
            A task document for the calculation.
        """
        logger.info(f"Getting task doc in: {dir_name}")
        print(f"Getting task doc in: {dir_name}".upper())

        if additional_fields is None:
            additional_fields = {}

        dir_name = Path(dir_name)
        task_files = _find_siesta_files(dir_name)
        print(f"{task_files=}")
        if len(task_files) == 0:
            raise FileNotFoundError("No Siesta files found!")

        calcs_reversed = []
        #all_abinit_objects = []
        for files_name, files in task_files.items():
            print(f"{files_name=} {files=}")
            calc_doc = Calculation.from_siesta_files(
                dir_name, siesta_output_file='siesta.out' #, **files 
            )
            calcs_reversed.append(calc_doc)
            #all_abinit_objects.append(abinit_objects)
        print(f"{calcs_reversed=}")
        print(f"in tags stage".upper())
        tags = additional_fields.get("tags")

        dir_name = get_uri(dir_name)  # convert to full uri path
        dir_name = strip_hostname(
            dir_name
        )  # VT: TODO to put here?necessary with laptop at least...

        # only store objects from last calculation
        # TODO: make this an option
        #abinit_objects = all_abinit_objects[-1]
        included_objects = None
        #if abinit_objects:
        #    included_objects = list(abinit_objects.keys())

        # rewrite the original structure save!

        if isinstance(calcs_reversed[-1].output.structure, Structure):
            attr = "from_structure"
            dat = {
                "structure": calcs_reversed[-1].output.structure,
                "meta_structure": calcs_reversed[-1].output.structure,
                "include_structure": True,
            }
        doc = getattr(cls, attr)(**dat)
        ddict = doc.dict()

        data = {
            #"abinit_objects": abinit_objects,
            "calcs_reversed": calcs_reversed,
            "completed_at": calcs_reversed[-1].completed_at,
            "dir_name": dir_name,
            "included_objects": included_objects,
            "input": InputDoc.from_siesta_calc_doc(calcs_reversed[0]),
            "meta_structure": calcs_reversed[-1].output.structure,
            "output": OutputDoc.from_siesta_calc_doc(calcs_reversed[-1]),
            "state": calcs_reversed[-1].has_siesta_completed,
            "structure": calcs_reversed[-1].output.structure,
            "tags": tags,
        }

        doc = cls(**ddict)
        doc = doc.model_copy(update=data)
        return doc.model_copy(update=additional_fields, deep=True)


def _find_siesta_files(
    path: Path | str,
) -> dict[str, Any]:
    """Find SIESTA files in a directory.

    Only files in folders with names matching a task name (or alternatively files
    with the task name as an extension, e.g., abinit.out) will be returned.

    Abinit files in the current directory will be given the task name "standard".

    Parameters
    ----------
    path: str or Path
        Path to a directory to search.

    Returns
    -------
    dict[str, Any]
        The filenames of the calculation outputs for each Abinit task,
        given as a ordered dictionary of::

            {
                task_name: {
                    "abinit_output_file": abinit_output_filename,
    """
    files =["siesta.out","MESSAGES"]
    path = Path(path)
    siesta_files = {}
    for file in files:
        siesta_files["siesta.out"] = Path(file)
        siesta_files["MESSAGES"] = Path(file)

    return siesta_files


