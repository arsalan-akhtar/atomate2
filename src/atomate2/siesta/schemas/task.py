"""A definition of a MSON document representing an SIESTA task."""

from __future__ import annotations

import json
import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np
from emmet.core.math import Matrix3D, Vector3D
from emmet.core.structure import MoleculeMetadata, StructureMetadata



class SiestaTaskDoc(StructureMetadata, MoleculeMetadata):
    """Definition of SIESTA task document.

    Parameters
    ----------
    dir_name: str
        The directory for this SIESTA task
    last_updated: str
        Timestamp for this task document was last updated
    completed_at: str
        Timestamp for when this task was completed
    input: .InputDoc
        The input to the first calculation
    output: .OutputDoc
        The output of the final calculation
    structure: Structure or Molecule
        Final output structure from the task
    state: .TaskState
        State of this task
    included_objects: List[.AimsObject]
        List of SIESTA objects included with this task document
    aims_objects: Dict[.AimsObject, Any]
        SIESTA objects associated with this task
    entry: ComputedEntry
        The ComputedEntry from the task doc
    analysis: .AnalysisDoc
        Summary of structural relaxation and forces
    task_label: str
        A description of the task
    tags: List[str]
        Metadata tags for this task document
    author: str
        Author extracted from transformations
    icsd_id: str
        International crystal structure database id of the structure
    calcs_reversed: List[.Calculation]
        The inputs and outputs for all SIESTA runs in this task.
    transformations: Dict[str, Any]
        Information on the structural transformations, parsed from a
        transformations.json file
    custodian: Any
        Information on the custodian settings used to run this
        calculation, parsed from a custodian.json file
    additional_json: Dict[str, Any]
        Additional json loaded from the calculation directory
    """
    pass 