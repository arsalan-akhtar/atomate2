from dataclasses import dataclass, field
from jobflow import Flow, Job, job
from pymatgen.core import Structure
from typing import ClassVar, Sequence
from pathlib import Path

from typing import TYPE_CHECKING
from jobflow import Flow, Maker
from abipy.flowtk.events import AbinitCriticalWarning , RelaxConvergenceWarning
from atomate2.abinit.jobs.base import BaseAbinitMaker
from atomate2.abinit.sets.base import AbinitInputGenerator
from atomate2.abinit.sets.core import RelaxSetGenerator
from atomate2.abinit.jobs.core import RelaxMaker,DeepRelaxMaker
if TYPE_CHECKING:
    from pathlib import Path

    from pymatgen.core.structure import Structure
    

from dataclasses import dataclass, field
from jobflow import Flow
from pymatgen.core import Structure
from pathlib import Path

@dataclass
class DeepRelaxFlowMaker(Maker):
    """
    Maker to generate a relaxation flow with DeepRelax and Abinit.

    Parameters
    ----------
    name : str
        Name of the flows produced by this maker.
    relaxation_makers : list[Maker]
        The maker or list of makers to use for the relaxation flow.
    """

    name: str = "DeepRelaxation Flow"
    
    # Makers for different relaxations
    #deep_relax_maker: DeepRelaxMaker = field(default_factory=lambda: DeepRelaxMaker(input_set_generator=RelaxSetGenerator()))
    deep_relax_maker: DeepRelaxMaker = field(default_factory=lambda: DeepRelaxMaker())
    deep_relax_ionmov_3_maker: RelaxMaker = field(default_factory=lambda: [RelaxMaker.ionic_relaxation()])
    relax_maker: RelaxMaker = field(default_factory=lambda: [RelaxMaker.ionic_relaxation()])
    #deep_relax_ml_32_delta_maker: DeepRelaxMaker = field(default_factory=lambda: [DeepRelaxMaker.ionic_relaxation_32_delta()])

    def make(
        self,
        structure: Structure | None = None,
        restart_from: str | Path | None = None,
    ) -> Flow:
        """Create a DeepRelax and Abinit relaxation flow.

        Parameters
        ----------
        structure : Structure
            A pymatgen structure object.
        restart_from : str or Path or None
            One previous directory to restart from.

        Returns
        -------
        Flow
            A relaxation flow.
        """
        jobs = []

        # Step 1:(DeepRelax structure prediction)
        if self.deep_relax_maker:
            deep_relax_job = self.deep_relax_maker.prediction(structure=structure)
            jobs.append(deep_relax_job)

        # Step 2: (Abinit relaxation with deeprelax strcutre prediction)
        if self.deep_relax_ionmov_3_maker:
            deep_relax_ionmov_3_job = self.deep_relax_ionmov_3_maker[0].make(structure=jobs[-1].output.structure ,restart_from=restart_from)
            jobs.append(deep_relax_ionmov_3_job)
        
        # Step 3: (Abinit relaxation without deeprelax strcutre prediction)    
        if self.relax_maker:
            relax_maker_job = self.relax_maker[0].make(structure=structure,restart_from=restart_from)
            jobs.append(relax_maker_job)


        return Flow(jobs, name=self.name)


