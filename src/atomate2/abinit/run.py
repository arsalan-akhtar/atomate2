"""Functions to run ABINIT."""

from __future__ import annotations

import logging
import subprocess
import time

from abipy.flowtk.qutils import time2slurm

from atomate2 import SETTINGS
from atomate2.abinit.utils.common import (
    INPUT_FILE_NAME,
    LOG_FILE_NAME,
    STDERR_FILE_NAME,
)

__all__ = ["run_abinit"]


SLEEP_TIME_STEP = 30


logger = logging.getLogger(__name__)


def run_abinit(
    abinit_cmd: str = None,
    mpirun_cmd: str = None,
    wall_time: int = None,
    start_time: float = None,
    mpirun_np:int=None, #AA
) -> None:
    """Run ABINIT."""
    abinit_cmd = abinit_cmd or SETTINGS.ABINIT_CMD
    mpirun_cmd = mpirun_cmd or SETTINGS.ABINIT_MPIRUN_CMD
    mpirun_np = mpirun_np or SETTINGS.ABINIT_MPIRUN_NP
    start_time = start_time or time.time()
    if mpirun_np is not None:
        command = [mpirun_cmd,"-np",str(mpirun_np), abinit_cmd] if mpirun_cmd is not None else [abinit_cmd]
    else:
        command = [mpirun_cmd, abinit_cmd] if mpirun_cmd is not None else [abinit_cmd]
    ##AA
    print(f"Debug: {command=}")
    print(f"Debug: {mpirun_np=}")

    max_end_time = 0.0
    if wall_time is not None:
        abinit_timelimit = wall_time
        if abinit_timelimit > 480:
            # TODO: allow tuning this timelimit buffer for abinit,
            #  e.g. using a config variable or possibly per job
            abinit_timelimit -= 240
        command.extend(["--timelimit", time2slurm(abinit_timelimit)])
        max_end_time = start_time + wall_time

    command.append(INPUT_FILE_NAME)
    ##AA
    print(f"Debug: {command=}")
    ##AA
    with open(LOG_FILE_NAME, "w") as stdout, open(STDERR_FILE_NAME, "w") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr)  # noqa: S603

        if wall_time is not None:
            while True:
                time.sleep(SLEEP_TIME_STEP)
                if process.poll() is not None:
                    break
                current_time = time.time()
                remaining_time = max_end_time - current_time
                if remaining_time < 5 * SLEEP_TIME_STEP:
                    process.terminate()

        process.wait()
    return
