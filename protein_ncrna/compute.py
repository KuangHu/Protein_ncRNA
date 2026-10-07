"""Refuse to run heavy work outside a SLURM allocation.

Every compute step in this project belongs in `sbatch/`. Login nodes here are
shared and heavily contended — measured at ~11% CPU per worker during an early
pilot, which made a login-node run roughly 10x slower as well as antisocial.

Call `require_slurm()` at the top of any script that calls prodigal, hmmsearch,
mmseqs or blast, or that fans out across genomes. `PN_ALLOW_LOGIN_NODE=1`
overrides it for a deliberate small test.
"""

from __future__ import annotations

import os
import socket
import sys


def in_slurm_job() -> bool:
    return bool(os.environ.get("SLURM_JOB_ID"))


def require_slurm(what: str, *, allow_env: str = "PN_ALLOW_LOGIN_NODE") -> None:
    """Exit unless inside a SLURM allocation.

    `what` names the step, so the message can point at the right sbatch script.
    """
    if in_slurm_job() or os.environ.get(allow_env):
        return
    sys.exit(
        f"\n{what} must run under SLURM, not on a login node ({socket.gethostname()}).\n"
        f"  submit:   sbatch sbatch/<script>.sbatch\n"
        f"  override: {allow_env}=1 <command>   # only for a small deliberate test\n"
    )


def describe() -> str:
    if in_slurm_job():
        return (f"SLURM job {os.environ['SLURM_JOB_ID']} on {socket.gethostname()}, "
                f"{os.environ.get('SLURM_CPUS_PER_TASK', '?')} cpus")
    return f"NOT in a SLURM job ({socket.gethostname()})"
