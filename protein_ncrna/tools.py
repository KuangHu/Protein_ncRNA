"""Resolve external binaries.

Everything this module shells out to (prodigal, hmmfetch, hmmsearch, and later
mafft/RNAalifold) ships in the same conda env as the interpreter running it, but
that env's bin/ is only on PATH when the env is activated. Jobs are usually
launched by absolute python path instead, so resolve against sys.executable
first and fall back to PATH.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

_ENV_BIN = Path(sys.executable).parent


def which(name: str) -> str:
    local = _ENV_BIN / name
    if local.is_file():
        return str(local)
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        f"{name} not found in {_ENV_BIN} or on PATH. "
        f"Install it into the env running this interpreter ({sys.executable})."
    )


def require(*names: str) -> dict[str, str]:
    """Resolve several binaries up front so a long job fails at startup."""
    missing, found = [], {}
    for n in names:
        try:
            found[n] = which(n)
        except FileNotFoundError:
            missing.append(n)
    if missing:
        raise FileNotFoundError(
            f"missing required tools: {', '.join(missing)} "
            f"(looked in {_ENV_BIN} and PATH)"
        )
    return found
