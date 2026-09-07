"""Shared contract for programmatically generated derived MAT artifacts."""

from __future__ import annotations

import os
from collections.abc import Mapping
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

from scipy.io import savemat, whosmat


def save_derived_mat(
    output_path: str | os.PathLike[str],
    *,
    readme: str,
    variables: Mapping[str, Any],
) -> Path:
    """Atomically save top-level variables plus a required README."""

    description = readme.strip()
    if len(description) < 20:
        raise ValueError(
            "README must provide a substantive description of the derived data."
        )
    if not variables:
        raise ValueError("At least one derived variable must be supplied.")
    if "README" in variables:
        raise ValueError("README is reserved for the description supplied separately.")

    path = Path(output_path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(
        dir=path.parent,
        prefix=f".{path.stem}-",
        suffix=".tmp.mat",
        delete=False,
    ) as temporary_file:
        temporary_path = Path(temporary_file.name)
    payload: dict[str, Any] = {"README": description, **variables}
    try:
        savemat(
            temporary_path,
            payload,
            do_compression=True,
            long_field_names=True,
            oned_as="column",
        )
        saved_variables = {name for name, _, _ in whosmat(temporary_path)}
        missing_variables = payload.keys() - saved_variables
        if missing_variables:
            missing = ", ".join(sorted(missing_variables))
            raise RuntimeError(
                f"Temporary MAT artifact is missing required variables: {missing}"
            )
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)
    return path
