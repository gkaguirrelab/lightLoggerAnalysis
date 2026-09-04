"""Define the temporal kernel that maps minispect illuminance to AGC response.

The historical AGC simulation output in ``data/agc_empirical_kernels.mat``
contains two step-response analyses and their mean empirical impulse response.
This script validates that source, normalizes the mean kernel exactly as the
convolution requires, and writes the standalone derived contract
``derived/MSIlluminanceToAGCKernel.mat``.
"""

from __future__ import annotations

import argparse
import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat

from derived_io import save_derived_mat


THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[2]
DEFAULT_SOURCE_PATH = PROJECT_ROOT / "data" / "agc_empirical_kernels.mat"
DEFAULT_OUTPUT_PATH = PROJECT_ROOT / "derived" / "MSIlluminanceToAGCKernel.mat"


@dataclass(frozen=True)
class EmpiricalKernelDefinition:
    """Exact temporal kernel and provenance used by the convolution."""

    source_file: str
    time_seconds: np.ndarray
    normalized_weights: np.ndarray
    raw_weight_sum: float
    sample_interval_seconds: float
    duration_seconds: float
    step_response_tau_fits_seconds: np.ndarray

    def as_matlab_struct(self) -> dict[str, Any]:
        """Return the standalone MATLAB kernel contract."""

        return {
            "sourceFile": self.source_file,
            "definition": (
                "normalized mean empirical impulse response from simulated "
                "dark-to-bright and bright-to-dark AGC steps"
            ),
            "timeSeconds": self.time_seconds,
            "normalizedWeights": self.normalized_weights,
            "rawMeanKernelSum": self.raw_weight_sum,
            "sampleIntervalSeconds": self.sample_interval_seconds,
            "durationSeconds": self.duration_seconds,
            "sampleCount": self.normalized_weights.size,
            "stepResponseTauFitsSeconds": self.step_response_tau_fits_seconds,
            "usesEmpiricalWeightsNotExponentialFit": True,
            "convolutionPadding": "edge",
            "convolutionMode": "valid",
        }


def _project_relative_path(path: Path) -> str:
    """Return stable project-relative provenance when possible."""

    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return os.fspath(path)


def _validate_kernel_vectors(
    time_seconds: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float]:
    """Validate a uniform kernel timebase and return interval and duration."""

    if time_seconds.size < 2 or time_seconds.shape != weights.shape:
        raise ValueError(
            "Kernel time and weight fields must be matching vectors with "
            "at least two samples."
        )
    if not np.all(np.isfinite(time_seconds)) or not np.all(np.isfinite(weights)):
        raise ValueError("The empirical kernel must contain only finite values.")

    time_steps = np.diff(time_seconds)
    sample_interval_seconds = float(np.median(time_steps))
    if sample_interval_seconds <= 0 or not np.allclose(
        time_steps,
        sample_interval_seconds,
    ):
        raise ValueError("The kernel timebase must be increasing and uniform.")
    return sample_interval_seconds, float(time_seconds[-1] - time_seconds[0])


def build_kernel_definition(
    source_path: str | os.PathLike[str] = DEFAULT_SOURCE_PATH,
) -> EmpiricalKernelDefinition:
    """Build the derived definition from the historical simulation output."""

    resolved_path = Path(source_path).expanduser().resolve()
    if not resolved_path.is_file():
        raise FileNotFoundError(f"AGC kernel source does not exist:\n{resolved_path}")

    source_data: dict[str, Any] = loadmat(resolved_path, simplify_cells=True)
    time_seconds = np.asarray(
        source_data.get("commonKernelTime"),
        dtype=float,
    ).reshape(-1)
    raw_kernel = np.asarray(source_data.get("meanKernel"), dtype=float).reshape(-1)
    sample_interval_seconds, duration_seconds = _validate_kernel_vectors(
        time_seconds,
        raw_kernel,
    )

    raw_weight_sum = float(np.sum(raw_kernel))
    if not np.isfinite(raw_weight_sum) or raw_weight_sum == 0:
        raise ValueError("meanKernel must have a finite, nonzero sum.")

    results = source_data.get("results", [])
    if isinstance(results, dict):
        results = [results]
    tau_fits = np.asarray(
        [result["tauFit"] for result in results if "tauFit" in result],
        dtype=float,
    )
    return EmpiricalKernelDefinition(
        source_file=_project_relative_path(resolved_path),
        time_seconds=time_seconds,
        normalized_weights=raw_kernel / raw_weight_sum,
        raw_weight_sum=raw_weight_sum,
        sample_interval_seconds=sample_interval_seconds,
        duration_seconds=duration_seconds,
        step_response_tau_fits_seconds=tau_fits,
    )


def define_agc_kernel(
    *,
    source_path: str | os.PathLike[str] = DEFAULT_SOURCE_PATH,
    output_path: str | os.PathLike[str] = DEFAULT_OUTPUT_PATH,
) -> Path:
    """Generate and return the standalone derived AGC-kernel artifact."""

    definition = build_kernel_definition(source_path)
    return save_derived_mat(
        output_path,
        readme=(
            "Created by defineMSIlluminanceToAGCKernel.py.\n"
            "MSIlluminanceToAGCKernel contains the exact normalized empirical "
            "kernel, timebase, source-response time constants, provenance, and "
            "boundary rules used to convolve minispect illuminance before AGC "
            "lag estimation."
        ),
        variables={"MSIlluminanceToAGCKernel": definition.as_matlab_struct()},
    )


def load_agc_kernel(
    kernel_path: str | os.PathLike[str] = DEFAULT_OUTPUT_PATH,
) -> EmpiricalKernelDefinition:
    """Load and validate the standalone derived AGC-kernel contract."""

    resolved_path = Path(kernel_path).expanduser().resolve()
    if not resolved_path.is_file():
        raise FileNotFoundError(
            "Derived AGC kernel does not exist. Run "
            "defineMSIlluminanceToAGCKernel.py first:\n"
            f"{resolved_path}"
        )

    artifact = loadmat(resolved_path, simplify_cells=True)
    readme = artifact.get("README")
    if not isinstance(readme, str) or len(readme.strip()) < 20:
        raise ValueError("Derived kernel file must contain a descriptive README.")
    kernel: Any = artifact.get("MSIlluminanceToAGCKernel")
    if not isinstance(kernel, dict):
        raise KeyError(
            "Derived kernel file must contain MSIlluminanceToAGCKernel: "
            f"{resolved_path}"
        )

    time_seconds = np.asarray(kernel.get("timeSeconds"), dtype=float).reshape(-1)
    normalized_weights = np.asarray(
        kernel.get("normalizedWeights"),
        dtype=float,
    ).reshape(-1)
    sample_interval_seconds, duration_seconds = _validate_kernel_vectors(
        time_seconds,
        normalized_weights,
    )
    if not np.isclose(np.sum(normalized_weights), 1.0):
        raise ValueError("Derived kernel weights must sum to one.")
    if int(kernel.get("sampleCount")) != normalized_weights.size:
        raise ValueError("Saved kernel sampleCount does not match its weights.")

    saved_interval = float(kernel.get("sampleIntervalSeconds"))
    saved_duration = float(kernel.get("durationSeconds"))
    if not np.isclose(saved_interval, sample_interval_seconds) or not np.isclose(
        saved_duration,
        duration_seconds,
    ):
        raise ValueError("Saved kernel sampling metadata does not match its vectors.")
    if (
        kernel.get("convolutionPadding") != "edge"
        or kernel.get("convolutionMode") != "valid"
    ):
        raise ValueError("Derived kernel has unsupported convolution rules.")

    return EmpiricalKernelDefinition(
        source_file=str(kernel.get("sourceFile")),
        time_seconds=time_seconds,
        normalized_weights=normalized_weights,
        raw_weight_sum=float(kernel.get("rawMeanKernelSum")),
        sample_interval_seconds=sample_interval_seconds,
        duration_seconds=duration_seconds,
        step_response_tau_fits_seconds=np.asarray(
            kernel.get("stepResponseTauFitsSeconds", []),
            dtype=float,
        ).reshape(-1),
    )


def _build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-path", type=Path, default=DEFAULT_SOURCE_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    return parser


def main(argv: Sequence[str] | None = None) -> Path:
    """Generate the standalone derived kernel from command-line arguments."""

    arguments = _build_argument_parser().parse_args(argv)
    output_path = define_agc_kernel(
        source_path=arguments.source_path,
        output_path=arguments.output_path,
    )
    print(f"Saved derived AGC kernel to {output_path}")
    return output_path


if __name__ == "__main__":
    main()
