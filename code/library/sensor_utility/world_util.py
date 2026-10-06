"""Utility functions and constants for the world camera sensor."""

from contextlib import contextmanager
from copy import deepcopy
from importlib import metadata
from functools import partial, lru_cache
import multiprocessing
from multiprocessing.pool import Pool
from multiprocessing.shared_memory import SharedMemory

import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
import sys
import warnings
import pathlib
from typing import Literal
import pandas as pd
import mat73
from natsort import natsorted 
from tqdm.auto import tqdm
import matlab
from typing import Iterable, Iterator
import dill
import pandas as pd
from numba import njit, prange
import scipy.io
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, QhullError, cKDTree
from scipy.special import ndtr
from scipy.sparse import csr_matrix

# Define constants for the project root and the derived calibration dir that we will consult later
_PROJECT_ROOT: pathlib.Path = pathlib.Path(__file__).resolve().parents[3]
_DERIVED_CALIBRATION_DIR: pathlib.Path = _PROJECT_ROOT / "derived"

 
def _load_pyagc() -> object | None:
    """Locate and import the optional ``PyAGC`` dependency.

    The world-camera utilities can run without ``PyAGC``, but when it is
    present we use it to populate AGC state tables and allowed setting
    ranges. This helper searches a small set of project-relative library
    locations, appends each candidate directory to ``sys.path`` if needed,
    and returns the first import that succeeds.

    Args:
        None.

    Returns:
        The imported ``PyAGC`` module, or ``None`` when the library is not
        available in any of the expected locations.
    """
    candidate_dirs = (
        pathlib.Path(__file__).resolve().parents[1] / "libraries_python" / "AGC_lib",
        pathlib.Path(__file__).resolve().parents[4] / "lightLogger" / "libraries_python" / "AGC_lib",
    )
    for candidate_dir in candidate_dirs:
        candidate_dir_str = str(candidate_dir)
        if candidate_dir.exists() and candidate_dir_str not in sys.path:
            sys.path.append(candidate_dir_str)
        try:
            import PyAGC  # type: ignore
            return PyAGC
        except ImportError:
            continue
    return None


PyAGC = _load_pyagc()


def _load_agc_lib() -> object | None:
    """Load the optional compiled AGC library through PyAGC.

    Args:
        None.

    Returns:
        The loaded library, or None if PyAGC is absent or loading fails.
    """
    if PyAGC is None:
        return None

    try:
        return PyAGC.import_AGC_lib()
    except Exception:
        return None


AGC_LIB = _load_agc_lib()

# Function to load in the fielding functions per dimension. Measured and saved in MATLAB.
def _import_fielding_functions() -> dict[tuple[int, int], np.ndarray]:
    """Load the measured spatial correction map from the derived calibration file.

    Args:
        None.

    Returns:
        A dictionary mapping (rows, cols) to the float64 fielding map.
    """
    fielding_functions: dict[tuple[int, int], np.ndarray] = {}
    fielding_function = scipy.io.loadmat(
        _DERIVED_CALIBRATION_DIR / "flatFieldingFunction.mat"
    )["correctionMap"].astype(np.float64, copy=False)

    fielding_functions[fielding_function.shape] = fielding_function

    return fielding_functions


def generate_RGB_mask(original_frame: np.ndarray, marker: Literal["str", "num"]="str", visualize_results: bool=False) -> tuple[np.ndarray, object] | np.ndarray:
    """Build a Bayer-pattern lookup mask for the frame geometry.

    The world camera is modeled here as blue on even/even coordinates, red
    on odd/odd coordinates, and green on the mixed-parity sites. This
    helper converts that parity rule into a reusable 2-D mask so later code
    can select raw Bayer samples by color without recomputing the coordinate
    sets.

    Args:
        original_frame: Example frame whose first two dimensions define the
            mask shape.
        marker: Output encoding for each pixel site. ``"str"`` stores the
            literal channel labels ``"R"``, ``"G"``, and ``"B"``; ``"num"``
            stores the RGB channel indices ``0``, ``1``, and ``2``.
        visualize_results: When ``True``, render a colorized depiction of
            the Bayer layout and return it alongside the mask.

    Returns:
        A 2-D mask array, or ``(mask, figure)`` when visualization is
        requested.
    """
    # Only the frame shape is needed. A BGGR tile repeats every two rows
    # and columns, so slices fill the entire mask without coordinate lists.
    image_shape: tuple[int, int] = original_frame.shape[:2]
    mask: np.ndarray = np.empty(
        image_shape, dtype="<U1" if marker == "str" else np.float64
    )
    red_marker: str | int = "R" if marker == "str" else 0
    green_marker: str | int = "G" if marker == "str" else 1
    blue_marker: str | int = "B" if marker == "str" else 2
    mask[1::2, 1::2] = red_marker
    mask[0::2, 1::2] = green_marker
    mask[1::2, 0::2] = green_marker
    mask[0::2, 0::2] = blue_marker

    if(visualize_results is True):
        # Give each Bayer site its display color using the same mask.
        colored_image: np.ndarray = np.zeros((*image_shape, 3), dtype=np.uint8)
        for channel_index, channel_marker in enumerate((red_marker, green_marker, blue_marker)):
            colored_image[..., channel_index][mask == channel_marker] = 255

        # Show the Bayer layout on a single axis.
        fig, ax = plt.subplots(1, 1)

        ax.imshow(colored_image)
        ax.set_title(f"{original_frame.shape} Bayer RGB Pattern")
        ax.axis("off")

        # Show the plot
        plt.show()

        return mask, fig

    return mask


# World temporal offset relative to other sensors. The world is the target, 
# so this is 0 (ms)
WORLD_TIME_OFFSET: float = 0

WORLD_CAMERA_DEFAULT_MODES: list[dict] = [
    {'format': "SRGGB10_CSI2P", 'unpacked': 'SRGGB10', 'bit_depth': 10, 'size': (640, 480), 'fps': 195.77, 'crop_limits': (360, 272, 2560, 1920), 'exposure_limits': (39, 6210271, None)},
    {'format': "SRGGB10_CSI2P", 'unpacked': 'SRGGB10', 'bit_depth': 10, 'size': (1600, 1200), 'fps': 42.94, 'crop_limits': (0, 0, 3200, 2400), 'exposure_limits': (75, 11766829, None)},
    {'format': "SRGGB10_CSI2P", 'unpacked': 'SRGGB10', 'bit_depth': 10, 'size': (1640, 1232), 'fps': 41.85, 'crop_limits': (0, 0, 3280, 2464), 'exposure_limits': (75, 11766829, None)},
    {'format': "SRGGB10_CSI2P", 'unpacked': 'SRGGB10', 'bit_depth': 10, 'size': (1920, 1080), 'fps': 47.57, 'crop_limits': (680, 692, 1920, 1080), 'exposure_limits': (75, 11766829, None)},
    {'format': "SRGGB10_CSI2P", 'unpacked': 'SRGGB10', 'bit_depth': 10, 'size': (3280, 2464), 'fps': 21.19, 'crop_limits': (0, 0, 3280, 2464), 'exposure_limits': (75, 11766829, None)},
    {'format': "SRGGB8", 'unpacked': 'SRGGB8', 'bit_depth': 8, 'size': (640, 480), 'fps': 195.77, 'crop_limits': (360, 272, 2560, 1920), 'exposure_limits': (39, 6210271, None)},
    {'format': "SRGGB8", 'unpacked': 'SRGGB8', 'bit_depth': 8, 'size': (1600, 1200), 'fps': 85.88, 'crop_limits': (0, 0, 3200, 2400), 'exposure_limits': (37, 5883414, None)},
    {'format': "SRGGB8", 'unpacked': 'SRGGB8', 'bit_depth': 8, 'size': (1640, 1232), 'fps': 83.7, 'crop_limits': (0, 0, 3280, 2464), 'exposure_limits': (37, 5883414, None)},
    {'format': "SRGGB8", 'unpacked': 'SRGGB8', 'bit_depth': 8, 'size': (1920, 1080), 'fps': 47.57, 'crop_limits': (680, 692, 1920, 1080), 'exposure_limits': (75, 11766829, None)},
    {'format': "SRGGB8", 'unpacked': 'SRGGB8', 'bit_depth': 8, 'size': (3296, 2464), 'fps': 21.19, 'crop_limits': (0, 0, 3280, 2464), 'exposure_limits': (75, 11766829, None)},
]

WORLD_CAMERA_CUSTOM_MODES: list[dict] = [
    {'format': "SRGGB8", 'unpacked': 'SRGGB8', 'bit_depth': 8, 'size': (640, 480), 'fps': 120, 'crop_limits': (360, 272, 2560, 1920), 'exposure_limits': (39, 6210271, None)},
    {'format': "SRGGB8", 'unpacked': 'SRGGB8', 'bit_depth': 8, 'size': (640, 480), 'fps': 180, 'crop_limits': (360, 272, 2560, 1920), 'exposure_limits': (39, 6210271, None)},
    {'format': "SRGGB8", 'unpacked': 'SRGGB8', 'bit_depth': 8, 'size': (3296, 2464), 'fps': 20, 'crop_limits': (0, 0, 3280, 2464), 'exposure_limits': (75, 11766829, None)},
]

WORLD_AGC_DEFAULT_TARGET: int = 127

# These contrast tables currently contain settings measured for an AGC
# target of 127. Additional targets can be added as sibling keys.
WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x9: dict[int | float, dict[float, tuple[float, float, float]]] = {
    WORLD_AGC_DEFAULT_TARGET: {float(NDF_level):  # Define the fixed settings for this camera per NDF filter
                               (1.0, 1.0, 468.0) if NDF_level == 0
                               else (1.0, 1.0, 4599.0) if NDF_level == 1
                               else (7.757575988769531, 1.0, 8290.0) if NDF_level == 2
                               else (10.666, 3.5039764011458963, 8290.0) if NDF_level == 3
                               else (10.666, 7.2417770421497565, 8290.0) if NDF_level == 4
                               else (10.666, 10.0, 8333)
                               for NDF_level in range(7)}
}

WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x75: dict[int | float, dict[float, tuple[float, float, float]]] = {
    WORLD_AGC_DEFAULT_TARGET: {float(NDF_level):  # Define the fixed settings for this camera per NDF filter
                               (1.0, 1.0, 508.0) if NDF_level == 0
                               else (1.000e+00, 1.000e+00, 5.537e+03) if NDF_level == 1
                               else (8.82758617e+00, 1.00000000e+00, 8333) if NDF_level == 2
                               else (10.666, 3.89023162e+00, 8333) if NDF_level == 3
                               else (10.666, 7.62447626e+00, 8333) if NDF_level == 4
                               else (10.666, 10.0, 8333)
                               for NDF_level in range(7)}
}

WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x5: dict[int | float, dict[float, tuple[float, float, float]]] = {
    WORLD_AGC_DEFAULT_TARGET: {
        0.0: (1.000e+00, 1.000e+00, 7045),
        0.1: (1.000e+00, 1.000e+00, 7045),
        0.2: (1.000e+00, 1.000e+00, 7045),
        0.3: (1.000e+00, 1.000e+00, 7045),
        0.4: (1.000e+00, 1.000e+00, 7045),
        0.5: (1.000e+00, 1.000e+00, 7045),
        0.6: (1.000e+00, 1.000e+00, 7045),
        0.7: (1.000e+00, 1.000e+00, 7045),
        0.8: (1.000e+00, 1.000e+00, 7045),
        0.9: (1.000e+00, 1.000e+00, 7045),
        1.0: (1.000e+00, 1.000e+00, 7045),
        1.1: (1.000e+00, 1.000e+00, 7045),
        1.2: (1.000e+00, 1.000e+00, 7045),
        1.3: (1.000e+00, 1.000e+00, 7045),
        1.4: (1.000e+00, 1.000e+00, 7045),
        1.5: (1.000e+00, 1.000e+00, 7045),
        1.6: (1.000e+00, 1.000e+00, 7045),
        1.7: (1.000e+00, 1.000e+00, 7045),
        1.8: (1.000e+00, 1.000e+00, 7045),
    },
}

WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x25: dict[int | float, dict[float, tuple[float, float, float]]] = {
    WORLD_AGC_DEFAULT_TARGET: {float(NDF_level):  # Define the fixed settings for this camera per NDF filter
                               (1.000e+00, 1.000e+00, 1.466e+03) if NDF_level == 0  # TODO: Fill this in
                               else (1.80281687, 1.0, 8333) if NDF_level == 1
                               else (10.666, 1.58156674e+00, 8333) if NDF_level == 2
                               else (10.666, 5.96331605e+00, 8333) if NDF_level == 3
                               else (10.666, 7.97561560e+00, 8333) if NDF_level == 4
                               else (10.666, 10.0, 8333)
                               for NDF_level in range(7)}
}

WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x1: dict[int | float, dict[float, tuple[float, float, float]]] = {
    WORLD_AGC_DEFAULT_TARGET: {float(NDF_level):  # Define the fixed settings for this camera per NDF filter
                               (1.0, 1.0, 3711.0) if NDF_level == 0
                               else (4.338983058929443, 1.0, 8290.0) if NDF_level == 1
                               else (10.666, 2.9977685360728445, 8290.0) if NDF_level == 2
                               else (10.666, 7.119079983538314, 8290.0) if NDF_level == 3
                               else (10.666, 8.039178161297238, 8290.0) if NDF_level == 4
                               else (10.666, 10.0, 8333)
                               for NDF_level in range(7)}
}
# Fractional NDFs we later measured
WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x1[WORLD_AGC_DEFAULT_TARGET] |= {
    0.2: (1.0, 1.0, 3771.0),
    0.3: (1.0, 1.0, 5860.0),
    0.4: (1.0, 1.0, 5926.0),
    0.5: (1.0, 1.0, 7752.0),
    0.6: (1.14798212, 1.0, 8333),
    0.7: (3.32467532, 1.0, 8333),
    0.8: (3.08433723, 1.0, 8333),
    0.9: (5, 1.0, 8333), # TODO: Manually changed this to 5 because it was not monotonic? What happened? 
    1.1: (7.31428576, 1.0, 8333),
    1.2: (10.23999998, 1.09120502, 8333),
    1.3: (10.23999998, 1.60061024, 8333),
    1.4: (10.23999998, 1.97600036, 8333),
    1.5: (10.23999998, 2.23682353, 8333),
    1.6: (10.23999998, 3.06064632, 8333),
    1.7: (10.23999998, 4.71412073, 8333),
    1.8: (10.23999998, 3.95815762, 8333),
}

WORLD_NDF_LEVEL_SETTINGS_CONTRAST_1x0: dict[int | float, dict[float, tuple[float, float, float]]] = {
    254: {
        0.0: (1.000e+00, 1.000e+00, 3.461e+03), 
        0.1: (1.000e+00, 1.000e+00, 5.148e+03), 
        0.2: (1.92481208e+00, 1.00000000e+00, 8333),
        0.3: (2.84444451e+00, 1.00000000e+00, 8333), 
        0.4: (3.1219511e+00, 1.0000000e+00, 8333),
        0.5: (3.93846154e+00, 1.00000000e+00, 8333), 
        0.6: (4.83018875e+00, 1.00000000e+00, 8333), 
        0.7: (7.75757599e+00, 1.00000000e+00, 8333),
        0.8: (9.14285755e+00, 1.00000000e+00, 8333),
        0.9: (9.14285755e+00, 1.00000000e+00, 8333),
        1.0: (10.333, 1.3084874, 8333),
        1.1: (10.333, 1.53624698, 8333),
        1.2: (10.333, 1.4731802, 8333),
        1.3: (10.333, 1.56616807, 8333),
        1.4: (10.333, 1.61156318, 8333),
        1.5: (10.333, 1.83363244, 8333),
        1.6: (10.333, 2.30311563, 8333),
        1.7: (10.333, 3.17007629, 8333),
        1.8: (10.333, 3.8178041, 8333),
    },
}

# NDF settings by AGC contrast target
WORLD_CONTRAST_LEVEL_NDF_SETTINGS: dict[float, dict[int | float, dict[float, tuple[float, float, float]]]] = {
    0.1: WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x1,
    0.25: WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x25,
    0.5: WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x5,
    0.75: WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x75,
    0.9: WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x9,
    1.0: WORLD_NDF_LEVEL_SETTINGS_CONTRAST_1x0,
}

# Maintain the historical default settings table used elsewhere in the MATLAB
# calibration pipeline. The 0.5 contrast table is the legacy baseline.
WORLD_NDF_LEVEL_SETTINGS: dict[float, tuple[float, float, float]] = WORLD_NDF_LEVEL_SETTINGS_CONTRAST_0x5[WORLD_AGC_DEFAULT_TARGET]
WORLD_CAM_FPS: int = 120
WORLD_FRAME_SHAPE: np.ndarray = np.array([480, 640], dtype=np.uint16)
WORLD_FRAME_DTYPE: object = np.uint8
WORLD_METADATA_DTYPE: object = np.float64
WORLD_USE_AGC: int = 1
WORLD_AGC_MODES_INT_STR: dict[int, str] = {0: "off", 1: "custom", 2: "built-in"}
WORLD_AGC_MODES_STR_INT: dict[str, int] = {val: key for key, val in WORLD_AGC_MODES_INT_STR.items()}
WORLD_SAVE_AGC_METADATA: bool = True
WORLD_AGC_CHANGE_INTERVAL: float = 0.250
WORLD_AGC_SPEED_SETTING: float = 0.95
WORLD_AGC_SETTINGS_RANGES: dict[str, np.ndarray] = PyAGC.retrieve_settings_ranges('W', AGC_LIB) if AGC_LIB is not None else {}
WORLD_AGC_DISCRETE_STATES: dict[str, dict[str, int | float]] = PyAGC.retrieve_discrete_states('W', AGC_LIB) if AGC_LIB is not None else {}


# The labels of the cols of the world AGC metdata 
# The world metadata files contain these columns with a 
# timestamp column before it
WORLD_AGC_METADATA_COLS: tuple = ("cameraAgain", "AGCDgain", "cameraExposure", "AGCAgain", "AGCExposure")

WORLD_RGB_MASK: np.ndarray = np.zeros(WORLD_FRAME_SHAPE, dtype=np.uint8)
WORLD_R_PIXELS: np.ndarray = np.array([(r, c)
                                       for r in range(WORLD_FRAME_SHAPE[0])
                                       for c in range(WORLD_FRAME_SHAPE[1])
                                       if(r % 2 != 0 and c % 2 != 0)],
                                      dtype=np.uint64)
WORLD_G_PIXELS: np.ndarray = np.array([(r, c)
                                       for r in range(WORLD_FRAME_SHAPE[0])
                                       for c in range(WORLD_FRAME_SHAPE[1])
                                       if((r % 2 == 0 and c % 2 != 0) or (r % 2 != 0 and c % 2 == 0))],
                                      dtype=np.uint64)
WORLD_B_PIXELS: np.ndarray = np.array([(r, c)
                                       for r in range(WORLD_FRAME_SHAPE[0])
                                       for c in range(WORLD_FRAME_SHAPE[1])
                                       if(r % 2 == 0 and c % 2 == 0)],
                                      dtype=np.uint64)
for idx, pixel_indices in enumerate((WORLD_R_PIXELS, WORLD_G_PIXELS, WORLD_B_PIXELS)):
    WORLD_RGB_MASK[pixel_indices[:, 0], pixel_indices[:, 1]] = idx


# Processing calibration: see code/defineWorldCameraCalibration/README.md#python-calibration-constants.

# Sensor offset and response curve (stages 1–2).
# Full-well response exponent; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_FULL_WELL_CLIPPING_EXPONENT: float = float(
    scipy.io.loadmat(_DERIVED_CALIBRATION_DIR / "nonLinearClippingExponent.mat")["clippingExponent"].item()
)

# Dark offset in 8-bit counts; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_DARK_SIGNAL: float = float(
    scipy.io.loadmat(_DERIVED_CALIBRATION_DIR / "darkSignal.mat")["darkSignal"].item()
)

# Maximum trusted inverse-response slope; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_MAX_ALLOWED_LINEARIZATION_DERIVATIVE: float = 4.0

# Spatial and Bayer-channel corrections (stages 4–5).
# Fielding maps indexed by frame shape; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_FIELDING_FUNCTIONS: dict[tuple[int, int], np.ndarray] = _import_fielding_functions()

# Shared RGB calibration data loaded once; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_RADIOMETRIC_CALIBRATION: dict[str, object] = scipy.io.loadmat(
    _DERIVED_CALIBRATION_DIR / "radiometricCorrectionRGB.mat"
)

# Channel weights in R, G, B order; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_RGB_SCALARS: np.ndarray = WORLD_RADIOMETRIC_CALIBRATION["radiometricCorrectionRGB"].astype(
    np.float64, copy=False
).reshape(-1)

# MATLAB's full Bayer correction map; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_RADIOMETRIC_CORRECTION_MAP: np.ndarray = WORLD_RADIOMETRIC_CALIBRATION["radiometricCorrectionMap"].astype(
    np.float64, copy=False
)


def generate_bayer_correction_matrix(image_shape: tuple[int, int]) -> np.ndarray:
    """Assign the calibrated RGB weight to each site in the existing Bayer mask.

    Args:
        image_shape: Frame dimensions as ``(rows, cols)``.

    Returns:
        A read-only float64 correction map with the requested shape. Its BGGR
        layout comes from ``generate_RGB_mask``, keeping one definition of
        which sensor sites measure red, green, and blue.
    """
    # We first create an empty image that is the shape we desire 
    frame_shape_template: np.ndarray = np.empty(image_shape, dtype=np.uint8)
    
    # Mark where the RGB pixels are in str format 
    rgb_mask: np.ndarray = generate_RGB_mask(frame_shape_template, marker="str")

    # Allocate a correction matrix that will store the associated weight 
    # for each of the pixels 
    correction: np.ndarray = np.empty(image_shape, dtype=np.float64)

    # Fill in the weights into the correction image
    for color, weight in zip("RGB", WORLD_RGB_SCALARS):
        correction[rgb_mask == color] = weight

    # Set this array to READ only so it cannot be further modified
    correction.setflags(write=False)
    return correction

# Cached BGGR weights for supported frames; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_BAYER_CORRECTION_MATRICES: dict[tuple[int, int], np.ndarray] = {
    (480, 640): generate_bayer_correction_matrix((480, 640)),
}

# Absolute radiance calibration (stage 6).
# Loaded camera-score fit and metadata; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_INTEGRATED_RADIANCE_CALIBRATION: dict[str, object] = scipy.io.loadmat(
    _DERIVED_CALIBRATION_DIR / "cameraScoreToIntegratedRadiance.mat"
)

# Log-space radiance slope and intercept; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_AGC_TO_RADIANCE_P: np.ndarray = np.asarray(
    WORLD_INTEGRATED_RADIANCE_CALIBRATION["agcToRadianceP"], dtype=np.float64
).reshape(-1)
if(WORLD_AGC_TO_RADIANCE_P.shape != (2,)
   or not np.all(np.isfinite(WORLD_AGC_TO_RADIANCE_P))):
    raise ValueError("agcToRadianceP must contain a finite slope and intercept")

# Product of harmonic spatial means; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_MEAN_SPATIAL_CORRECTIONS: dict[tuple[int, int], float] = {
    shape: float(1 / np.nanmean(1 / fielding)
                 / np.nanmean(1 / WORLD_RADIOMETRIC_CORRECTION_MAP))
    for shape, fielding in WORLD_FIELDING_FUNCTIONS.items()
}

# Six imputation workers are the measured baseline for 120-frame buffers.
# See tmp/benchmark_imputation_workers.json for the worker-count comparison.
WORLD_IMPUTATION_WORKERS: int = 6

# Best measured demosaicing count for sequential 3,600-frame buffers.
# See tmp/benchmark_demosaic_large_buffers.json for startup and warmed timings.
WORLD_DEMOSAIC_WORKERS: int = 16

# Interpolation used for imputation and demosaicing.
# Coordinate shear to resolve Bayer-grid ties; see README.md#python-calibration-constants in code/defineWorldCameraCalibration.
WORLD_TRIANGULATION_SHEAR: float = 1e-4



def get_world_contrast_level_ndf_settings(
    contrast_level: float,
    agc_target: int | float = WORLD_AGC_DEFAULT_TARGET,
) -> dict[float, tuple[float, float, float]]:
    """Retrieve NDF settings for a given contrast level and AGC target.

    Args:
        contrast_level: The contrast level to look up.
        agc_target: The AGC target value to look up within the contrast
            table.

    Returns:
        A dictionary mapping NDF values to tuples of calibration parameters
        for the requested contrast level and AGC target.

    Raises:
        ValueError: If ``contrast_level`` is not among the supported keys in
            ``WORLD_CONTRAST_LEVEL_NDF_SETTINGS``, or if ``agc_target`` is
            not available for that contrast level.
    """
    contrast_level = float(contrast_level)
    if contrast_level not in WORLD_CONTRAST_LEVEL_NDF_SETTINGS:
        raise ValueError(
            f"Unsupported world-camera contrast target {contrast_level}. "
            f"Expected one of {tuple(WORLD_CONTRAST_LEVEL_NDF_SETTINGS.keys())}."
        )

    agc_target_settings = WORLD_CONTRAST_LEVEL_NDF_SETTINGS[contrast_level]
    if agc_target not in agc_target_settings:
        raise ValueError(
            f"Unsupported world-camera AGC target {agc_target} for contrast {contrast_level}. "
            f"Expected one of {tuple(agc_target_settings.keys())}."
        )

    return agc_target_settings[agc_target]


def restore_settings_dict_types(sensor_mode: dict) -> None:
    """Restore the types of a sensor-mode dict after deserialization.

    Converts fields of ``sensor_mode`` back to their proper Python types
    in place (e.g. strings, ints, tuples, floats).

    Args:
        sensor_mode: A dictionary describing a camera sensor mode whose
            values may have been coerced to generic types during
            deserialization.

    Returns:
        None. The supplied sensor-mode dictionary is updated in place.
    """
    sensor_mode['format'] = str(sensor_mode['format'])
    sensor_mode['unpacked'] = str(sensor_mode['unpacked'])
    sensor_mode['bit_depth'] = int(sensor_mode['bit_depth'])
    for field in ('size', 'crop_limits'):
        sensor_mode[field] = tuple(int(val) for val in sensor_mode[field])
    sensor_mode['fps'] = float(sensor_mode['fps'])
    sensor_mode['exposure_limits'] = tuple([int(val) for val in sensor_mode['exposure_limits'][:2]] + [None])

@njit(parallel=True)
def _generate_debayered_provenance_map_numba(rows: int, cols: int) -> np.ndarray:
    """Compute raw Bayer-frame contributor coordinates for each debayered pixel.

    This is a Numba-compiled helper used by
    ``generate_debayered_provenance_map``. It fills a dense array with the
    ``(row, col)`` coordinates from the raw Bayer frame that contribute to
    each output pixel after bilinear demosaicing.

    Args:
        rows: Number of rows in the image.
        cols: Number of columns in the image.

    Returns:
        A ``np.ndarray`` of shape ``(rows, cols, 9, 2)`` with dtype
        ``uint16``, where the last two dimensions list up to 9
        ``(row, col)`` contributor coordinates per output pixel.
    """
    # Initialize the provenance map. We store a fixed-size list of 9 (row, col)
    # locations per debayered output pixel so downstream code can index the result
    # with a regular dense NumPy array.
    provenance_map: np.ndarray = np.empty((rows, cols, 9, 2), dtype=np.uint16)

    def clamp_row(row_num: int) -> int:
        """Clamp a candidate row index into valid image bounds.

        Args:
            row_num: Requested row index from a local interpolation
                neighborhood.

        Returns:
            The nearest valid row in the closed interval
            ``[0, rows - 1]``.
        """
        return min(max(row_num, 0), rows - 1)

    def clamp_col(col_num: int) -> int:
        """Clamp a candidate column index into valid image bounds.

        Args:
            col_num: Requested column index from a local interpolation
                neighborhood.

        Returns:
            The nearest valid column in the closed interval
            ``[0, cols - 1]``.
        """
        return min(max(col_num, 0), cols - 1)

    for r in prange(rows):
        for c in range(cols):
            # OpenCV's bilinear Bayer path computes interior pixels and then fills the
            # image borders by copying adjacent output pixels. To match that behavior,
            # map border output pixels back onto the interior output pixel whose
            # already-debayered value OpenCV would copy from.
            source_r: int = r
            source_c: int = c

            if(rows > 1):
                if(source_r == 0):
                    source_r = 1
                elif(source_r == rows - 1):
                    source_r = rows - 2

            if(cols > 1):
                if(source_c == 0):
                    source_c = 1
                elif(source_c == cols - 1):
                    source_c = cols - 2

            # Determine the Bayer site type according to the measured layout used
            # throughout this module:
            #   blue  at (even row, even col)
            #   green at mixed parity coordinates
            #   red   at (odd row, odd col)
            row_even: bool = (source_r % 2 == 0)
            col_even: bool = (source_c % 2 == 0)
            is_green_site: bool = row_even != col_even

            # Mirror OpenCV's bilinear demosaic support, which differs by Bayer
            # site:
            #   - Red/blue sites draw on the full local 3x3: the measured center,
            #     the four cross neighbors (green estimates), and the four diagonal
            #     neighbors (opposite-corner color estimate).
            #   - Green sites use only the cross-shaped support: the measured green
            #     center plus the two vertical and two horizontal neighbors, which
            #     carry the red and blue estimates. A green site has no diagonal
            #     contributors.
            if(is_green_site):
                # Green site: 5-tap cross support (center + up/down/left/right).
                provenance_map[r, c, 0, 0] = clamp_row(source_r)
                provenance_map[r, c, 0, 1] = clamp_col(source_c)

                provenance_map[r, c, 1, 0] = clamp_row(source_r - 1)
                provenance_map[r, c, 1, 1] = clamp_col(source_c)

                provenance_map[r, c, 2, 0] = clamp_row(source_r + 1)
                provenance_map[r, c, 2, 1] = clamp_col(source_c)

                provenance_map[r, c, 3, 0] = clamp_row(source_r)
                provenance_map[r, c, 3, 1] = clamp_col(source_c - 1)

                provenance_map[r, c, 4, 0] = clamp_row(source_r)
                provenance_map[r, c, 4, 1] = clamp_col(source_c + 1)

                # The provenance map stores a fixed 9 contributors per pixel, but a
                # green site has only 5. Pad the unused slots with the center
                # coordinate (already contributor 0) so any scan over all 9 slots
                # sees no taps beyond the true cross support.
                provenance_map[r, c, 5, 0] = clamp_row(source_r)
                provenance_map[r, c, 5, 1] = clamp_col(source_c)

                provenance_map[r, c, 6, 0] = clamp_row(source_r)
                provenance_map[r, c, 6, 1] = clamp_col(source_c)

                provenance_map[r, c, 7, 0] = clamp_row(source_r)
                provenance_map[r, c, 7, 1] = clamp_col(source_c)

                provenance_map[r, c, 8, 0] = clamp_row(source_r)
                provenance_map[r, c, 8, 1] = clamp_col(source_c)

            else:
                # Red/blue site: full 3x3 support (center + 4 cross + 4 diagonal).
                provenance_map[r, c, 0, 0] = clamp_row(source_r)
                provenance_map[r, c, 0, 1] = clamp_col(source_c)

                provenance_map[r, c, 1, 0] = clamp_row(source_r - 1)
                provenance_map[r, c, 1, 1] = clamp_col(source_c)

                provenance_map[r, c, 2, 0] = clamp_row(source_r + 1)
                provenance_map[r, c, 2, 1] = clamp_col(source_c)

                provenance_map[r, c, 3, 0] = clamp_row(source_r)
                provenance_map[r, c, 3, 1] = clamp_col(source_c - 1)

                provenance_map[r, c, 4, 0] = clamp_row(source_r)
                provenance_map[r, c, 4, 1] = clamp_col(source_c + 1)

                provenance_map[r, c, 5, 0] = clamp_row(source_r - 1)
                provenance_map[r, c, 5, 1] = clamp_col(source_c - 1)

                provenance_map[r, c, 6, 0] = clamp_row(source_r - 1)
                provenance_map[r, c, 6, 1] = clamp_col(source_c + 1)

                provenance_map[r, c, 7, 0] = clamp_row(source_r + 1)
                provenance_map[r, c, 7, 1] = clamp_col(source_c - 1)

                provenance_map[r, c, 8, 0] = clamp_row(source_r + 1)
                provenance_map[r, c, 8, 1] = clamp_col(source_c + 1)

    return provenance_map


def generate_debayered_provenance_map(debayered_image: np.ndarray,
                                      pattern: Literal["RGGB"] = "RGGB"
                                    ) -> np.ndarray:
    """Generate a provenance map for a debayered frame.

    Returns an array of shape ``(rows, cols, 9, 2)`` whose last two
    dimensions store the raw Bayer-frame ``(row, col)`` coordinates that
    contributed to each debayered output pixel.

    This implementation is designed to match the contributor layout implied
    by ``cv2.cvtColor(..., cv2.COLOR_BayerRG2RGB)`` more closely than a
    simple clamped 3x3 neighborhood model:

        - Border output pixels are mapped to the adjacent interior output
          pixel whose debayered value OpenCV copies onto the border.
        - Red and blue sites use the full 3x3 neighborhood (center + 4 cross
          + 4 diagonal). Green sites use only OpenCV's cross-shaped support
          (center + the 2 vertical and 2 horizontal neighbors) and have no
          diagonal contributors; the unused slots of the fixed-size 9-entry
          contributor list are padded with the center coordinate.

    For ``pattern="RGGB"``, this function intentionally follows the measured
    layout used elsewhere in this module:

        - blue  at ``(even row, even col)``
        - green at mixed parity coordinates
        - red   at ``(odd row, odd col)``

    Args:
        debayered_image: A debayered image array whose first two dimensions
            ``(rows, cols)`` define the output size.
        pattern: The Bayer pattern to use. Currently only ``"RGGB"`` is
            supported.

    Returns:
        A ``np.ndarray`` of shape ``(rows, cols, 9, 2)`` with dtype
        ``uint16`` mapping each output pixel to its raw-frame contributors.

    Raises:
        ValueError: If ``pattern`` is not ``"RGGB"``.
    """
    if(pattern != "RGGB"):
        raise ValueError(f"Unsupported Bayer pattern for provenance mapping: {pattern}")

    # Get the image size 
    rows, cols = debayered_image.shape[:2]

    # Generate the contributor map with a compiled helper. We cast to uint16 at the
    # boundary to preserve the previous compact storage format for typical world
    # camera frame sizes.
    provenance_map: np.ndarray = _generate_debayered_provenance_map_numba(rows, cols)

    return provenance_map.astype(np.uint16, copy=False)


def debayer(image_or_video: np.ndarray,
            visualize_results: bool=False
        ) -> np.ndarray | tuple[np.ndarray, object]:
    """Debayer raw world-camera Bayer data into RGB.

    Args:
        image_or_video: A single raw Bayer frame with shape ``(rows, cols)``
            or a raw Bayer frame buffer with shape ``(frames, rows, cols)``.
            The input must already be in a dtype accepted by OpenCV's Bayer
            conversion path.
        visualize_results: When ``True``, display a before/after figure and
            return it with the debayered result. Visualization supports only
            a single ``(rows, cols)`` frame and asserts otherwise.

    Returns:
        The RGB image or frame buffer. When ``visualize_results`` is
        ``True``, returns ``(debayered, figure)``.
    """
    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert image_or_video.ndim == 2, (
            f"Debayer visualization only supports a single frame with shape "
            f"(rows, cols). Got shape {image_or_video.shape}."
        )
        unmodified_image_or_video = image_or_video.copy() 

    # Allocate the output buffer 
    # for the debayering. This is the same shame as the 
    # input buffer, but with RGB channels dimension        
    debayered: np.ndarray = np.empty(tuple(list(image_or_video.shape) + [3]), dtype=image_or_video.dtype)

    # If we passed in a single image, just generate that single image, 
    # otherwise, populate the buffer 
    if(image_or_video.ndim == 2):
        debayered = cv2.cvtColor(image_or_video, cv2.COLOR_BayerRG2RGB) 
    elif(image_or_video.ndim == 3):
        for frame_num, frame in enumerate(image_or_video):
           cv2.cvtColor(frame, cv2.COLOR_BayerRG2RGB, dst=debayered[frame_num])
    else:
        raise Exception(f"Unsupported N Dimensions: {image_or_video.ndim}. N dimensions must be 2 or 3")

    # Visualize the results if desired
    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Debayering (Before / After)", fontweight='bold', fontsize=18)

        axes[0].imshow(unmodified_image_or_video, cmap="gray")
        axes[0].set_title("Before")
        axes[0].axis("off")

        axes[1].imshow(debayered)
        axes[1].set_title("After")
        axes[1].axis("off")

        plt.tight_layout()
        plt.show()

        return debayered, fig
        
    return debayered


def debayer_image(image: np.ndarray,
                  visualize_results: bool=False,
                  dst: np.ndarray | None=None
                  ) -> np.ndarray | tuple[np.ndarray, object] | None:
    """Backward-compatible wrapper for :func:`debayer`.

    Args:
        image: A single raw Bayer frame with shape ``(rows, cols)``.
        visualize_results: When ``True``, display a before/after figure and
            return it with the debayered image. Visualization supports only a
            single ``(rows, cols)`` frame and asserts otherwise.
        dst: Optional pre-allocated RGB output array. When provided, the
            debayered image is copied into ``dst`` and ``None`` is returned.

    Returns:
        The debayered image, ``(debayered, figure)`` when visualization is
        requested, or ``None`` when ``dst`` is provided and visualization is
        disabled.
    """
    if(dst is not None):
        assert visualize_results is False, "dst cannot be used with visualize_results=True"
        dst[:] = debayer(image)
        return None

    return debayer(image, visualize_results=visualize_results)


def calculate_world_saturation_threshold(
    original_bit_depth: int = 8,
    dark_noise: float = WORLD_DARK_SIGNAL,
    clipping_exponent: float = WORLD_FULL_WELL_CLIPPING_EXPONENT,
    max_allowed_derivative: float = WORLD_MAX_ALLOWED_LINEARIZATION_DERIVATIVE,
) -> int:
    """Return Geoff's raw-count threshold for unreliable response inversion.

    The threshold is the raw sensor value at which the derivative of the
    inverse full-well response reaches ``max_allowed_derivative``. Values at
    or above it are represented as ``Inf`` and subsequently imputed.

    Args:
        original_bit_depth: Bit depth used to determine the largest raw count.
        dark_noise: Dark offset in the same count units as the raw input.
        clipping_exponent: Positive exponent of the fitted full-well response.
        max_allowed_derivative: Largest trusted inverse-response slope; must be greater than one.

    Returns:
        Integer raw-count threshold at which samples are marked as saturated.
    """
    if(original_bit_depth <= 0):
        raise ValueError("original_bit_depth must be positive")
    if(clipping_exponent <= 0):
        raise ValueError("clipping_exponent must be positive")
    if(max_allowed_derivative <= 1):
        raise ValueError("max_allowed_derivative must be greater than one")

    sensor_max: float = float(2 ** original_bit_depth - 1)
    smax: float = sensor_max - float(dark_noise)
    if(smax <= 0):
        raise ValueError(
            f"dark_noise={dark_noise} leaves no usable range for "
            f"original_bit_depth={original_bit_depth}."
        )

    exponent_ratio: float = clipping_exponent / (clipping_exponent + 1)
    threshold_above_dark: float = smax * (
        1 - max_allowed_derivative ** (-exponent_ratio)
    ) ** (1 / clipping_exponent)
    return int(np.floor(threshold_above_dark + dark_noise))


def linearize_camera_responsivity(image_or_video: np.ndarray,
                                  original_bit_depth: int = 8,
                                  dark_noise: float = WORLD_DARK_SIGNAL,
                                  clipping_exponent: float = WORLD_FULL_WELL_CLIPPING_EXPONENT,
                                  max_allowed_derivative: float = WORLD_MAX_ALLOWED_LINEARIZATION_DERIVATIVE,
                                  visualize_results: bool=False
                                  ) -> None:
    """Linearize world-camera values using the fitted full-well model.

    This uses the same equation as Stage 2 of ``reconstructionPipeline.m``.
    Negative dark-subtracted values retain their sign. Only the nonlinear
    gain uses a non-negative signal. Quantization correction belongs to Stage 1.

    Args:
        image_or_video: Writable floating-point array of camera counts. The
            values are modified in place. Pass a copy if the original counts
            need to be retained, and convert integer arrays to float first.
        original_bit_depth: Bit depth of the input image values.
        dark_noise: The measured dark offset to remove before inversion.
        clipping_exponent: The fitted soft-clipping exponent from the
            full-well calibration.
        max_allowed_derivative: Maximum inverse-response derivative. Raw
            values at or above the resulting threshold are marked ``Inf``.
        visualize_results: When ``True``, display a before/after figure. Visualization supports only
            a single ``(rows, cols)`` frame and asserts otherwise.

    Returns:
        None. The supplied array is modified in place. When visualization is
        enabled, the before/after figure is displayed without returning it.
    Raises:
        TypeError: If the input does not have a floating-point dtype.
    """

    # Linearization always writes into the supplied array. Check its dtype
    # before changing any values; integer arrays cannot hold the result.
    if not np.issubdtype(image_or_video.dtype, np.floating):
        raise TypeError(
            f"In-place linearization requires a floating-point array; got {image_or_video.dtype}."
        )

    # Keep the original values only when a before/after figure is requested.
    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert image_or_video.ndim == 2, (
            f"Camera responsivity linearization visualization only supports "
            f"a single frame with shape (rows, cols). Got shape "
            f"{image_or_video.shape}."
        )
        unmodified_image_or_video = image_or_video.copy()

    # First, let's find the max sensor value based on the bit depth
    # and also correct this for the dark noise observed in the signal
    max_sensor_value: float = float(2 ** original_bit_depth - 1)
    dark_signal: float = float(dark_noise)
    smax: float = max_sensor_value - dark_signal

    # If the dark noise is somehow so great it's larger than the signal
    # something has really gone wrong
    if(smax <= 0):
        raise ValueError(
            f"dark_noise={dark_noise} leaves no usable range for "
            f"original_bit_depth={original_bit_depth}."
        )

    saturation_threshold: int = calculate_world_saturation_threshold(
        original_bit_depth=original_bit_depth,
        dark_noise=dark_noise,
        clipping_exponent=clipping_exponent,
        max_allowed_derivative=max_allowed_derivative,
    )
    saturation_mask: np.ndarray = image_or_video >= saturation_threshold

    # Preserve signed dark-subtracted values, as in MATLAB.
    image_or_video -= dark_signal

    # Match MATLAB's yPrime ./ (1 - (yPrime ./ Smax).^n).^(1./n).
    with np.errstate(divide="ignore", invalid="ignore"):
        image_or_video[:] = image_or_video / (
            1 - (np.maximum(image_or_video, 0) / smax) ** clipping_exponent
        ) ** (1 / clipping_exponent)

    # Match reconstructionPipeline.m: values in the unstable upper part of
    # the inverse response are treated as ceiling samples even if they have
    # not reached the integer sensor maximum.
    image_or_video[saturation_mask] = np.inf

    # If visualize results is true, we will print an output of what the image looks like
    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Full Well Correction (Before / After)", fontweight='bold', fontsize=18)

        axes[0].imshow(unmodified_image_or_video, cmap="gray")
        axes[0].set_title("Before")
        axes[0].axis("off")

        axes[1].imshow(image_or_video, cmap="gray")
        axes[1].set_title("After")
        axes[1].axis("off")

        plt.tight_layout()
        plt.show()

        return None

    return None


def _scattered_interpolant(sample_x: np.ndarray,
                           sample_y: np.ndarray,
                           sample_values: np.ndarray,
                           image_shape: tuple[int, int]
                          ) -> np.ndarray:
    """Interpolate Bayer samples using linear interpolation and nearest extrapolation.

    This follows MATLAB ``scatteredInterpolant(x, y, v, 'linear', 'nearest')``,
    but does not guarantee identical results around missing-sample regions:
    MATLAB and SciPy can choose different valid Delaunay triangulations there.

    MATLAB builds the interpolant and then evaluates it over the full pixel
    grid with ``F(X, Y)``. Both steps happen here, so this returns the grid
    directly rather than a callable.

    What the MATLAB call does, and therefore what this reproduces:

    * It triangulates the scattered sample points, splitting the plane into
      triangles whose corners are samples (a Delaunay triangulation).
    * For a query inside a triangle, the value is a weighted blend of that
      triangle's three corners. That is the ``'linear'`` method.
    * A query outside the outer boundary of all the samples sits in no
      triangle, so it instead copies its single closest sample. That is the
      ``'nearest'`` extrapolation method.

    Args:
        sample_x: One-based integer column coordinate of each Bayer sample.
        sample_y: One-based integer row coordinate of each Bayer sample.
            Samples must belong to one color channel. Input order may vary;
            nearest-neighbor ties use the MATLAB Bayer-coordinate rank.
        sample_values: Value at each sample.
        image_shape: The ``(rows, cols)`` grid to evaluate over.

    Returns:
        The interpolated grid, shaped ``image_shape``.
    """
    rows, cols = image_shape

    # MATLAB queries every pixel via [X, Y] = meshgrid(1:W, 1:H).
    query_x, query_y = np.meshgrid(np.arange(1, cols + 1, dtype=np.float64),
                                   np.arange(1, rows + 1, dtype=np.float64))
    query_x = query_x.ravel()
    query_y = query_y.ravel()

    sample_x = np.asarray(sample_x, dtype=np.float64)
    sample_y = np.asarray(sample_y, dtype=np.float64)
    sample_values = np.asarray(sample_values, dtype=np.float64)

    # Shear y by a multiple of x before triangulating. Every four neighbouring
    # Bayer samples lie on a common circle, which makes a Delaunay
    # triangulation of them non-unique; the shear breaks that tie the same way
    # MATLAB breaks it on an intact lattice. It is affine, so barycentric
    # weights for a given triangle are unchanged. Missing-sample regions can
    # still select different triangles and therefore different values.
    sheared_samples: np.ndarray = np.column_stack(
        (sample_x, sample_y + WORLD_TRIANGULATION_SHEAR * sample_x)
    )
    sheared_queries: np.ndarray = np.column_stack(
        (query_x, query_y + WORLD_TRIANGULATION_SHEAR * query_x)
    )

    # 'linear' interpolates inside the convex hull of the samples.
    try:
        interpolated: np.ndarray = np.asarray(
            LinearNDInterpolator(Delaunay(sheared_samples), sample_values, fill_value=np.nan)(sheared_queries),
            dtype=np.float64,
        )
    except QhullError:
        # Collinear or otherwise degenerate samples cannot be triangulated, so
        # every query falls through to the nearest-neighbour step below.
        interpolated = np.full(query_x.size, np.nan, dtype=np.float64)

    # 'nearest' extrapolates every query that fell outside the convex hull.
    outside_hull: np.ndarray = np.isnan(interpolated)
    if(np.any(outside_hull)):
        tree: cKDTree = cKDTree(np.column_stack((sample_x, sample_y)))
        outside_points: np.ndarray = np.column_stack((query_x[outside_hull], query_y[outside_hull]))

        # Gather equally close samples, then select by MATLAB's Bayer ordering:
        # cell-position block first, followed by column, then row. This rank
        # preserves MATLAB's last-sample tie rule regardless of input order.
        sample_rows = sample_y.astype(np.int64) - 1
        sample_cols = sample_x.astype(np.int64) - 1
        cell_blocks = (sample_rows % 2) * 2 + sample_cols % 2
        matlab_rank = cell_blocks * (rows * cols) + sample_cols * rows + sample_rows
        nearest_distance, _ = tree.query(outside_points)
        tied_samples: list[list[int]] = tree.query_ball_point(
            outside_points, nearest_distance * (1 + 1e-9) + 1e-12
        )
        interpolated[outside_hull] = sample_values[[tied[np.argmax(matlab_rank[tied])] for tied in tied_samples]]

    return interpolated.reshape(rows, cols)


# ---------------------------------------------------------------------------
# How the imputation works, in plain terms
#
# The problem. A pixel that blew out sits at the sensor ceiling and a pixel
# that saw nothing sits on the floor. Either way the true radiance was lost:
# all we know is that it was "at least this bright" or "at most this dim". The
# linearization stage marks those pixels Inf and 0 respectively. This function
# replaces them with a best estimate of what they would have read.
#
# The idea. Colour channels are correlated. If a red pixel blew out but the
# green and blue around it are still valid, those neighbours say a lot about
# how bright red must have been. So the estimate is built by conditioning on
# whichever other channels survived at that pixel.
#
# The four steps below:
#
#   1. Every pixel physically measures only ONE colour, because of the Bayer
#      filter. Interpolate each channel across the whole frame so that every
#      pixel carries an estimate of all three colours. That gives us something
#      to condition on. This grid is rgb_map.
#
#   2. Fit a 3-D Gaussian (a mean vector and a 3x3 covariance) to log RGB over
#      the pixels whose three channels are all valid. This is the prior: it
#      captures what colours this particular scene tends to contain, and how
#      the channels move together. Logs are used because radiance spans orders
#      of magnitude and is far closer to Gaussian once logged.
#
#   3. Record the brightest and dimmest value actually observed in each
#      channel. A ceiling pixel must be at least as bright as the brightest
#      thing we did manage to measure; a floor pixel at most as dim as the
#      dimmest. These become the bounds in step 4.
#
#   4. For each ruined pixel, condition the prior on its surviving channels.
#      That yields a mean and variance for the missing channel. But we also
#      know the answer lies beyond the step-3 bound, so the estimate is the
#      mean of that Gaussian restricted to the far side of the bound, which
#      has a closed form. Exponentiate to get back to linear sensor units.
# ---------------------------------------------------------------------------
def _impute_pixel_values_single(radiance_map: np.ndarray,
                                bayer_pattern: str,
                                dst: np.ndarray | None=None
                               ) -> np.ndarray:
    """Replace floor and ceiling samples in one Bayer frame.

    The estimate uses the positive RGB samples in this frame to fit a
    log-space color distribution, then conditions that distribution on the
    usable neighboring colors. Pixels with valid measurements pass through.

    Args:
        radiance_map: Two-dimensional linearized Bayer frame; non-positive and infinite values are
            imputation targets.
        bayer_pattern: Bayer layout: BGGR, RGGB, GRBG, or GBRG.
        dst: Optional float64 output with the same shape. It must not share memory with the input.

    Returns:
        An independent float64 frame, or the supplied dst filled with the result.
        The input frame is not modified.
    """
    rows, cols = radiance_map.shape

    pattern = str(bayer_pattern).upper()
    if(pattern not in ("BGGR", "RGGB", "GRBG", "GBRG")):
        raise ValueError(f"Unknown Bayer pattern: {bayer_pattern}")

    # Initialize the output directly in its final storage. Buffered callers
    # supply a slice, avoiding a separate frame allocation and later stack copy.
    if(dst is None):
        dst = np.empty(radiance_map.shape, dtype=np.float64)
    elif(dst.shape != radiance_map.shape or dst.dtype != np.float64):
        raise ValueError("Imputation destination must match the frame shape and have dtype float64.")
    elif(np.shares_memory(dst, radiance_map)):
        raise ValueError("Imputation destination must not overlap the input frame.")
    np.copyto(dst, radiance_map, casting="unsafe")

    # A finite, positive frame has no samples to impute. Skip interpolation
    # and prior fitting, but preserve the function's independent float64 output.
    if(not np.any(~np.isfinite(radiance_map) | (radiance_map <= 0))):
        return dst

    if((rows, cols) == tuple(WORLD_FRAME_SHAPE) and pattern == "BGGR"):
        # Reuse row-ordered coordinates directly for sequential image access.
        # The interpolant handles MATLAB nearest-neighbor ties explicitly.
        bayer_idx = (WORLD_R_PIXELS, WORLD_G_PIXELS,
                     WORLD_B_PIXELS)
    else:
        # Preserve support for other Bayer patterns and frame dimensions.
        cell_positions = ((0, 0), (0, 1), (1, 0), (1, 1))
        bayer_idx = [
            np.array([(r, c)
                      for position, (row_parity, col_parity) in enumerate(cell_positions)
                      if(pattern[position] == channel)
                      for c in range(col_parity, cols, 2)
                      for r in range(row_parity, rows, 2)], dtype=np.uint64)
            for channel in "RGB"
        ]

    # Interpolate to get cross-channel conditioning data. After this loop every
    # pixel carries an estimate of all three colours, not only the one its own
    # Bayer position measured.
    rgb_map: np.ndarray = np.zeros((rows, cols, 3), dtype=np.float64)
    for channel in range(3):
        channel_rows: np.ndarray = bayer_idx[channel][:, 0]
        channel_cols: np.ndarray = bayer_idx[channel][:, 1]

        # Read this channel's own samples out of the raw map.
        sub_val: np.ndarray = radiance_map[channel_rows, channel_cols]

        # A ceiling sample is Inf and a floor sample is non-positive.
        inf_mask: np.ndarray = np.isinf(sub_val)
        floor_mask: np.ndarray = (sub_val <= 0)

        # Neither carries usable information, so blank both out for the fit.
        work_sub: np.ndarray = sub_val.copy()
        work_sub[inf_mask | floor_mask] = np.nan
        valid_idx: np.ndarray = ~np.isnan(work_sub)

        # MATLAB works in one-based (x=column, y=row) coordinates.
        sub_x: np.ndarray = channel_cols.astype(np.float64) + 1
        sub_y: np.ndarray = channel_rows.astype(np.float64) + 1

        # With no usable sample the channel plane is left at zero.
        if(np.any(valid_idx)):
            rgb_map[:, :, channel] = _scattered_interpolant(
                sub_x[valid_idx], sub_y[valid_idx], work_sub[valid_idx], (rows, cols)
            )

        # Restore Inf and 0 at sub-grid locations so the imputation step can
        # find them again after the interpolation has smoothed over them.
        channel_grid: np.ndarray = rgb_map[:, :, channel]
        channel_grid[channel_rows[inf_mask], channel_cols[inf_mask]] = np.inf
        channel_grid[channel_rows[floor_mask], channel_cols[floor_mask]] = 0.0
        rgb_map[:, :, channel] = channel_grid

    # Extract prior statistics. Only pixels whose three channels are all
    # present and positive describe the scene's colour distribution.
    pixels: np.ndarray = rgb_map.reshape(-1, 3)
    valid_mask: np.ndarray = np.all(np.isfinite(pixels) & (pixels > 0), axis=1)
    # A completely dark frame has no positive RGB triples. Match Geoff's
    # fallback by fitting finite triples after lifting non-positive values
    # to a small positive floor, so their logarithms remain defined.
    if not np.any(valid_mask):
        valid_mask = np.all(np.isfinite(pixels), axis=1)
        pixels[pixels <= 0] = 1e-6
    valid_pixels: np.ndarray = pixels[valid_mask]
    if not valid_pixels.size:
        raise ValueError("No finite RGB samples available for Bayesian imputation.")

    # The model is Gaussian in log space, so the prior is fit to log radiance.
    log_valid: np.ndarray = np.log(valid_pixels)
    mu: np.ndarray = np.mean(log_valid, axis=0)
    covariance: np.ndarray = (np.cov(log_valid, rowvar=False, ddof=1)
                              if len(log_valid) > 1 else np.zeros((3, 3)))

    # Record how bright a ceiling sample must be, and how dim a floor sample
    # must be, from each channel's observed range.
    s_log: np.ndarray = np.zeros(3, dtype=np.float64)
    f_log: np.ndarray = np.zeros(3, dtype=np.float64)
    for channel in range(3):
        channel_values: np.ndarray = pixels[:, channel]
        valid_c: np.ndarray = channel_values[(channel_values > 0) & ~np.isinf(channel_values)]
        if(valid_c.size == 0):
            # Fall back to a fixed range when a channel has no valid sample.
            s_log[channel] = np.log(1.0)
            f_log[channel] = np.log(1e-4)
        else:
            s_log[channel] = np.log(np.max(valid_c))
            f_log[channel] = np.log(max(np.min(valid_c), 1e-6))

    # The destination already contains the original samples; only ceiling
    # and floor samples are replaced below.
    raw_fixed: np.ndarray = dst

    # Pre-compute the log of every interpolated pixel once. log(0) is -Inf and
    # log(Inf) is Inf, so the finiteness test below rejects both.
    with np.errstate(divide="ignore", invalid="ignore"):
        log_fixed_pixels: np.ndarray = np.log(pixels.reshape(rows, cols, 3))

    # Constant factor of the Gaussian density, hoisted out of the pixel loop.
    sqrt_two_pi: float = np.sqrt(2 * np.pi)

    # Identify the replacement targets once and reuse them for all channels.
    target_mask: np.ndarray = np.isinf(radiance_map) | (radiance_map <= 0)
    for c_target in range(3):

        # Restricted to the Bayer positions belonging to this channel.
        bayer_target_mask: np.ndarray = np.zeros((rows, cols), dtype=bool)
        bayer_target_mask[bayer_idx[c_target][:, 0], bayer_idx[c_target][:, 1]] = True

        # The pixels this pass will replace. MATLAB walks these in column-major
        # order; each pixel is independent, so the order does not matter.
        active_impute_indices: np.ndarray = np.argwhere(target_mask & bayer_target_mask)

        if(active_impute_indices.size == 0):
            continue

        # Within a frame, conditioning depends only on the target color and
        # which other channels are known. Precompute the three nonempty cases
        # once per target, rather than repeating matrix inversions per pixel.
        # This requires at most nine inversions for the entire frame.
        other_channels: list[int] = [channel for channel in range(3) if channel != c_target]
        known_channel_cases: tuple[tuple[int, ...], ...] = (
            (other_channels[0],),
            (other_channels[1],),
            tuple(other_channels),
        )
        conditioning: dict[tuple[int, ...], tuple[np.ndarray, np.ndarray, float]] = {}
        for known_channels in known_channel_cases:
            known_cols: np.ndarray = np.array(known_channels, dtype=np.intp)
            mu_k: np.ndarray = mu[known_cols]
            s_k: np.ndarray = covariance[np.ix_(known_cols, known_cols)]
            s_sk: np.ndarray = covariance[c_target, known_cols]
            # Preserve the existing ridge regularization and operation order.
            s_k_inv: np.ndarray = np.linalg.inv(s_k + 1e-6 * np.eye(known_cols.size))
            coefficients: np.ndarray = s_sk @ s_k_inv
            conditional_variance: float = float(covariance[c_target, c_target] - coefficients @ s_sk.T)
            conditional_std: float = np.sqrt(max(conditional_variance, 1e-8))
            conditioning[known_channels] = (mu_k, coefficients, conditional_std)
        prior_std: float = np.sqrt(max(float(covariance[c_target, c_target]), 1e-8))

        # Group target pixels by which of the other two channels survived.
        # Cases 0, 1, 2, and 3 mean neither, first, second, or both channels.
        # Every pixel in a group shares the same conditional covariance, so
        # NumPy can calculate all conditional means in one operation.
        target_rows: np.ndarray = active_impute_indices[:, 0]
        target_cols: np.ndarray = active_impute_indices[:, 1]
        evidence: np.ndarray = log_fixed_pixels[target_rows, target_cols]
        cases: np.ndarray = (
            np.isfinite(evidence[:, other_channels[0]]).astype(np.uint8)
            + 2 * np.isfinite(evidence[:, other_channels[1]])
        )
        for case in range(4):
            selected: np.ndarray = cases == case
            if not np.any(selected):
                continue
            pixel_rows: np.ndarray = target_rows[selected]
            pixel_cols: np.ndarray = target_cols[selected]
            conditional_mean: np.ndarray
            conditional_std: float
            if case:
                known: tuple[int, ...] = tuple(
                    other_channels[i] for i in range(2) if case & (1 << i)
                )
                mu_k, coefficients, conditional_std = conditioning[known]
                conditional_mean = mu[c_target] + (evidence[selected][:, known] - mu_k) @ coefficients
            else:
                # With no usable neighboring color, use the scene-wide prior.
                conditional_std = prior_std
                conditional_mean = np.full(pixel_rows.size, mu[c_target])

            # A ceiling sample constrains the answer above the observed range;
            # a floor sample constrains it below. Apply the corresponding
            # truncated-normal expectation in log space to the whole group.
            ceiling: np.ndarray = np.isinf(radiance_map[pixel_rows, pixel_cols])
            bound: np.ndarray = np.where(ceiling, s_log[c_target], f_log[c_target])
            z_score: np.ndarray = (bound - conditional_mean) / conditional_std
            z_score = np.where(np.isnan(z_score), 0, z_score)
            tail: np.ndarray = np.where(ceiling, 1 - ndtr(z_score), ndtr(z_score))
            numerator: np.ndarray = (conditional_std / sqrt_two_pi) * np.exp(-z_score * z_score / 2)

            # Extremely unlikely tails make the division unstable. Geoff pins
            # those estimates to the bound; compute the ratio only elsewhere.
            adjustment: np.ndarray = np.zeros_like(tail)
            np.divide(numerator, tail, out=adjustment, where=tail >= 1e-15)
            expected_log: np.ndarray = np.where(
                tail < 1e-15, bound,
                conditional_mean + np.where(ceiling, adjustment, -adjustment),
            )
            raw_fixed[pixel_rows, pixel_cols] = np.exp(expected_log)

    return raw_fixed


def _impute_pixel_values_shared_worker(frame_index: int,
                                      bayer_pattern: str,
                                      shared_name: str,
                                      output_shape: tuple[int, ...]
                                     ) -> None:
    """Impute one frame in the shared buffer owned by the parent process.

    Args:
        frame_index: Index of the frame assigned to this worker.
        bayer_pattern: Bayer layout used to interpret the frame.
        shared_name: Name of the existing shared-memory allocation.
        output_shape: Shape of the shared (frames, rows, cols) float64 buffer.

    Returns:
        None. The worker replaces its frame slice in place and closes its handle;
        the parent remains responsible for releasing the shared allocation.
    """
    shared_memory = SharedMemory(name=shared_name)
    try:
        output = np.ndarray(output_shape, dtype=np.float64, buffer=shared_memory.buf)
        try:
            output[frame_index] = _impute_pixel_values_single(output[frame_index], bayer_pattern)
        finally:
            del output
    finally:
        # The parent owns the allocation and unlinks it after the pool finishes.
        shared_memory.close()


def impute_pixel_values(linearized_image_or_buffer: np.ndarray,
                        bayer_pattern: str="BGGR",
                        visualize_results: bool=False,
                        n_workers: int=WORLD_IMPUTATION_WORKERS
                       ) -> np.ndarray | tuple[np.ndarray, object]:
    """Impute floor and ceiling Bayer samples using Geoff's Bayesian model.

    This is the Python equivalent of MATLAB ``imputePixelValues``, which
    implements a Bayesian estimate of the linearized sensor value at pixels
    that sit at the ceiling (``Inf``) or on the floor (non-positive), following:

        Zhang X, Brainard DH. Estimation of saturated pixel values in digital
        color imaging. Journal of the Optical Society of America A. 2004 Dec
        1;21(12):2301-10.

    Modified to add imputation of floor values, and to consider the
    distribution of pixel values in the log transformed space.

    MATLAB handles one frame. A ``(frames, rows, cols)`` buffer is accepted
    here and each frame is modeled independently, matching repeated calls to
    the MATLAB function.

    Args:
        linearized_image_or_buffer: Nonempty float64 NumPy frame or frame buffer whose
            ceiling samples are ``Inf`` and whose floor samples are non-positive.
            Callers supply this dtype; the function does not convert the input.
        bayer_pattern: Bayer layout of the sensor.
        visualize_results: When ``True``, display a before/after figure and
            return it alongside the result. Supported for one frame only.
        n_workers: Positive integer number of processes for a frame buffer,
            default WORLD_IMPUTATION_WORKERS (6). Callers supply a valid value. Use 1
            for serial processing. Single frames are processed directly.
            Input and output use shared memory; only frame indices are sent.
            Script callers must guard their entry point with
            ``if __name__ == "__main__":`` for multiprocessing spawn.

    Returns:
        The imputed frame or buffer, or ``(imputed, figure)`` when
        visualization is requested.

    Notes:
        Benchmarked on 2026-10-05 on a Mac Studio (Mac15,14), Apple M3 Ultra,
        28-core CPU (20 performance + 8 efficiency), 256 GB unified memory,
        macOS 15.3 (arm64).
        Timed complete calls on 120 real 480 x 640 indoor frames after
        quantization correction and linearization; all frames needed imputation
        (0.068% floor pixels, no ceiling pixels). Compared 1/2/4/6/8/12 workers,
        including spawned-process startup, shared-memory copies, and cleanup.
        Six took 48.2 s, eight 47.3 s, and twelve 59.7 s in the initial screen.
        Six was selected as the baseline despite eight being slightly faster.
        Every tested result matched serial exactly. This imputation benchmark
        used 120-frame buffers, not the 3,600-frame demosaicing workload.
        Measurements: tmp/benchmark_imputation_workers.json.
    """
    # Input must either be a single image or a buffer of images.
    assert linearized_image_or_buffer.ndim in (2, 3), (
        "Imputation requires a single frame with shape (rows, cols) "
        "or a frame buffer with shape (frames, rows, cols). Got shape "
        f"{linearized_image_or_buffer.shape}."
    )

    # Visualization supports one frame. Keep its original values for the
    # before/after figure only when visualization is requested.
    unmodified_image_or_buffer: np.ndarray | None = None
    if(visualize_results is True):
        assert linearized_image_or_buffer.ndim == 2, (
            "Imputation visualization only supports a single frame "
            f"with shape (rows, cols). Got shape {linearized_image_or_buffer.shape}."
        )
        unmodified_image_or_buffer = linearized_image_or_buffer.copy()

    # Case 1: a single (rows, cols) frame. Process it directly in this
    # process, regardless of n_workers, and return a separate output array.
    if(linearized_image_or_buffer.ndim == 2):
        result: np.ndarray = _impute_pixel_values_single(linearized_image_or_buffer, bayer_pattern)

    # Case 2: a (frames, rows, cols) buffer with one worker.
    # Allocate the output once and process frames one at a time.
    elif(n_workers == 1):
        result = np.empty(linearized_image_or_buffer.shape, dtype=np.float64)
        for frame_index, frame in enumerate(linearized_image_or_buffer):
            # Write directly into this frame's output slice instead of making
            # a separate result for every frame and stacking them afterward.
            _impute_pixel_values_single(frame, bayer_pattern, dst=result[frame_index])

    # Case 3: a frame buffer with multiple workers requested.
    else:
        # Step 1: find which frames need imputation.
        # The comparison produces one True/False value per pixel. np.any over
        # rows and columns reduces that to one flag per frame. np.flatnonzero
        # then returns the indices of the flagged frames, such as [0, 3, 7].
        frames_to_impute: np.ndarray = np.flatnonzero(
            np.any(~np.isfinite(linearized_image_or_buffer) | (linearized_image_or_buffer <= 0), axis=(1, 2))
        )
        if not frames_to_impute.size:
            # There are frames in the buffer, but none need any values replaced.
            # Return a copy to keep the output independent of the caller's input.
            # No worker processes or shared memory are needed in this case.
            return linearized_image_or_buffer.copy()

        # Step 2: reserve a block of memory that all worker processes can access.
        # Ordinary NumPy arrays are not automatically shared between spawned
        # processes. SharedMemory provides the shared bytes; it does not yet
        # contain an image or know the shape of our frame buffer.
        # Both input and output are float64, so the input's byte count fits.
        shared_memory: SharedMemory = SharedMemory(create=True, size=linearized_image_or_buffer.nbytes)
        try:
            # Give those shared bytes a NumPy shape and dtype. This creates a
            # view of the shared memory, not another allocation of image data.
            shared_output: np.ndarray = np.ndarray(
                linearized_image_or_buffer.shape, dtype=np.float64, buffer=shared_memory.buf
            )
            try:
                # Step 3: initialize the shared buffer with every input frame.
                # Workers will replace only the frames that need imputation.
                # All other frames are already correct and stay as copied here.
                # The caller's original array is separate and remains unchanged.
                np.copyto(shared_output, linearized_image_or_buffer)

                # Step 4: start the workers and assign one task per flagged frame.
                # For example, three flagged frames need at most three workers,
                # even if the caller requested twelve.
                worker_count: int = min(n_workers, len(frames_to_impute))
                with multiprocessing.get_context("spawn").Pool(processes=worker_count) as pool:
                    # Each tuple below supplies the arguments for one worker call.
                    # The worker uses the shared-memory name to open this buffer,
                    # reads its assigned frame, and writes the imputed frame back
                    # into that same slice. Workers use different frame indices,
                    # so they never write over one another's results.
                    # Only these small arguments are sent through the pool;
                    # the image data stay in shared memory. chunksize=1 assigns
                    # one frame per task, and starmap waits for all tasks to finish.
                    pool.starmap(
                        _impute_pixel_values_shared_worker,
                        ((index, bayer_pattern, shared_memory.name, linearized_image_or_buffer.shape)
                         for index in frames_to_impute),
                        chunksize=1,
                    )

                # Step 5: all frames are now ready. Copy the shared result into
                # an ordinary NumPy array that owns its memory. The caller can
                # keep this result after we release the temporary shared buffer.
                result = shared_output.copy()
            finally:
                # Step 6a: discard our NumPy view before closing the shared bytes
                # it refers to. The independent result copy is not affected.
                del shared_output
        finally:
            # Step 6b: release the temporary shared memory. These finally blocks
            # also run if a worker fails, so an error does not skip cleanup.
            # close() releases this process's handle; unlink() requests removal
            # of the shared allocation. Workers close their own handles.
            try:
                shared_memory.close()
            finally:
                # Attempt removal even if closing our handle raises an error.
                shared_memory.unlink()

    # If visualize results is true, we will print an output of what the image looks like
    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Floor / Ceiling Imputation (Before / After)", fontweight="bold", fontsize=18)

        axes[0].imshow(unmodified_image_or_buffer, cmap="gray")
        axes[0].set_title("Before")
        axes[0].axis("off")

        axes[1].imshow(result, cmap="gray")
        axes[1].set_title("After")
        axes[1].axis("off")

        plt.tight_layout()
        plt.show()

        return result, fig

    return result


@lru_cache(maxsize=12)
def _demosaic_interpolation_weights(image_shape: tuple[int, int],
                                     coordinate_bytes: bytes
                                    ) -> csr_matrix:
    """Cache linear interpolation and nearest-edge weights for an intact channel.

    Args:
        image_shape: Output dimensions as (rows, cols).
        coordinate_bytes: Contiguous int64 (row, col) Bayer coordinates encoded
            as bytes, so channel geometry forms an immutable cache key.

    Returns:
        A sparse matrix mapping measured channel samples to all image pixels.
        At most twelve channel geometries are retained across calls.
    """
    rows, cols = image_shape
    coordinates: np.ndarray = np.frombuffer(coordinate_bytes, dtype=np.int64).reshape(-1, 2)
    sample_points: np.ndarray = coordinates[:, ::-1].astype(np.float64) + 1
    query_y, query_x = np.indices(image_shape, dtype=np.float64)
    query_points: np.ndarray = np.column_stack((query_x.ravel() + 1, query_y.ravel() + 1))

    # Use the same tie-breaking shear as the scattered interpolation fallback.
    sheared_samples: np.ndarray = sample_points.copy()
    sheared_queries: np.ndarray = query_points.copy()
    sheared_samples[:, 1] += WORLD_TRIANGULATION_SHEAR * sheared_samples[:, 0]
    sheared_queries[:, 1] += WORLD_TRIANGULATION_SHEAR * sheared_queries[:, 0]
    triangulation: Delaunay = Delaunay(sheared_samples)
    simplex: np.ndarray = triangulation.find_simplex(sheared_queries)
    inside: np.ndarray = np.flatnonzero(simplex >= 0)

    # Store three barycentric weights per interior pixel; reuse them every frame.
    transforms: np.ndarray = triangulation.transform[simplex[inside]]
    weights: np.ndarray = np.einsum(
        "nij,nj->ni", transforms[:, :2], sheared_queries[inside] - transforms[:, 2]
    )
    weights = np.column_stack((weights, 1 - weights.sum(axis=1)))
    vertices: np.ndarray = triangulation.simplices[simplex[inside]]

    # Edge pixels copy the nearest sample, using the existing MATLAB tie rule.
    outside: np.ndarray = np.flatnonzero(simplex < 0)
    tree: cKDTree = cKDTree(sample_points)
    distances, _ = tree.query(query_points[outside])
    ties: list[list[int]] = tree.query_ball_point(
        query_points[outside], distances * (1 + 1e-9) + 1e-12
    )
    sample_rows, sample_cols = coordinates.T
    rank: np.ndarray = ((sample_rows % 2) * 2 + sample_cols % 2) * rows * cols + sample_cols * rows + sample_rows
    nearest: np.ndarray = np.array([indices[np.argmax(rank[indices])] for indices in ties], dtype=np.int64)
    return csr_matrix(
        (np.concatenate((weights.ravel(), np.ones(outside.size))),
         (np.concatenate((np.repeat(inside, 3), outside)),
          np.concatenate((vertices.ravel(), nearest)))),
        shape=(rows * cols, len(coordinates)),
    )


def _interpolate_demosaic_samples(channel_idx: np.ndarray,
                                  sample_values: np.ndarray,
                                  image_shape: tuple[int, int]
                                 ) -> np.ndarray:
    """Apply cached Bayer weights, falling back for missing or degenerate samples.

    Args:
        channel_idx: Measured channel coordinates as (row, col) pairs.
        sample_values: Channel values or color-to-green ratios at those sites.
        image_shape: Output dimensions as (rows, cols).

    Returns:
        A float64 interpolated channel plane with nearest extrapolation.
    """
    # Normal reconstructed frames are finite: triangulate only on the first call.
    if(np.all(np.isfinite(sample_values))):
        coordinate_bytes: bytes = np.asarray(channel_idx, dtype=np.int64).tobytes()
        try:
            weights: csr_matrix = _demosaic_interpolation_weights(image_shape, coordinate_bytes)
            return (weights @ sample_values).reshape(image_shape)
        except QhullError:
            pass

    # Missing samples change the triangulation, so retain Geoff's scattered fit.
    valid: np.ndarray = ~np.isnan(sample_values)
    return _scattered_interpolant(
        channel_idx[valid, 1] + 1, channel_idx[valid, 0] + 1,
        sample_values[valid], image_shape,
    )


def _interpolate_channel(radiance_map: np.ndarray,
                         channel_idx: np.ndarray,
                         image_shape: tuple[int, int]
                        ) -> np.ndarray:
    """Fill a channel plane from its Bayer samples using linear interpolation.

    Outside the sample region, use the nearest sample. Restore infinite
    values at their original sensor sites after interpolation.

    Args:
        radiance_map: Two-dimensional Bayer radiance frame.
        channel_idx: Array of (row, col) positions belonging to this color channel.
        image_shape: Output dimensions as (rows, cols).

    Returns:
        A float64 channel plane. With no usable samples, the plane contains NaN.
    """
    rows, cols = image_shape
    channel_rows: np.ndarray = channel_idx[:, 0]
    channel_cols: np.ndarray = channel_idx[:, 1]

    # Extract values for this channel from the raw radiance map.
    sub_val: np.ndarray = radiance_map[channel_rows, channel_cols]

    # Identify which specific sub-pixel locations are Inf in the raw map.
    inf_mask_sub: np.ndarray = np.isinf(sub_val)

    # Prepare work values by turning Inf and NaN into NaN for the interpolant.
    work_sub: np.ndarray = sub_val.copy()
    work_sub[np.isinf(work_sub)] = np.nan

    # Filter out NaN values for fitting.
    valid_idx: np.ndarray = ~np.isnan(work_sub)

    # Nothing to fit, so the whole plane is undefined.
    if(not np.any(valid_idx)):
        return np.full((rows, cols), np.nan, dtype=np.float64)

    # Perform piecewise-linear interpolation with nearest extrapolation for edge robustness.
    full_grid: np.ndarray = _interpolate_demosaic_samples(channel_idx, work_sub, (rows, cols))

    # Restore Inf values at their original sub-grid locations.
    if(np.any(inf_mask_sub)):
        full_grid[channel_rows[inf_mask_sub], channel_cols[inf_mask_sub]] = np.inf

    return full_grid


def _interpolate_ratio_channel(radiance_map: np.ndarray,
                               channel_idx: np.ndarray,
                               full_grid_g: np.ndarray,
                               image_shape: tuple[int, int]
                              ) -> np.ndarray:
    """Interpolate a red or blue channel using green as a guide.

    Interpolate the measured color-to-green ratios, then multiply by the
    full green plane. This lets the color channel follow detail resolved
    by green while preserving its measured color ratios.

    Args:
        radiance_map: Two-dimensional Bayer radiance frame.
        channel_idx: Array of (row, col) positions for the red or blue channel.
        full_grid_g: Interpolated green plane, with the same spatial shape.
        image_shape: Output dimensions as (rows, cols).

    Returns:
        The reconstructed color plane, with original infinite sites restored.
        If no ratios are usable, the plane contains NaN.
    """
    rows, cols = image_shape
    channel_rows: np.ndarray = channel_idx[:, 0]
    channel_cols: np.ndarray = channel_idx[:, 1]

    # Extract values for the target colour channel (R or B).
    sub_val: np.ndarray = radiance_map[channel_rows, channel_cols]

    # Identify which specific sub-pixel locations are Inf in the raw map.
    inf_mask_sub: np.ndarray = np.isinf(sub_val)

    # Prepare work values by turning Inf into NaN for the interpolant.
    work_sub: np.ndarray = sub_val.copy()
    work_sub[np.isinf(work_sub)] = np.nan

    # Extract corresponding Green guide values at these specific locations.
    guide_val: np.ndarray = full_grid_g[channel_rows, channel_cols]

    # Compute the colour ratio (R/G or B/G) using the NaN-masked working values.
    # Add a small epsilon to the denominator to prevent division by zero.
    epsilon: float = 1e-6
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_sub: np.ndarray = work_sub / (guide_val + epsilon)

    # Filter out NaN values for fitting.
    valid_idx: np.ndarray = ~np.isnan(ratio_sub)

    # Nothing to fit, so the whole plane is undefined.
    if(not np.any(valid_idx)):
        return np.full((rows, cols), np.nan, dtype=np.float64)

    # Perform piecewise-linear interpolation of the colour ratio.
    full_ratio_grid: np.ndarray = _interpolate_demosaic_samples(channel_idx, ratio_sub, (rows, cols))

    # Calculate the final grid by multiplying the interpolated ratio by the full Green guide.
    full_grid: np.ndarray = full_ratio_grid * full_grid_g

    # Restore Inf values at their original sub-grid locations.
    if(np.any(inf_mask_sub)):
        full_grid[channel_rows[inf_mask_sub], channel_cols[inf_mask_sub]] = np.inf

    return full_grid


# ---------------------------------------------------------------------------
# How the ratio-corrected demosaicing works, in plain terms
#
# The problem. Because of the Bayer filter each pixel measured only one colour.
# Demosaicing fills in the other two so that every pixel has a full RGB triple.
#
# The naive approach interpolates each channel on its own. That blurs edges,
# because red and blue are sampled at only a quarter of the pixels each and so
# carry very little spatial detail on their own.
#
# The RCD idea. Green is sampled at half of all pixels, twice as densely as red
# or blue, so green is the sharpest record of where the edges in the scene are.
# Meanwhile the RATIO of red to green varies slowly and smoothly across a
# scene, even across an edge, because an edge usually changes brightness rather
# than hue. So instead of interpolating red directly:
#
#   1. Interpolate green across the whole frame. This is the guide.
#   2. At each red sample, form the ratio red / green.
#   3. Interpolate that smooth ratio field across the whole frame.
#   4. Multiply back by the full green guide to recover red everywhere.
#
# Red inherits its detail from green, which actually resolved it, while the
# ratio supplies the colour. Blue is handled identically.
# ---------------------------------------------------------------------------
def _demosaic_radiance_map_rcd_single(radiance_map: np.ndarray,
                                      bayer_pattern: str
                                     ) -> np.ndarray:
    """Reconstruct RGB radiance from one Bayer frame using green-guided ratios.

    Args:
        radiance_map: Two-dimensional Bayer radiance frame.
        bayer_pattern: Bayer layout: BGGR, RGGB, GRBG, or GBRG.

    Returns:
        A float64 array shaped (rows, cols, 3), in red, green, blue order.
    """
    rows, cols = radiance_map.shape

    pattern = str(bayer_pattern).upper()
    if(pattern not in ("BGGR", "RGGB", "GRBG", "GBRG")):
        raise ValueError(f"Unknown Bayer pattern: {bayer_pattern}")

    if((rows, cols) == tuple(WORLD_FRAME_SHAPE) and pattern == "BGGR"):
        # Reuse row-ordered coordinates directly for sequential image access.
        # The interpolant handles MATLAB nearest-neighbor ties explicitly.
        rgb_idx = (WORLD_R_PIXELS, WORLD_G_PIXELS,
                     WORLD_B_PIXELS)
    else:
        # Preserve support for other Bayer patterns and frame dimensions.
        cell_positions = ((0, 0), (0, 1), (1, 0), (1, 1))
        rgb_idx = [
            np.array([(r, c)
                      for position, (row_parity, col_parity) in enumerate(cell_positions)
                      if(pattern[position] == channel)
                      for c in range(col_parity, cols, 2)
                      for r in range(row_parity, rows, 2)], dtype=np.uint64)
            for channel in "RGB"
        ]

    # 1. Interpolate Green channel first to act as the spatial luminance guide.
    # Green is sampled twice as densely as red or blue, so it carries the most
    # spatial detail.
    full_grid_g: np.ndarray = _interpolate_channel(radiance_map, rgb_idx[1], (rows, cols))

    # 2. Interpolate Red and Blue channels using the Ratio-Corrected approach.
    # This interpolates the R/G and B/G ratios, then multiplies by the full Green channel.
    full_grid_r: np.ndarray = _interpolate_ratio_channel(
        radiance_map, rgb_idx[0], full_grid_g, (rows, cols)
    )
    full_grid_b: np.ndarray = _interpolate_ratio_channel(
        radiance_map, rgb_idx[2], full_grid_g, (rows, cols)
    )

    # Assemble into a rows x cols x 3 matrix (R, G, B).
    return np.stack((full_grid_r, full_grid_g, full_grid_b), axis=-1)


@contextmanager
def demosaic_worker_pool(n_workers: int=WORLD_DEMOSAIC_WORKERS) -> Iterator[Pool | None]:
    """Reuse demosaicing workers and their interpolation caches across buffers.

    Args:
        n_workers: Positive process count, default WORLD_DEMOSAIC_WORKERS (16).
            Selected from sequential 3,600-frame benchmarks including startup
            amortized over ten buffers. One yields None for serial processing.

    Yields:
        A spawned pool to pass as worker_pool to demosaic_radiance_map_rcd or
        as demosaic_pool to world_transformation_pipeline. The pool is closed
        on context exit, including when processing fails. Guard script entry
        points with ``if __name__ == "__main__":`` for multiprocessing spawn.
        Use one context around the entire chunk loop to amortize startup.

    Notes:
        Benchmarked on 2026-10-05 on a Mac Studio (Mac15,14), Apple M3 Ultra,
        28-core CPU (20 performance + 8 efficiency), 256 GB unified memory,
        macOS 15.3 (arm64).
        Timed one 3,600 x 480 x 640 float64 Bayer buffer at a time, using a
        120-frame seeded positive-radiance sequence repeated 30 times. Compared
        8/12/16/24/28 spawned workers with a persistent pool per configuration.
        Timed startup separately, then three warmed calls including shared-memory
        allocation, input/output copies, and cleanup. Sixteen was fastest tested:
        22.77 s for the first call and 10.39 s warmed median, versus 11.09 s
        warmed with twelve. Every configuration matched the serial reference
        exactly. One first call plus nine warmed medians estimates 116.3 s
        for ten buffers; this is not an end-to-end recording measurement and
        excludes file I/O and other stages. Direct calls without a pool stay
        serial to avoid repeatedly paying startup and cache-construction costs.
        Measurements: tmp/benchmark_demosaic_large_buffers.json.
    """
    if(not isinstance(n_workers, int) or isinstance(n_workers, bool) or n_workers < 1):
        raise ValueError("n_workers must be a positive integer.")
    if(n_workers == 1):
        yield None
    else:
        with multiprocessing.get_context("spawn").Pool(processes=n_workers) as pool:
            yield pool


def _demosaic_shared_worker(frame_index: int,
                             bayer_pattern: str,
                             input_name: str,
                             output_name: str,
                             input_shape: tuple[int, ...]
                            ) -> None:
    """Demosaic one frame into its exclusive shared RGB output slice.

    Args:
        frame_index: Frame assigned to this worker.
        bayer_pattern: Bayer layout used to interpret the input.
        input_name: Parent-owned shared float64 Bayer allocation.
        output_name: Parent-owned shared float64 RGB allocation.
        input_shape: Shape of the input (frames, rows, cols) buffer.

    Returns:
        None. Only this frame's output slice is modified. The parent unlinks
        both allocations after all tasks finish; workers retain cached weights.
    """
    input_memory: SharedMemory = SharedMemory(name=input_name)
    try:
        output_memory: SharedMemory = SharedMemory(name=output_name)
        try:
            source: np.ndarray = np.ndarray(input_shape, dtype=np.float64, buffer=input_memory.buf)
            output: np.ndarray = np.ndarray((*input_shape, 3), dtype=np.float64, buffer=output_memory.buf)
            try:
                output[frame_index] = _demosaic_radiance_map_rcd_single(source[frame_index], bayer_pattern)
            finally:
                del source, output
        finally:
            output_memory.close()
    finally:
        input_memory.close()


def _demosaic_shared_buffer(values: np.ndarray,
                             bayer_pattern: str,
                             worker_pool: Pool
                            ) -> np.ndarray:
    """Process a buffer with existing workers, transferring only task descriptors.

    Args:
        values: Nonempty float64 Bayer buffer shaped (frames, rows, cols).
        bayer_pattern: Bayer layout used to interpret each frame.
        worker_pool: Persistent spawned pool owned by the caller.

    Returns:
        An independent float64 RGB buffer. Shared allocations are released on
        success or failure, and no shared-memory handles escape to the caller.
    """
    input_memory: SharedMemory = SharedMemory(create=True, size=values.nbytes)
    try:
        output_memory: SharedMemory = SharedMemory(create=True, size=values.nbytes * 3)
        try:
            shared_input: np.ndarray = np.ndarray(values.shape, dtype=np.float64, buffer=input_memory.buf)
            shared_output: np.ndarray = np.ndarray((*values.shape, 3), dtype=np.float64, buffer=output_memory.buf)
            try:
                # Copy input once. Each worker writes a disjoint RGB frame.
                np.copyto(shared_input, values)
                worker_pool.starmap(
                    _demosaic_shared_worker,
                    ((index, bayer_pattern, input_memory.name, output_memory.name, values.shape)
                     for index in range(len(values))),
                    chunksize=1,
                )
                # Return ordinary NumPy storage before releasing shared memory.
                return shared_output.copy()
            finally:
                del shared_input, shared_output
        finally:
            output_memory.close()
            output_memory.unlink()
    finally:
        input_memory.close()
        input_memory.unlink()


def demosaic_radiance_map_rcd(radiance_map_or_buffer: np.ndarray,
                              bayer_pattern: str="BGGR",
                              visualize_results: bool=False,
                              worker_pool: Pool | None = None
                             ) -> np.ndarray | tuple[np.ndarray, object]:
    """Demosaic a Bayer-pattern radiance map into a 3-D RGB image.

    This is the Python equivalent of MATLAB ``demosaicRadianceMap``. It
    takes a 2-D radiance map with a specified Bayer pattern and returns a
    ``rows x cols x 3`` array containing the interpolated Red, Green, and Blue
    channels using a Ratio-Corrected Demosaicing (RCD) algorithm.

    Pixels assigned an ``Inf`` value in the input retain their ``Inf`` value in
    their respective output channel after demosaicing when usable samples exist.
    A channel with no usable samples is NaN, matching MATLAB's early return.
    Input values are preserved and output uses float64 in the same radiance units.
    Green is interpolated first; red and blue use color / (green + 1e-6).
    Intact channel geometry is cached to avoid triangulating every video frame.
    Missing samples use the scattered fallback, where SciPy and MATLAB may
    choose different valid triangles.

    MATLAB handles one map. A ``(frames, rows, cols)`` buffer is accepted here
    and each frame is demosaiced independently.

    Args:
        radiance_map_or_buffer: Bayer radiance map shaped ``(rows, cols)`` or a
            buffer shaped ``(frames, rows, cols)``.
        bayer_pattern: Bayer layout of the sensor.
        visualize_results: When ``True``, display a before/after figure and
            return it alongside the result. Supported for one frame only.
        worker_pool: Optional persistent pool from demosaic_worker_pool(n_workers).
            Buffers use shared memory when supplied; None runs serially. Single
            frames run serially. Reuse the pool across chunks to amortize startup
            and interpolation-cache construction. Input and output are copied
            once across the shared-memory boundary.

    Returns:
        ``(rows, cols, 3)`` for one map or ``(frames, rows, cols, 3)`` for a
        buffer, or that result paired with a figure when visualization is
        requested.

    Notes:
        Benchmarked on 2026-10-05 on a Mac Studio (Mac15,14), Apple M3 Ultra,
        28-core CPU (20 performance + 8 efficiency), 256 GB unified memory,
        macOS 15.3 (arm64).
        Timed one 3,600 x 480 x 640 float64 Bayer buffer at a time, using a
        120-frame seeded positive-radiance sequence repeated 30 times. Compared
        8/12/16/24/28 spawned workers with a persistent pool per configuration.
        Timed startup separately, then three warmed calls including shared-memory
        allocation, input/output copies, and cleanup. Sixteen was fastest tested:
        22.77 s for the first call and 10.39 s warmed median, versus 11.09 s
        warmed with twelve. Every configuration matched the serial reference
        exactly. One first call plus nine warmed medians estimates 116.3 s
        for ten buffers; this is not an end-to-end recording measurement and
        excludes file I/O and other stages. Direct calls without a pool stay
        serial to avoid repeatedly paying startup and cache-construction costs.
        Measurements: tmp/benchmark_demosaic_large_buffers.json.
    """
    values: np.ndarray = np.asarray(radiance_map_or_buffer, dtype=np.float64)

    # Accept one radiance map or a buffer of them, and nothing else.
    if(values.ndim not in (2, 3)):
        raise ValueError(
            "RCD demosaicing requires shape (rows, cols) or "
            f"(frames, rows, cols). Got {values.shape}."
        )
    if(visualize_results is True and values.ndim != 2):
        raise AssertionError("RCD demosaicing visualization supports a single frame only")

    # Demosaic one map directly, or every frame of a buffer independently.
    if(values.ndim == 2):
        result: np.ndarray = _demosaic_radiance_map_rcd_single(values, bayer_pattern)
    elif(worker_pool is not None and len(values) > 1):
        result = _demosaic_shared_buffer(values, bayer_pattern, worker_pool)
    else:
        # Allocate the buffer once instead of retaining frames and stacking a copy.
        result = np.empty((*values.shape, 3), dtype=np.float64)
        for frame_index, frame in enumerate(values):
            result[frame_index] = _demosaic_radiance_map_rcd_single(frame, bayer_pattern)

    # If visualize results is true, we will print an output of what the image looks like
    if(visualize_results is True):
        # Radiance is unbounded, so stretch the middle 98% into the display range.
        finite_values: np.ndarray = result[np.isfinite(result)]
        display_result: np.ndarray = result.copy()
        if(finite_values.size):
            lower, upper = np.percentile(finite_values, (1, 99))
            if(upper > lower):
                display_result = np.clip((display_result - lower) / (upper - lower), 0, 1)

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Radiance Demosaicing (Before / After)", fontweight="bold", fontsize=18)

        axes[0].imshow(values, cmap="gray")
        axes[0].set_title("Bayer radiance")
        axes[0].axis("off")

        axes[1].imshow(display_result)
        axes[1].set_title("RCD RGB radiance")
        axes[1].axis("off")

        plt.tight_layout()
        plt.show()

        return result, fig

    return result


def apply_fielding_function(image_or_video: np.ndarray, 
                            visualize_results: bool=False,
                            fielding_function: np.ndarray | None = None
                           ) -> None:
    """Apply the registered fielding correction in place.

    Args:
        image_or_video: A floating-point raw Bayer frame with shape
            ``(rows, cols)`` or frame buffer with shape
            ``(frames, rows, cols)``. The spatial frame size must exist in
            ``WORLD_FIELDING_FUNCTIONS``.
        visualize_results: When ``True``, display a before/after figure. Visualization supports
            only a single ``(rows, cols)`` frame and asserts otherwise.
        fielding_function: Optional two-dimensional correction map. If omitted, use the calibrated
            map for the frame shape.

    Returns:
        None. The supplied array is modified in place. When visualization is
        enabled, the before/after figure is displayed without returning it.
    """

    # Input must either be a single image or a buffer of images 
    assert image_or_video.ndim in (2, 3), (
        f"Fielding function requires a single frame with shape (rows, cols) "
        f"or a frame buffer with shape (frames, rows, cols). Got shape "
        f"{image_or_video.shape}."
    )

    # If we want to visualize we need to copy an unmodified version of the input
    # as this is an in-place operation
    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert image_or_video.ndim == 2, (
            f"Fielding function visualization only supports a single frame "
            f"with shape (rows, cols). Got shape {image_or_video.shape}."
        )
        unmodified_image_or_video = image_or_video.copy() 

    # If the fielding function was not passed in,
    # load it manually.
    image_shape: tuple[int, int] = image_or_video.shape[-2:]
    if(fielding_function is None):
        fielding_function: np.ndarray = WORLD_FIELDING_FUNCTIONS[image_shape]

    # Element wise multiply the two matrices together in place
    image_or_video *= fielding_function

    # Visualize the results if desired 
    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Fielding Function (Before / After)", fontweight='bold', fontsize=18)

        axes[0].imshow(unmodified_image_or_video, cmap="gray")
        axes[0].set_title("Before")
        axes[0].axis("off")

        axes[1].imshow(image_or_video, cmap="gray")
        axes[1].set_title("After")
        axes[1].axis("off")

        plt.tight_layout()
        plt.show()

        return None

    return None

# Apply the per-color weights to the color pixels of a frame.
def apply_color_correction(image_or_video: np.ndarray, 
                           visualize_results: bool=False,
                           ) -> None:
    """Apply calibrated Bayer-site RGB scaling factors in place.

    The correction constants in ``WORLD_RGB_SCALARS`` are applied to raw
    Bayer samples using a cached BGGR weight map for each image size. The
    same 2-D map broadcasts across every frame of a buffer.

    Args:
        image_or_video: A floating-point raw Bayer frame with shape
            ``(rows, cols)`` or frame buffer with shape
            ``(frames, rows, cols)``.
        visualize_results: When ``True``, display a before/after figure. Visualization supports
            only a single ``(rows, cols)`` frame and asserts otherwise.

    Returns:
        None. The supplied array is modified in place. When visualization is
        enabled, the before/after figure is displayed without returning it.
    Raises:
        ValueError: If no cached correction map exists for the frame shape.
    """

    # Input must either be a single image or a buffer of images 
    assert image_or_video.ndim in (2, 3), (
        f"Color correction requires a single raw Bayer frame with shape "
        f"(rows, cols) or a raw Bayer frame buffer with shape "
        f"(frames, rows, cols). Got shape {image_or_video.shape}."
    )

    # If we want to do visualization we need do make a copy of the original image 
    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert image_or_video.ndim == 2, (
            f"Color correction visualization only supports a single frame "
            f"with shape (rows, cols). Got shape {image_or_video.shape}."
        )
        unmodified_image_or_video = image_or_video.copy() 

    # Use the spatial shape for either a single image or a frame buffer.
    image_shape: tuple[int, int] = image_or_video.shape[-2:]

    # Each supported shape has a correction map cached when the module loads.
    if(image_shape not in WORLD_BAYER_CORRECTION_MATRICES):
        raise ValueError(f"Unsupported image shape for Bayer color correction: {image_shape}.")
    bayer_correction_matrix: np.ndarray = WORLD_BAYER_CORRECTION_MATRICES[image_shape]

    # Elementwise multiplication writes directly to the input. Broadcasting
    # avoids creating a separate weight map for each frame in a buffer.
    image_or_video *= bayer_correction_matrix

    # Visualize the results if desired 
    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Color Correction (Before / After)", fontweight='bold', fontsize=18)

        axes[0].imshow(unmodified_image_or_video, cmap="gray")
        axes[0].set_title("Before")
        axes[0].axis("off")

        axes[1].imshow(image_or_video, cmap="gray")
        axes[1].set_title("After")
        axes[1].axis("off")

        plt.tight_layout()
        plt.show()

        return None

    # Return the modified frame
    return None


def apply_color_weights(original_frame: np.ndarray,
                        visualize_results: bool=False,
                        ) -> None:
    """Backward-compatible alias for :func:`apply_color_correction`.

    Older code in this project referred to the radiometric RGB correction
    factors as "color weights." The actual implementation lives in
    ``apply_color_correction``; this wrapper preserves the older public API
    and forwards the arguments unchanged.

    Args:
        original_frame: Raw Bayer frame to correct.
        visualize_results: Whether to request the optional before/after
            display from ``apply_color_correction``.

    Returns:
        None. The supplied array is modified in place. When visualization is
        enabled, the before/after figure is displayed without returning it.
    """
    apply_color_correction(original_frame, visualize_results=visualize_results)

def apply_digital_gain(image_or_video: np.ndarray, 
                       dgain: int | float | np.ndarray, 
                       visualize_results: bool=False
                    ) -> None | tuple[np.ndarray, object]:
    """Apply world-camera digital gain in place.

    Args:
        image_or_video: A floating-point raw Bayer frame with shape
            ``(rows, cols)`` or frame buffer with shape
            ``(frames, rows, cols)``.
        dgain: A scalar gain for a single frame, or a one-dimensional array
            with one scalar gain per frame for a frame buffer.
        visualize_results: When ``True``, display a before/after figure and
            return it with the gain-corrected input array. Visualization
            supports only a single ``(rows, cols)`` frame and asserts
            otherwise.

    Returns:
        ``None`` after in-place correction, or ``(image_or_video, figure)``
        when visualization is requested.
    """
    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert image_or_video.ndim == 2, (
            f"Digital gain visualization only supports a single frame with "
            f"shape (rows, cols). Got shape {image_or_video.shape}."
        )
        unmodified_image_or_video = image_or_video.copy()

    # Multiply the frames by the dgains 
    # (Note: this does not do any allocation, so all of this stuff is very fast)
    if(image_or_video.ndim == 2):
        # When a single 2D RAW frame is passed, we expect a single 
        # numeric value for the digital gain 
        assert np.isscalar(dgain), f"When passing in a 2D image, Dgain must be a scalar"
        
        image_or_video *= dgain
    
    elif(image_or_video.ndim == 3):
        # when a frame buffer is passed we expect the dgain to 
        # be a 1D array of numeric values 
        dgain = np.asarray(dgain)
        assert dgain.ndim == 1, f"When passing in a frame buffer, Dgain must be a 1D array of numeric values, one per each frame"
        assert dgain.shape[0] == image_or_video.shape[0], (
            f"Dgain length {dgain.shape[0]} must equal the number of frames "
            f"{image_or_video.shape[0]}"
        )

        # Expand the dgain 
        dgain_expanded: np.ndarray = dgain[:, np.newaxis, np.newaxis]
        image_or_video *= dgain_expanded

    else:
        raise Exception(f"Unsupported N dimensions: {image_or_video.ndim}")

    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Digital Gain (Before / After)", fontweight='bold', fontsize=18)

        axes[0].imshow(unmodified_image_or_video, cmap="gray")
        axes[0].set_title("Before")
        axes[0].axis("off")

        axes[1].imshow(image_or_video, cmap="gray")
        axes[1].set_title("After")
        axes[1].axis("off")

        plt.tight_layout()
        plt.show()

        return image_or_video, fig

    return None

@njit(parallel=True)
def _apply_floor_ceiling_compiled(debayered_frame_buffer: np.ndarray,
                                  raw_frame_buffer: np.ndarray,
                                  debayered_provenance_map: np.ndarray,
                                  floor: int | float,
                                  ceiling: int | float,
                                  clip_floor: int | float,
                                  clip_ceiling: int | float,
                                  ) -> None:
    """Propagate raw floor and ceiling samples into RGB frames in place.

    Each output pixel is checked against its raw contributors. A floor
    contributor takes precedence when a pixel has both floor and ceiling
    contributors. Numba processes independent output rows in parallel.

    Args:
        debayered_frame_buffer: RGB output buffer shaped (frames, rows, cols, 3), modified in place.
        raw_frame_buffer: Corresponding raw Bayer buffer shaped (frames, rows, cols).
        debayered_provenance_map: Raw (row, col) contributors for each output pixel.
        floor: Raw values at or below this level count as floor samples.
        ceiling: Raw values at or above this level count as ceiling samples.
        clip_floor: Value assigned to all RGB channels when a floor contributor is found.
        clip_ceiling: Value assigned when a ceiling contributor is found without a floor
            contributor.

    Returns:
        None. The supplied RGB buffer holds the corrected values.
    """
    debayered_height: int = debayered_provenance_map.shape[0]
    debayered_width: int = debayered_provenance_map.shape[1]

    n_frames: int = raw_frame_buffer.shape[0]
    n_contributors: int = debayered_provenance_map.shape[2]

    # Parallelize over rows since rows are independent
    for output_row in prange(debayered_height):
        for output_col in range(debayered_width):
            # Iterate over every frame for this output pixel
            for frame_index in range(n_frames):
                contains_saturated_contributor: bool = False
                contains_dark_contributor: bool = False

                # Iterate over RAW Bayer contributors for this debayered pixel
                for contributor_index in range(n_contributors):

                    contributor_row: int = debayered_provenance_map[output_row, output_col, contributor_index, 0]
                    contributor_col: int = debayered_provenance_map[output_row, output_col, contributor_index, 1]

                    # Directly read RAW Bayer value without advanced indexing allocation
                    raw_value = raw_frame_buffer[frame_index, contributor_row, contributor_col]

                    # Track whether any contributor is saturated/dark
                    if(raw_value >= ceiling):
                        contains_saturated_contributor = True

                    if(raw_value <= floor):
                        contains_dark_contributor = True

                # This ended up faster than the previous NumPy vectorization approach
                # because it avoids repeatedly allocating temporary arrays from:
                #
                # raw_frame_buffer[:, contributor_rows, contributor_cols]
                #
                # which uses advanced indexing and therefore creates copies.

                # Dark clipping takes precedence over saturation clipping
                if(contains_dark_contributor):
                    debayered_frame_buffer[frame_index, output_row, output_col, 0] = clip_floor
                    debayered_frame_buffer[frame_index, output_row, output_col, 1] = clip_floor
                    debayered_frame_buffer[frame_index, output_row, output_col, 2] = clip_floor

                elif(contains_saturated_contributor):
                    debayered_frame_buffer[frame_index, output_row, output_col, 0] = clip_ceiling
                    debayered_frame_buffer[frame_index, output_row, output_col, 1] = clip_ceiling
                    debayered_frame_buffer[frame_index, output_row, output_col, 2] = clip_ceiling

    return


def apply_floor_ceiling(debayered_image_or_video: np.ndarray,
                        raw_image_or_video: np.ndarray,
                        debayered_provenance_map: np.ndarray | None=None,
                        floor_ceiling: tuple[int | float, int | float] = (0, 255), 
                        clip_values: tuple[int | float, int | float] = (0, 255), 
                        visualize_results: bool=False,
                        ) -> None | tuple[np.ndarray, object]:
    """Propagate raw Bayer floor/ceiling hits into debayered RGB pixels.

    ``debayered_provenance_map`` records which raw Bayer samples contribute
    to each debayered output pixel. This function scans those raw
    contributors for every debayered pixel. If any contributor is at or below
    ``floor``, the entire RGB output pixel is set to ``clip_values[0]``.
    Otherwise, if any contributor is at or above ``ceiling``, the entire RGB
    output pixel is set to ``clip_values[1]``. The debayered input is
    modified in place.

    Args:
        debayered_image_or_video: Debayered RGB image with shape
            ``(height, width, 3)`` or debayered RGB frame buffer with shape
            ``(frames, height, width, 3)``.
        raw_image_or_video: Raw Bayer image with shape ``(height, width)``
            or raw Bayer frame buffer with shape ``(frames, height, width)``.
        debayered_provenance_map: Optional contributor map with shape
            ``(height, width, contributors, 2)``. When omitted, it is
            generated from ``debayered_image_or_video``.
        floor_ceiling: ``(floor, ceiling)`` raw-value sentinels. Floor
            propagation takes precedence if a pixel has both dark and
            saturated contributors.
        clip_values: Output values to write for floor and ceiling
            contributors, respectively. This can include ``np.nan`` when the
            debayered input is floating point.
        visualize_results: When ``True``, display a before/after figure and
            return it with the modified debayered image. Visualization
            supports only a single debayered image.

    Returns:
        ``None`` after in-place correction, or ``(debayered_image_or_video,
        figure)`` when visualization is requested.
    """

    assert len(floor_ceiling) == 2, f"floor_ceiling must contain exactly 2 values. Received: {floor_ceiling}"
    floor, ceiling = floor_ceiling
    assert floor <= ceiling, f"floor must be <= ceiling. Received floor={floor}, ceiling={ceiling}"
    assert len(clip_values) == 2, f"clip_values must contain exactly 2 values. Received: {clip_values}"
    clip_floor, clip_ceiling = clip_values

    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert debayered_image_or_video.ndim == 3, (
            "Floor/ceiling visualization only supports a single debayered "
            f"image with shape (rows, cols, 3). Got shape {debayered_image_or_video.shape}."
        )
        unmodified_image_or_video = debayered_image_or_video.copy()


    # Normalize single-image and video inputs to the compiled buffer shapes:
    # raw:       (frames, rows, cols)
    # debayered: (frames, rows, cols, 3)
    if(debayered_image_or_video.ndim == 3):
        assert raw_image_or_video.ndim == 2, (
            "A single debayered image must be paired with a single raw Bayer "
            f"image. Got raw shape {raw_image_or_video.shape} and debayered "
            f"shape {debayered_image_or_video.shape}."
        )
        debayered_frame_buffer: np.ndarray = debayered_image_or_video[np.newaxis, ...]
        raw_frame_buffer: np.ndarray = raw_image_or_video[np.newaxis, ...]

    elif(debayered_image_or_video.ndim == 4):
        assert raw_image_or_video.ndim == 3, (
            "A debayered frame buffer must be paired with a raw Bayer frame "
            f"buffer. Got raw shape {raw_image_or_video.shape} and debayered "
            f"shape {debayered_image_or_video.shape}."
        )
        debayered_frame_buffer = debayered_image_or_video
        raw_frame_buffer = raw_image_or_video

    else:
        raise Exception(f"Unsupported debayered N dimensions: {debayered_image_or_video.ndim}")

    assert debayered_frame_buffer.shape[-1] == 3, (
        f"Debayered input must have 3 RGB channels. Got shape {debayered_image_or_video.shape}."
    )
    assert raw_frame_buffer.shape[0] == debayered_frame_buffer.shape[0], (
        f"Raw and debayered frame counts differ: {raw_frame_buffer.shape[0]} | "
        f"{debayered_frame_buffer.shape[0]}."
    )
    assert raw_frame_buffer.shape[1:3] == debayered_frame_buffer.shape[1:3], (
        f"Raw and debayered spatial shapes differ: {raw_frame_buffer.shape[1:3]} | "
        f"{debayered_frame_buffer.shape[1:3]}."
    )

    if(debayered_provenance_map is None):
        debayered_provenance_map = generate_debayered_provenance_map(debayered_frame_buffer[0])

    assert debayered_provenance_map.shape[:2] == debayered_frame_buffer.shape[1:3], (
        f"Debayered provenance map shape {debayered_provenance_map.shape[:2]} "
        f"does not match debayered frame shape {debayered_frame_buffer.shape[1:3]}."
    )
    assert debayered_provenance_map.ndim == 4 and debayered_provenance_map.shape[-1] == 2, (
        "debayered_provenance_map must have shape (rows, cols, contributors, 2). "
        f"Got shape {debayered_provenance_map.shape}."
    )

    clip_values_array: np.ndarray = np.array(clip_values, dtype=np.float64)
    if(np.any(np.isnan(clip_values_array))):
        assert np.issubdtype(debayered_frame_buffer.dtype, np.floating), (
            "clip_values contains NaN, so debayered_image_or_video must have "
            f"a floating dtype. Got {debayered_frame_buffer.dtype}."
        )

    _apply_floor_ceiling_compiled(
        debayered_frame_buffer,
        raw_frame_buffer,
        debayered_provenance_map,
        floor,
        ceiling,
        clip_floor,
        clip_ceiling,
    )

    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Floor / Ceiling Propagation (Before / After)", fontweight='bold', fontsize=18)

        axes[0].imshow(unmodified_image_or_video)
        axes[0].set_title("Before")
        axes[0].axis("off")

        axes[1].imshow(debayered_image_or_video)
        axes[1].set_title("After")
        axes[1].axis("off")

        plt.tight_layout()
        plt.show()

        return debayered_image_or_video, fig

    return None


# Embed a world-frame timestamp into the 8-bit image itself.
def embed_timestamp(original_frame: np.ndarray, timestamp: np.float64, visualize_results: bool=False) -> tuple[np.ndarray, object] | np.ndarray:
    """Embed a ``float64`` timestamp into the first 8 bytes of a frame.

    The frame is flattened, the timestamp is converted into its little-endian
    8-byte representation, and those bytes overwrite the first 8 pixels of
    the flattened buffer. The result is then reshaped back to the original
    image shape. This preserves the frame container while intentionally
    sacrificing a small number of pixel values.

    Args:
        original_frame: ``uint8`` frame that will carry the embedded
            timestamp bytes.
        timestamp: ``np.float64`` timestamp value to encode.
        visualize_results: When ``True``, display the image before and after
            the byte overwrite and return the figure as well.

    Returns:
        The modified frame, or ``(frame, figure)`` when visualization is
        requested.
    """
    # Initialize a variable for the figure handle that 
    # will be used to visualize (if desired)
    fig: object | None = None
    assert original_frame.dtype == np.uint8, f"To do proper embedding, the frame must be uint8"
    assert type(timestamp) == np.float64, f"To do proper embedding, the timestamp must be float64"

    # Allocate a copy of the image we will edit and flatten it so we can directly 
    # edit the pixels 
    embedded_frame: np.ndarray = original_frame.copy().flatten() 

    # Let's convert the timestamp to its bytes representation 
    timestamp_as_u8: np.ndarray = np.frombuffer(np.array(timestamp, dtype='<f8').tobytes(), dtype=np.uint8) # little endian bytes
    embedded_frame[:len(timestamp_as_u8)] = timestamp_as_u8
    embedded_frame = embedded_frame.reshape(original_frame.shape)

    # Visualize the results if desired 
    if(visualize_results is True):
        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Timestamp Embedding (Before / After)", fontweight='bold', fontsize=18)
        
        # Left axis will be the "before"
        # image 
        axes[0].imshow(original_frame, cmap="gray")
        axes[0].set_title('Before')

        axes[1].imshow(embedded_frame.reshape(original_frame.shape))
        axes[1].set_title("After")

        # Show the plot 
        plt.show() 

        return embedded_frame, fig

    return embedded_frame 


# Extract a world-frame timestamp from an embedded 8-bit image.
def extract_timestamp(embedded_frame: np.ndarray) -> np.float64:
    """Decode a timestamp that was embedded by :func:`embed_timestamp`.

    Args:
        embedded_frame: ``uint8`` frame whose first 8 flattened bytes store
            a little-endian ``float64`` timestamp.

    Returns:
        The recovered timestamp as ``np.float64``.
    """
    assert embedded_frame.dtype == np.uint8, f"To extract a timestamp, the frame must be uint8"
    
    # First, let's flatten the image to get the bytes as a sequence 
    flattened_image: np.ndarray = embedded_frame.flatten() 

    # Next, we will take the first 8 bytes and view them as np.float64, 
    # interpreting them as little endian 
    timestamp_bytes: np.ndarray = flattened_image[:8]

    # Convert the timestamp back to np.float64, knowing they were assigned as little endian 
    timestamp: np.float64 = timestamp_bytes.view('<f8')[0]

    return timestamp

# Convert a Bayer image to LMS color space.
def bayer_to_lms(original_image: np.ndarray, 
                 camera: Literal["IMX219", "standard"]="IMX219",
                 path_to_spectral_sensitivities: str="", 
                 visualize_results: bool=False,
                 matlab_engine: object | None=None
                ) -> tuple[object, np.ndarray] | np.ndarray:  
    # We need a Psychtoobox function to complete this (sadge)
    # namely  WlsToS 
    """Prototype for mapping a raw Bayer image into LMS cone space.

    The current implementation performs only the setup stages of that
    pipeline. It optionally starts MATLAB, loads the camera spectral
    sensitivity table, converts the wavelength axis into Psychtoolbox's
    sampling format, and asks MATLAB for human cone sensitivities. The final
    numerical image projection has not been completed yet, so the function
    currently returns ``None`` after those preparations.

    Args:
        original_image: Raw Bayer image that will eventually be converted.
        camera: Camera model label reserved for the full conversion
            pipeline.
        path_to_spectral_sensitivities: Path to a ``.mat`` or spreadsheet
            file containing camera spectral sensitivity data.
        visualize_results: Reserved for future visualization of the LMS
            output.
        matlab_engine: Optional existing MATLAB engine session. If omitted,
            one is started and the ``lightLoggerAnalysis`` project is
            activated.

    Returns:
        Currently ``None`` because the Bayer-to-LMS projection is not yet
        implemented.
    """
    if(matlab_engine is None):
        import matlab.engine
        matlab_engine = matlab.engine.start_matlab()
        matlab_engine.tbUseProject('lightLoggerAnalysis', nargout=0)

    # First, let's make a copy of the original image
    modified_image: np.ndarray = original_image.copy() 

    # Load in the dataframe and then convert to numpy and not evil pandas (because i am less fluent in it ;-;)
    table: np.ndarray | None = None
    if(path_to_spectral_sensitivities.endswith(".mat")):
        table = pd.DataFrame(mat73.loadmat(path_to_spectral_sensitivities)["T"]).to_numpy()
    else:
        table = pd.read_excel(path_to_spectral_sensitivities).to_numpy() 

    # Extract the wavelengths from the table 
    wavelengths: np.ndarray = table[:, 0]
    rgb_values: np.ndarray = table[:, (-3, -2, -1)]

    # Convert wavelengths to sampling format
    samples: np.ndarray = np.array(matlab_engine.WlsToS(matlab.double(wavelengths), nargout=1), dtype=np.float64)
    
    # Get the sensitivities for the foveal cone classes
    field_size_degrees: int = 30
    observer_age_in_years: int = 30
    pupil_diameter_mm: int = 2
    t_receptors: np.ndarray = np.array(matlab_engine.GetHumanPhotoreceptorSS(matlab.double(samples),
                                                                    {'LConeTabulatedAbsorbance2Deg', 'MConeTabulatedAbsorbance2Deg','SConeTabulatedAbsorbance2Deg'},
                                                                    *[matlab.double(item) for item in (field_size_degrees, observer_age_in_years, pupil_diameter_mm)], [], [], [], [], 
                                                                    nargout=1
                                                                   ), dtype=np.float64)

    #Create the spectrum implied by the rgb camera weights, and then project
    # % that on the receptors

    # Splice the RGB pixels from the image into another matrix  [nPixels, [R, G, B] ]

    # Multiply the result below 
   #  cone_vec: np.ndarray = np.transpose(( t_receptors * np.transpose((rgbVec @ rgb_values)) ))

    # Then put them back in their place

    """
    % Convert RGB --> LMS contrast relative to background
    background_RGB = mean(rgbSignal, 1);
    modulation_RGB = rgbSignal - background_RGB;
    modulation_LMS = cameraToCones(modulation_RGB, options.camera);
    background_LMS = cameraToCones(background_RGB, options.camera);
    
    lmsSignal = modulation_LMS ./ background_LMS;
    
    % Select a post-receptoral channel
    switch options.postreceptoralChannel
        case {'LM'}
            signal = (lmsSignal(:,1)+lmsSignal(:,2))/2;
        case {'L-M'}
            signal = (lmsSignal(:,1)-lmsSignal(:,2));
        case {'S'}
            signal = ((lmsSignal(:,3)-lmsSignal(:,1))+lmsSignal(:,2))/2;
    end
    """


    return 

# Convert an RGB image to LMS color space.
def rgb_to_lms(original_image: np.ndarray, visualize_results: bool=False) -> tuple[object, np.ndarray] | np.ndarray:


    """Placeholder for converting an RGB image into LMS space.

    Args:
        original_image: Debayered RGB image to convert.
        visualize_results: Reserved for future visualization hooks.

    Returns:
        Currently ``None`` because the conversion has not been implemented.
    """
    return 


# Calculate per-color statistics used to derive color weights.
def calculate_color_weights(sorted_calibration_measurements: dict, visualize_results: bool=False) -> np.ndarray:
    # Let's generate the bayer pattern for a 480, 640 frame 
    """Measure raw Bayer-channel means from calibration recordings.

    Despite the historical name, this function does not normalize the
    outputs into multiplicative weights. Instead, it traverses the
    ``contrast_gamma`` calibration structure, averages a selected subset of
    repeated world-camera recordings for each contrast/frequency condition,
    extracts a square ROI centered in the frame, and computes the temporal
    mean intensity of the ROI separately for the Bayer ``R``, ``G``, and
    ``B`` pixel classes.

    Args:
        sorted_calibration_measurements: Nested calibration structure whose
            inner elements contain world-camera frame stacks at
            ``measurement['W']['v']``.
        visualize_results: When ``True``, plot the ROI time-series and the
            per-channel temporal means for each condition.

    Returns:
        Array of shape ``(n_contrast_levels, n_frequencies, 3)`` containing
        the raw temporal means for ``R``, ``G``, and ``B``.
    """
    bayer_pattern: np.ndarray = generate_RGB_mask(np.zeros((480, 640), dtype=np.uint8))
    R_pixel_locations: set[tuple] = set([(y, x) for (y, x) in zip(*np.where(bayer_pattern == 'R')) ])
    G_pixel_locations: set[tuple] = set([(y, x) for (y, x) in zip(*np.where(bayer_pattern == 'G')) ])
    B_pixel_locations: set[tuple] = set([(y, x) for (y, x) in zip(*np.where(bayer_pattern == 'B')) ])

    # Let's splice out only contrast gamma NDF 0 
    contrast_gamma_NDF0: list = sorted_calibration_measurements["contrast_gamma"][0]
    
    # Then, we will iterate over the contrast levels
    RGB_weights: np.ndarray | None = None
        
    num_contrast_levels: int = len(contrast_gamma_NDF0)
    for contrast_idx in range(num_contrast_levels):
        # Then, we will extract just the first frequency (there is only 1)
        num_frequencies: int = len(contrast_gamma_NDF0[contrast_idx])
        for frequency_idx in range(num_frequencies):
            frequency_measurement = contrast_gamma_NDF0[contrast_idx][frequency_idx]

            # Then, we will go over the 3 measurements 
            avg_world_camera_v: None | np.ndarray = None
            
            # Find the min length of the measurements, because they may not all be the same length
            min_world_v_length: float = float("inf")
            num_measurements: int = len(frequency_measurement)
            for measurement_idx in range(num_measurements):
                measurement = frequency_measurement[measurement_idx][0]
                
                # Next, let's take just the world camera V values 
                world_camera_v = measurement['W']['v']
                min_world_v_length = min(min_world_v_length, len(world_camera_v))

            # TODO: In a previous measurement, a measurement was 
            #       weird, so I am just using a single out of the 3 measurements 
            #       as this weird measurement messed with the averages 
            measurements_to_consider: tuple[int] = (1, 2) # exclusive range
            start, end = measurements_to_consider
            for measurement_idx in range(start, end):
                measurement = frequency_measurement[measurement_idx][0]
                
                # Next, let's take just the world camera V values 
                world_camera_v = measurement['W']['v']

                if(avg_world_camera_v is None):
                    avg_world_camera_v = world_camera_v[:min_world_v_length, :, :].astype(np.float64)
                else:
                    avg_world_camera_v += world_camera_v[:min_world_v_length, :, :].astype(np.float64)

            avg_world_camera_v /= (end - start)

            # Next, let's splice out the target region 
            [midpt_y, midpt_x] = np.array(avg_world_camera_v.shape[1:]) // 2

            # Let's record the ROI coords wrt to entire frame 
            roi_frame_coords = [ (midpt_y + dy, midpt_x + dx)
                                for dy in range(-20, 20)
                                for dx in range(-20, 20)
                                ]
            
            roi_R_pixels = np.array([coord for coord in roi_frame_coords 
                            if coord in R_pixel_locations
                        ])


            roi_B_pixels = np.array([coord for coord in roi_frame_coords 
                            if coord in B_pixel_locations
                        ])

            roi_G_pixels = np.array([coord for coord in roi_frame_coords 
                            if coord in G_pixel_locations
                        ])

            # Now, we calculate the temporal means 
            # for each of RGB 
            temporal_means: list = []
            for means_idx, (colorname, pixels) in enumerate(zip("RGB", (roi_R_pixels, roi_G_pixels, roi_B_pixels))):
                rows: np.ndarray = pixels[:, 0]
                cols: np.ndarray = pixels[:, 1]

                roi_color_intensity = np.mean(avg_world_camera_v[:, rows, cols], axis=1) 

                roi_color_temporal_mean = float(np.mean(roi_color_intensity))
                temporal_means.append((colorname, roi_color_intensity, roi_color_temporal_mean))

            # Save the weights for this combination of contrast + frequency 
            if(RGB_weights is None):
                RGB_weights = np.zeros((num_contrast_levels, num_frequencies, 3), dtype=np.float64)
            
            RGB_weights[contrast_idx][frequency_idx][:] = [m for (c, v, m) in temporal_means]


            # If we do not want to visualize, just continue
            if(visualize_results is False):
                continue

            # Now, plot them over time for this measurement
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))
            ax_ts, ax_tm = axes  # time-series, temporal-mean

            frame_numbers: np.ndarray = np.arange(min_world_v_length)
            for (colorname, spatial_mean, temporal_mean) in temporal_means:
                ax_ts.plot(
                    frame_numbers,
                    spatial_mean,
                    color=colorname.lower(),
                    linewidth=2,
                    label=f"{colorname} ROI"
                )

            # ---- Left: time-series plot ----
            ax_ts.set_title(f"ContrastIdx: {contrast_idx} | ROI Pixel Intensities by Frame Number", fontsize=14)
            ax_ts.set_xlabel("Frame Number", fontsize=12)
            ax_ts.set_ylabel("Avg Intensity", fontsize=12)
            ax_ts.set_ylim(0, 255)
            ax_ts.set_xlim(0, min(100, min_world_v_length - 1))
            ax_ts.legend()
            ax_ts.grid(alpha=0.3)

            # ---- Right: temporal mean per channel ----
            x = np.arange(len(temporal_means))

            # make bars match RGB colors
            heights = [m for (c, v, m) in temporal_means]
            labels  = [c for (c, v, m) in temporal_means]
            bar_colors = [c.lower() for c in labels]

            bars = ax_tm.bar(x,
                            heights,
                            color=bar_colors,
                            edgecolor="black",
                            linewidth=1.2
                            )

            ax_tm.set_xticks(x)
            ax_tm.set_xticklabels(labels)

            ax_tm.set_title(f"ContrastIdx: {contrast_idx} | Temporal Mean (ROI)", fontsize=14)
            ax_tm.set_xlabel("Channel", fontsize=12)
            ax_tm.set_ylabel("Temporal Mean Intensity", fontsize=12)
            ax_tm.set_xticks(x)
            ax_tm.set_xticklabels([c.upper() for c in bar_colors])
            ax_tm.set_ylim(0, 255)
            ax_tm.grid(alpha=0.3, axis="y")

            # Optional: annotate bar values
            for xi, val in zip(x, heights):
                ax_tm.text(xi, val + 3, f"{val:.1f}", ha="center", va="bottom", fontsize=10)


            fig.tight_layout()
            plt.show()
            plt.close(fig)

    return RGB_weights

"""
Calculate the per-color weights to apply to each color in order to equalize them
to the R channel, but for a SINGLE measurement.

Expected `measurement` shape/format (same as your inner-loop usage):
    measurement['W']['v'] is a numpy array with shape (T, H, W)
"""

def calculate_color_weights_single_measurement(measurement: dict,
                                               visualize_results: bool = False,
                                               roi_half_size: int = 20 
                                               ) -> np.ndarray:
    # Generate the Bayer pattern for a 480x640 frame
    """Compute raw ROI means for one world-camera measurement.

    The function selects a square ROI around the frame center, partitions
    those ROI pixels by their Bayer color class, computes a per-frame mean
    trace for each class, and then collapses each trace to a single temporal
    mean. As with ``calculate_color_weights``, the return value is the raw
    ``[R, G, B]`` mean vector rather than a normalized set of gains.

    Args:
        measurement: Measurement dictionary containing
            ``measurement['W']['v']`` with shape ``(time, height, width)``.
        visualize_results: When ``True``, display the ROI time-series and
            temporal mean summary plots.
        roi_half_size: Half-width of the square ROI in pixels.

    Returns:
        Length-3 ``float64`` vector containing the ROI temporal means for
        the Bayer ``R``, ``G``, and ``B`` samples.
    """
    bayer_pattern: np.ndarray = generate_RGB_mask(np.zeros((480, 640), dtype=np.uint8))

    R_pixel_locations: set[tuple[int, int]] = set(zip(*np.where(bayer_pattern == "R")))
    G_pixel_locations: set[tuple[int, int]] = set(zip(*np.where(bayer_pattern == "G")))
    B_pixel_locations: set[tuple[int, int]] = set(zip(*np.where(bayer_pattern == "B")))

    # Extract world camera V values
    world_camera_v: np.ndarray = measurement["W"]["v"]

    # We'll treat this as the averaged signal 
    avg_world_camera_v: np.ndarray = world_camera_v.astype(np.float64)
    min_world_v_length: int = avg_world_camera_v.shape[0]

    # Splice out the target region (center ROI)
    midpt_y, midpt_x = np.array(avg_world_camera_v.shape[1:]) // 2

    roi_frame_coords: list[tuple] = [(midpt_y + dy, midpt_x + dx)
                                    for dy in range(-roi_half_size, roi_half_size)
                                    for dx in range(-roi_half_size, roi_half_size)
                                    ]

    roi_R_pixels: np.ndarray = np.array([coord for coord in roi_frame_coords if coord in R_pixel_locations])
    roi_G_pixels: np.ndarray = np.array([coord for coord in roi_frame_coords if coord in G_pixel_locations])
    roi_B_pixels: np.ndarray = np.array([coord for coord in roi_frame_coords if coord in B_pixel_locations])

    # Compute per-frame ROI mean (time-series) + temporal mean for each channel
    temporal_means: list[tuple[str, np.ndarray, float]] = []
    for colorname, pixels in zip("RGB", (roi_R_pixels, roi_G_pixels, roi_B_pixels)):
        rows: np.ndarray = pixels[:, 0]
        cols: np.ndarray = pixels[:, 1]

        # Per-frame mean across all ROI pixels in this channel -> shape (T,)
        roi_mean_by_frame: np.ndarray = np.mean(avg_world_camera_v[:, rows, cols], axis=1)
        roi_temporal_mean: float = float(np.mean(roi_mean_by_frame))

        temporal_means.append((colorname, roi_mean_by_frame, roi_temporal_mean))

    # Weights output: same idea as your RGB_weights[..., :] but just a single (3,) vector
    rgb_temporal_means: np.ndarray = np.array([m for (c, v, m) in temporal_means], dtype=np.float64)  # (3,)

    # If your intent is "weights to equalize to R", you typically want:
    #   wR = 1
    #   wG = meanR / meanG
    #   wB = meanR / meanB
    # But your original function returns raw per-channel means, so we keep that EXACT behavior.
    RGB_weights: np.ndarray = rgb_temporal_means

    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        ax_ts, ax_tm = axes

        frame_numbers = np.arange(min_world_v_length)

        # Left: time-series
        for (colorname, series, tmean) in temporal_means:
            ax_ts.plot(
                frame_numbers,
                series,
                color=colorname.lower(),
                linewidth=2,
                label=f"{colorname} ROI",
            )

        ax_ts.set_title("ROI Pixel Intensities by Frame Number", fontsize=14)
        ax_ts.set_xlabel("Frame Number", fontsize=12)
        ax_ts.set_ylabel("Avg Intensity", fontsize=12)
        ax_ts.set_ylim(0, 255)
        ax_ts.set_xlim(0, min(100, min_world_v_length - 1))
        ax_ts.legend()
        ax_ts.grid(alpha=0.3)

        # Right: temporal mean bars
        x = np.arange(3)
        heights = [m for (c, v, m) in temporal_means]
        labels = [c for (c, v, m) in temporal_means]
        bar_colors = [c.lower() for c in labels]

        ax_tm.bar(x, heights, color=bar_colors, edgecolor="black", linewidth=1.2)
        ax_tm.set_title("Temporal Mean (ROI)", fontsize=14)
        ax_tm.set_xlabel("Channel", fontsize=12)
        ax_tm.set_ylabel("Temporal Mean Intensity", fontsize=12)
        ax_tm.set_xticks(x)
        ax_tm.set_xticklabels(labels)
        ax_tm.set_ylim(0, 255)
        ax_tm.grid(alpha=0.3, axis="y")

        for xi, val in zip(x, heights):
            ax_tm.text(xi, val + 3, f"{val:.1f}", ha="center", va="bottom", fontsize=10)

        fig.tight_layout()
        plt.show()
        plt.close(fig)

    return RGB_weights

# Given a recording path, return all frame timestamps.
def world_timestamps_from_chunks(raw_chunks_path: str,
                                 convert_to_seconds: bool=True,
                                 verbose: bool=False,
                                 fill_missing_frames: bool=True
                                ) -> np.ndarray:
    """Return the concatenated world-camera timestamp vector.

    This is a convenience wrapper around ``world_metadata_from_chunks`` that
    extracts only the timestamp column after chunk concatenation and any gap
    interpolation.

    Args:
        raw_chunks_path: Directory containing world metadata chunks.
        convert_to_seconds: Whether to convert the stored nanosecond
            timestamps into seconds.
        verbose: Whether to display chunk-loading progress.
        fill_missing_frames: Whether to insert timestamps for gaps between
            chunks. Defaults to True to preserve compatibility with existing
            callers that relied on gaps always being filled. Set False to
            return only timestamps of physically stored frames.

    Returns:
        One-dimensional NumPy array of world-frame timestamps.
    """
    return world_metadata_from_chunks(
        raw_chunks_path,
        convert_to_seconds=convert_to_seconds,
        verbose=verbose,
        fill_missing_frames=fill_missing_frames,
    )["timestamp"].to_numpy()


def generate_real_dummy_frame_distribution(path_to_video: str) -> np.ndarray:
    """Identify the real and missing frames in a chunked world recording.

    The returned array describes the frame sequence that a constant-frame-rate video needs in order
    to preserve the timing of the captured world frames. Every captured frame is represented by
    ``True``. Timestamp gaps large enough to imply dropped frames insert one or more ``False``
    entries before the next captured frame.

    Args:
        path_to_video: Recording directory containing ``config.pkl`` and the world metadata chunks.

    Returns:
        One-dimensional Boolean array where ``True`` denotes a captured frame and ``False`` denotes
        a dummy frame needed to fill a timestamp gap.
    """
    # Accept either the GKA recording directory itself or its parent activity directory.
    gka_path: str = os.path.join(path_to_video, "GKA")
    if(not os.path.exists(os.path.join(path_to_video, "config.pkl")) and os.path.isdir(gka_path)):
        path_to_video = gka_path

    # Read the nominal frame rate from the recording config, or use the standard 120 FPS fallback when it is unavailable.
    config_filepath: str = os.path.join(path_to_video, "config.pkl")
    if(os.path.exists(config_filepath)):
        with open(config_filepath, "rb") as config_file:
            config_data: dict = dill.load(config_file)
        recording_fps: float = config_data["sensors"]["W"]["sensor_mode"]["fps"]
    else:
        warnings.warn(f"Config filepath does not exist: {config_filepath}. Assuming world-camera FPS is {WORLD_CAM_FPS}.", RuntimeWarning)
        recording_fps = WORLD_CAM_FPS

    assert recording_fps > 0, f"World camera FPS must be positive. Got: {recording_fps}"
    frame_period_seconds: float = 1 / recording_fps

    # Natural sorting keeps the timestamps in recording order when chunk numbers reach multiple digits.
    metadata_paths: list[str] = natsorted([os.path.join(path_to_video, filename) for filename in os.listdir(path_to_video) if filename.startswith("world") and "metadata" in filename])
    assert len(metadata_paths) > 0, f"0 world metadata chunks found @ {path_to_video}"

    # Load only the timestamp column from each nonempty chunk. The source timestamps are nanoseconds, so convert them to seconds to match the writer's original calculation.
    timestamp_chunks: list[np.ndarray] = []
    for metadata_path in metadata_paths:
        metadata: np.ndarray = np.load(metadata_path, mmap_mode="r")
        assert metadata.ndim == 2 and metadata.shape[1] > 0, f"World metadata must be a 2D array with a timestamp column. Path: {metadata_path}, shape: {metadata.shape}"
        if len(metadata) > 0:
            timestamp_chunks.append(np.asarray(metadata[:, 0], dtype=np.float64) / (10 ** 9)) # convert to seconds

    # An empty recording has no real or dummy frames to describe.
    if len(timestamp_chunks) == 0:
        return np.empty(0, dtype=bool)

    # Concatenating the timestamps lets the loop detect gaps both within chunks and across chunk boundaries.
    timestamps: np.ndarray = np.concatenate(timestamp_chunks)

    # Build the writing plan in timestamp order. False entries are inserted before the real frame whose timestamp revealed the gap.
    frame_distribution: list[bool] = []
    previous_timestamp: float | None = None
    for timestamp in timestamps:
        if previous_timestamp is not None:
            # Work out how many nominal frame periods fit between the two captured frames. Subtract one because the current timestamp belongs to a real frame that is appended below.
            time_between_frames: float = timestamp - previous_timestamp
            assert time_between_frames >= 0, f"World frame timestamps must be nondecreasing. Previous: {previous_timestamp}, current: {timestamp}"
            missing_frame_count: int = max(0, int((time_between_frames / frame_period_seconds) - 1))

            if missing_frame_count > 5 * recording_fps:
                warnings.warn(f"Missed {missing_frame_count} frames between world frames. This may be unnaturally large")

            # Every missing frame becomes a False entry immediately before the next real frame.
            frame_distribution.extend([False] * missing_frame_count)

        # Every timestamp came from captured metadata, so it always contributes one real frame.
        frame_distribution.append(True)
        previous_timestamp = timestamp

    return np.array(frame_distribution, dtype=bool)

def world_frame_to_visual_angle(world_frame: np.ndarray, matlab_engine: object | None=None) -> np.ndarray:
    """Map each raw world-camera pixel to azimuth and elevation in degrees.

    The frame values are not used; its spatial dimensions define the pixel
    coordinates evaluated with the calibrated fisheye model. The returned
    array has shape ``(rows, columns, 2)``, where the final axis contains
    ``[azimuth_degrees, elevation_degrees]`` for the corresponding input
    pixel.

    Args:
        world_frame: Raw Bayer ``(rows, columns)`` or spatially equivalent
            ``(rows, columns, channels)`` world-camera frame.
        matlab_engine: Optional running MATLAB Engine session. If omitted, the function starts its
            own session.

    Returns:
        A float64 array with shape ``(rows, columns, 2)`` containing the
        calibrated visual-angle coordinate at every spatial pixel.
    """
    # The calibration applies only to full-size raw world-camera coordinates.
    world_frame = np.asarray(world_frame)
    if(world_frame.ndim not in (2, 3)):
        raise ValueError(f"world_frame must have shape (rows, columns) or (rows, columns, channels). Got {world_frame.shape}.")

    rows: int = int(world_frame.shape[0])
    columns: int = int(world_frame.shape[1])
    expected_shape: tuple[int, int] = (480, 640)
    if((rows, columns) != expected_shape):
        raise ValueError(f"world_frame spatial shape must match the calibrated raw camera shape {expected_shape}. Got {(rows, columns)}.")

    # Match MATLAB's one-based [x, y] image coordinates at every pixel center.
    x_coordinates: np.ndarray | None = None
    y_coordinates: np.ndarray | None = None
    x_coordinates, y_coordinates = np.meshgrid(np.arange(1, columns + 1, dtype=np.float64), np.arange(1, rows + 1, dtype=np.float64))
    sensor_points: np.ndarray = np.column_stack((x_coordinates.reshape(-1), y_coordinates.reshape(-1)))

    # Resolve and validate the MATLAB function and calibration paths before starting MATLAB.
    project_root: pathlib.Path = pathlib.Path(__file__).resolve().parents[3]
    intrinsics_path: pathlib.Path = project_root / "derived" / "arducamB0392cameraIntrinsics.mat"
    if(not intrinsics_path.is_file()):
        raise FileNotFoundError(f"World-camera intrinsics calibration does not exist: {intrinsics_path}")

    # Load the calibrated fisheye object and call the shared MATLAB conversion.
    need_to_initialize_engine: bool = matlab_engine is None
    if(need_to_initialize_engine):
        import matlab.engine as matlab_engine_module
        matlab_engine = matlab_engine_module.start_matlab()
        matlab_engine.tbUseProject('lightLoggerAnalysis', nargout=0)

    # Generate the visual field points.
    try:
        calibration_data: dict = matlab_engine.load(os.fspath(intrinsics_path), "arducamB0392cameraIntrinsics", nargout=1)
        calibration_results: object = calibration_data["arducamB0392cameraIntrinsics"]["results"]
        fisheye_intrinsics: object = matlab_engine.getfield(calibration_results, "Intrinsics", nargout=1)
        visual_field_points: object = matlab_engine.anglesFromIntrinsics(matlab.double(sensor_points.tolist()), fisheye_intrinsics, nargout=1)
    finally:
        if(need_to_initialize_engine):
            matlab_engine.quit()

    # Reshape so they are the same shape as the world frame
    return np.asarray(visual_field_points, dtype=np.float64).reshape(rows, columns, 2)

def world_frame_visual_angle_to_steradians(world_frame_visual_angle: np.ndarray) -> np.ndarray:
    """Convert a visual-angle coordinate image into steradians per pixel.

    Args:
        world_frame_visual_angle: Float array shaped ``(rows, columns, 2)``.
            The final axis must contain ``[azimuth_degrees,
            elevation_degrees]`` as returned by
            :func:`world_frame_to_visual_angle`.

    Returns:
        A float64 array shaped ``(rows, columns)`` containing the local solid
        angle subtended by each spatial pixel in steradians.

    Notes:
        The input gives angular coordinates at pixel centers rather than
        pixel corners. Solid angle is therefore evaluated from finite
        differences of the corresponding unit-direction field. Summing the
        result approximates the full calibrated camera field of view.
    """
    # Validate the azimuth/elevation image before calculating spatial derivatives.
    visual_angle: np.ndarray = np.asarray(world_frame_visual_angle, dtype=np.float64)
    if(visual_angle.ndim != 3 or visual_angle.shape[2] != 2):
        raise ValueError(f"world_frame_visual_angle must have shape (rows, columns, 2). Got {visual_angle.shape}.")

    # Convert each azimuth/elevation coordinate into its 3-D unit viewing direction.
    azimuth: np.ndarray = np.deg2rad(visual_angle[:, :, 0])
    elevation: np.ndarray = np.deg2rad(visual_angle[:, :, 1])
    cos_elevation: np.ndarray = np.cos(elevation)
    unit_directions: np.ndarray = np.stack((cos_elevation * np.sin(azimuth), -np.sin(elevation), cos_elevation * np.cos(azimuth)), axis=-1)

    # The cross product of the row and column direction derivatives is the local spherical-area Jacobian.
    edge_order: int = 2 if visual_angle.shape[0] >= 3 and visual_angle.shape[1] >= 3 else 1
    direction_change_per_row: np.ndarray = np.gradient(unit_directions, axis=0, edge_order=edge_order)
    direction_change_per_column: np.ndarray = np.gradient(unit_directions, axis=1, edge_order=edge_order)
    steradians_per_pixel: np.ndarray = np.linalg.norm(np.cross(direction_change_per_column, direction_change_per_row), axis=-1)

    return steradians_per_pixel


def world_camera_field_of_view_steradians(matlab_engine: object | None=None) -> float:
    """Calculate the calibrated world camera's total field of view in steradians.

    This independently integrates the fisheye model over the rectangular
    sensor boundary. It does not use the per-pixel visual-angle or steradian
    maps, so it can validate their summed solid angle.

    Args:
        matlab_engine: Optional caller-owned MATLAB engine to reuse.

    Returns:
        Total calibrated camera field of view in steradians.
    """
    project_root: pathlib.Path = pathlib.Path(__file__).resolve().parents[3]
    intrinsics_path: pathlib.Path = project_root / "derived" / "arducamB0392cameraIntrinsics.mat"
    if(not intrinsics_path.is_file()):
        raise FileNotFoundError(f"World-camera intrinsics calibration does not exist: {intrinsics_path}")

    need_to_initialize_engine: bool = matlab_engine is None
    if(need_to_initialize_engine):
        import matlab.engine as matlab_engine_module
        matlab_engine = matlab_engine_module.start_matlab()
        matlab_engine.tbUseProject('lightLoggerAnalysis', nargout=0)

    try:
        calibration_data: dict = matlab_engine.load(os.fspath(intrinsics_path), "arducamB0392cameraIntrinsics", nargout=1)
        calibration_results: object = calibration_data["arducamB0392cameraIntrinsics"]["results"]
        fisheye_intrinsics: object = matlab_engine.getfield(calibration_results, "Intrinsics", nargout=1)
        solid_angle: object = matlab_engine.calculateFisheyeSolidAngle(fisheye_intrinsics, nargout=1)
    finally:
        if(need_to_initialize_engine):
            matlab_engine.quit()

    return float(solid_angle)


# Given a recording path, return all frame metadata.
def world_metadata_from_chunks(raw_chunks_path: str,
                                 convert_to_seconds: bool=True,
                                 verbose: bool=False,
                                 fill_missing_frames: bool=True
                                ) -> pd.DataFrame:
    # Find the config file 
    # of the recording. This will tell us about the FPS 
    # and how to interpolate between chunks 
    """Assemble timestamps and AGC metadata across all world chunks.

    The function reads ``config.pkl`` to recover the nominal world-camera
    frame rate, loads each naturally sorted metadata chunk, and concatenates
    them into a single table. When a timestamp gap appears between adjacent
    chunks, it estimates how many frames are missing from the nominal frame
    period and inserts synthetic rows whose timestamps span the gap while
    the AGC-setting columns are filled with ``NaN``. Set
    ``fill_missing_frames=False`` to concatenate only stored rows, preserving
    physical frame indices even when recorded camera settings contain NaNs.

    Args:
        raw_chunks_path: Directory containing ``config.pkl`` and the world
            metadata chunk files.
        convert_to_seconds: Whether to convert timestamps from nanoseconds
            since boot into seconds.
        verbose: Whether to show progress while loading the metadata chunks.
        fill_missing_frames: Whether to insert timestamped NaN rows for gaps
            between chunks. Defaults to True to preserve compatibility with
            existing callers that relied on gaps always being filled.

    Returns:
        ``pandas.DataFrame`` in frame order with a zero-based index and a
        timestamp column followed by legacy ``Again``, ``Dgain``, ``exposure``
        or the modern ``WORLD_AGC_METADATA_COLS`` camera settings.
    """
    chunks_path: str = os.path.abspath(os.path.expanduser(raw_chunks_path))
    if(not os.path.isdir(chunks_path)):
        raise FileNotFoundError(f"Raw chunks path does not exist: {chunks_path}")

    # Only gap filling needs the configured frame rate.
    if fill_missing_frames:
        # Read the configured frame rate, or use the standard 120 FPS fallback when the config file is unavailable.
        config_filepath: str = os.path.join(chunks_path, "config.pkl")
        if(os.path.exists(config_filepath)):
            with open(config_filepath, 'rb') as config_file:
                config_data: dict = dill.load(config_file)
            recording_fps: float = config_data['sensors']['W']['sensor_mode']['fps']
        else:
            warnings.warn(f"Config filepath does not exist: {config_filepath}. Assuming world-camera FPS is {WORLD_CAM_FPS}.", RuntimeWarning)
            recording_fps = WORLD_CAM_FPS

        assert recording_fps > 0, f"World camera FPS must be positive. Got: {recording_fps}"
        frame_period_ns: float = (10 ** 9) / recording_fps

    # First, let's find the world metadata chunks
    world_metadata_chunks: list[str] = natsorted([os.path.join(chunks_path, filename)
                                                  for filename in os.listdir(chunks_path)
                                                  if filename.startswith("world")
                                                  and "metadata" in filename
                                                 ]
                                            )
    assert len(world_metadata_chunks) > 0, f"0 world metadata chunks found @ {chunks_path}"
    
    # Once we have them, let's iterate over the paths 
    metadata: list[np.ndarray] | np.ndarray = []

    # Once we have the paths to the metadata chunks, we will simply read them in 
    path_iterator: Iterable = range(len(world_metadata_chunks)) if verbose is False else tqdm(range(len(world_metadata_chunks)), desc="Loading world metadata chunks")

    previous_chunk_end: float | None = None
    current_chunk_start: float | None = None
    for chunk_num in path_iterator:
        # Retrieve the path to this chunk 
        metadata_chunk_path: str = world_metadata_chunks[chunk_num]

        # Load in the metadata
        world_metadata: np.ndarray = np.load(metadata_chunk_path)
        # Sometimes the last chunk is empty. If this is true, skip it
        if(len(world_metadata) == 0):
            break

        current_chunk_start = world_metadata[0, 0]

        # If this is the first chunk, we can save both the previous end
        # and current start from it 
        if(chunk_num == 0):
            previous_chunk_end = world_metadata[-1, 0]    

        # Otherwise, we need to interpolate 
        # the timestamps in between the chunks 
        # that comes BEFORE world_metadata 
        elif fill_missing_frames:
            # Find the missing time in nano seconds
            gap_ns: float = current_chunk_start - previous_chunk_end

            # Estimate how many missing frames are between chunks
            # We subtract 1 because the boundary timestamps already exist
            num_missing_frames: int = max(0, gap_ns // frame_period_ns) - 1

            if(num_missing_frames > 0):
                # Calculate the timestamps for this downtime
                downtime_timestamps: np.ndarray = (previous_chunk_end + frame_period_ns * np.arange(1, num_missing_frames + 1))

                # Allocate an array for the downtime metadata, which will be filled with NaNs except for the timestamps 
                downtime_metadata: np.ndarray = np.full( (len(downtime_timestamps), world_metadata.shape[1]), np.nan, dtype=world_metadata.dtype )
                downtime_metadata[:, 0] = downtime_timestamps

                # Combine them into the world metadata BEFORE the existing world metadata
                world_metadata = np.concatenate([downtime_metadata, world_metadata])

            # Now, the previous chunk end is the last timestamp in this buffer
            previous_chunk_end = world_metadata[-1, 0]

        # Save it to the running list 
        metadata.append(world_metadata)

    # convert to standardized np array 
    metadata = np.vstack(metadata)

    # World timestamps are in nanoseconds by default,
    # so convert to seconds if desired
    if(convert_to_seconds is True):
        metadata[:, 0] /= ( 10 ** 9) 

    # Make a dataframe so that the columns are clearly labeled
    # NOTE: Legacy columns are timestamp, AG, DG, EXP 
    #       Modern columns are WORLD_AGC_METADATA_COLS
    metadata: pd.DataFrame = pd.DataFrame(metadata, columns=["timestamp"] + (list(WORLD_AGC_METADATA_COLS) if metadata.shape[-1] == 6 else ["Again", "Dgain", "exposure"]))

    return metadata


def plot_world_camera_settings(
    world_metadata: pd.DataFrame,
    events: Iterable[dict] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, tuple[plt.Axes, plt.Axes]]:
    """Plot world-camera gain and exposure settings over time and frame number.

    The modern metadata layout plots the camera and AGC-requested analog
    gains, the camera digital gain, and both camera and AGC-requested
    exposures. The legacy ``Again``, ``Dgain``, and ``exposure`` layout is
    also supported. Gain traces use shades of blue on the left axis, while
    exposure traces use shades of orange on the right axis so the two axis
    color families never overlap. The bottom X axis shows elapsed seconds and
    the top X axis spans zero to the last metadata row on a linear frame scale.

    Args:
        world_metadata: DataFrame returned by ``world_metadata_from_chunks``
            with timestamps converted to seconds (the default behavior).
        events: Iterable of dictionaries with either ``timestamp`` (absolute
            seconds, on the metadata clock) or ``frame_num`` (zero-based row).
            If both are supplied, timestamp takes precedence. Optional ``label``
            overrides the default one-based "event i" annotation. Markers use
            green/purple/red colors and are excluded from the settings legend.
        ax: Optional gain axis to draw on; otherwise create a larger figure.
            The caller controls layout when supplying an axis.

    Returns:
        The figure and a ``(gain_axis, exposure_axis)`` tuple.

    Raises:
        TypeError: If ``world_metadata`` is not a pandas DataFrame.
        ValueError: If the DataFrame is empty, has no finite timestamps, or
            does not contain either a complete modern or legacy metadata
            layout.
    """
    # Fail early with a clear message instead of allowing plotting or column
    # lookup errors to surface later in the function.
    if(not isinstance(world_metadata, pd.DataFrame)):
        raise TypeError("world_metadata must be a pandas DataFrame")
    if(world_metadata.empty):
        raise ValueError("world_metadata must contain at least one row")
    if("timestamp" not in world_metadata.columns):
        raise ValueError("world_metadata must contain a 'timestamp' column")

    # Recordings can use one of two metadata layouts. Modern recordings store
    # the settings applied by the camera alongside the settings requested by
    # the custom AGC. Legacy recordings store only one value per setting.
    modern_gain_columns: tuple[str, ...] = ("cameraAgain", "AGCAgain", "AGCDgain")
    modern_exposure_columns: tuple[str, ...] = ("cameraExposure", "AGCExposure")
    legacy_gain_columns: tuple[str, ...] = ("Again", "Dgain")
    legacy_exposure_columns: tuple[str, ...] = ("exposure",)

    # Select a complete layout rather than plotting a misleading partial set
    # of camera settings when one or more expected columns are absent.
    modern_columns: set[str] = set(modern_gain_columns + modern_exposure_columns)
    legacy_columns: set[str] = set(legacy_gain_columns + legacy_exposure_columns)
    available_columns: set[str] = set(world_metadata.columns)
    if(modern_columns.issubset(available_columns)):
        gain_columns = modern_gain_columns
        exposure_columns = modern_exposure_columns
    elif(legacy_columns.issubset(available_columns)):
        gain_columns = legacy_gain_columns
        exposure_columns = legacy_exposure_columns
    else:
        raise ValueError(
            "world_metadata must contain either the modern setting columns "
            f"{sorted(modern_columns)} or legacy setting columns {sorted(legacy_columns)}"
        )

    # world_metadata_from_chunks returns timestamps in seconds by default.
    # Subtract the first valid timestamp so the plot begins at zero while
    # retaining NaNs in their original positions.
    timestamps: np.ndarray = world_metadata["timestamp"].to_numpy(dtype=np.float64)
    finite_timestamps: np.ndarray = timestamps[np.isfinite(timestamps)]
    if(finite_timestamps.size == 0):
        raise ValueError("world_metadata must contain at least one finite timestamp")
    elapsed_seconds: np.ndarray = timestamps - finite_timestamps[0]

    # Keep every gain trace in the blue family and every exposure trace in the
    # orange family. Shades distinguish related camera and AGC values, while
    # the separate families make the two Y axes visually unambiguous.
    gain_colors: dict[str, str] = {
        "cameraAgain": "#08519c",
        "AGCAgain": "#6baed6",
        "AGCDgain": "#08306b",
        "Again": "#2171b5",
        "Dgain": "#6baed6",
    }
    exposure_colors: dict[str, str] = {
        "cameraExposure": "#d94801",
        "AGCExposure": "#fd8d3c",
        "exposure": "#e6550d",
    }
    display_names: dict[str, str] = {
        "cameraAgain": "Camera AGain",
        "AGCAgain": "AGC AGain",
        "AGCDgain": "Camera DGain",
        "cameraExposure": "Camera exposure",
        "AGCExposure": "AGC exposure",
        "Again": "AGain",
        "Dgain": "DGain",
        "exposure": "Exposure",
    }

    # twinx shares the elapsed-time X axis while allowing exposure, whose
    # numeric range is much larger than gain, to retain its own Y scale.
    if ax is None:
        figure, gain_axis = plt.subplots(figsize=(16, 7))
    else:
        gain_axis = ax
        figure = ax.figure
    exposure_axis: plt.Axes = gain_axis.twinx()

    # Collect lines from both axes so they can appear in one combined legend.
    lines: list = []

    # Draw gain values against the left-hand scale.
    for column_name in gain_columns:
        lines.extend(gain_axis.plot(
            elapsed_seconds,
            world_metadata[column_name].to_numpy(dtype=np.float64),
            color=gain_colors[column_name],
            linewidth=1.5,
            label=display_names[column_name],
        ))

    # Draw exposure values against the independent right-hand scale.
    for column_name in exposure_columns:
        lines.extend(exposure_axis.plot(
            elapsed_seconds,
            world_metadata[column_name].to_numpy(dtype=np.float64),
            color=exposure_colors[column_name],
            linewidth=1.5,
            label=display_names[column_name],
        ))

    # Use a simple, evenly spaced frame scale across the recording.
    frame_axis = gain_axis.twiny()
    frame_axis.set_xlim(0, max(len(world_metadata) - 1, 1))
    frame_axis.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, integer=True))
    if len(world_metadata) == 1:
        frame_axis.set_xticks([0])
    frame_axis.set_xlabel("Frame number")

    # Resolve events against the same metadata clock/physical-row coordinates.
    event_colors = ("#238b45", "#7a0177", "#cb181d", "#525252")
    for event_number, event in enumerate(events if events is not None else (), 1):
        if not isinstance(event, dict):
            raise TypeError(f"Event {event_number} must be a dict")
        if "timestamp" in event:
            event_x = float(event["timestamp"]) - finite_timestamps[0]
        elif "frame_num" in event:
            frame_num = float(event["frame_num"])
            if not np.isfinite(frame_num) or not frame_num.is_integer() or not 0 <= frame_num < len(timestamps):
                raise ValueError(f"Event {event_number} frame_num must be an in-range integer")
            event_x = elapsed_seconds[int(frame_num)]
        else:
            raise ValueError(f"Event {event_number} requires timestamp or frame_num")
        if not np.isfinite(event_x):
            raise ValueError(f"Event {event_number} must have a finite timestamp")
        color = event_colors[(event_number - 1) % len(event_colors)]
        gain_axis.axvline(event_x, color=color, linestyle="--", linewidth=1,
                         label="_nolegend_")
        gain_axis.annotate(
            str(event.get("label", f"event {event_number}")),
            xy=(event_x, 0.98 - 0.08 * ((event_number - 1) % 3)),
            xycoords=gain_axis.get_xaxis_transform(), xytext=(3, 0),
            textcoords="offset points", rotation=90, fontsize=8,
            color=color, va="top", ha="left",
        )

    # Match each Y axis's labels, ticks, and visible spine to its line family.
    # This reinforces which scale belongs to which group of traces.
    gain_axis.set_xlabel("Elapsed time (s)")
    gain_axis.set_ylabel("Gain", color="#08306b")
    exposure_axis.set_ylabel("Exposure", color="#d94801")
    gain_axis.set_ylim(bottom=1)
    exposure_axis.set_ylim(bottom=0)
    gain_axis.tick_params(axis="y", colors="#08306b")
    exposure_axis.tick_params(axis="y", colors="#d94801")
    gain_axis.spines["left"].set_color("#08306b")
    exposure_axis.spines["right"].set_color("#d94801")
    gain_axis.grid(True, alpha=0.25)
    gain_axis.margins(x=0)
    gain_axis.set_title("World camera settings")

    # Matplotlib otherwise creates one legend per axis, so explicitly pass all
    # lines to the left axis to produce a single complete legend.
    gain_axis.legend(lines, [line.get_label() for line in lines],
                     loc="upper left", bbox_to_anchor=(1.14, 1.0), borderaxespad=0)
    if ax is None:
        figure.tight_layout()

    return figure, (gain_axis, exposure_axis)


def world_raw_frames_from_chunks(path_to_recording: str, 
                                 use_mean_frame: bool=False,
                                 verbose: bool = False
                                ) -> np.ndarray:
    """Load raw world-camera frames in natural chunk order.

    Args:
        path_to_recording: Directory containing world frame .npy chunks and their metadata.
        use_mean_frame: If True, return the spatial mean of each frame instead of the full frame.
        verbose: Whether to show a progress bar while reading chunks.

    Returns:
        An array shaped (frames, rows, cols), or a one-dimensional array of
        per-frame means when use_mean_frame is True.
    """
    frame_chunk_paths: list[str] = [os.path.join(path_to_recording, folder_name) for folder_name in natsorted(os.listdir(path_to_recording)) if folder_name.startswith("world") and not folder_name.endswith("metadata.npy")]
    assert len(frame_chunk_paths) > 0, f"No frame chunks found in {path_to_recording}"
    
    frames: list[np.ndarray] = []
    chunk_iterator: Iterable = range(len(frame_chunk_paths)) if verbose is False else tqdm(range(len(frame_chunk_paths)), desc="Processing chunks")
    for chunk_num in chunk_iterator:
        chunk_path: str = frame_chunk_paths[chunk_num]

        frame_chunk: np.ndarray = np.load(chunk_path)
        if(use_mean_frame is False):
            frames.extend(frame_chunk)
            continue 
        frames.extend(np.mean(frame_chunk, axis=(1, 2)))

    return np.array(frames)

def world_counts_to_radiance(image_or_video: np.ndarray,
                             agc_settings: dict[str, float | np.ndarray],
                             visualize_results: bool=False
                            ) -> None:
    """Convert corrected world-camera counts to absolute radiance units.

    Args:
        image_or_video: Flat-fielded and color-corrected world-camera counts
            with shape ``(rows, cols)`` or ``(frames, rows, cols)``. The input
            is scaled in place.
        agc_settings: Camera settings containing either the legacy keys
            ``exposure``, ``Again``, and ``Dgain`` or the modern keys
            ``cameraExposure``, ``cameraAgain``, and ``AGCDgain``. Values may
            be scalars or one-dimensional arrays containing one value per
            frame.
        visualize_results: When ``True``, display the counts before and after
            conversion. Visualization supports a single frame only.

    Returns:
        None. The supplied array is modified in place. When visualization is
        enabled, the before/after figure is displayed without returning it.
    """
    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert image_or_video.ndim == 2, "Radiance visualization only supports a single frame"
        unmodified_image_or_video = image_or_video.copy()

    # Accept both generations of world-camera metadata names.
    exposure: float | np.ndarray = agc_settings["exposure"] if "exposure" in agc_settings else agc_settings["cameraExposure"]
    analog_gain: float | np.ndarray = agc_settings["Again"] if "Again" in agc_settings else agc_settings["cameraAgain"]
    digital_gain: float | np.ndarray = agc_settings["Dgain"] if "Dgain" in agc_settings else agc_settings["AGCDgain"]

    # Match reconstructionPipeline.m: undo digital gain at the AGC target,
    # subtract dark signal, then invert the sensor response for each frame.
    set_point: np.ndarray = 127.0 / np.asarray(digital_gain, dtype=np.float64) - WORLD_DARK_SIGNAL
    smax: float = 255.0 - WORLD_DARK_SIGNAL
    exponent: float = WORLD_FULL_WELL_CLIPPING_EXPONENT
    linearized_set_point: np.ndarray = set_point / (1 - (set_point / smax) ** exponent) ** (1 / exponent)
    effective_set_point: np.ndarray = (
        linearized_set_point
        * WORLD_MEAN_SPATIAL_CORRECTIONS[image_or_video.shape[-2:]]
    )

    # Match MATLAB: 10.^polyval(agcToRadianceP, log10(cameraScore)).
    # The fitted coefficients directly predict mean integrated radiance.
    this_camera_score: np.ndarray = np.asarray(np.asarray(exposure) * np.asarray(analog_gain) / effective_set_point, dtype=np.float64)
    mean_integrated_radiance: np.ndarray = np.power(
        10.0, np.polyval(WORLD_AGC_TO_RADIANCE_P, np.log10(this_camera_score))
    )

    # Give each buffered frame its own broadcastable radiance scale. Scalar
    # settings naturally apply the same scale to every frame.
    radiance_scale: np.ndarray = mean_integrated_radiance / effective_set_point
    if(image_or_video.ndim == 3 and np.ndim(radiance_scale) > 0):
        radiance_scale: np.ndarray = radiance_scale.reshape(-1, 1, 1)

    # Scale corrected counts by the fitted mean integrated radiance.
    image_or_video *= radiance_scale

    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("World Counts to Radiance (Before / After)", fontweight='bold', fontsize=18)
        axes[0].imshow(unmodified_image_or_video, cmap="gray")
        axes[0].set_title("Before")
        axes[0].axis("off")
        axes[1].imshow(image_or_video, cmap="gray")
        axes[1].set_title("After")
        axes[1].axis("off")
        plt.tight_layout()
        plt.show()
        return None

    return None


def world_transformation_pipeline(raw_frame_or_buffer: np.ndarray,
                                  agc_settings: dict[str, float | np.ndarray],
                                  gaze_angles_xy: np.ndarray | None = None,
                                  age: int | None = None,
                                  n_workers: int=WORLD_IMPUTATION_WORKERS,
                                  stagewise_output: list[np.ndarray] | None = None,
                                  demosaic_pool: Pool | None = None
                                 ) -> dict[str, np.ndarray | dict[str, object]]:
    """Reconstruct and demosaic world-camera radiance using Geoff's seven stages.

    The result is RGB radiance, matching MATLAB ``reconstructionPipeline``
    followed by ``demosaicRadianceMap``. Each processing stage calls its standalone
    helper, whose ``visualize_results`` option can be used independently.

    Raw input is preserved. After imputation, the fielding, color, and radiance
    corrections reuse one output array to avoid copying an entire buffer at
    every stage. ``n_workers`` controls frame-wise imputation; use 1 for serial
    processing. Digital gain enters only through the radiance AGC set point.

    Args:
        raw_frame_or_buffer: Raw 8-bit counts shaped (rows, cols) or (frames, rows, cols). The input
            is preserved.
        agc_settings: Exposure, analog gain, and digital gain as scalars or one value per frame.
            Accepts legacy or modern metadata keys.
        gaze_angles_xy: Reserved for later processing; unused here.
        age: Reserved for later processing; unused here.
        n_workers: Maximum number of processes used for frame-wise imputation. Use 1 for serial
            execution.
        stagewise_output: Optional empty list that receives eight independent
            float64 snapshots: raw counts, quantization correction, linearization,
            imputation, fielding, RGB correction, Bayer radiance, and demosaiced
            RGB radiance. The final snapshot adds a trailing RGB dimension. Omit this argument to avoid the snapshot copies.

        demosaic_pool: Optional persistent pool from demosaic_worker_pool(n_workers).
            Reuse it across recording chunks for shared-memory demosaicing.
            None keeps demosaicing serial, avoiding startup for one-off calls.

    Returns:
        Dictionary with "data", the float64 RGB radiance in W/m²/sr shaped
        (rows, cols, 3) or (frames, rows, cols, 3), and "metadata", a dictionary
        containing the processing calibration constants under their module names,
        source calibration records, input AGC settings, and processing settings.
        Shape-indexed constants contain the selected frame's value. Metadata is
        copied once per call so editing it cannot change shared calibrations.

    Notes:
        Benchmarked on 2026-10-05 on a Mac Studio (Mac15,14), Apple M3 Ultra,
        28-core CPU (20 performance + 8 efficiency), 256 GB unified memory,
        macOS 15.3 (arm64).
        Imputation compared 1/2/4/6/8/12 workers on 120 real 480 x 640 frames,
        including process startup and shared-memory transfers. Six was chosen
        as the baseline (48.2 s; eight was slightly faster at 47.3 s).
        Demosaicing compared 8/12/16/24/28 workers on sequential 3,600-frame
        float64 buffers using repeated seeded radiance data and persistent pools.
        Three warmed calls included allocation, copies, and cleanup; sixteen
        was fastest tested at 10.39 s/buffer, with a 22.77 s first call.
        Estimated demosaicing time for ten buffers is 116.3 s, excluding file I/O
        and other stages. All configurations matched their serial references.
        See tmp/benchmark_imputation_workers.json and
        tmp/benchmark_demosaic_large_buffers.json for the measurements.
    """

    # First, we need to ensure that the transformation pipeline
    # can accept either a single raw frame or a buffer of frames,
    # that's it
    if(raw_frame_or_buffer.ndim not in (2, 3)):
        raise ValueError(
            "World reconstruction requires shape (rows, cols) or "
            f"(frames, rows, cols). Got {raw_frame_or_buffer.shape}."
        )

    # Get the image shape (r, c). We can get this from .shape[-2]
    # because if it's a buffer [f, r, c] we would extract r, c
    # and if it's a single frame [r, c] then it's just the entire thing
    image_shape: tuple[int, int] = raw_frame_or_buffer.shape[-2:]

    # We need to ensure we support this size of image
    # for the fielding functions and the bayer correction matrices.
    # These calibrations were made for specific image shapes
    if(image_shape not in WORLD_FIELDING_FUNCTIONS
       or image_shape not in WORLD_BAYER_CORRECTION_MATRICES):
        raise ValueError(f"Unsupported calibrated image shape: {image_shape}.")

    # Save a boolean of whether we want to store stagewise output or not
    store_stagewise_output: bool = stagewise_output is not None
    if(store_stagewise_output is True):
        if not isinstance(stagewise_output, list):
            raise TypeError("stagewise_output must be an empty list or None.")
        if stagewise_output:
            raise ValueError("stagewise_output must be empty before reconstruction.")

        # Stage 0 of this would be to supply the raw image / buffer
        stagewise_output.append(raw_frame_or_buffer.astype(np.float64, copy=True))

    # Stage 1: convert raw sensor counts to double precision and correct the
    # quantization bias. The camera's 10-to-8-bit shift discards an average of
    # 0.375 counts in 8-bit units.
    # We make a single copy at the start to not edit the input array
    quantization_corrected: np.ndarray = raw_frame_or_buffer.astype(np.float64, copy=True)
    quantization_corrected += 0.375

    # Store this stage in the stagewise output if desired
    if(store_stagewise_output is True):
        stagewise_output.append(quantization_corrected.copy())

    # Stage 2: invert the fitted sensor response in place. Negative counts
    # remain signed after dark subtraction; only the nonlinear gain is clamped
    # to a non-negative input. Unreliable high-count samples become Inf.
    linearize_camera_responsivity(
        quantization_corrected,
        original_bit_depth=8,
        dark_noise=WORLD_DARK_SIGNAL,
        clipping_exponent=WORLD_FULL_WELL_CLIPPING_EXPONENT,
        visualize_results=False,
    )
    # Rename the same array after the in-place operation; no data are copied.
    linearized: np.ndarray = quantization_corrected
    if(store_stagewise_output is True):
        stagewise_output.append(linearized.copy())

    # Stage 3: replace non-positive floor samples and infinite ceiling samples
    # using each frame's cross-channel Bayesian model. Geoff's current pipeline
    # does not multiply pixel values by digital gain before this step.
    imputed: np.ndarray = impute_pixel_values(
        linearized,
        bayer_pattern="BGGR",
        visualize_results=False,
        n_workers=n_workers,
    )
    # Imputation allocates its output, so this stage uses the returned array.
    if(store_stagewise_output is True):
        stagewise_output.append(imputed.copy())

    # Stage 4: correct spatial sensitivity with the calibrated fielding map.
    # The helper modifies the array in place. The new name below refers to the
    # same array, making the next stage clear without allocating another copy.
    apply_fielding_function(imputed, visualize_results=False)
    # Rename the same array after fielding; this alias does not copy data.
    fielding_corrected: np.ndarray = imputed
    if(store_stagewise_output is True):
        stagewise_output.append(fielding_corrected.copy())

    # Stage 5: equalize the Bayer RGB channels with their calibrated weights.
    # This helper also operates in place and uses the cached correction map.
    apply_color_correction(
        fielding_corrected,
        visualize_results=False,
    )
    # Rename the same array after RGB correction; this alias does not copy data.
    color_corrected: np.ndarray = fielding_corrected
    if(store_stagewise_output is True):
        stagewise_output.append(color_corrected.copy())

    # Stage 6: convert corrected counts to absolute integrated radiance.
    # The helper accounts for digital gain in the AGC set point, uses harmonic
    # spatial means, and applies Geoff's updated camera-score calibration.
    world_counts_to_radiance(
        color_corrected,
        agc_settings=agc_settings,
        visualize_results=False,
    )
    # Rename the same array after radiance scaling; this alias does not copy data.
    radiance_map: np.ndarray = color_corrected
    if(store_stagewise_output is True):
        stagewise_output.append(radiance_map.copy())

    # Stage 7: use green-guided color ratios to reconstruct full RGB radiance.
    # The standalone helper reuses interpolation weights across buffered frames.
    demosaiced: np.ndarray = demosaic_radiance_map_rcd(
        radiance_map, bayer_pattern="BGGR", worker_pool=demosaic_pool,
    )
    if(store_stagewise_output is True):
        stagewise_output.append(demosaiced.copy())

    # Record every processing calibration, including the source fit metadata.
    # Store only the applied geometry's maps, once per buffer rather than per frame.
    calibration_metadata: dict[str, object] = deepcopy({
        "WORLD_FULL_WELL_CLIPPING_EXPONENT": WORLD_FULL_WELL_CLIPPING_EXPONENT,
        "WORLD_DARK_SIGNAL": WORLD_DARK_SIGNAL,
        "WORLD_MAX_ALLOWED_LINEARIZATION_DERIVATIVE": WORLD_MAX_ALLOWED_LINEARIZATION_DERIVATIVE,
        "WORLD_FIELDING_FUNCTIONS": WORLD_FIELDING_FUNCTIONS[image_shape],
        "WORLD_RADIOMETRIC_CALIBRATION": WORLD_RADIOMETRIC_CALIBRATION,
        "WORLD_RGB_SCALARS": WORLD_RGB_SCALARS,
        "WORLD_RADIOMETRIC_CORRECTION_MAP": WORLD_RADIOMETRIC_CORRECTION_MAP,
        "WORLD_BAYER_CORRECTION_MATRICES": WORLD_BAYER_CORRECTION_MATRICES[image_shape],
        "WORLD_INTEGRATED_RADIANCE_CALIBRATION": WORLD_INTEGRATED_RADIANCE_CALIBRATION,
        "WORLD_AGC_TO_RADIANCE_P": WORLD_AGC_TO_RADIANCE_P,
        "WORLD_MEAN_SPATIAL_CORRECTIONS": WORLD_MEAN_SPATIAL_CORRECTIONS[image_shape],
        "WORLD_TRIANGULATION_SHEAR": WORLD_TRIANGULATION_SHEAR,
        "agc_settings": agc_settings,
        "image_shape": image_shape,
        "bayer_pattern": "BGGR",
        "original_bit_depth": 8,
        "quantization_bias": 0.375,
        "radiance_agc_target": 127.0,
        "demosaic_ratio_epsilon": 1e-6,
        "units": "W/m²/sr",
    })
    return {"data": demosaiced, "metadata": calibration_metadata}



def plot_reconstruction_stages(
    stagewise_output: list[np.ndarray],
    title: str = "World-camera reconstruction",
) -> tuple[plt.Figure, np.ndarray]:
    """Display the saved Bayer and RGB stages using Geoff's MATLAB figure layout.

    Each column shows one stage above its RGB-channel histogram. Images use
    their finite 85th percentile as a display scale, with ceiling samples in
    red and non-positive floor samples in blue. Histograms show the percentage
    of channel samples versus the percentage of that stage's finite maximum.
    Hover coordinates report the original values rather than display-scaled ones.

    Args:
        stagewise_output: Eight snapshots (seven Bayer planes and one RGB image) collected by
            world_transformation_pipeline for a single frame. For a buffer,
            select the same frame index from each snapshot before plotting.
        title: Title for the complete figure.

    Returns:
        The displayed figure and its (2, 8) axes array. Saved stages are not
        modified; display scaling operates on separate arrays.

    Raises:
        ValueError: If the seven Bayer planes and final RGB image have incompatible shapes.
    """
    stage_names: tuple[str, ...] = (
        "Raw", "Quantization", "Linearization", "Imputation",
        "Flat field", "RGB correction", "Radiance", "Demosaiced",
    )
    if len(stagewise_output) != len(stage_names):
        raise ValueError("Provide all eight snapshots from world_transformation_pipeline.")
    image_shape: tuple[int, ...] = stagewise_output[0].shape
    if (len(image_shape) != 2
        or any(stage.shape != image_shape for stage in stagewise_output[:-1])
        or stagewise_output[-1].shape != (*image_shape, 3)):
        raise ValueError("Provide seven same-shaped Bayer planes followed by their RGB image.")

    def original_pixel_coordinates(x: float, y: float, stage: np.ndarray) -> str:
        """Format the unscaled stage value under the mouse pointer.

        Args:
            x: Horizontal image coordinate in pixels.
            y: Vertical image coordinate in pixels.
            stage: Original Bayer values for the selected axes.

        Returns:
            A coordinate/value label, or an empty string outside the image.
        """
        if not np.isfinite(x) or not np.isfinite(y):
            return ""
        row: int = int(np.floor(y + 0.5))
        col: int = int(np.floor(x + 0.5))
        if 0 <= row < stage.shape[0] and 0 <= col < stage.shape[1]:
            value: str = (f"{stage[row, col]:.8g}" if stage.ndim == 2
                          else np.array2string(stage[row, col], precision=8))
            return f"row={row}, col={col}, value={value}"
        return ""

    figure: plt.Figure
    axes: np.ndarray
    figure, axes = plt.subplots(2, len(stage_names), figsize=(25, 7), constrained_layout=True)
    figure.suptitle(title, fontsize=15, fontweight="bold")
    edges: np.ndarray = np.linspace(0, 1, 256)

    for stage_index, (stage_name, stage) in enumerate(zip(stage_names, stagewise_output)):
        # Display scaling never changes the saved numerical stage. Empty or
        # non-positive ranges use a scale of one so all-dark images still plot.
        finite_values: np.ndarray = stage[np.isfinite(stage)]
        display_scale: float = float(np.percentile(finite_values, 85)) if finite_values.size else 1.0
        if display_scale <= 0:
            display_scale = 1.0
        normalized: np.ndarray = np.clip(
            np.nan_to_num(stage / display_scale, nan=0.0, posinf=1.0, neginf=0.0), 0, 1
        )
        # The final stage already has RGB channels; Bayer stages need replication.
        display_rgb: np.ndarray = (np.repeat(normalized[..., None], 3, axis=-1)
                                   if stage.ndim == 2 else normalized.copy())
        ceiling: np.ndarray = np.isinf(stage)
        floor: np.ndarray = stage <= 0
        if(stage.ndim == 3):
            ceiling = np.any(ceiling, axis=-1)
            floor = np.any(floor, axis=-1)
        display_rgb[ceiling] = (1, 0, 0)
        display_rgb[floor] = (0, 0, 1)

        image_axis: plt.Axes = axes[0, stage_index]
        image_axis.imshow(display_rgb)
        image_axis.set_title(f"{stage_index}. {stage_name}", fontsize=11)
        image_axis.axis("off")
        # Bind this stage now so every axes reports its own original values.
        image_axis.format_coord = partial(original_pixel_coordinates, stage=stage)

        # Bayer samples belong to exactly one measured color. Combine both
        # green subgrids, keeping the channel denominators separate as in MATLAB.
        channel_values: tuple[np.ndarray, ...] = (
            stage[1::2, 1::2].ravel(),
            np.concatenate((stage[0::2, 1::2].ravel(), stage[1::2, 0::2].ravel())),
            stage[0::2, 0::2].ravel(),
        ) if stage.ndim == 2 else tuple(stage[..., channel].ravel() for channel in range(3))
        histogram_scale: float = float(np.max(finite_values)) if finite_values.size else 1.0
        if histogram_scale <= 0:
            histogram_scale = 1.0
        minimum_labels: list[str] = []
        maximum_labels: list[str] = []
        histogram_axis: plt.Axes = axes[1, stage_index]
        for channel_name, color, values in zip("RGB", ("red", "green", "blue"), channel_values):
            finite_channel: np.ndarray = values[np.isfinite(values)]
            minimum_labels.append(f"{np.min(finite_channel):.3g}" if finite_channel.size else "n/a")
            maximum_labels.append(f"{np.max(finite_channel):.3g}" if finite_channel.size else "n/a")
            # MATLAB puts non-finite samples in the highest bin. Negative
            # finite samples fall outside [0, 1] but remain in the denominator.
            histogram_values: np.ndarray = np.where(np.isfinite(values), values, histogram_scale)
            counts: np.ndarray = np.histogram(histogram_values / histogram_scale, bins=edges)[0]
            percentages: np.ndarray = counts * (100.0 / values.size) if values.size else counts.astype(float)
            histogram_axis.plot(edges[:-1] * 100, percentages, color=color, linewidth=1, label=channel_name)

        histogram_axis.text(
            0.02, 0.98,
            "RGB min: " + ", ".join(minimum_labels)
            + "\nRGB max: " + ", ".join(maximum_labels)
            + f"\nCeiling: {np.count_nonzero(np.isinf(stage)):,}"
            + f"  Floor: {np.count_nonzero(stage <= 0):,}",
            transform=histogram_axis.transAxes, va="top", fontsize=7,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )
        histogram_axis.set_xlim(0, 100)
        histogram_axis.set_ylim(bottom=0)
        histogram_axis.set_xlabel("% of stage maximum", fontsize=9)
        histogram_axis.tick_params(labelsize=8)
        histogram_axis.spines[["top", "right"]].set_visible(False)
        if stage_index == 0:
            histogram_axis.set_ylabel("% of channel pixels", fontsize=9)
            histogram_axis.legend(loc="lower right", fontsize=8)

    plt.show()
    return figure, axes


def visualize_pipeline(
    recording_path: str,
    chunk_range: tuple[int | None, int | None] = (0, None),
) -> list[tuple[plt.Figure, np.ndarray]]:
    """Run the reconstruction once per selected chunk and plot its saved stages.

    The middle frame of each chunk is passed to world_transformation_pipeline
    with an empty stagewise_output list. The resulting snapshots are displayed
    by plot_reconstruction_stages, following Geoff's MATLAB image/histogram
    layout. This function does not implement a separate processing pipeline.

    Args:
        recording_path: Directory containing world-frame .npy chunks and their
            matching <chunk stem>_metadata.npy files. Metadata may use either
            timestamp/Again/Dgain/exposure or timestamp plus WORLD_AGC_METADATA_COLS.
        chunk_range: Half-open (start, end) slice of the naturally sorted chunk
            list. None uses the corresponding beginning or end of the list.

    Returns:
        One (figure, axes) tuple per selected chunk. Each axes array has shape
        (2, 8): stage images on top and RGB histograms underneath.

    Raises:
        FileNotFoundError: If no world chunks or a matching metadata file exists.
        ValueError: If a chunk is empty, has the wrong shape, or its metadata
            does not match the frame count or a supported column layout.
    """
    recording_directory: pathlib.Path = pathlib.Path(recording_path).expanduser()
    chunk_paths: list[pathlib.Path] = natsorted(
        [path for path in recording_directory.glob("*world*.npy") if "metadata" not in path.name],
        key=lambda path: path.name,
    )
    if not chunk_paths:
        raise FileNotFoundError(f"No world chunk files found in: {recording_directory}")

    figures: list[tuple[plt.Figure, np.ndarray]] = []
    start, end = chunk_range
    for chunk_path in chunk_paths[start:end]:
        # Match metadata by filename rather than zipping two independently
        # sorted lists, which could silently pair unrelated files.
        metadata_path: pathlib.Path = chunk_path.with_name(f"{chunk_path.stem}_metadata.npy")
        if not metadata_path.is_file():
            raise FileNotFoundError(f"Missing metadata for {chunk_path.name}: {metadata_path}")

        # Memory mapping reads only the chosen frame instead of loading the
        # entire recording chunk into memory for a single diagnostic image.
        chunk: np.ndarray = np.load(chunk_path, mmap_mode="r")
        metadata: np.ndarray = np.load(metadata_path, mmap_mode="r")
        if chunk.ndim != 3 or chunk.shape[0] == 0:
            raise ValueError(f"Expected a nonempty (frames, rows, cols) chunk: {chunk_path}")
        if metadata.ndim != 2 or metadata.shape[0] != chunk.shape[0]:
            raise ValueError(f"Metadata rows must match the frame count: {metadata_path}")

        selected_index: int = len(chunk) // 2
        settings_values: np.ndarray = metadata[selected_index, 1:]
        setting_names: tuple[str, ...]
        if settings_values.size == len(WORLD_AGC_METADATA_COLS):
            setting_names = WORLD_AGC_METADATA_COLS
        elif settings_values.size == 3:
            setting_names = ("Again", "Dgain", "exposure")
        else:
            raise ValueError(f"Unsupported world metadata columns: {metadata_path}")
        agc_settings: dict[str, float] = {
            name: float(value) for name, value in zip(setting_names, settings_values)
        }

        # Capture independent copies through the actual reconstruction path.
        # All display logic consumes these snapshots and cannot alter the result.
        stagewise_output: list[np.ndarray] = []
        world_transformation_pipeline(
            chunk[selected_index], agc_settings, n_workers=1,
            stagewise_output=stagewise_output,
        )
        title: str = f"{recording_directory.name} | {chunk_path.name} | frame {selected_index}"
        figures.append(plot_reconstruction_stages(stagewise_output, title=title))

    return figures


def main() -> None:
    """Reserved entry point; no command-line processing is implemented yet.

    Args:
        None.

    Returns:
        None.
    """
    pass


if(__name__ == '__main__'):
    pass
