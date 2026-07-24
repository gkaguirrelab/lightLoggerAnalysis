"""Utility functions and constants for the world camera sensor."""

import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
import sys
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
 
def _load_pyagc():
    """Locate and import the optional ``PyAGC`` dependency.

    The world-camera utilities can run without ``PyAGC``, but when it is
    present we use it to populate AGC state tables and allowed setting
    ranges. This helper searches a small set of project-relative library
    locations, appends each candidate directory to ``sys.path`` if needed,
    and returns the first import that succeeds.

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


def _load_agc_lib():
    """Import the optional compiled AGC library through ``PyAGC``."""
    if PyAGC is None:
        return None

    try:
        return PyAGC.import_AGC_lib()
    except Exception:
        return None


AGC_LIB = _load_agc_lib()

# Function to load in the fielding functions in per dimensions. Measured and saved in MATLAB 
def _import_fielding_functions() -> dict[tuple[int], np.ndarray]:
    # First find the path to the fielding function directory 
    fielding_function_folder: str = os.path.join(pathlib.Path(__file__).parents[3], "code", "fieldingFunction") 
    assert os.path.exists(fielding_function_folder), f"Fielding function path: {fielding_function_folder} does not exist"

    # Find the correction maps for different image shapes 
    correction_maps_folder: str = os.path.join(fielding_function_folder, "correctionMaps")
    assert os.path.exists(correction_maps_folder), f"Correction maps folder: {correction_maps_folder} does not exist"

    # Generate a dictionary based off of these correction maps 
    correction_filepaths: list[str] = [os.path.join(correction_maps_folder, file) for file in os.listdir(correction_maps_folder) if file.endswith(".mat")]
    assert len(correction_filepaths) > 0, f"No correction files found at path: {correction_maps_folder}"

    fielding_functions: dict[tuple[int], np.ndarray] = {}
    for filepath in correction_filepaths:
        fielding_function: dict = scipy.io.loadmat(filepath)["correctionMap"].astype(np.float64, copy=False)
        
        fielding_functions[fielding_function.shape] = fielding_function

    return fielding_functions


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

# Store the scalar multipliers for all of the different colors of pixel's in an image 
# We calculated this by making a measurement with the light logger device 
# on the roof of Goddard on a cloudy day. We took the asymptotic RGB values of 
# the center region of this recording, scaled them relative to the blue channel, 
# and found individual multiplers a, b, c such that the original RGB values 
# times these numbers equaled the relative weights according to the PR670, 
# which Geoff also used to take a measurement of this cloudy day. We then 
# Solved the following formula to make sure the mean of these new weights 
# was equal to 1 
"""
Ra Ga Ba
Example:
1. 1.1. 1.2

Solve for:
Rb Gb Bb

where:
Mean[Rb Gb Bb] = 1
Gb/Rb = Ga/Ra Gb/Bb = Ga/Ba Rb/Bb = Ra/Ba
"""
WORLD_RGB_SCALARS: np.ndarray = np.array([1.032, 0.803, 1.164], dtype=np.float64) 

# Define a mapping between frame sizes and fielding functions of the camera
WORLD_FIELDING_FUNCTIONS: dict[tuple[int], np.ndarray] = _import_fielding_functions()

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

# This the dark noise of the camera. That is, we measured a recording 
# from the camera when it is entirely wrapped in black cloth 
# and this was the result. This is 
WORLD_DARK_NOISE: float = 16
WORLD_FULL_WELL_CLIPPING_EXPONENT: float = 5.3918


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


def linearize_camera_responsivity(image_or_video: np.ndarray,
                                  original_bit_depth: int = 8,
                                  dark_noise: float = WORLD_DARK_NOISE,
                                  clipping_exponent: float = WORLD_FULL_WELL_CLIPPING_EXPONENT,
                                  dst: np.ndarray | None = None,
                                  visualize_results: bool=False
                                  ) -> np.ndarray | tuple[np.ndarray, object]:
    """Linearize world-camera values using the fitted full-well model.

    This is the Python translation of ``linearizeY`` in
    ``fitFullWellCapacityEffect.m``. Inputs at or below ``dark_noise`` are
    mapped to zero, the top sensor value is treated as saturated, and the
    valid range ``dark_noise`` through ``2 ** original_bit_depth - 2`` is
    expanded onto ``0`` through ``2 ** original_bit_depth - 2``.

    Args:
        image_or_video: Raw camera frame or frame buffer to linearize.
        original_bit_depth: Bit depth of the input image values.
        dark_noise: The measured dark offset to remove before inversion.
        clipping_exponent: The fitted soft-clipping exponent from the
            full-well calibration.
        dst: Optional output array with the same shape as ``image_or_video``.
            When provided, the result is written into this array and the same
            array is returned.
        visualize_results: When ``True``, display a before/after figure and
            return it with the linearized result. Visualization supports only
            a single ``(rows, cols)`` frame and asserts otherwise.

    Returns:
        A rounded linearized array, or ``(linearized, figure)`` when
        visualization is requested.
    """
    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert image_or_video.ndim == 2, (
            f"Camera responsivity linearization visualization only supports "
            f"a single frame with shape (rows, cols). Got shape "
            f"{image_or_video.shape}."
        )
        unmodified_image_or_video = image_or_video.copy() 

    max_sensor_value: float = float(2 ** original_bit_depth - 1)
    max_linearized_value: float = float(2 ** original_bit_depth - 2)
    smin: float = float(dark_noise)
    smax: float = max_sensor_value - smin

    if(smax <= 1):
        raise ValueError(
            f"dark_noise={dark_noise} leaves no usable range for "
            f"original_bit_depth={original_bit_depth}."
        )

    y_prime: np.ndarray = np.clip(image_or_video.astype(np.float64, copy=False) - smin, 0, smax)
    y_max: float = smax - 1
    a_max: float = y_max / (1 - (y_max / smax) ** clipping_exponent) ** (1 / clipping_exponent)
    saturated_mask: np.ndarray = y_prime >= smax
    positive_mask: np.ndarray = (y_prime > 0) & ~saturated_mask

    if(dst is None):
        dst = np.empty_like(y_prime, dtype=np.float64)

    dst[~positive_mask & ~saturated_mask] = 0
    dst[saturated_mask] = max_sensor_value
    dst[positive_mask] = (
        (
            y_prime[positive_mask]
            / (1 - (y_prime[positive_mask] / smax) ** clipping_exponent) ** (1 / clipping_exponent)
        )
        / a_max
        * max_linearized_value
    )

    np.round(dst, out=dst)
    np.clip(dst, 0, max_sensor_value, out=dst)

    # If visualize results is true, we will print an output of what the image looks like 
    if(visualize_results is True):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Full Well Correction (Before / After)", fontweight='bold', fontsize=18)

        axes[0].imshow(unmodified_image_or_video, cmap="gray")
        axes[0].set_title("Before")
        axes[0].axis("off")

        axes[1].imshow(dst, cmap="gray")
        axes[1].set_title("After")
        axes[1].axis("off")

        plt.tight_layout()
        plt.show()

        return dst, fig

    return dst

def apply_fielding_function(image_or_video: np.ndarray, 
                            visualize_results: bool=False
                           ) -> None | tuple[np.ndarray, object]:
    """Apply the registered fielding correction in place.

    Args:
        image_or_video: A floating-point raw Bayer frame with shape
            ``(rows, cols)`` or frame buffer with shape
            ``(frames, rows, cols)``. The spatial frame size must exist in
            ``WORLD_FIELDING_FUNCTIONS``.
        visualize_results: When ``True``, display a before/after figure and
            return it with the corrected input array. Visualization supports
            only a single ``(rows, cols)`` frame and asserts otherwise.

    Returns:
        ``None`` after in-place correction, or ``(image_or_video, figure)``
        when visualization is requested.
    """
    assert image_or_video.ndim in (2, 3), (
        f"Fielding function requires a single frame with shape (rows, cols) "
        f"or a frame buffer with shape (frames, rows, cols). Got shape "
        f"{image_or_video.shape}."
    )
    
    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert image_or_video.ndim == 2, (
            f"Fielding function visualization only supports a single frame "
            f"with shape (rows, cols). Got shape {image_or_video.shape}."
        )
        unmodified_image_or_video = image_or_video.copy() 

    # Get the fielding function for this frame size (depending on if we passed in a frame buffer)
    # or a single image
    image_size: tuple = image_or_video.shape[1:3] if image_or_video.ndim == 3 else image_or_video.shape[:2]
    fielding_function: np.ndarray = WORLD_FIELDING_FUNCTIONS[image_size]
    assert image_size == fielding_function.shape, (
        f"Fielding function shape: {fielding_function.shape} is not equal to "
        f"frame shape: {image_size}"
    )

    # Apply the fielding function 
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

        return image_or_video, fig

    return None

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
    # Initialize a variable for the figure handle that
    # will be used to visualize (if desired)
    fig: object | None = None

    # Initialize an array of characters. RGB will represent 
    # the indices where there are the respective colors in 
    # bayer pattern 
    mask: np.ndarray = np.full(original_frame.shape[:2], 'x', dtype='<U1') if marker == "str" else np.full(original_frame.shape[:2], -1, dtype=np.float64)
    
    # Extract the dimensions of the frame 
    height, width = original_frame.shape[:2]

    # Create a list for only the red pixels in a frame
    world_r_pixels: np.ndarray = np.array( [ (r, c) 
                                            for r in range(height) 
                                            for c in range(width)
                                            if(r % 2 != 0 and c % 2 != 0)
                                        ], 
                                        dtype=np.uint64
                                        )

    # Create a list for only the green pixels in a frame 
    world_g_pixels: np.ndarray = np.array( [ (r, c) 
                                            for r in range(height)
                                            for c in range(width)
                                            if(r % 2 == 0 and c % 2 != 0)
                                            or(r % 2 != 0 and c % 2 == 0) 
                                        ],  
                                        dtype=np.uint64
                                        )

    # Create a list for only the blue pixels in a frame
    world_b_pixels: np.ndarray = np.array([ (r, c)
                                            for r in range(height)
                                            for c in range(width)
                                            if(r % 2 == 0 and c % 2 == 0)
                                          ], 
                                          dtype=np.uint64
                                        )
    
    # Set the values in the mask 
    for idx, (color, pixel_coords) in enumerate(zip("RGB", (world_r_pixels, world_g_pixels, world_b_pixels))):
        rows: np.ndarray = pixel_coords[:, 0]
        cols: np.ndarray = pixel_coords[:, 1]

        mask[rows, cols] = color if marker == "str" else idx

    # Visualize the results if desired
    if(visualize_results is True):
        # Convert the mask to have color 
        colored_image: np.ndarray = np.empty(tuple(list(original_frame.shape[:2]) + [3]), dtype=np.uint8)
        for idx, (color, pixel_coords) in enumerate(zip("RGB", (world_r_pixels, world_g_pixels, world_b_pixels))):
            rows: np.ndarray = pixel_coords[:, 0]
            cols: np.ndarray = pixel_coords[:, 1]
            color_vector: np.ndarray = np.zeros((3,), dtype=np.uint8)
            color_vector[idx] = 255
            colored_image[rows, cols] = color_vector

        # Initialize a figure with a single axis 
        fig, ax = plt.subplots(1, 1)

        # Middle axis will be the after but without 
        # color correction 
        ax.imshow(colored_image)
        ax.set_title(f"{original_frame.shape} Bayer RGB Pattern")
        ax.axis("off")

        # Show the plot 
        plt.show() 

        return mask, fig
    
    return mask

# Apply the per-color weights to the color pixels of a frame.
def apply_color_correction(image_or_video: np.ndarray, 
                           bayer_pixel_locations: list[np.ndarray] | None=None, 
                           visualize_results: bool=False,
                           ) -> None | tuple[np.ndarray, object]:
    """Apply calibrated Bayer-site RGB scaling factors in place.

    The correction constants in ``WORLD_RGB_SCALARS`` are applied to raw
    Bayer samples by color class. For a frame buffer, ``bayer_pixel_locations``
    can be precomputed once with :func:`generate_RGB_mask` and reused.

    Args:
        image_or_video: A floating-point raw Bayer frame with shape
            ``(rows, cols)`` or frame buffer with shape
            ``(frames, rows, cols)``.
        bayer_pixel_locations: Optional list of three ``(n, 2)`` coordinate
            arrays for red, green, and blue Bayer sites, respectively.
        visualize_results: When ``True``, display a before/after figure and
            return it with the corrected input array. Visualization supports
            only a single ``(rows, cols)`` frame and asserts otherwise.

    Returns:
        ``None`` after in-place correction, or ``(image_or_video, figure)``
        when visualization is requested.
    """
    assert image_or_video.ndim in (2, 3), (
        f"Color correction requires a single raw Bayer frame with shape "
        f"(rows, cols) or a raw Bayer frame buffer with shape "
        f"(frames, rows, cols). Got shape {image_or_video.shape}."
    )

    unmodified_image_or_video: np.ndarray | None = None
    if(visualize_results is True):
        assert image_or_video.ndim == 2, (
            f"Color correction visualization only supports a single frame "
            f"with shape (rows, cols). Got shape {image_or_video.shape}."
        )
        unmodified_image_or_video = image_or_video.copy() 
    
    # Generate tbe bayer mask if not already passedi n
    if(bayer_pixel_locations is None): 
        bayer_RGB_mask = generate_RGB_mask(image_or_video if image_or_video.ndim == 2 else image_or_video[0])
        bayer_pixel_locations = [ np.argwhere(bayer_RGB_mask == color) for color in "RGB" ]

    # Apply the color correction
    for (pixels, weight) in zip(bayer_pixel_locations, WORLD_RGB_SCALARS):
        rows: np.ndarray = pixels[:, 0]
        cols: np.ndarray = pixels[:, 1]

        # Apply the weight to the specified pixels 
        # on either this single frame or entire video 
        if(image_or_video.ndim == 3):
            image_or_video[:, rows, cols] *= weight
        else:
            image_or_video[rows, cols] *= weight

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

        return image_or_video, fig

    # Return the modified frame
    return None


def apply_color_weights(original_frame: np.ndarray,
                        visualize_results: bool=False,
                        ) -> None | tuple[np.ndarray, object]:
    """Backward-compatible alias for :func:`apply_color_correction`.

    Older code in this project referred to the radiometric RGB correction
    factors as "color weights." The actual implementation lives in
    ``apply_color_correction``; this wrapper preserves the older public API
    and forwards the arguments unchanged.

    Args:
        original_frame: Raw Bayer frame to correct.
        visualize_results: Whether to request the optional before/after
            figure from ``apply_color_correction``.

    Returns:
        ``None`` after in-place correction, or ``(original_frame, figure)``
        when visualization is enabled.
    """
    return apply_color_correction(original_frame, visualize_results=visualize_results)

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
    """Compiled floor/ceiling propagation over normalized frame buffers."""
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

        # Middle axis will be the after but without 
        # color correction 
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
def world_timestamps_from_chunks(recording_path: str, 
                                 convert_to_seconds: bool=True,
                                 verbose: bool=False
                                ) -> np.ndarray:
    """Return the concatenated world-camera timestamp vector.

    This is a convenience wrapper around ``world_metadata_from_chunks`` that
    extracts only the timestamp column after chunk concatenation and any gap
    interpolation.

    Args:
        recording_path: Recording directory containing world metadata chunks.
        convert_to_seconds: Whether to convert the stored nanosecond
            timestamps into seconds.
        verbose: Whether to display chunk-loading progress.

    Returns:
        One-dimensional NumPy array of world-frame timestamps.
    """
    return world_metadata_from_chunks(recording_path, convert_to_seconds, verbose)["timestamp"].to_numpy()

# Given a recording path, return all frame metadata.
def world_metadata_from_chunks(recording_path: str, 
                                 convert_to_seconds: bool=True,
                                 verbose: bool=False
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
    the AGC-setting columns are filled with ``NaN``.

    Args:
        recording_path: Recording directory containing ``config.pkl`` and
            ``world*_metadata*.npy`` files.
        convert_to_seconds: Whether to convert timestamps from nanoseconds
            since boot into seconds.
        verbose: Whether to show progress while loading the metadata chunks.

    Returns:
        ``pandas.DataFrame`` with columns ``["timestamp", "Again",
        "Dgain", "exposure"]`` in frame order.
    """
    config_filepath: str = os.path.join(recording_path, "config.pkl")
    assert os.path.exists(config_filepath), f"Config filepath does not exist: {config_filepath}"
    config_data: dict | None = None
    with open(config_filepath, 'rb') as f:
        config_data = dill.load(f)
    recording_fps: float = config_data['sensors']['W']['sensor_mode']['fps']
    frame_period_ns: float = (10 ** 9) / recording_fps

    # First, let's find the world metadata chunks
    world_metadata_chunks: list[str] = natsorted([os.path.join(recording_path, filename)
                                                  for filename in os.listdir(recording_path)
                                                  if filename.startswith("world")
                                                  and "metadata" in filename
                                                 ]
                                            )
    assert len(world_metadata_chunks) > 0, f"0 world metadata chunks found @ {recording_path}"
    
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
        else:
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
    metadata: pd.DataFrame = pd.DataFrame(metadata, columns=["timestamp"] + list(WORLD_AGC_METADATA_COLS))

    return metadata


def world_raw_frames_from_chunks(path_to_recording: str, 
                                 use_mean_frame: bool=False,
                                 verbose: bool = False
                                ) -> np.ndarray:
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

def main():
    """Run the command-line entry point."""
    pass


if(__name__ == '__main__'):
    pass
