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
from typing import Iterable
import dill
import pandas as pd
from numba import njit, prange
 
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
    WORLD_AGC_DEFAULT_TARGET: {float(NDF_level):  # Define the fixed settings for this camera per NDF filter
                               (1.0, 1.0, 747.0) if NDF_level == 0
                               else (1.0, 1.0, 7085.0) if NDF_level == 1
                               else (10.666, 1.14, 8333) if NDF_level == 2
                               else (10.666, 4.297, 8333) if NDF_level == 3
                               else (10.666, 10.0, 8333)
                               for NDF_level in range(7)}
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
        0.7: (10.333, 1.28238666, 8333),
        0.8: (10.333, 1.28885415, 8333),
        0.9: (10.333, 1.29904864, 8333),
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
WORLD_AGC_METADATA_COLS: tuple = ("Again", "Dgain", "exposure")

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
WORLD_FIELDING_FUNCTIONS: dict[tuple[int], np.ndarray] = {(480, 640): np.ones((480, 640), dtype=np.float64)}
WORLD_FIELDING_FUNCTION: np.ndarray = WORLD_FIELDING_FUNCTIONS[(480, 640)]
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

            # For red and blue sites, OpenCV bilinear demosaicing draws from the full
            # local 3x3 support when considering the union of all raw contributors to
            # that output pixel's RGB value.
            #
            # For green sites, the exact OpenCV interpolation is channel-dependent:
            # one missing color comes from the horizontal pair and the other missing
            # color comes from the vertical pair. The union of contributors is thus
            # the 5-pixel cross centered on the site. We repeat the center pixel to
            # keep the fixed output size of 9 coordinates.
            if(is_green_site):
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

                provenance_map[r, c, 5, 0] = clamp_row(source_r)
                provenance_map[r, c, 5, 1] = clamp_col(source_c)

                provenance_map[r, c, 6, 0] = clamp_row(source_r)
                provenance_map[r, c, 6, 1] = clamp_col(source_c)

                provenance_map[r, c, 7, 0] = clamp_row(source_r)
                provenance_map[r, c, 7, 1] = clamp_col(source_c)

                provenance_map[r, c, 8, 0] = clamp_row(source_r)
                provenance_map[r, c, 8, 1] = clamp_col(source_c)

            else:
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
        - Red and blue sites use the full 3x3 bilinear support.
        - Green sites use the 5-pixel cross that is the union of the
          horizontal and vertical bilinear contributors used for the two
          missing color channels.

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

def debayer_image(image: np.ndarray,
                  visualize_results: bool=False,
                  dst: np.ndarray | None=None
                ) -> np.ndarray | tuple[np.ndarray, object]:
    """Debayer a world camera image into an RGB image.

    Args:
        image: A single-channel Bayer-pattern image.
        visualize_results: If ``True``, display a before/after plot and
            return the figure handle alongside the debayered image.
        dst: Optional pre-allocated output array for in-place debayering.

    Returns:
        The debayered RGB image as a ``np.ndarray``. When
        ``visualize_results`` is ``True``, returns a tuple of
        ``(debayered_image, figure)``. When ``dst`` is provided and
        ``visualize_results`` is ``False``, returns ``None`` (the result
        is written into ``dst``).
    """
    # Initialize a variable for the figure handle that
    # will be used to visualize (if desired)
    fig: object | None = None

    # Ensure the image is correctly rounded 
    guarded_image: np.ndarray = image if image.dtype == np.uint8 else np.clip(np.round(image), 0, 255).astype(np.uint8)

    # Debayer the image according to the sensor's bayer pattern
    # either in place, or not in place depending on input args
    if(dst is None):
        debayered_image: np.ndarray = cv2.cvtColor(guarded_image, cv2.COLOR_BayerRG2RGB)
    else:
        cv2.cvtColor(guarded_image, cv2.COLOR_BayerRG2RGB, dst=dst)

    # Visualize the results if desired
    if(visualize_results is True):
        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Debayering (Before / After)", fontweight='bold', fontsize=18)
        
        # Left axis will be the "before"
        # image 
        axes[0].imshow(guarded_image, cmap="gray")
        axes[0].set_title('Before')

        # Middle axis will be the after but without 
        # color correction 
        axes[1].imshow( (debayered_image if dst is None else dst) )
        axes[1].set_title("After")

        # Show the plot 
        plt.show() 

        return debayered_image, fig
        
    return debayered_image if dst is None else None

def generate_fielding_function(image: np.ndarray) -> np.ndarray:
    """Construct a fielding correction image for a world-camera frame.

    This function is still a stub. Its intended role is to estimate a
    spatial gain map that compensates for vignetting or pixel-to-pixel
    sensitivity differences, but that calibration has not been implemented
    here yet.

    Args:
        image: Example frame whose shape determines the desired correction
            map shape.

    Returns:
        A zero-valued array with the same shape and dtype as ``image``.
        Callers should treat this as a placeholder, not a calibrated
        correction surface.
    """

    return np.zeros_like(image)

def apply_fielding_function(original_frame: np.ndarray, visualize_results: bool=False) -> tuple[np.ndarray, object] | np.ndarray:
    """Apply the registered fielding correction for this frame shape.

    The function looks up a multiplicative gain image from
    ``WORLD_FIELDING_FUNCTIONS`` using ``original_frame.shape`` as the key,
    multiplies the frame by that gain image in ``float64``, and then clips
    and rounds the result back to ``uint8``.

    Args:
        original_frame: Raw frame to correct. Its shape must match one of
            the fielding maps stored in ``WORLD_FIELDING_FUNCTIONS``.
        visualize_results: When ``True``, display a before/after figure and
            return that figure together with the corrected frame.

    Returns:
        The corrected ``uint8`` frame, or ``(frame, figure)`` when
        ``visualize_results`` is ``True``.
    """
    # Initialize a variable for the figure handle that
    # will be used to visualize (if desired)
    fig: object | None = None

    # Retrieve the fielding function for the current camera configuration
    fielding_function: np.ndarray = WORLD_FIELDING_FUNCTIONS[original_frame.shape]
    assert original_frame.shape == fielding_function.shape, f"Fielding function shape: {fielding_function.shape} is not equal to frame shape: {original_frame.shape}"

    # Apply the fielding function 
    modified_image: np.ndarray = np.clip(np.round(fielding_function.astype(np.float64) * original_frame.astype(np.float64)), 0, 255).astype(np.uint8)

    # Visualize the results if desired 
    if(visualize_results is True):
        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Fielding function (Before / After)", fontweight='bold', fontsize=18)
        
        # Left axis will be the "before"
        # image 
        axes[0].imshow(original_frame, cmap="gray")
        axes[0].set_title('Before')

        # Middle axis will be the after but without 
        # color correction 
        axes[1].imshow(modified_image)
        axes[1].set_title("After")

        # Show the plot 
        plt.show() 

        return modified_image, fig

    return modified_image

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
        colored_image: np.ndarray = np.empty(original_frame.shape, dtype=np.uint8)
        for idx, (color, pixel_coords) in enumerate(zip("RGB", (world_r_pixels, world_g_pixels, world_b_pixels))):
            rows: np.ndarray = pixel_coords[:, 0]
            cols: np.ndarray = pixel_coords[:, 1]
            color_vector: np.ndarray = np.zeros((3,), dtype=np.uint8)
            color_vector[idx] = 255
            colored_image[rows, cols] = color_vector

        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 1)

        # Middle axis will be the after but without 
        # color correction 
        axes[0].imshow(colored_image)
        axes[0].set_title(f"{original_frame.shape} Bayer RGB Pattern")

        # Show the plot 
        plt.show() 

        return mask, fig
    
    return mask

# Apply the per-color weights to the color pixels of a frame.
def apply_color_correction(original_frame: np.ndarray, 
                           visualize_results: bool=False,
                           ) -> tuple[np.ndarray, object] | np.ndarray:
    """Apply the calibrated world-camera RGB scaling factors.

    The correction constants in ``WORLD_RGB_SCALARS`` are applied in one of
    two ways depending on the input layout:

    1. For a 2-D raw Bayer frame, the function first builds the Bayer color
       mask and then scales only the pixels belonging to each color class.
    2. For a 3-D RGB image, it scales each channel plane directly.

    The computation is done in ``float64`` and then rounded and clipped back
    into the 8-bit image range.

    Args:
        original_frame: Either a raw Bayer image or a debayered RGB image.
        visualize_results: When ``True``, show a before/after figure and
            return it with the corrected frame.

    Returns:
        The corrected frame as ``uint8``, or ``(frame, figure)`` when
        visualization is requested.
    """
    # Initialize a variable for the figure handle that 
    # will be used to visualize (if desired)
    fig: object | None = None
    
    # First, we must cast the original frame to a float array 
    # to apply float scalars 
    # If it is already in float form, then no need to copy
    frame_as_float: np.ndarray = original_frame.astype(np.float64)

    # Next, we need to generate a bayer pattern for this size of image 
    # if the pixel locations have not been passed in
    is_grayscale: bool = not (len(original_frame.shape) == 3)
    mask: np.ndarray | None = None
    if(is_grayscale is True):
        mask = generate_RGB_mask(original_frame)
        assert mask.shape[:2] == original_frame.shape[:2], f"Mask: {mask.shape[:2]} and original frame shape {original_frame.shape[:2]} are unequal"

    # Next, we will apply the weights 
    for channel_idx, (color, weight) in enumerate(zip("RGB", WORLD_RGB_SCALARS)):
        # Find the pixels that match this color 
        if(is_grayscale is True):
            pixels: np.ndarray = np.argwhere(mask == color)
            rows: np.ndarray = pixels[:, 0]
            cols: np.ndarray = pixels[:, 1]

            # Apply the weight to the specified pixels 
            frame_as_float[rows, cols] *= weight

        # If debayered image, just apply to the channel number
        else:
            frame_as_float[:, :, channel_idx] *= weight    

    # Round and clip values in the 255 range and cast back to uint8
    modified_frame: np.ndarray = np.clip(np.round(frame_as_float), 0, 255).astype(np.uint8)

    # Visualize the results if desired 
    if(visualize_results is True):
        # Initialize a figure with two axes 
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle("Color correction (Before / After)", fontweight='bold', fontsize=18)
        
        # Left axis will be the "before"
        # image 
        axes[0].imshow(original_frame)
        axes[0].set_title('Before')

        # Middle axis will be the after but without 
        # color correction 
        axes[1].imshow(modified_frame)
        axes[1].set_title("After")

        # Show the plot 
        plt.show() 

        return modified_frame, fig

    # Return the modified frame
    return modified_frame


def apply_color_weights(original_frame: np.ndarray,
                        visualize_results: bool=False,
                        ) -> tuple[np.ndarray, object] | np.ndarray:
    """Backward-compatible alias for :func:`apply_color_correction`.

    Older code in this project referred to the radiometric RGB correction
    factors as "color weights." The actual implementation lives in
    ``apply_color_correction``; this wrapper preserves the older public API
    and forwards the arguments unchanged.

    Args:
        original_frame: Bayer or RGB frame to correct.
        visualize_results: Whether to request the optional before/after
            figure from ``apply_color_correction``.

    Returns:
        The corrected frame, or ``(frame, figure)`` when visualization is
        enabled.
    """
    return apply_color_correction(original_frame, visualize_results=visualize_results)

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
    metadata: pd.DataFrame = pd.DataFrame(metadata, columns=["timestamp", "Again", "Dgain", "exposure"])

    return metadata


def main():
    """Run the command-line entry point."""
    pass


if(__name__ == '__main__'):
    pass
