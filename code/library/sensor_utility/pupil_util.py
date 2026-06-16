import numpy as np
import pathlib
import sys


def _load_pyagc():
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


AGC_lib_path: str = str(pathlib.Path(__file__).resolve().parents[1] / "libraries_python" / "AGC_lib")
PyAGC = _load_pyagc()

"""Define constant values associated with the Pupil Camera"""
PUPIL_CAM_FPS: int = 120 # Capture at 120 FPS
PUPIL_FRAME_SHAPE: np.ndarray = np.array([400, 400], dtype=np.uint16) # Frames from the camera are delivered at (400, 400)
PUPIL_FRAME_DTYPE: object = np.uint8 # Define the datatype of the frames captured from the pupil camera 
PUPIL_METADATA_DTYPE: object = np.float64 # Define the datatype of the metadata captured per frame (timestamps)
PUPIL_USE_AGC: int = 0 # Which AGC method to use (0:off, 1: custom)
PUPIL_AGC_MODES_INT_STR: dict[int, str] = {0: "off", 1: "custom", 2:"built-in"} # Mapping between enum types for AGC
PUPIL_AGC_MODES_STR_INT: dict[str, int] = {val: key for key, val in PUPIL_AGC_MODES_INT_STR.items()} # Mapping between enum types for AGC
PUPIL_SAVE_AGC_METADATA: bool = True # Whether or not to save all metadata from the AGC as well as timestamps
PUPIL_AGC_METADATA_COLS: tuple = ("Again", "Dgain", "exposure") # The columns of the AGC metadata
PUPIL_AGC_SETTINGS_RANGES: dict[str, np.ndarray] = PyAGC.retrieve_settings_ranges('P') if PyAGC is not None else {}
PUPIL_AGC_DISCRETE_STATES: dict[str, dict[str, int | float]] = PyAGC.retrieve_discrete_states('P') if PyAGC is not None else {}
PUPIL_AGC_ROI: tuple[tuple[int]] = ( (50, 50), (350, 350)) # Define the ROI of the pupil camera that will be used to calculate frame mean. 
                                                             # This attempts to only mean where the eye is in the image. This tuple is the 
                                                             # top left and the bottom right pixels of the box 
PUPIL_NDF_LEVEL_SETTINGS: dict[int, tuple[int, int]] = {NDF_level: [1, 1.0, 500] # Define the fixed settings for this camera per integer NDF filter
                                                        for NDF_level in range(7)
                                                       }
PUPIL_FOCAL_LENGTH: float = 561.5 # Fx = 561.471804, Fy = 562.494105

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


def main():
    pass


if(__name__ == '__main__'):
    pass
