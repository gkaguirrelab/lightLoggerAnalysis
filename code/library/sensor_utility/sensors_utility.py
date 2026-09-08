"""Utilities that operate on complete light-logger recordings."""

import os
import re
import shutil
import tempfile
import numpy as np
from natsort import natsorted

import world_util
import ms_util


def group_sensors_files(recording_path: str) -> dict[str, list[tuple]]:
    """Pair per-sensor chunk metadata files with their data payloads.

    Each sensor chunk is stored as two sibling files: a metadata array whose
    filename contains ``"_metadata"`` and a value file with the same stem
    minus that marker. This helper performs that filename bookkeeping once
    and returns naturally sorted pairs for the world, pupil, and MS
    sensors.

    Args:
        recording_path: Directory containing the chunk files for one
            recording.

    Returns:
        Dictionary keyed by ``"W"``, ``"P"``, and ``"M"`` whose values are
        ordered ``(metadata_path, value_path)`` tuples.
    """
    def group_sensor_files(sensor_name: str) -> list:
        """Collect the ordered chunk-file pairs for one sensor namespace.

        Args:
            sensor_name: Filename fragment identifying the sensor, such as
                ``"world"``, ``"pupil"``, or ``"ms"``.

        Returns:
            Naturally sorted list of ``(metadata_path, value_path)`` tuples
            for that sensor.
        """
        return [ ( os.path.join(recording_path, file), os.path.join(recording_path, file.replace("_metadata", "") ) ) 
                   for file in natsorted(os.listdir(recording_path)) 
                   if sensor_name in file and "metadata" in file
               ]

    return {sensor[0].upper(): group_sensor_files(sensor) # n chunks = [ (metadata_path, frame_buffer_path), ...  ]
            for sensor in ("world", "pupil", "ms")
           }


def _process_raw_world_helper(files: tuple[str, str]) -> None:
        # Iterate over the files
        for (metadata_path, data_path) in files:
            # Load in metadata and associated frame buffer 
            metadata_buffer, frame_buffer = np.load(metadata_path), np.load(data_path) 

            # Extract just the metadata cols and format into a Python dictionary
            # Note: we have to handle old vs new column format here 
            agc_settings: dict[str, np.ndarray] = {"Again", "Dgain", "exposure"}

            # Transform the recordings with the pipeline
            assert frame_buffer.dtype == np.float64, f"Frame buffer was dtype: {frame_buffer.dtype}, expected np.float64"
            processed_buffer = world_util.world_transformation_pipeline(frame_buffer, agc_settings)

            # Output to the destination # TODO: Think about how this should look


        return
             

def _process_raw_ms_helper(files: tuple[str, str]) -> None:
        # Iterate over the files
        for (metadata_path, data_path) in files:
            # Load in metadata and associated frame buffer 
            metadata_buffer, frame_buffer = np.load(metadata_path), np.load(data_path) 

            # Transform the recordings with the pipeline
            radiance, S, f_val = ms_util.estimate_radiance_spectrum_form_ms(frame_buffer, visualize_results=False)

            # Output to the destination # TODO: Think about how this should look

def process_raw_recording(path_to_raw: str,
                          output_path: str
                        ) -> None:

    # Ensure the raw recording
    assert os.path.exists(path_to_raw) and os.path.isdir(path_to_raw) and len(os.listdir(path_to_raw)) > 0, f"{path_to_raw} must be an existing, non-empty directory"

    # Step 1: Find the chunk files for all of the different sensors 
    files_by_sensor: dict = group_sensors_files(path_to_raw)

    # Step 2: Now, let's process each of the sensors data 
    helper_map: dict[str, object] = {"W": _process_raw_world_helper, 
                                     "M": _process_raw_ms_helper
                                    }
    for sensor_name, sensor_files in files_by_sensor.items():
        helper_map[sensor_name](sensor_files, output_path)


    return 


def main() -> None:
    """Entry point placeholder."""
    pass


if(__name__ == "__main__"):
    main()
