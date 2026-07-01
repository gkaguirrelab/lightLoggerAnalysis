"""Utilities for loading, sorting, and parsing radiometric calibration measurement files."""

import gc
import re
import os
from natsort import natsorted
import copy
import pathlib
import sys
from typing import Iterable, Literal
from tqdm.auto import tqdm

light_logger_analysis_path: str = pathlib.Path(__file__).parents[2]
chunk_io_path: str = os.path.join(light_logger_analysis_path, "library", "matlabIO", "python_libraries")

for path in (chunk_io_path,):
    sys.path.append(path)

import chunk_io

def _sort_by_setting_measurement(flattened_readings: list[str]) -> list[list[str]]:
    """Sort MS linearity measurements by setting index, then measurement index.

    Args:
        flattened_readings: Flat list of folder path strings, each containing
            'settingsIdx' and 'measurementIdx' markers.

    Returns:
        A 2D list where rows correspond to settings levels and columns
        correspond to measurement indices.
    """
    # Define the strings that represent each number to sort by
    rows_sort_string: str = "settingsIdx"
    cols_sort_string: str = "measurementIdx"

    # First, let's find the max MATLAB idx of the settings levels. This will
    # inform us how many rows we will have
    num_settings_levels: int = max([ int(re.search(rf"\d+{rows_sort_string}", reading).group()[:-len(rows_sort_string)])
                                     for reading in flattened_readings
                                   ]
                                  )

    # Then, let's find the max MATLAB idx of the measurement numbers. This will
    # inform us how many cols we will have
    num_measurements: int = max([ int(re.search(rf"\d+{cols_sort_string}", reading).group()[:-len(cols_sort_string)])
                                      for reading in flattened_readings
                                ]
                               )

    assert(all(num > 0 for num in (num_settings_levels, num_measurements)))

    # Initialize the readings matrix
    readings_matrix: list[list[str]] = []

    # Iterate over the settings levels. These will serve as the rows
    # of the matrix
    for settings_level in range(1, num_settings_levels+1):
        # Initialize a row of measurements
        measurements_row: list[str] = []

        # Iterate over the measurements and append them to the row
        for measurement_number in range(1, num_measurements+1):
            # Construct the labels for this video
            labels: set[str] = set([f"_{settings_level}{rows_sort_string}",
                                    f"_{measurement_number}{cols_sort_string}"
                                    ]
                                    )
            # Find the measurements at this contrast level and frequency
            for folder in flattened_readings:
                if(all(label in folder for label in labels)):
                    measurements_row.append(folder)

        # Append the row of measurements to the matrix
        readings_matrix.append(measurements_row)

    return readings_matrix


def _sort_by_contrast_frequency_measurement(flattened_readings: list[str]) -> list[list[list[str]]]:
    """Sort temporal sensitivity measurements by contrast, frequency, and measurement index.

    Args:
        flattened_readings: Flat list of folder path strings, each containing
            'contrastIdx', 'freqIdx', and 'measurementIdx' markers.

    Returns:
        A 3D list where the first dimension is contrast level, the second is
        frequency level, and the third is measurement index.
    """
    # Define the strings that represent each number to sort by
    rows_sort_string: str = "contrastIdx"
    cols_sort_string: str = "freqIdx"
    subcols_sort_string: str = "measurementIdx"

    # First, let's find the max MATLAB idx of the contrast levels. This will
    # inform us how many rows we will have
    num_contrast_levels: int = max([ int(re.search(rf"\d+{rows_sort_string}", reading).group()[:-len(rows_sort_string)])
                                     for reading in flattened_readings
                                   ]
                                  )

    # Then, let's find the max MATLAB idx of the frequency levels. This will
    # inform us how many cols we will have
    num_frequency_levels: int = max([ int(re.search(rf"\d+{cols_sort_string}", reading).group()[:-len(cols_sort_string)])
                                      for reading in flattened_readings
                                    ]
                                   )

    # Then, let's find the max MATLAb idx of the number of measurements at this combination of
    # contrast and frequency
    num_measurements: int = max([ int(re.search(rf"\d+{subcols_sort_string}", reading).group()[:-len(subcols_sort_string)])
                                  for reading in flattened_readings
                                ]
                               )

    assert(all(num > 0 for num in (num_contrast_levels, num_contrast_levels, num_measurements)))

    # Initialize the readings matrix
    readings_matrix: list[list[list[str]]] = []

    # Build the matrix of readings and initialize a set
    for contrast_level in range(1, num_contrast_levels+1):
        # Initialize a row for this contrast level
        # and all the frequencies recorded at it
        contrast_row: list[list[str]] = []

        # Find the frequencies at this contrast level
        for frequency_level in range(1, num_frequency_levels+1):
            # Initialize a sub matrix for this
            frequency_mat: list[str] = []

            # Iterate over the measurements at this contrast and frequency level
            for measurement_number in range(1, num_measurements+1):
                # Construct the labels for this video
                labels: set[str] = set([f"_{contrast_level}{rows_sort_string}",
                                        f"_{frequency_level}{cols_sort_string}",
                                        f"_{measurement_number}{subcols_sort_string}"
                                       ]
                                      )

                # Find the measurements at this contrast level and frequency
                for folder in flattened_readings:
                    if(all(label in folder for label in labels)):
                        frequency_mat.append(folder)

            # Append the frequency mat to the contrast row
            contrast_row.append(frequency_mat)

        # Append the row of matrices to the outer matrix
        readings_matrix.append(contrast_row)

    return readings_matrix


def _parse_ms_linearity_readings(folders: list[list[str]],
                                 apply_digital_gain: bool=False,
                                 apply_radiometric_correction: bool=False,
                                 use_mean_frame: bool=False,
                                 convert_time_units: bool=False,
                                 convert_to_float: bool=False,
                                 mean_axes: dict[str, tuple]= {'W': (1, 2), 'P': (1, 2), 'M': (0,)},
                                 verbose: Literal["tqdm", "text", "off"] = "off",
                                 differentiate_color: bool=False
                                 ) -> list[list[dict]]:
    """Parse MS linearity readings into Python from their folder paths.

    Args:
        folders: 2D list of folder path strings, sorted by settings level
            and measurement index.
        apply_digital_gain: Whether to apply digital gain correction.
        apply_radiometric_correction: Whether to apply RGB radiometric
            correction.
        use_mean_frame: Whether to average frames.
        convert_time_units: Whether to convert time units.
        convert_to_float: Whether to convert data to float.
        mean_axes: Dictionary mapping sensor keys to axes over which to
            compute the mean.
        verbose: Verbosity level. 'tqdm' for progress bars, 'text' for
            printed progress, 'off' for silent.
        differentiate_color: Whether to differentiate color channels.

    Returns:
        A 2D list of parsed measurement dictionaries with the same
        dimensional structure as the input folders.
    """
    # Initialize a return list
    parsed_folders = copy.deepcopy(folders)

    n_settings: int = len(folders)

    # Iterate over the settings levels (the rows)
    settings_iterator: Iterable = tqdm(range(n_settings), desc="Parsing MS linearity settings levels", leave=False) if verbose == "tqdm" else range(n_settings)
    for settings_level in settings_iterator:
        if verbose == "text":
            print(f"  Settings level {settings_level + 1}/{n_settings}")

        n_measurements: int = len(folders[settings_level])
        # Iterate over the cols (measurement number)
        measurement_iterator: Iterable = tqdm(range(n_measurements), desc="Parsing measurements", leave=False) if verbose == "tqdm" else range(n_measurements)
        for measurement_number in measurement_iterator:
            if verbose == "text":
                print(f"    Measurement {measurement_number + 1}/{n_measurements}")
            parsed_folders[settings_level][measurement_number] = chunk_io.parse_chunks(folders[settings_level][measurement_number],
                                                                              apply_RGB_correction=apply_radiometric_correction,
                                                                              apply_digital_gain=apply_digital_gain,
                                                                              use_mean_frame=use_mean_frame,
                                                                              convert_time_units=convert_time_units,
                                                                              convert_to_float=convert_to_float,
                                                                              mean_axes=mean_axes,
                                                                              differentiate_color=differentiate_color
                                                                             )
            gc.collect()

    return parsed_folders


def _sort_by_contrast_target_settings_measurement(folders: list[str]) -> list[list[list]]:
    """Sort camera linearity measurement folder paths by contrast target, settings, and measurement.

    Args:
        folders: Flat list of folder path strings, each containing
            'contrastTargetIdx', 'settingsIdx', and 'measurementIdx' markers.

    Returns:
        A 3D list where the first dimension is contrast target, the second
        is settings index, and the third is measurement index.
    """
    # Define the output 
    sorted_folders: list[list[list]] = []

    # Let's first find the dimensions of the measurement
    agc_contrast_target_pattern: str = r"(\d+)contrastTargetIdx"
    num_contrast_targets: int = max([int(re.search(agc_contrast_target_pattern, foldername).group(1)) for foldername in folders]) 

    settings_idx_pattern: str = r"_(\d+)settingsIdx"
    num_settings_idx: int = max([int(re.search(settings_idx_pattern, foldername).group(1)) for foldername in folders]) 

    measurement_idx_pattern: str = r"_(\d+)measurementIdx"
    num_measurements: int = max([int(re.search(measurement_idx_pattern, foldername).group(1)) for foldername in folders]) 

    # Iterate over the measurements 
    for contrast_target_num in range(1, num_contrast_targets + 1):
        contrast_target_matrix: list[list] = []
        for settings_num in range(1, num_settings_idx + 1):
            settings_row: list[str] = []
            for measurement_num in range(1, num_measurements + 1):

                # Find the measurement that has this contrast target num, 
                # settingsIdx and measurement num 
                idx_markers: tuple[str] = (f"{contrast_target_num}contrastTargetIdx", f"{settings_num}settingsIdx", f"{measurement_num}measurementIdx")
                measurement_path: str = [folderpath for folderpath in folders if all(marker in folderpath for marker in idx_markers)][0]

                # Save this measurement
                settings_row.append(measurement_path)

            # Save this row of measurements at this settings level 
            contrast_target_matrix.append(settings_row)

        # Save this contrast target matrix to the sorted folders 
        sorted_folders.append(contrast_target_matrix)

    return sorted_folders

def _parse_camera_linearity_readings(folders: list[list[list[str]]],
                                     verbose: Literal["tqdm", "text", "off"] = "off",
                                     apply_radiometric_correction: bool=False,
                                     apply_digital_gain: bool=False,
                                     convert_time_units: bool=False,
                                     differentiate_color: bool=False
                                    ) -> None:
    """Parse camera linearity readings into Python from their folder paths.

    Args:
        folders: 3D list of folder path strings, sorted by contrast target,
            settings index, and measurement index.
        verbose: Verbosity level. 'tqdm' for progress bars, 'text' for
            printed progress, 'off' for silent.
        apply_radiometric_correction: Whether to apply RGB radiometric
            correction.
        apply_digital_gain: Whether to apply digital gain correction.
        convert_time_units: Whether to convert time units.
        differentiate_color: Whether to differentiate color channels.

    Returns:
        A 3D list of parsed measurement dictionaries with the same
        dimensional structure as the input folders.
    """
    # Deifne the output
    parsed_folders: list[list[list[str]]] | list[list[list[dict]]] = copy.deepcopy(folders)

    n_contrast_targets: int = len(parsed_folders)

    # Now, let's iterate over the folders
    contrast_target_iterator: Iterable = tqdm(range(n_contrast_targets), desc="Parsing camera linearity contrast targets", leave=False) if verbose == "tqdm" else range(n_contrast_targets)
    for contrast_target_idx in contrast_target_iterator:
        if verbose == "text":
            print(f"  Contrast target {contrast_target_idx + 1}/{n_contrast_targets}")

        contrast_target_matrix = parsed_folders[contrast_target_idx]
        n_settings: int = len(contrast_target_matrix)
        settings_iterator: Iterable = tqdm(range(n_settings), desc="Parsing camera linearity settings levels", leave=False) if verbose == "tqdm" else range(n_settings)
        for settings_idx in settings_iterator:
            if verbose == "text":
                print(f"    Settings level {settings_idx + 1}/{n_settings}")

            settings_row = contrast_target_matrix[settings_idx]
            n_measurements: int = len(settings_row)
            measurement_iterator: Iterable = tqdm(range(n_measurements), desc="Parsing measurements", leave=False) if verbose == "tqdm" else range(n_measurements)
            for measurement_num in measurement_iterator:
                if verbose == "text":
                    print(f"      Measurement {measurement_num + 1}/{n_measurements}")
                measurement_path = settings_row[measurement_num]
                parsed_folders[contrast_target_idx][settings_idx][measurement_num] = chunk_io.parse_chunks(measurement_path,
                                                                                                           apply_RGB_correction=apply_radiometric_correction,
                                                                                                           apply_digital_gain=apply_digital_gain,
                                                                                                           convert_time_units=convert_time_units,
                                                                                                           differentiate_color=differentiate_color
                                                                                                        )
                gc.collect()

    return parsed_folders

def _parse_temporal_sensitivity_readings(folders: list[list[str]],
                                         apply_digital_gain: bool=False,
                                         apply_radiometric_correction: bool=False,
                                         use_mean_frame: bool=False,
                                         convert_time_units: bool=False,
                                         convert_to_float: bool=False,
                                         mean_axes: dict[str, tuple]= {'W': (1, 2), 'P': (1, 2), 'M': (0,)},
                                         contains_agc_metadata_dict: dict[str, bool]={'W': True, 'P': False, 'M': False},
                                         verbose: Literal["tqdm", "text", "off"] = "off",
                                         differentiate_color: bool=False
                                        ) -> list[list[dict]]:
    """Parse temporal sensitivity readings into Python from their folder paths.

    Args:
        folders: 3D list of folder path strings, sorted by contrast level,
            frequency level, and measurement index.
        apply_digital_gain: Whether to apply digital gain correction.
        apply_radiometric_correction: Whether to apply RGB radiometric
            correction.
        use_mean_frame: Whether to average frames.
        convert_time_units: Whether to convert time units.
        convert_to_float: Whether to convert data to float.
        mean_axes: Dictionary mapping sensor keys to axes over which to
            compute the mean.
        contains_agc_metadata_dict: Dictionary mapping sensor keys to
            whether they contain AGC metadata.
        verbose: Verbosity level. 'tqdm' for progress bars, 'text' for
            printed progress, 'off' for silent.
        differentiate_color: Whether to differentiate color channels.

    Returns:
        A 3D list of parsed measurement dictionaries with the same
        dimensional structure as the input folders.
    """
    # Initialize a return list
    parsed_folders: list[list[str]] | list[list[dict]] = copy.deepcopy(folders)

    n_contrasts: int = len(folders)

    # Iterate over the folders and load them in
    contrast_iterator: Iterable = tqdm(range(n_contrasts), desc="Parsing contrast levels", leave=False) if verbose == "tqdm" else range(n_contrasts)
    for contrast_level in contrast_iterator:
        if verbose == "text":
            print(f"  Contrast level {contrast_level + 1}/{n_contrasts}")

        n_frequencies: int = len(folders[contrast_level])
        frequency_iterator: Iterable = tqdm(range(n_frequencies), desc="Parsing frequency levels", leave=False) if verbose == "tqdm" else range(n_frequencies)
        for frequency_level in frequency_iterator:
            if verbose == "text":
                print(f"    Frequency level {frequency_level + 1}/{n_frequencies}")

            n_measurements: int = len(parsed_folders[contrast_level][frequency_level])
            measurement_results: list = []
            for measurement_num in range(n_measurements):
                if verbose == "text":
                    print(f"      Measurement {measurement_num + 1}/{n_measurements}")
                measurement_results.append(
                    chunk_io.parse_chunks(folders[contrast_level][frequency_level][measurement_num],
                                          apply_RGB_correction=apply_radiometric_correction,
                                          apply_digital_gain=apply_digital_gain,
                                          use_mean_frame=use_mean_frame,
                                          convert_time_units=convert_time_units,
                                          convert_to_float=convert_to_float,
                                          mean_axes=mean_axes,
                                          contains_agc_metadata_dict=contains_agc_metadata_dict,
                                          differentiate_color=differentiate_color
                                         )
                )
                gc.collect()
            parsed_folders[contrast_level][frequency_level] = measurement_results

    return parsed_folders

def load_sorted_calibration_files(experiment_path: str,
                                  apply_digital_gain: bool=False,
                                  use_mean_frame: bool=False,
                                  apply_radiometric_correction: bool=False,
                                  convert_time_units: bool=False,
                                  convert_to_float: bool=False,
                                  mean_axes: dict[str, tuple]= {'W': (1, 2), 'P': (1, 2), 'M': (0,)},
                                  verbose: Literal["tqdm", "text", "off"] = "off",
                                  differentiate_color: bool=False,
                                  parse_files: bool=True
                                  ) -> dict[str, list]:
    """Load and return sorted calibration measurements after decompressing and decrypting.

    Given a directory to a measurement (e.g. "W1P2M3/NDF_0/"), returns sorted
    lists of the measurements per calibration operation. Note: this sorts them
    by the index of their occurrence. If the stimulus was not exposed in the
    same order per measurement, the frequency and corresponding contrast will
    not match.

    Args:
        experiment_path: Path to the experiment directory containing NDF
            subfolders.
        apply_digital_gain: Whether to apply digital gain correction.
        use_mean_frame: Whether to average frames.
        apply_radiometric_correction: Whether to apply RGB radiometric
            correction.
        convert_time_units: Whether to convert time units.
        convert_to_float: Whether to convert data to float.
        mean_axes: Dictionary mapping sensor keys to axes over which to
            compute the mean.
        verbose: Verbosity level. 'tqdm' for progress bars, 'text' for
            printed progress, 'off' for silent.
        differentiate_color: Whether to differentiate color channels.
        parse_files: If True, parse the files into Python dictionaries.
            If False, return only the sorted folder paths.

    Returns:
        A dictionary with keys 'ms_linearity', 'temporal_sensitivity',
        'phase_fitting', 'contrast_gamma', and 'world_linearity'. Values
        are nested lists of parsed measurement dictionaries (if parse_files
        is True) or sorted folder path strings (if parse_files is False).
    """    
    
    # First, let's find the NDFs for this experiment 
    NDF_folders: list[str] = [os.path.join(experiment_path, folder) 
                              for folder in natsorted(os.listdir(experiment_path))
                              if 'NDF' in folder
                             ]

    assert len(NDF_folders) > 0, "Must be at least a single NDF folder"

    n_NDFs: int = len(NDF_folders)

    # First, we will group folders together by the calibration operation and the NDF level it occured at,
    # then sort them into their proper dimensional structure

    ##################
    # MS LINEARITY
    ##################

    ms_linearity_folders: list[str] = [ [os.path.join(NDF_folder, folder)
                                         for folder in os.listdir(NDF_folder)
                                         if "mslinearity" in folder.lower()
                                        ]
                                        for NDF_folder in NDF_folders
                                      ]

    ms_linearity_folders_sorted: list[list[str]] = [ _sort_by_setting_measurement(NDF_measurement)
                                                     if len(NDF_measurement) > 0
                                                     else
                                                     []
                                                     for NDF_measurement in ms_linearity_folders
                                                   ]

    ##################
    # WORLD CAMERA LINEARITY
    ##################

    camera_linearity_folders: list[str] = [ [os.path.join(NDF_folder, folder)
                                             for folder in os.listdir(NDF_folder)
                                             if folder.lower().startswith("worldcameralinearity")
                                            ]
                                            for NDF_folder in NDF_folders
                                          ]

    camera_linearity_folders_sorted: list[list[list[str]]] = [_sort_by_contrast_target_settings_measurement(NDF_folder)
                                                              if len(NDF_folder) > 0
                                                              else []
                                                              for NDF_folder in camera_linearity_folders
                                                            ]

    ##################
    # TEMPORAL SENSITIVITY
    ##################

    temporal_sensitivity_folders: list[str] = [ [os.path.join(NDF_folder, folder)
                                                 for folder in os.listdir(NDF_folder)
                                                 if "temporalsensitivity" in folder.lower()
                                                 and "phasefitting" not in folder.lower()
                                                 and "contrastgamma" not in folder.lower()
                                                ]
                                                for NDF_folder in NDF_folders
                                              ]

    temporal_sensitivity_folders_sorted: list[list[list[str]]] = [ _sort_by_contrast_frequency_measurement(NDF_measurement)
                                                                   if len(NDF_measurement) > 0
                                                                   else
                                                                   []
                                                                   for NDF_measurement in temporal_sensitivity_folders
                                                                 ]

    ##################
    # PHASE FITTING
    ##################

    phase_fitting_folders: list[str] = [ [ os.path.join(NDF_folder, folder)
                                           for folder in os.listdir(NDF_folder)
                                           if "phasefitting" in folder.lower()
                                         ]
                                         for NDF_folder in NDF_folders
                                       ]

    phase_fitting_folders_sorted: list[list[list[str]]] = [ _sort_by_contrast_frequency_measurement(NDF_measurement)
                                                            if len(NDF_measurement) > 0
                                                            else
                                                            []
                                                            for NDF_measurement in phase_fitting_folders
                                                          ]

    ##################
    # CONTRAST GAMMA
    ##################

    contrast_gamma_folders: list[str] = [ [os.path.join(NDF_folder, folder)
                                           for folder in os.listdir(NDF_folder)
                                           if "contrastgamma" in folder.lower()
                                          ]
                                          for NDF_folder in NDF_folders
                                        ]

    contrast_gamma_folders_sorted: list[list[list[str]]] = [ _sort_by_contrast_frequency_measurement(NDF_measurement)
                                                             if len(NDF_measurement) > 0
                                                             else
                                                             []
                                                             for NDF_measurement in contrast_gamma_folders
                                                           ]

    ##################
    # RETURN SORTED PATHS OR PARSED DATA
    ##################

    if not parse_files:
        return {"ms_linearity": ms_linearity_folders_sorted,
                "temporal_sensitivity": temporal_sensitivity_folders_sorted,
                "phase_fitting": phase_fitting_folders_sorted,
                "contrast_gamma": contrast_gamma_folders_sorted,
                "world_linearity": camera_linearity_folders_sorted
               }

    # If parse_files is True, parse all the data in Python and return the parsed readings

    if verbose == "text":
        print(f"[MS Linearity]")

    ms_linearity_readings: list[list[dict]] = []
    for ndf_idx, NDF_measurement in enumerate(ms_linearity_folders_sorted):
        if verbose == "text":
            print(f"NDF {ndf_idx + 1}/{n_NDFs}")
        if len(NDF_measurement) > 0:
            ms_linearity_readings.append(
                _parse_ms_linearity_readings(NDF_measurement,
                                             apply_digital_gain=apply_digital_gain,
                                             apply_radiometric_correction=apply_radiometric_correction,
                                             use_mean_frame=use_mean_frame,
                                             convert_time_units=convert_time_units,
                                             convert_to_float=convert_to_float,
                                             mean_axes=mean_axes,
                                             verbose=verbose,
                                             differentiate_color=differentiate_color
                                            )
            )
        else:
            ms_linearity_readings.append([])
    gc.collect()

    if verbose == "text":
        print(f"[World Camera Linearity]")

    camera_linearity_readings: list[list[list[dict]]] = []
    for ndf_idx, NDF_folder in enumerate(camera_linearity_folders_sorted):
        if verbose == "text":
            print(f"NDF {ndf_idx + 1}/{n_NDFs}")
        camera_linearity_readings.append(
            _parse_camera_linearity_readings(NDF_folder,
                                             verbose=verbose,
                                             apply_radiometric_correction=apply_radiometric_correction,
                                             apply_digital_gain=apply_digital_gain,
                                             convert_time_units=convert_time_units,
                                             differentiate_color=differentiate_color
                                            )
        )
    gc.collect()

    if verbose == "text":
        print(f"[Temporal Sensitivity]")

    temporal_sensitivity_readings: list[list[list[dict]]] = []
    for ndf_idx, NDF_measurement in enumerate(temporal_sensitivity_folders_sorted):
        if verbose == "text":
            print(f"NDF {ndf_idx + 1}/{n_NDFs}")
        if len(NDF_measurement) > 0:
            temporal_sensitivity_readings.append(
                _parse_temporal_sensitivity_readings(NDF_measurement,
                                                     apply_digital_gain=apply_digital_gain,
                                                     apply_radiometric_correction=apply_radiometric_correction,
                                                     use_mean_frame=use_mean_frame,
                                                     convert_time_units=convert_time_units,
                                                     convert_to_float=convert_to_float,
                                                     mean_axes=mean_axes,
                                                     contains_agc_metadata_dict={'W': True, 'P': False, 'M': False},
                                                     verbose=verbose,
                                                     differentiate_color=differentiate_color
                                                    )
            )
        else:
            temporal_sensitivity_readings.append([])
    gc.collect()

    if verbose == "text":
        print(f"[Phase Fitting]")

    phase_fitting_readings: list[list[list[dict]]] = []
    for ndf_idx, NDF_measurement in enumerate(phase_fitting_folders_sorted):
        if verbose == "text":
            print(f"NDF {ndf_idx + 1}/{n_NDFs}")
        if len(NDF_measurement) > 0:
            phase_fitting_readings.append(
                _parse_temporal_sensitivity_readings(NDF_measurement,
                                                     apply_digital_gain=apply_digital_gain,
                                                     apply_radiometric_correction=apply_radiometric_correction,
                                                     use_mean_frame=use_mean_frame,
                                                     convert_time_units=convert_time_units,
                                                     convert_to_float=convert_to_float,
                                                     mean_axes=mean_axes,
                                                     contains_agc_metadata_dict={'W': True, 'P': False, 'M': False},
                                                     verbose=verbose,
                                                     differentiate_color=differentiate_color
                                                    )
            )
        else:
            phase_fitting_readings.append([])
    gc.collect()

    if verbose == "text":
        print(f"[Contrast Gamma]")

    contrast_gamma_readings: list[list[list[dict]]] = []
    for ndf_idx, NDF_measurement in enumerate(contrast_gamma_folders_sorted):
        if verbose == "text":
            print(f"NDF {ndf_idx + 1}/{n_NDFs}")
        if len(NDF_measurement) > 0:
            contrast_gamma_readings.append(
                _parse_temporal_sensitivity_readings(NDF_measurement,
                                                     apply_digital_gain=apply_digital_gain,
                                                     apply_radiometric_correction=apply_radiometric_correction,
                                                     use_mean_frame=use_mean_frame,
                                                     convert_time_units=convert_time_units,
                                                     convert_to_float=convert_to_float,
                                                     mean_axes=mean_axes,
                                                     contains_agc_metadata_dict={'W': True, 'P': False, 'M': False},
                                                     verbose=verbose,
                                                     differentiate_color=differentiate_color
                                                    )
            )
        else:
            contrast_gamma_readings.append([])
    gc.collect()

    return {"ms_linearity": ms_linearity_readings,
            "temporal_sensitivity": temporal_sensitivity_readings,
            "phase_fitting": phase_fitting_readings,
            "contrast_gamma": contrast_gamma_readings,
            "world_linearity": camera_linearity_readings
           }


if(__name__ == "__main__"):
    experiment_path: str = "/Volumes/T7 Shield/5cameraLinearity" 
    results = load_sorted_calibration_files(experiment_path, apply_digital_gain=True, use_mean_frame=True, convert_time_units=True, convert_to_float=True, verbose="tqdm")
