import re 
import os
from natsort import natsorted
import copy
import pathlib
import sys
from typing import Iterable
from tqdm.auto import tqdm

light_logger_analysis_path: str = pathlib.Path(__file__).parents[2]
chunk_io_path: str = os.path.join(light_logger_analysis_path, "library", "matlabIO", "python_libraries")

for path in (chunk_io_path,):
    sys.path.append(path)

import chunk_io

"""Define a subfunction to sort MS linearity measurements
   first by the setting idx of the recording and then
   the measurement idx of that setting
"""
def _sort_by_setting_measurement(flattened_readings: list[str]) -> list[list[str]]:
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


"""Define a subfunction to sort any temporal sensitivity
   measurement first by the contrast level and then
   by the frequency. Out of this, you get a
   2D list where the rows are the contrast
   and the cols are the frequencies
"""
def _sort_by_contrast_frequency_measurement(flattened_readings: list[str]) -> list[list[list[str]]]:
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


"""Define a subfunction to parse the MS linearity readings
   into Python from their folder paths
"""
def _parse_ms_linearity_readings(folders: list[list[str]],
                                 apply_digital_gain: bool=False,
                                 apply_radiometric_correction: bool=False,
                                 use_mean_frame: bool=False,
                                 convert_time_units: bool=False,
                                 convert_to_float: bool=False,
                                 mean_axes: dict[str, tuple]= {'W': (1, 2), 'P': (1, 2), 'M': (0,)},
                                 verbose: bool = False
                                 ) -> list[list[dict]]:
    # Initialize a return list
    parsed_folders = copy.deepcopy(folders)

    # Iterate over the settings levels (the rows)
    settings_iterator: Iterable = range(len(folders)) if verbose is False else tqdm(range(len(folders)), desc="Parsing MS linearity settings levels", leave=False)
    for settings_level in range(len(folders)):
        # Iterate over the cols (measurement number)
        measurement_iterator: Iterable = range(len(folders[settings_level])) if verbose is False else tqdm(range(len(folders)), desc="Parsing measurements", leave=False)
        for measurement_number in measurement_iterator:
            parsed_folders[settings_level][measurement_number] = chunk_io.parse_chunks(folders[settings_level][measurement_number],
                                                                              apply_RGB_correction=apply_radiometric_correction,
                                                                              apply_digital_gain=apply_digital_gain,
                                                                              use_mean_frame=use_mean_frame,
                                                                              convert_time_units=convert_time_units,
                                                                              convert_to_float=convert_to_float,
                                                                              mean_axes=mean_axes
                                                                             )

    return parsed_folders


"""Define a helper function to sort the camera linearity measurements folderpaths"""
def _sort_by_contrast_target_settings_measurement(folders: list[str]) -> list[list[list]]:
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

"""Define a subfunction to parse the camera linearity readings into Python from their folder paths"""
def _parse_camera_linearity_readings(folders: list[list[list[str]]],
                                     verbose: bool = False, 
                                     apply_radiometric_correction: bool=False,
                                     apply_digital_gain: bool=False, 
                                     convert_time_units: bool=False
                                    ) -> None:
    # Deifne the output 
    parsed_folders: list[list[list[str]]] | list[list[list[dict]]] = copy.deepcopy(folders)

    # Now, let's iterate over the folders 
    contrast_target_iterator: Iterable = range(len(parsed_folders)) if verbose is False else tqdm(range(len(parsed_folders)), desc="Parsing camera linearity contrast targets", leave=False)
    for contrast_target_idx in contrast_target_iterator:
        contrast_target_matrix = parsed_folders[contrast_target_idx]
        settings_iterator: Iterable = range(len(contrast_target_matrix)) if verbose is False else tqdm(range(len(contrast_target_matrix)), desc="Parsing camera linearity settings levels", leave=False)
        for settings_idx in settings_iterator:
            settings_row = contrast_target_matrix[settings_idx]
            measurement_iterator: Iterable = range(len(settings_row)) if verbose is False else tqdm(range(len(settings_row)), desc="Parsing measurements", leave=False)
            for measurement_num in measurement_iterator:
                measurement_path = settings_row[measurement_num]
                parsed_folders[contrast_target_idx][settings_idx][measurement_num] = chunk_io.parse_chunks(measurement_path, 
                                                                                                           apply_RGB_correction=apply_radiometric_correction, 
                                                                                                           apply_digital_gain=apply_digital_gain,
                                                                                                           convert_time_units=convert_time_units
                                                                                                        )

    return parsed_folders

"""Define a subfunction to parse the temporal sensitivity
   readings into Python from their folder paths
"""
def _parse_temporal_sensitivity_readings(folders: list[list[str]],
                                         apply_digital_gain: bool=False,
                                         apply_radiometric_correction: bool=False,
                                         use_mean_frame: bool=False,
                                         convert_time_units: bool=False,
                                         convert_to_float: bool=False,
                                         mean_axes: dict[str, tuple]= {'W': (1, 2), 'P': (1, 2), 'M': (0,)},
                                         contains_agc_metadata_dict: dict[str, bool]={'W': True, 'P': False, 'M': False},
                                         verbose: bool = False
                                        ) -> list[list[dict]]:
    # Initialize a return list
    parsed_folders: list[list[str]] | list[list[dict]] = copy.deepcopy(folders)

    # Iterate over the folders and load them in
    contrast_iterator: Iterable = range(len(folders)) if verbose is False else tqdm(range(len(folders)), desc="Parsing contrast levels", leave=False)
    for contrast_level in contrast_iterator:
        frequency_iterator: Iterable = range(len(folders[contrast_level])) if verbose is False else tqdm(range(len(folders[contrast_level])), desc="Parsing frequency levels", leave=False)
        for frequency_level in frequency_iterator:
            parsed_folders[contrast_level][frequency_level] = [ chunk_io.parse_chunks(folders[contrast_level][frequency_level][measurement_num],
                                                                           apply_RGB_correction=apply_radiometric_correction,
                                                                           apply_digital_gain=apply_digital_gain,
                                                                           use_mean_frame=use_mean_frame,
                                                                           convert_time_units=convert_time_units,
                                                                           convert_to_float=convert_to_float,
                                                                           mean_axes=mean_axes,
                                                                           contains_agc_metadata_dict=contains_agc_metadata_dict
                                                                          )
                                                                for measurement_num in range(len(parsed_folders[contrast_level][frequency_level]))
                                                              ]

    return parsed_folders

"""Given a directory to a measurement (e.g. "W1P2M3/NDF_0/),
   return sorted lists of the measurements per 
   calibration operation after decompressing 
   and decrypting. Note: this sorts them by 
   the index of their occurence. If the stimulate was not 
   exposed in the same order per measurement for instance, 
   the frequency and corresponding contrast will not be the same
"""
def load_sorted_calibration_files(experiment_path: str, 
                                  apply_digital_gain: bool=False, 
                                  use_mean_frame: bool=False, 
                                  apply_radiometric_correction: bool=False,
                                  convert_time_units: bool=False, 
                                  convert_to_float: bool=False, 
                                  mean_axes: dict[str, tuple]= {'W': (1, 2), 'P': (1, 2), 'M': (0,)},
                                  verbose: bool = False
                                  ) -> dict[str, list[dict] | list[list[dict]]]:    
    
    # First, let's find the NDFs for this experiment 
    NDF_folders: list[str] = [os.path.join(experiment_path, folder) 
                              for folder in natsorted(os.listdir(experiment_path))
                              if 'NDF' in folder
                             ]

    assert len(NDF_folders) > 0, "Must be at least a single NDF folder"


    # First, we will group folders together by the calibration operation and the NDF level it occured at 

    ##################

    # MS LINEARITY

    ##################

    # First, we will group the folders together by the calibration 
    # operation that was done to create them, then parse 
    ms_linearity_folders: list[str] = [ [os.path.join(NDF_folder, folder)  
                                         for folder in os.listdir(NDF_folder)
                                         if "mslinearity" in folder.lower()
                                        ]
                                        for NDF_folder in NDF_folders
                                      ]

    # We will sort the MS linearity folders into a matrix of the settings level and the measurement numbers per NDF
    ms_linearity_folders_sorted: list[list[str]] = [ _sort_by_setting_measurement(NDF_measurement)
                                                     if len(NDF_measurement) > 0 
                                                     else 
                                                     []
                                                     for NDF_measurement in ms_linearity_folders
                                                   ]

    # Parse the MS linearity readings
    ms_linearity_readings: list[list[dict]] = [ _parse_ms_linearity_readings(NDF_measurement,
                                                                              apply_digital_gain=apply_digital_gain,
                                                                              apply_radiometric_correction=apply_radiometric_correction,
                                                                              use_mean_frame=use_mean_frame,
                                                                              convert_time_units=convert_time_units,
                                                                              convert_to_float=convert_to_float,
                                                                              mean_axes=mean_axes,
                                                                              verbose=verbose
                                                                             )
                                                if len(NDF_measurement) > 0 
                                                else 
                                                []
                                                for NDF_measurement in ms_linearity_folders_sorted
                                              ] 
    


    ##################

    # WORLD CAMERA LINEARITY

    ##################
    # First, gather all of the calibration measurements
    # for this type of calibration
    camera_linearity_folders: list[str] = [ [os.path.join(NDF_folder, folder)
                                             for folder in os.listdir(NDF_folder)
                                             if folder.lower().startswith("worldcameralinearity")
                                            ]
                                            for NDF_folder in NDF_folders
                                          ]

    # Next, we will sort these measurements per NDF into contrastAGCTarget, settingsIdx, measurementNum 
    camera_linearity_folders_sorted: list[list[list[str]]] = [_sort_by_contrast_target_settings_measurement(NDF_folder)
                                                              if len(NDF_folder) > 0 
                                                              else []
                                                              for NDF_folder in camera_linearity_folders
                                                            ] 

    # Next, parse the readings 
    camera_linearity_readings: list[list[list[dict]]] = [_parse_camera_linearity_readings(NDF_folder,
                                                                                           verbose=verbose,
                                                                                           apply_radiometric_correction=apply_radiometric_correction, 
                                                                                           apply_digital_gain=apply_digital_gain,
                                                                                           convert_time_units=convert_time_units
                                                                                          )
                                                         for NDF_folder in camera_linearity_folders_sorted
                                                        ]


    ##################

    # TEMPORAL SENSITIVITY

    ##################

    # Gather the paths to the temporal sensitivity measurements
    temporal_sensitivity_folders: list[str] = [ [os.path.join(NDF_folder, folder)
                                                 for folder in os.listdir(NDF_folder)
                                                 if "temporalsensitivity" in folder.lower()
                                                 and "phasefitting" not in folder.lower()
                                                 and "contrastgamma" not in folder.lower() 
                                                ] 
                                                for NDF_folder in NDF_folders
                                              ]

    # Sort each NDF measurment into a matrix of Contrast, Frequency, and Measurement for each NDF level 
    temporal_sensitivity_folders_sorted: list[list[list[str]]] =  [ _sort_by_contrast_frequency_measurement(NDF_measurement)
                                                                    if len(NDF_measurement) > 0 
                                                                    else 
                                                                    []
                                                                    for NDF_measurement in temporal_sensitivity_folders
                                                                  ] 
    

    # Parse the readings
    temporal_sensitivity_readings: list[list[list[dict]]] = [ _parse_temporal_sensitivity_readings(NDF_measurement,
                                                                                                     apply_digital_gain=apply_digital_gain,
                                                                                                     apply_radiometric_correction=apply_radiometric_correction,
                                                                                                     use_mean_frame=use_mean_frame,
                                                                                                     convert_time_units=convert_time_units,
                                                                                                     convert_to_float=convert_to_float,
                                                                                                     mean_axes=mean_axes,
                                                                                                     contains_agc_metadata_dict={'W': True, 'P': False, 'M': False},
                                                                                                     verbose=verbose
                                                                                                    )
                                                              if len(NDF_measurement) > 0 
                                                              else 
                                                              []
                                                              for NDF_measurement in temporal_sensitivity_folders_sorted
                                                            ]
    


    ##################

    # PHASE FITTING

    ##################

    # Find the phase fitting folders
    phase_fitting_folders: list[str] = [ [ os.path.join(NDF_folder, folder)
                                           for folder in os.listdir(NDF_folder)
                                           if "phasefitting" in folder.lower()
                                         ] 
                                         for NDF_folder in NDF_folders
                                       ]
    
    # Sort the folders per NDF level 
    phase_fitting_folders_sorted: list[list[list[str]]] = [ _sort_by_contrast_frequency_measurement(NDF_measurement)
                                                            if len(NDF_measurement) > 0 
                                                            else 
                                                            []
                                                            for NDF_measurement in phase_fitting_folders
                                                          ]

    # Parse the the phase fitting readings
    # TODO: check the types of these are correct
    phase_fitting_readings: list[list[list[dict]]] = [ _parse_temporal_sensitivity_readings(NDF_measurement,
                                                                                             apply_digital_gain=apply_digital_gain,
                                                                                             apply_radiometric_correction=apply_radiometric_correction,
                                                                                             use_mean_frame=use_mean_frame,
                                                                                             convert_time_units=convert_time_units,
                                                                                             convert_to_float=convert_to_float,
                                                                                             mean_axes=mean_axes,
                                                                                             contains_agc_metadata_dict={'W': True, 'P': False, 'M': False},
                                                                                             verbose=verbose
                                                                                            )
                                                       if len(NDF_measurement) > 0
                                                       else 
                                                       []
                                                       for NDF_measurement in phase_fitting_folders_sorted
                                                    ]


    ##################

    # CONTRAST GAMMA

    ##################


    # Find the contrast gamma folders                                                                                
    contrast_gamma_folders: list[str] =  [ [os.path.join(NDF_folder, folder)
                                            for folder in os.listdir(NDF_folder)
                                            if "contrastgamma" in folder.lower()
                                           ] 
                                            for NDF_folder in NDF_folders
                                         ] 

    # Sort the folder paths per NDF level
    contrast_gamma_folders_sorted: list[list[list[str]]] = [ _sort_by_contrast_frequency_measurement(NDF_measurement)
                                                             if len(NDF_measurement) > 0
                                                             else 
                                                             []
                                                             for NDF_measurement in contrast_gamma_folders 
                                                           ]
    
    # Parse the contrast gamma readings
    contrast_gamma_readings: list[list[list[dict]]] = [ _parse_temporal_sensitivity_readings(NDF_measurement,
                                                                                              apply_digital_gain=apply_digital_gain,
                                                                                              apply_radiometric_correction=apply_radiometric_correction,
                                                                                              use_mean_frame=use_mean_frame,
                                                                                              convert_time_units=convert_time_units,
                                                                                              convert_to_float=convert_to_float,
                                                                                              mean_axes=mean_axes,
                                                                                              contains_agc_metadata_dict={'W': True, 'P': False, 'M': False},
                                                                                              verbose=verbose
                                                                                             )
                                                        if len(NDF_measurement) > 0 
                                                        else 
                                                        []
                                                        for NDF_measurement in contrast_gamma_folders_sorted
                                                      ]



    ##################

    # ASSEMBLE RETURN VALUE

    ##################

    # Package the readings together into a dictionary
    parsed_readings_dict: dict[str, list[dict] | list[list[list[dict]]]] = {"ms_linearity": ms_linearity_readings, 
                                                                            "temporal_sensitivity": temporal_sensitivity_readings, 
                                                                            "phase_fitting": phase_fitting_readings,
                                                                            "contrast_gamma": contrast_gamma_readings, 
                                                                            "world_linearity": camera_linearity_readings
                                                                           }
    
    return parsed_readings_dict


if(__name__ == "__main__"):
    experiment_path: str = "/Volumes/T7 Shield/5cameraLinearity" 
    results = load_sorted_calibration_files(experiment_path, apply_digital_gain=True, use_mean_frame=True, convert_time_units=True, convert_to_float=True, verbose=True)
