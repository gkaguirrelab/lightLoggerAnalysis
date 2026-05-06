import os 
import re
from typing import Iterable, Any
from natsort import natsorted
from scipy.io import loadmat, savemat


def gather_spds(src_dir: str,
                subjects_to_skip: Iterable=set(),
                subjects_to_process: Iterable=set(), 
                activities_to_skip: Iterable=set(), 
                activities_to_process: Iterable=set(), 
                color_modes_to_skip: Iterable=set(), 
                color_modes_to_process: Iterable=set(), 
                projection_types_to_skip: Iterable=set(), 
                projection_types_to_process: Iterable=set(), 
                include_best_fit: bool=False, 
                output_as_mat: str = ""
                ) -> dict | None:
    
    # Convert everything to set if not set already for fast lookups
    # TODO: 

    # Define the output dictionary
    output_dict: dict = {}

    # First, let's gather all the color modes we want to explore 
    color_mode_paths: list[str] = [os.path.join(src_dir, filename)
                                   for filename in natsorted(os.listdir(src_dir))
                                   if os.path.isdir(os.path.join(src_dir, filename))
                                   and not filename.startswith(".")
                                   and not filename.startswith("actigraphy")
                                   and _is_desired(filename, color_modes_to_process, color_modes_to_skip)
                                ]
    assert len(color_modes_to_process) == 0 or ( len(color_modes_to_process) > 0 and set(color_mode_paths) == color_modes_to_process), f"Color modes to process requested were: {color_modes_to_process} but color modes found were: {color_mode_paths}"

    # Then, let's gather all of the subjects in the src_dir
    for color_mode_path in color_mode_paths:
        color_mode: str = os.path.basename(color_mode_path)

        # Initialize this field in the output dict
        output_dict[color_mode] = {}

        subject_paths: list[str] = [ os.path.join(color_mode_path, filename)
                                    for filename in natsorted(os.listdir(color_mode_path))
                                    if os.path.isdir(os.path.join(color_mode_path, filename))
                                    and not filename.startswith(".")
                                    and _is_desired(filename, subjects_to_process, subjects_to_skip)
                                ]
        assert len(subjects_to_process) == 0 or ( len(subjects_to_process) > 0 and set(subject_paths) == subjects_to_process), f"Subjects to process requested were: {subjects_to_process} but color modes found were: {subject_paths}"

        # Now, let's gather all of the activities we want to load in
        for subject_path in subject_paths:
            subject_id: str = os.path.basename(subject_path)

            # Initialize this field in the colormode of this output dict 
            output_dict[color_mode][subject_id] = {}

            activity_paths: list[str] = [ os.path.join(subject_path, filename)
                                          for filename in natsorted(os.listdir(subject_path))
                                          if os.path.isdir(os.path.join(subject_path, filename))
                                          and not filename.startswith(".")
                                          and _is_desired(filename, activities_to_process, activities_to_skip)
                                        ]
            assert len(activities_to_process) == 0 or ( len(activities_to_process) > 0 and set(subject_paths) == activities_to_process), f"Activities to process requested were: {activities_to_process} but color modes found were: {activity_paths}"

            # Now, let's gather the projection paths 
            for activity_path in activity_paths:
                activity_name: str = os.path.basenane(activity_path)

                # Initialize this field in the subject of this colormode 
                output_dict[color_mode][subject_id][activity_name] = {}

                valid_projection_types: list[str] = natsorted([projection_type for projection_type in 
                                                             ("justProjection", "virtualFoveation") 
                                                              if _is_desired(projection_type, projection_types_to_process, projection_types_to_skip)

                                                              ]
                                                            )
                assert len(valid_projection_types) > 0

                # Now, let's load in the data from this projection type
                for projection_type in valid_projection_types:
                    projection_spd_path: str = ""
                    projection_best_fit_path: str = ""

                    # Load in this data 
                    output_dict[color_mode][subject_id][activity_name][projection_type] = {"spd": loadmat(projection_spd_path),
                                                                                           "best_fit": loadmat(projection_best_fit_path) if include_best_fit else projection_best_fit_path
                                                                                         }

    # If we want to output as a mat file (does a bunch of auto conversion for us, do it now)
    if(output_as_mat != ""):
        savemat(output_as_mat, {"spds": output_dict})
        return 

    return output_dict

def _is_desired(item: Any, to_process: Iterable, to_skip: Iterable) -> bool:
    return True



if(__name__ == "__main__"):
    pass 