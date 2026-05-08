import os 
import re
from typing import Iterable, Literal, Any
from natsort import natsorted
from scipy.io import loadmat, savemat


def load_spds(src_dir: str,
                subjects_to_skip: Iterable=set(),
                subjects_to_process: Iterable=set(), 
                activities_to_skip: Iterable=set(), 
                activities_to_process: Iterable=set(), 
                color_modes_to_skip: Iterable[Literal["a", "c_lm", "c_s"]]=set(), 
                color_modes_to_process: Iterable[Literal["a", "c_lm", "c_s"]]=set(), 
                projection_types_to_skip: Iterable=set(), 
                projection_types_to_process: Iterable=set(), 
                paths_only: bool=False, 
                output_as_mat: str = ""
                ) -> dict | None:
    
    # Convert everything to set if not set already for fast lookups
    subjects_to_skip = set(subjects_to_skip)
    subjects_to_process = set(subjects_to_process)
    activities_to_skip = set(activities_to_skip)
    activities_to_process = set(activities_to_process)
    color_modes_to_skip = set(color_modes_to_skip)
    color_modes_to_process = set(color_modes_to_process)
    projection_types_to_skip = set(projection_types_to_skip)
    projection_types_to_process = set(projection_types_to_process)

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
    assert len(color_modes_to_process) == 0 or ( len(color_modes_to_process) > 0 and set(os.path.basename(path) for path in color_mode_paths) == color_modes_to_process), f"Color modes to process requested were: {color_modes_to_process} but color modes found were: {color_mode_paths}"

    # Then, let's gather all of the subjects in the src_dir
    for color_mode_path in color_mode_paths:
        color_mode: str = os.path.basename(color_mode_path)

        # Initialize this field in the output dict
        output_dict[color_mode] = {}


        assert len(os.listdir(color_mode_path)) > 0, f"No files found in path: {color_mode_path}"
        subject_paths: list[str] = [ os.path.join(color_mode_path, filename)
                                    for filename in natsorted(os.listdir(color_mode_path))
                                    if os.path.isdir(os.path.join(color_mode_path, filename))
                                    and not filename.startswith(".")
                                    and filename.startswith("FLIC")
                                    and _is_desired(_extract_num_from_id(os.path.basename(filename)), subjects_to_process, subjects_to_skip)
                                ]
        assert len(subjects_to_process) == 0 or ( len(subjects_to_process) > 0 and set( _extract_num_from_id(os.path.basename(path)) for path in subject_paths) == subjects_to_process), f"Subjects to process requested were: {subjects_to_process} but subjects found were: {subject_paths}"

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
            assert len(activities_to_process) == 0 or ( len(activities_to_process) > 0 and set(os.path.basename(path) for path in activity_paths) == activities_to_process), f"Activities to process requested were: {activities_to_process} but color modes found were: {activity_paths}"

            # Now, let's gather the projection paths 
            for activity_path in activity_paths:
                activity_name: str = os.path.basename(activity_path)

                # Initialize this field in the subject of this colormode 
                output_dict[color_mode][subject_id][activity_name] = {}

                valid_projection_types: list[str] = natsorted([projection_type for projection_type in 
                                                             ("justProjection", "virtuallyFoveated") 
                                                              if _is_desired(projection_type, projection_types_to_process, projection_types_to_skip)

                                                              ]
                                                            )
                assert len(valid_projection_types) > 0

                # Now, let's load in the data from this projection type
                for projection_type in valid_projection_types:
                    projection_spd_path: str = os.path.join(activity_path, f"{subject_id}_{activity_name}_{projection_type}_SPDResults.mat")
                    assert os.path.exists(projection_spd_path), f"Path does not exist: {projection_spd_path}"

                    projection_best_fit_path: str = os.path.join(activity_path, f"{color_mode}_{projection_type}_bestFit.mat")
                    assert os.path.exists(projection_best_fit_path), f"Path does not exist: {projection_best_fit_path}"

                    # Load in this data (either as path or loading in the mat file)
                    output_dict[color_mode][subject_id][activity_name][projection_type] = {"spd": loadmat(projection_spd_path) if not paths_only else projection_spd_path,
                                                                                           "best_fit": loadmat(projection_best_fit_path) if not paths_only else projection_best_fit_path
                                                                                          }
                    
    # If we want to output as a mat file (does a bunch of auto conversion for us, do it now)
    if(output_as_mat != ""):
        savemat(output_as_mat, {"spds": output_dict})
        return 

    return output_dict

def _is_desired(item: Any, to_process: Iterable, to_skip: Iterable) -> bool:
    if(item in to_process):
        return True 
    
    return len(to_process) == 0 and item not in to_skip

def _extract_num_from_id(subject_id: str) -> int:
    assert re.fullmatch("FLIC_\d+", subject_id), f"{subject_id} does not fit the format FLIC_[NUM]"
    return int(re.search(r"\d+", subject_id).group())



if(__name__ == "__main__"):
    pass 
