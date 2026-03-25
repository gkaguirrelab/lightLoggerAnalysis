import os
import shutil
from natsort import natsorted
import re
from tqdm.auto import tqdm
from typing import Iterable, Literal 
import pathlib
import sys
import zipfile
import warnings
import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import copy


# Construct the paths to our custom utility libraries 
light_loger_analysis_dir: str = str(pathlib.Path(__file__).parents[2]) 
video_util_path: str = os.path.join(light_loger_analysis_dir, "code", "library", "matlabIO", "python_libraries")
virtual_foveation_util_path: str = os.path.join(light_loger_analysis_dir, "code", "applyVirtualFoveation", "pythonCode")

custom_library_paths: list[str] = (light_loger_analysis_dir, video_util_path, virtual_foveation_util_path)
assert all(os.path.exists(path) for path in custom_library_paths)

for path in custom_library_paths:
    sys.path.append(path)

# Import the custom libraries 
import video_io 
import virtual_foveation

# -----------------------------------------------------------------------------
# generate_world_videos
#
# Build playable world-camera videos from raw chunk recordings.
# Iterates over all subjects and activities, combines chunk files from the
# GKA directory, optionally applies preprocessing (debayer, color weights,
# digital gain, frame filling), and writes a single W.avi per recording.
#
# Inputs:
#   src_dir            root directory containing raw FLIC recordings
#   dst_dir            destination directory for processed videos
#   overwrite_existing regenerate videos even if output exists
#   apply_color_weights apply camera color calibration
#   debayer_images     perform Bayer demosaicing
#   apply_digital_gain apply gain correction
#   fill_missing_frames interpolate missing frames
#   verbose            print progress information
# -----------------------------------------------------------------------------
def generate_world_videos(src_dir: str="/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026", 
                          dst_dir: str="/Volumes/FLIC_processing/NEWscriptedIndoorOutdoorVideos2026", 
                          overwrite_existing: bool=False,
                          apply_color_weights: bool=True, 
                          debayer_images: bool=True, 
                          apply_digital_gain: bool=True, 
                          fill_missing_frames: bool=True, 
                          verbose: bool=False
                         ) -> None:
    
    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 


    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Construct the path to the chunks that form the world video 
            world_video_in: str = os.path.join(activity_path, "GKA")
            assert os.path.exists(world_video_in) and os.path.isdir(world_video_in) and len(os.listdir(world_video_in)) > 0, f"Problem with: {world_video_in}"

            # Construct the output path to where this video will be output 
            world_video_out_dir: str = os.path.join(dst_dir, subject_id, activity_name, "GKA")
            world_video_out: str = os.path.join(world_video_out_dir, "W.avi")
            
            # Skip existing videos if so desired 
            if(os.path.exists(world_video_out) and overwrite_existing is False):
                continue

            # Otherwise, we will generate the video 
            # to the output location 
            os.makedirs(world_video_out_dir, exist_ok=True)

            if(verbose is True):
                print(f"Input: {world_video_in}")
                print(f"Output: {world_video_out}")

            video_io.world_chunks_to_video(world_video_in, 
                                           world_video_out,
                                           apply_color_weights=apply_color_weights, 
                                           debayer_images=debayer_images,
                                           apply_digital_gain=apply_digital_gain, 
                                           fill_missing_frames=fill_missing_frames,
                                           convert_to_seconds=True,
                                           verbose=verbose
                                          )
            
    return 

# -----------------------------------------------------------------------------
# generate_egocentric_mapper_results
#
# Runs the egocentric video mapper to align Neon gaze data with the world
# camera video. Produces gaze/fixation mappings and optional visualization
# outputs for each subject/activity recording.
#
# Inputs:
#   src_dir                    raw dataset directory
#   dst_dir                    processing output directory
#   mapping_choice             Fixations / Gaze / Both
#   refresh_time_threshold_sec mapping refresh threshold
#   render_video               render mapping visualization
#   render_video_comparison    render comparison video
#   optic_flow_algorithm       optical flow algorithm
#   image_matcher              feature matcher used for alignment
#   show_video_preview         show live preview
#   overwrite_existing         overwrite existing results
#   verbose                    print progress
# -----------------------------------------------------------------------------
def generate_egocentric_mapper_results(src_dir: str="/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026", 
                                       dst_dir: str="/Volumes/FLIC_processing/NEWscriptedIndoorOutdoorVideos2026",
                                       mapping_choice: Literal["Fixations", "Gaze", "Both"]='Both',
                                       refresh_time_threshold_sec: float=0.5 ,
                                       render_video: bool= False,
                                       render_video_comparison: bool= False,
                                       optic_flow_algorithm: Literal["Lucas-Kanade", "Gunnar Farneback"] = "Lucas-Kanade", 
                                       image_matcher: Literal["Efficient_LOFTR", "LOFTR_indoor"] = "Efficient_LOFTR",
                                       show_video_preview: bool=False, 
                                       overwrite_existing: bool=False,
                                       verbose: bool=False
                                      ) -> None:
    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 


    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Next, let's find the paths to both the Neon and the world videos 
            neon_dir: str = os.path.join(activity_path, "Neon")
            neon_video_path: str = os.path.join(neon_dir, [filename for filename in os.listdir(neon_dir) if "." not in filename][0])
            world_video_path: str = os.path.join(dst_dir, subject_id, activity_name, "GKA", "W.avi")

            # Assert these paths exist 
            assert os.path.exists(neon_video_path) and len(os.listdir(neon_video_path)) > 0, f"Problem with: {neon_video_path}"
            assert os.path.exists(world_video_path), f"Problem with: {world_video_path}"

            # Define the output directory 
            neon_output_dir: str = os.path.join(dst_dir, subject_id, activity_name, "Neon", "egocentric_mapper_results")

            # If the output already exsits adn we do not want to overwrite, then just skip 
            if(os.path.exists(neon_output_dir) and overwrite_existing is False):
                continue

            # Otherwise, run the egocentric video mapper 
            if(verbose is True):
                print(f"Input:")
                print(f"\tNeon: {neon_video_path}")
                print(f"\tWorld: {world_video_path}")
                print(f"Output: {neon_output_dir}")

            virtual_foveation.run_egocentric_video_mapper(neon_timeseries_dir=neon_video_path, 
                                                          alternative_vid_path=world_video_path, 
                                                          output_dir=neon_output_dir,
                                                          mapping_choice=mapping_choice, 
                                                          refresh_time_threshold_sec=refresh_time_threshold_sec, 
                                                          render_video=render_video,
                                                          render_video_comparison=render_video_comparison, 
                                                          optic_flow_algorithm=optic_flow_algorithm, 
                                                          image_matcher=image_matcher, 
                                                          show_video_preview=show_video_preview
                                                        )

    return 

# -----------------------------------------------------------------------------
# generate_virtually_foveated_videos
#
# Generates retinally-centered (virtually foveated) videos using MATLAB
# routines. For each subject/activity, launches MATLAB, runs the virtual
# foveation pipeline, and moves the resulting video into the processing
# directory.
#
# Inputs:
#   src_dir            raw dataset directory
#   dst_dir            processed dataset directory
#   overwrite_existing regenerate existing outputs
#   verbose            print progress
# -----------------------------------------------------------------------------
def generate_virtually_foveated_videos(src_dir: str="/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026", 
                                       dst_dir: str="/Volumes/FLIC_processing/NEWscriptedIndoorOutdoorVideos2026",
                                       overwrite_existing: bool=False,
                                       verbose: bool=False,
                                       video_types: Iterable[Literal["tag", "task"]]= ("tag", "task"),
                                       subjects_to_skip: Iterable=set(), 
                                       activities_to_skip: Iterable= set(),
                                       projection_types: Iterable[Literal["virtuallyFoveated", "justProjection"]] = set(["virtuallyFoveated", "justProjection"])
                                      ) -> None:
    
    import matlab.engine

    # Initialize the MATLAB engine to utilize the MATLAB function we have developed for this purpose 
    eng: object = matlab.engine.start_matlab()  
    eng.pyenv('Version', '~/Documents/MATLAB/projects/lightLoggerAnalysis/analysis_env/bin/python', nargout=0)
    eng.tbUseProject('lightLoggerAnalysis', nargout=0)
    

    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 


    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)
        subject_id_number: int = int(re.search("\d+", subject_id).group())
        
        # Skip desired subjects 
        if(subject_id_number in subjects_to_skip):
            continue 

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Skip any desired activities
            if(activity_name in activities_to_skip):
                continue 

            # Define a temporary output location (to guard against permission issues)
            temp_output_dir: str = os.path.join(os.path.expanduser("~/Desktop"), "temp_output_dir_foveation")
            output_dir: str = os.path.join(dst_dir, subject_id, activity_name)
            os.makedirs(output_dir, exist_ok=True) # Okay for this to exist 

            # Generate april tag and task for this subjecft/video
            for video_type in video_types:    
                # Determine if we should generate both with/without projection or not
                # We JUST virtually foveate the tag videos
                for projection_type in projection_types:
                    projection_only_flag: bool = projection_type == "justProjection"

                    # We JUST want to foveate the tags
                    if(video_type == "tag" and projection_type == "justProjection"):
                        continue 

                    # This needs to be hardcoded here because the matlab routine will receive
                    # just the temp output, so if structure changes we will also need to change this 
                    output_filepath: str = os.path.join(dst_dir, subject_id, activity_name, f"{subject_id}_{activity_name}_{video_type}_{projection_type}.avi")
                    if(os.path.exists(output_filepath) and overwrite_existing is False):
                        continue 

                    # Make the temp output dir if necessary 
                    if(os.path.exists(temp_output_dir)):
                        shutil.rmtree(temp_output_dir)
                    os.makedirs(temp_output_dir)
                                
                    if(verbose is True):
                        print("Input: ")
                        print(f"\t Subject id: {subject_id}")
                        print(f"\t Subject id number: {subject_id_number}")
                        print(f"\t Activity: {activity_name}")
                        print(f"\t Video type: {video_type}")
                        print(f"\t Projection only: {projection_only_flag}")

                        print("Output: ")
                        print(f"\t Output dir: {output_dir}")

                    # Generate for jsut projection only 
                    eng.generateVirtuallyFoveatedVideos([subject_id_number], 
                                                        "output_dir", temp_output_dir, 
                                                        "activity", activity_name, 
                                                        "video_type", video_type, 
                                                        "overwrite_existing", overwrite_existing,
                                                        "verbose", verbose,
                                                        "just_projection", projection_only_flag,
                                                        nargout=0
                                                    )

                    # Move the temp output to the target 
                    temp_output_filenames: list[str] = [ filename for filename in os.listdir(temp_output_dir) 
                                                        if filename.endswith(".avi")
                                                    ] 
                
                    assert len(temp_output_filenames) == 1, f"Found {len(temp_output_filenames)} @ {temp_output_filenames} temp output files. There should only be 1"
                    temp_output_filepath: str = os.path.join(temp_output_dir, temp_output_filenames[0])

                    if(os.path.exists(output_filepath)):
                        os.remove(output_filepath)
                    shutil.move(temp_output_filepath, output_filepath)

                    # Delete the temporary dir 
                    shutil.rmtree(temp_output_dir)

    # Close the MATLAB engine 
    eng.quit() 
    
    return 

# -----------------------------------------------------------------------------
# generate_spds
#
# Computes temporal spatial power spectra (SPD) statistics from processed
# videos using MATLAB analysis functions. Generates exponent maps, variance
# maps, and spectral summaries for each subject/activity.
#
# Inputs:
#   src_dir             directory containing processed videos
#   dst_dir             directory for SPD outputs and figures
#   overwrite_existing  recompute existing results
#   activities_to_skip  set of activity names to ignore
#   verbose             print progress
# -----------------------------------------------------------------------------
def generate_spds(src_dir: str="/Volumes/FLIC_processing/NEWscriptedIndoorOutdoorVideos2026", 
                  dst_dir: str="/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/NEWscriptedIndoorOutdoorVideos2026",
                  overwrite_existing: bool=False,
                  subjects_to_skip: Iterable=set(), 
                  activities_to_skip: Iterable=set(["lunch", "phone"]), 
                  projection_types: Iterable[Literal["virtuallyFoveated", "justProjection"]] = set(["virtuallyFoveated", "justProjection"]), 
                  common_axes: bool=False, 
                  verbose: bool=False) -> None:
    
    import matlab.engine

    # Initialize the MATLAB engine to utilize the MATLAB function we have developed for this purpose 
    eng: object = matlab.engine.start_matlab()  
    eng.pyenv('Version', '~/Documents/MATLAB/projects/lightLoggerAnalysis/analysis_env/bin/python', nargout=0)
    eng.tbUseProject('lightLoggerAnalysis', nargout=0)
    

     # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 

    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)
        subject_id_number: int = int(re.search("\d+", subject_id).group())

        # Skip subjects we dont want to process 
        if(subject_id_number in subjects_to_skip):
            continue

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            if(activity_name in activities_to_skip):
                continue

            # Generate the output path
            output_dir: str = os.path.join(dst_dir, subject_id, activity_name)
            os.makedirs(output_dir, exist_ok=True)

            # Generate the temp dir wher we will
            # temporarily copy the src files to. This is to avoid network interrupts 
            # reading large files over e.g. the NAS
            temp_output_dir: str = os.path.join(os.path.expanduser("~/Desktop"), "temp_output_dir_SPDs")
            if(os.path.exists(temp_output_dir)):
                shutil.rmtree(temp_output_dir)
            os.makedirs(temp_output_dir)

            # Copy the activity to the temporary directory 
            temp_activity_path: str = os.path.join(temp_output_dir, subject_id, activity_name)
            shutil.copytree(activity_path, temp_activity_path)

            # Generate SPDs for desired projection types (e.g. justProjection and virtuallyFoveated)
            for projection_type in projection_types:
                if(verbose is True):
                        print("Input: ")
                        print(f"\t Subject id: {subject_id}")
                        print(f"\t Subject id number: {subject_id_number}")
                        print(f"\t Activity: {activity_name}")
                        print(f"\t Projection Type: {projection_type}")

                        print("Output: ")
                        print(f"\t Output dir: {output_dir}")

                # Generate the SPDs and save at targeted output path 
                eng.processSPDs(temp_output_dir, 
                                output_dir, 
                                "subjects", [subject_id_number], 
                                "activities", [activity_name],
                                "verbose", verbose,
                                "overwrite_existing", overwrite_existing,
                                "save_figures", True, 
                                "projection_type", projection_type, 
                                nargout=0
                            )
            
            # Remove the temporary directory for this activity      
            shutil.rmtree(temp_output_dir)

    # If we want the output plots to share a common axis, 
    if(common_axes is True):
        adjust_spd_axes(src_dir=dst_dir,
                        dst_dir=dst_dir,
                        overwrite_existing=True, 
                        verbose=verbose, 
                        subjects_to_skip=subjects_to_skip,
                        projection_types=projection_types, 
                        activites_to_skip=activities_to_skip
                        )

    
    # Close the MATLAB engine
    eng.quit() 

    return 

def group_spds(src_dir: str="/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/NEWscriptedIndoorOutdoorVideos2026", 
               dst_dir: str="/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/NEWscriptedIndoorOutdoorVideos2026",
               groups: dict[str, dict[str, str | None]] = {"indoor": {"work": None, "chat": None, "walkIndoor": None}, 
                                                           "outdoor": {"walkOutdoor": None, "walkBiopond": None, "sitBiopond": None}
                                                          }, 
               overwrite_existing: bool=False,
               subjects_to_skip: Iterable=set(), 
               activities_to_skip: Iterable=set(["lunch", "phone"]), 
               projection_types: Iterable[Literal["virtuallyFoveated", "justProjection"]] = set(["virtuallyFoveated", "justProjection"]), 
               sort_by_experiment_ordering: bool=True, 
               verbose: bool=False
            ) -> None:
    
    import matlab.engine

    # Now we will make an inverse mapping of activities to groups rather than groups to activities 
    activites_to_groups: dict = {activity_name: group_name 
                                 for group_name, activity_dict in groups.items()
                                 for activity_name, _ in activity_dict.items()  
                                }



    # Initialize the MATLAB engine to utilize the MATLAB function we have developed for this purpose 
    eng: object = matlab.engine.start_matlab()  
    eng.pyenv('Version', '~/Documents/MATLAB/projects/lightLoggerAnalysis/analysis_env/bin/python', nargout=0)
    eng.tbUseProject('lightLoggerAnalysis', nargout=0)

    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 

    # First, we will find the min/max axes across all activites. 
    # This will let us build the power by freq graph 
    across_all_axes_min_maxes: dict[str, np.ndarray[float]] = _find_spd_axes_across_all(subject_paths, 
                                                                        subjects_to_skip,
                                                                        activities_to_skip,
                                                                        projection_types,
                                                                        verbose=False
                                                                        )

    # Next we will find the colorbar min/maxs by subject 
    per_subject_axes_min_maxes: dict[str, np.ndarray[float]] = _find_spd_axes_per_subject(subject_paths, 
                                                                        subjects_to_skip,
                                                                        activities_to_skip,
                                                                        projection_types,
                                                                        verbose=False
                                                                    )

    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)
        subject_id_number: int = int(re.search("\d+", subject_id).group())

        # Skip subjects we dont want to process 
        if(subject_id_number in subjects_to_skip):
            continue
        
         # Initialize per subject copy of the activity grouping 
        per_subject_activity_grouping: dict = {group: 
                                                    {activity_name: 
                                                        {projection_type: ""
                                                         for projection_type in projection_types
                                                        }
                                                    for activity_name, _ in activity_dict.items()
                                                    if activity_name not in activities_to_skip
                                                    }
                                               for group, activity_dict in groups.items()
                                              }
        

        # Construct the min maxes for this subject 
        axes_min_maxes: dict[str, np.ndarray] = per_subject_axes_min_maxes[subject_id_number]
        axes_min_maxes["spdByRegion"] = across_all_axes_min_maxes["spdByRegion"]
        axes_min_maxes["spdByRegion"][0] = max(axes_min_maxes["spdByRegion"][0], 10e-9)
        axes_min_maxes["frq"] = across_all_axes_min_maxes["frq"]

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            if(activity_name in activities_to_skip or activity_name not in activites_to_groups):
                continue

            # Now, we will gather the path to the SPD 
            for projection_type in projection_types:
                spd_results_mat_path: str = os.path.join(activity_path, f"{subject_id}_{activity_name}_{projection_type}_SPDResults.mat")
                assert os.path.exists(spd_results_mat_path), f"Problem with: {spd_results_mat_path}"    
                group_name: str = activites_to_groups[activity_name]
                per_subject_activity_grouping[group_name][activity_name][projection_type] = spd_results_mat_path
        
        # Assert there is no empty paths before sending to matlab 
        for group, activities_dict in per_subject_activity_grouping.items():
            for activity, projection_dict in activities_dict.items():
                for projection_type in projection_dict:
                    assert projection_dict[projection_type] != "", f"{group} | {activity} | {projection_type} is not assigned a path"

    

                
        # Now that we have the activities grouped together for this subject, pass it over to MATLAB 
        # to do the plotting 
        # First, output to a temp .mat file to make transfer easier between the two languages 
        temp_output_filepath: str = os.path.expanduser("~/Desktop/group_spds_temp.mat")
        scipy.io.savemat(temp_output_filepath, {"groupedActivityData": per_subject_activity_grouping})

        # Construct the output directory 
        output_dir: str = os.path.join(dst_dir, subject_id)
        eng.groupSPDs(temp_output_filepath, 
                      "exponent_clim", axes_min_maxes["exponentMap"], 
                      "variance_clim", axes_min_maxes["varianceMap"], 
                      "spd_xlim", axes_min_maxes["frq"], 
                      "spd_ylim", axes_min_maxes["spdByRegion"],
                      "title", str(subject_id_number),
                      "output_dir", output_dir, 
                      "overwrite_existing", overwrite_existing,
                      "sort_by_preferred_order", sort_by_experiment_ordering,
                      nargout=0
                    )


        # Remove the figure after it is finisehd 
        os.remove(temp_output_filepath)

    # Close the matlab engine 
    eng.quit()

    return 


def adjust_spd_axes(src_dir: str="/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/NEWscriptedIndoorOutdoorVideos2026",
                    dst_dir: str="/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/NEWscriptedIndoorOutdoorVideos2026",
                    subjects_to_skip: Iterable= set(),
                    activities_to_skip: Iterable= set(),
                    projection_types: Iterable[Literal["virtuallyFoveated", "justProjection"]] = set(["virtuallyFoveated", "justProjection"]), 
                    combine_figures: bool=True,
                    overwrite_existing: bool=False, 
                    verbose: bool=False
                   ) -> None:
    import matlab.engine

    # Initialize the MATLAB engine to utilize the MATLAB function we have developed for this purpose 
    eng: object = matlab.engine.start_matlab()  
    eng.pyenv('Version', '~/Documents/MATLAB/projects/lightLoggerAnalysis/analysis_env/bin/python', nargout=0)
    eng.tbUseProject('lightLoggerAnalysis', nargout=0)
    

     # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 

    # First, we will find the min/max axes across all activites. 
    # This will let us build the power by freq graph 
    across_all_axes_min_maxes: dict[str, np.ndarray[float]] = _find_spd_axes_across_all(subject_paths, 
                                                                        subjects_to_skip,
                                                                        activities_to_skip,
                                                                        projection_types,
                                                                        verbose=False
                                                                        )

    # Next we will find the colorbar min/maxs by subject 
    per_subject_axes_min_maxes: dict[str, np.ndarray[float]] = _find_spd_axes_per_subject(subject_paths, 
                                                                        subjects_to_skip,
                                                                        activities_to_skip,
                                                                        projection_types,
                                                                        verbose=False
                                                                    )

    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=False)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)
        subject_id_number: int = int(re.search("\d+", subject_id).group())

        # Skip unwatned subjects 
        if(subject_id_number in subjects_to_skip):
            continue 

        # Construct the min maxes for this subject 
        axes_min_maxes: dict[str, np.ndarray] = per_subject_axes_min_maxes[subject_id_number]
        axes_min_maxes["spdByRegion"] = across_all_axes_min_maxes["spdByRegion"]
        axes_min_maxes["spdByRegion"][0] = max(axes_min_maxes["spdByRegion"][0], 10e-9)
        axes_min_maxes["frq"] = across_all_axes_min_maxes["frq"]

        # We need to ensure x and y of the loglog spd plot are finite and > 0 
        if( not (all(axes_min_maxes["spdByRegion"] > 0) and all(axes_min_maxes["frq"] > 0))):
            raise Exception(f"Zero or negative X or Y axis in Log SPD plot: {axes_min_maxes}")

        if( not (all(np.isfinite(axes_min_maxes["spdByRegion"]))) and not all(np.isfinite(axes_min_maxes["frq"]))):
            raise Exception(f"Infinite X or Y axis in Log SPD plot: {axes_min_maxes}")

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                    if os.path.isdir(os.path.join(subject_path, filename))
                                    ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Skip unwanted activites 
            if(activity_name in activities_to_skip):
                continue 


            # We will put both projection types onto the same graph
            # so make a temporary .mat file that contains the path to do this 
            if(combine_figures is True):
                # Generate the path to the SPD results file
                temp_combined_dict: dict = {activity_name: {projection_type: os.path.join(activity_path, f"{subject_id}_{activity_name}_{projection_type}_SPDResults.mat")
                                                            for projection_type in projection_types}
                                            }
                
                # Assert that these paths exist
                assert all(os.path.exists(path) for projection_type, path in temp_combined_dict[activity_name].items())
                
            

                output_dir: str = os.path.join(dst_dir, subject_id, activity_name)
                title: str = f"{subject_id}_{activity_name}_modesCombined"
                eng.plotSPDs(temp_combined_dict[activity_name]["virtuallyFoveated"],
                                "exponent_clim", axes_min_maxes["exponentMap"], 
                                "variance_clim", axes_min_maxes["varianceMap"], 
                                "spd_xlim", axes_min_maxes["frq"], 
                                "spd_ylim", axes_min_maxes["spdByRegion"],
                                "overwrite_existing", overwrite_existing, 
                                "output_dir", output_dir,
                                "title", title,
                                "justProjectionActivityData", temp_combined_dict[activity_name]["justProjection"], 
                                nargout=0
                            ) 
                
            # Otherwise output separate figures
            else:
                # Iterate over the projection types 
                for projection_type in projection_types:
                    # Generate the path to the SPD results file
                    spd_results_filepath: str = os.path.join(activity_path, f"{subject_id}_{activity_name}_{projection_type}_SPDResults.mat")
                    assert os.path.exists(spd_results_filepath), f"SPD Results path does not exist: {spd_results_filepath}"
                    
                    output_dir: str = os.path.join(dst_dir, subject_id, activity_name)
                    title: str = f"{subject_id}_{activity_name}_{projection_type}"
                    eng.plotSPDs(spd_results_filepath,
                                 "exponent_clim", axes_min_maxes["exponentMap"], 
                                 "variance_clim", axes_min_maxes["varianceMap"], 
                                 "spd_xlim", axes_min_maxes["frq"], 
                                 "spd_ylim", axes_min_maxes["spdByRegion"],
                                 "overwrite_existing", overwrite_existing, 
                                 "output_dir", output_dir,
                                 "title", title,
                                 nargout=0
                                ) 


    # Close the matlab engine 
    eng.quit() 

    return 


def _find_spd_axes_across_all(subject_paths: list[str],
                                   subjects_to_skip: Iterable[str] = set(),
                                   activities_to_skip: Iterable[str] = set(),
                                   projection_types: Iterable[Literal["virtuallyFoveated", "justProjection"]] = set(["virtuallyFoveated", "justProjection"]), 
                                   verbose: bool=False
                                 ) -> None:
    # Initialize min max per type of graph 
    min_maxes: dict[str, list[float]] = {
        "exponentMap": [float("inf"), float("-inf")],
        "varianceMap": [float("inf"), float("-inf")],
        "spdByRegion": [float("inf"), float("-inf")],
        "frq": [float("inf"), float("-inf")],
    }
    

    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=False)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)
        subject_id_number: int = int(re.search("\d+", subject_id).group())

        # Skip unwatned subjects 
        if(subject_id_number in subjects_to_skip):
            continue 

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                    if os.path.isdir(os.path.join(subject_path, filename))
                                    ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Skip unwanted activites 
            if(activity_name in activities_to_skip):
                continue 

            # Iterate over the projection types 
            for projection_type in projection_types:
                # load the SPD results from this activity
                spd_results_filepath: str = os.path.join(activity_path, f"{subject_id}_{activity_name}_{projection_type}_SPDResults.mat")
                assert os.path.exists(spd_results_filepath), f"SPD Results path does not exist: {spd_results_filepath}"
                spd_results: object = scipy.io.loadmat(spd_results_filepath)['activityData'][0, 0][activity_name]

                # Extract the per graph info 
                for graph_type in min_maxes:
                    graph_info: np.ndarray = spd_results[graph_type][0, 0].astype(np.float64)
                    
                    # Find the min and max of this graph info 
                    graph_min: float = np.nanmin(graph_info)
                    graph_max: float = np.nanmax(graph_info)

                    # Compare to the global min/max and update if needed
                    global_min, global_max = min_maxes[graph_type]
                    if(graph_min < global_min):
                        min_maxes[graph_type][0] = graph_min
                    if(graph_max > global_max):
                        min_maxes[graph_type][1] = graph_max

    # convert min maxes to np.ndarray 
    for graph_type in min_maxes:
        min_maxes[graph_type] = np.array( min_maxes[graph_type], dtype=np.float64)

    return min_maxes


def _find_spd_axes_per_subject(subject_paths: list[str],
                                   subjects_to_skip: Iterable[str] = set(),
                                   activities_to_skip: Iterable[str] = set(),
                                   projection_types: Iterable[Literal["virtuallyFoveated", "justProjection"]] = set(["virtuallyFoveated", "justProjection"]),
                                   verbose: bool=False
                              ) -> None:
    
    # Initialize min max per type of graph 
    min_maxes: dict[int, list[float]] = {}


    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=False)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)
        subject_id_number: int = int(re.search("\d+", subject_id).group())

        # Skip unwatned subjects 
        if(subject_id_number in subjects_to_skip):
            continue 

        # Insert this subject into the min-maxes dict 
        min_maxes[subject_id_number] = {"exponentMap": [float("inf"), float("-inf")],
                                        "varianceMap": [float("inf"), float("-inf")],
                                        "spdByRegion": [float("inf"), float("-inf")],
                                        "frq": [float("inf"), float("-inf")]
                                       }

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                    if os.path.isdir(os.path.join(subject_path, filename))
                                    ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Skip unwanted activites 
            if(activity_name in activities_to_skip):
                continue 

            # Iterate over the projection types 
            for projection_type in projection_types:
                # load the SPD results from this activity
                spd_results_filepath: str = os.path.join(activity_path, f"{subject_id}_{activity_name}_{projection_type}_SPDResults.mat")
                assert os.path.exists(spd_results_filepath), f"SPD Results path does not exist: {spd_results_filepath}"
                spd_results: object = scipy.io.loadmat(spd_results_filepath)['activityData'][0, 0][activity_name]

                # Extract the per graph info 
                for graph_type in min_maxes[subject_id_number]:
                    graph_info: np.ndarray = spd_results[graph_type][0, 0].astype(np.float64)
                    
                    # Find the min and max of this graph info 
                    graph_min: float = np.nanmin(graph_info)
                    graph_max: float = np.nanmax(graph_info)

                    # Compare to the global min/max and update if needed
                    global_min, global_max = min_maxes[subject_id_number][graph_type]
                    if(graph_min < global_min):
                        min_maxes[subject_id_number][graph_type][0] = graph_min
                    if(graph_max > global_max):
                        min_maxes[subject_id_number][graph_type][1] = graph_max

    # convert min maxes to np.ndarray 
    for subject_id_number in min_maxes:
        for graph_type in min_maxes[subject_id_number]:
            min_maxes[subject_id_number][graph_type] = np.array( min_maxes[subject_id_number][graph_type], dtype=np.float64)

    return min_maxes


def _find_spd_axes_per_activity(subject_paths: list[str],
                                activities_list: list[str],
                                activities_to_skip: list[str], 
                                subjects_to_skip: list[str], 
                                projection_types: Iterable[Literal["virtuallyFoveated", "justProjection"]] = set(["virtuallyFoveated", "justProjection"]),
                                verbose: bool=False 
                              ) -> None:

    # Initialize min max per type of graph 
    min_maxes: dict[int, list[float]] = {}

    # First, let's go over all of the activities 
    activities_iterator: Iterable = range(len(activities_list)) if verbose is False else tqdm(range(len(activities_list)), desc="Processing Activities", leave=False)
    for activity_num in activities_iterator:
        # Retrieve the activity path and activity name
        activity_name: str = activities_list[activity_num]
        if(activity_name in activities_to_skip):
            continue
        
        # Insert this subject into the min-maxes dict 
        min_maxes[activity_name] = {"exponentMap": [float("inf"), float("-inf")],
                                        "varianceMap": [float("inf"), float("-inf")],
                                        "spdByRegion": [float("inf"), float("-inf")],
                                        "frq": [float("inf"), float("-inf")]
                                       }


        # Now, let's iterate over all the subject paths 
        subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=False)
        for subject_num in subject_iterator:
            # Retrieve the subject path and subject name
            subject_path: str = subject_paths[subject_num]
            subject_id: str = os.path.basename(subject_path)
            subject_id_number: int = int(re.search("\d+", subject_id).group())

            # Skip unwatned subjects 
            if(subject_id_number in subjects_to_skip):
                continue 
            
            # Construct the path to this subject's activity
            activity_path: str = os.path.join(subject_path, activity_name)
            assert os.path.exists(activity_path), f"Problem with: {activity_path}"

            # Iterate over the projection types 
            for projection_type in projection_types:
                # load the SPD results from this activity
                spd_results_filepath: str = os.path.join(activity_path, f"{subject_id}_{activity_name}_{projection_type}_SPDResults.mat")
                assert os.path.exists(spd_results_filepath), f"SPD Results path does not exist: {spd_results_filepath}"
                spd_results: object = scipy.io.loadmat(spd_results_filepath)['activityData'][0, 0][activity_name]

                # Extract the per graph info 
                for graph_type in min_maxes[activity_name]:
                    graph_info: np.ndarray = spd_results[graph_type][0, 0].astype(np.float64)
                    
                    # Find the min and max of this graph info 
                    graph_min: float = np.nanmin(graph_info)
                    graph_max: float = np.nanmax(graph_info)

                    # Compare to the global min/max and update if needed
                    global_min, global_max = min_maxes[activity_name][graph_type]
                    if(graph_min < global_min):
                        min_maxes[activity_name][graph_type][0] = graph_min
                    if(graph_max > global_max):
                        min_maxes[activity_name][graph_type][1] = graph_max

    return min_maxes


def generate_spds_across_subject(src_dir: str="/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/NEWscriptedIndoorOutdoorVideos2026/", 
                                 dst_dir: str="/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/NEWscriptedIndoorOutdoorVideos2026/acrossSubjects",
                                 overwrite_existing: bool=False,
                                 subjects_to_skip: Iterable=set(), 
                                 activities_to_skip: Iterable=set(["lunch", "phone"]), 
                                 projection_types: Iterable[Literal["virtuallyFoveated", "justProjection"]] = ["virtuallyFoveated", "justProjection"], 
                                 common_axes: bool=False, 
                                 combine_figures: bool=False, 
                                 verbose: bool=False
                                ) -> None:
    import matlab.engine

     # Initialize the MATLAB engine to utilize the MATLAB function we have developed for this purpose 
    eng: object = matlab.engine.start_matlab()  
    eng.pyenv('Version', '~/Documents/MATLAB/projects/lightLoggerAnalysis/analysis_env/bin/python', nargout=0)
    eng.tbUseProject('lightLoggerAnalysis', nargout=0)
    

    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    subject_id_numbers: list[int] = [ int(re.search(r"\d+", os.path.basename(path)).group()) 
                                        for path in subject_paths 
                                    ]
    subject_id_numbers = [num for num in subject_id_numbers if num not in subjects_to_skip]

    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 

    # First, let's get all of the activities  we will go over 
    activities_list: list[str] = sorted(  set([filename
                                        for subject_path in subject_paths
                                        for filename in os.listdir(subject_path)
                                        if os.path.isdir(os.path.join(subject_path, filename))
                                        and filename not in activities_to_skip
                                              ]
                                            )
                                        )


    # Define the axes limits (by default just automatically calculated, false)
    # or find the axes limits per activity for all graph types
    default_axes_min_maxes: dict[str, bool | list[float]] = {"exponentMap": False,
                                                         "varianceMap": False,
                                                         "spdByRegion": False,
                                                         "frq": False
                                                        }
    if(common_axes is True):
        all_axes_min_maxes = _find_spd_axes_per_activity(subject_paths=subject_paths,
                                                        activities_list=activities_list, 
                                                        activities_to_skip=activities_to_skip,
                                                        subjects_to_skip=subjects_to_skip,
                                                        projection_types=projection_types,
                                                        verbose=False
                                                        )
        # Loglog  plots cannot have x or y <= 0     
        # so fix that here 
        for activity_name in activities_list:
            all_axes_min_maxes[activity_name]["frq"][0] = max(all_axes_min_maxes[activity_name]["frq"][0], 10e-9)
            all_axes_min_maxes[activity_name]["spdByRegion"][0] = max(all_axes_min_maxes[activity_name]["spdByRegion"][0], 10e-9)

            # Make all of them numpy arrays for easy matlab conversion too 
            for field in all_axes_min_maxes[activity_name]:
                all_axes_min_maxes[activity_name][field] = np.array(all_axes_min_maxes[activity_name][field])

    # First, let's go over all of the activities 
    activities_iterator: Iterable = range(len(activities_list)) if verbose is False else tqdm(range(len(activities_list)), desc="Processing Activities", leave=False)
    for activity_num in activities_iterator:
        # Retrieve the activity path and activity name
        activity_name: str = activities_list[activity_num]

        # Gather the path to this activity for all subjects (and ensure it exists)
        activities_paths: list[str] = [ os.path.join(subject_path, activity_name) 
                                        for subject_path, subject_id_num in zip(subject_paths, subject_id_numbers)
                                        if subject_id_num not in subjects_to_skip
                                      ]
        for path in activities_paths:
            if(not os.path.exists(path)):
                raise Exception(f"Path: {path} does not exist")


        # Gather the list of numerical subject IDS 
        subject_id_numbers: list[int] = [ int(re.search(r"\d+", os.path.basename(path)).group()) 
                                          for path in subject_paths 
                                        ]
        subject_id_numbers = [num for num in subject_id_numbers if num not in subjects_to_skip]

        # If common axes is true, we will find the limits of the plots 
        # across all subjects for this activity
        
        # Initialize min max per type of graph 
        axes_min_maxes = default_axes_min_maxes if common_axes is False else all_axes_min_maxes[activity_name]
        
        # Generate the output dir
        output_dir: str = os.path.join(dst_dir, activity_name)
        os.makedirs(output_dir, exist_ok=True)

        # iterate over desired projection types 
        for projection_type in projection_types if combine_figures is False else (list(projection_types)[0],):
            # Call the MATLAB function to do the processing
            eng.processSPDsAcrossSubjects(src_dir, 
                                          output_dir, 
                                          "subjects", subject_id_numbers,
                                          "activities", [activity_name],
                                          "verbose", False, 
                                          "save_figures", True, 
                                          "overwrite_existing", overwrite_existing, 
                                          "projection_type", projection_type,
                                          "exponent_clim", axes_min_maxes["exponentMap"],
                                          "variance_clim", axes_min_maxes["varianceMap"], 
                                          "spd_xlim", axes_min_maxes["frq"],
                                          "spd_ylim", axes_min_maxes["spdByRegion"], 
                                          "combine_figures", combine_figures,
                                          nargout=0
                                        )

    # Close the matlab engine 
    eng.quit()

    return 

# -----------------------------------------------------------------------------
# unpack_neon_recordings
#
# Extracts Neon eye-tracking recordings from zipped archives inside each
# activity folder and places them into a "Neon" directory. The original
# zip archive is removed after extraction.
#
# Inputs:
#   src_dir            raw dataset directory
#   overwrite_existing re-extract if Neon folder already exists
#   verbose            print progress
# -----------------------------------------------------------------------------
def unpack_neon_recordings(src_dir: str="/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026",
                           overwrite_existing: bool=False,
                           verbose: bool=False
                          ) -> None:
    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 


    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)
        
            # Define the output path 
            neon_recording_output_dir: str = os.path.join(activity_path, "Neon")

            # Skip files we do not want to overwrite
            if(os.path.exists(neon_recording_output_dir) and len(os.listdir(neon_recording_output_dir)) > 0 and overwrite_existing is False):
               continue

            # Find the .zip file containing the neon recording 
            try:
                neon_recording_filename: str = [filename for filename in os.listdir(activity_path)
                                            if "Timeseries Data + Scene Video" in filename
                                           ][0]
            except: 
                raise Exception(f"No Neon timeseries .zip file in {activity_path}")    
            
            neon_recording_zip: str = os.path.join(activity_path, neon_recording_filename)

            if(verbose is True):
                print(f"Input: {neon_recording_zip}")
                print(f"Output: {neon_recording_output_dir}")

            # Unzip the file
            try:
                with zipfile.ZipFile(neon_recording_zip, 'r') as zip_ref:
                    zip_ref.extractall(neon_recording_output_dir)
            except Exception as e:
                shutil.rmtree(neon_recording_output_dir)
                raise Exception(f"Error when unzipping: {e}")


            # Remove the original file 
            os.remove(neon_recording_zip)

    return 


# -----------------------------------------------------------------------------
# rename_world_recordings
#
# Normalizes naming of raw world-camera recordings by renaming the original
# recording directory to "GKA". Ensures the dataset follows the expected
# structure used by downstream processing scripts.
#
# Inputs:
#   src_dir            raw dataset directory
#   overwrite_existing rename even if destination exists
#   verbose            print progress
# -----------------------------------------------------------------------------
def rename_world_recordings(src_dir: str="/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026",
                            overwrite_existing: bool=False,
                            verbose: bool=False
                           ) -> None:
    
    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 


    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Define the output location of this world video 
            world_recording_output_dir: str = os.path.join(activity_path, "GKA")

            # Skip files we do not want to overwrite
            if(os.path.exists(world_recording_output_dir) and len(os.listdir(world_recording_output_dir)) > 0 and overwrite_existing is False):
               continue

            # Try to find the world recording dir
            try:
                world_recording_filename: str = [filename for filename in os.listdir(activity_path)
                                                 if subject_id in filename and activity_name in filename
                                           ][0]
            except: 
                raise Exception(f"No world recording file in {activity_path}")    

            # Rename the file            
            os.rename(os.path.join(activity_path, world_recording_filename), world_recording_output_dir)

    return 

# -----------------------------------------------------------------------------
# verify_neon_integrity
#
# Performs sanity checks on extracted Neon eye-tracking recordings. Ensures
# required files exist (gaze, blink, fixation, timestamps, etc.) and that
# the recording directory contains the expected video and metadata files.
# Missing items generate warnings.
#
# Inputs:
#   src_dir   raw dataset directory
#   verbose   print additional diagnostics
# -----------------------------------------------------------------------------
def verify_neon_integrity(src_dir: str="/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026",
                          verbose: bool=False
                         ) -> None:
    
    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 


    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)
        
            # Construct the path to the Neon folder
            neon_folder_path: str = os.path.join(activity_path, "Neon")
            folder_exists: bool = os.path.exists(neon_folder_path)
            if(not folder_exists or len(os.listdir(neon_folder_path)) == 0):
                warnings.warn(f"{neon_folder_path} does not exist or is empty")
                continue 
                
            # If the folder exists, ensure it has the desired content 
            for filename in ("enrichment_info.txt", "sections.csv"):
                filepath: str = os.path.join(neon_folder_path, filename)
                if(not os.path.exists(filepath)):
                    warnings.warn(f"{filepath} does not exist")

            # Then, there should be a single directory in this file containing other 
            # information 
            neon_recording_folder: str | None = None
            try:
                neon_recording_folder = [os.path.join(neon_folder_path, filename)
                                            for filename in os.listdir(neon_folder_path)
                                            if os.path.isdir(os.path.join(neon_folder_path, filename))
                                            ][0]
            except:
                warnings.warn(f"Subfolder does not exist in {neon_folder_path}")

            # Next, make sure the required files exist in this folder 
            for filename in ("3d_eye_states.csv", "blinks.csv", "events.csv", "fixations.csv", "gaze.csv", "world_timestamps.csv", "saccades.csv", "template.csv"):
                filepath: str = os.path.join(neon_recording_folder, filename)
                if(not os.path.exists(filepath)):
                    warnings.warn(f"{filepath} does not exist")
            
            # Make sure there is a .mp4 video in this folder 
            try:
                mp4_video_name: str = [filename for filename in os.listdir(neon_recording_folder)
                                       if filename.endswith(".mp4")
                                      ][0]
            except:
                warnings.warn(f"{neon_recording_folder} does not have an .mp4 video")


    return 

def verify_world_neon_pairing(raw_dir: str="/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026",
                              processing_dir: str="/Volumes/FLIC_processing/NEWscriptedIndoorOutdoorVideos2026",
                              subjects_to_skip: Iterable= set(),
                              activities_to_skip: Iterable= set(),
                              verbose: bool=False 
                             ) -> dict[str, str]:
    
    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(raw_dir, subject_name) 
                                          for subject_name in os.listdir(raw_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(raw_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {raw_dir}" 

    # Initialize return dict 
    analyzed_data: dict[str, list] = {}

    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)
        subject_id_number: int = int(re.search("\d+", subject_id).group())

        # Skip subjects we are not interested in examining 
        if(subject_id_number in subjects_to_skip):
            continue

        analyzed_data[subject_id_number] = {}

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Skip desired activites 
            if(activity_name in activities_to_skip):
                continue 
                
            analyzed_data[subject_id_number][activity_name] = True

            # Now, we will find the Neon recording in the RAW dir 
            # and the world recording in the processing dir 
            neon_dir: str = os.path.join(activity_path, "Neon")
            assert os.path.exists(neon_dir) and len(os.listdir(neon_dir)) > 0, f"Problem with: {neon_dir}"
            neon_recording_subdir_name: str | None = None
            try:
                neon_recording_subdir_name = [filename for filename in os.listdir(neon_dir) 
                                              if os.path.isdir(os.path.join(neon_dir, filename))
                                             ][0]
            except:
                raise Exception(f"No neon recording subdir @: {neon_dir}")
            neon_recording_subdir: str = os.path.join(neon_dir, neon_recording_subdir_name)

            neon_recording_filename: str | None
            try:   
                neon_recording_filename = [filename for filename in os.listdir(neon_recording_subdir)
                                           if filename.endswith(".mp4")
                                          ][0]
            except:
                raise Exception(f"No neon recording @: {neon_recording_subdir}")
            neon_recording_path: str = os.path.join(neon_recording_subdir, neon_recording_filename)
            assert os.path.exists(neon_recording_path), f"Problem with: {neon_recording_path}"

            world_recording_dir: str = os.path.join(processing_dir, subject_id, activity_name, "GKA")
            assert os.path.exists(world_recording_dir) and len(os.listdir(world_recording_dir)) > 0, f"Problem with: {world_recording_dir}"
            world_recording_path: str = os.path.join(world_recording_dir, "W.avi")
            assert os.path.exists(world_recording_path), f"Problem with: {world_recording_path}"

            # Now, let's do a comparison between the first frames of the videos. 
            # We can just take the first world frame. That is always non null 
            world_frame_number: int = 0 
            world_frame: np.ndarray = video_io.destruct_video(world_recording_path, start_frame=world_frame_number, end_frame=world_frame_number+1)[0]

            # Neon has a period at the start full of gray frames. We will iterate through the frames and find the 
            # first frame that is not all gray 
            neon_num_frames: int = video_io.inspect_video_frame_count(neon_recording_path)
            neon_frame: np.ndarray | None = None
            for neon_frame_number in range(neon_num_frames):
                neon_frame = video_io.destruct_video(neon_recording_path, start_frame=neon_frame_number, end_frame=neon_frame_number+1)[0]

                # If not all pixels of the 3D neon frame are equal, then it's the first 
                # frame of video
                if(not np.all(neon_frame == neon_frame.flat[0])):
                    break

            fig, axes = plt.subplots(1, 2)
            fig.suptitle(f"{subject_id} | {activity_name}", y=0.8)
            for ax, frame, title, frame_number in zip(axes, (world_frame, neon_frame), "WN", (world_frame_number, neon_frame_number)):
                ax.set_title(f"{title} | Frame {frame_number}")
                ax.imshow(frame)
            plt.show()

    return analyzed_data

def generate_tag_task_start_ends(raw_dir: str="/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026",
                                 processing_dir: str="/Volumes/FLIC_processing/NEWscriptedIndoorOutdoorVideos2026",
                                 overwrite_existing: bool=False, 
                                 verbose: bool=False, 
                                 tag_end_offset: int= -120 * 2,
                                 task_start_offset: int=120 * 5, 
                                 task_length_seconds: int = 4 * 60, 
                                 min_confidence: float=30.0,
                                 frame_step: int=10
                                ) -> dict[str, dict[str, bool]]:

    raise NotImplementedError()

    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(raw_dir, subject_name) 
                                          for subject_name in os.listdir(raw_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(raw_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {raw_dir}" 

    # Initialize a dict to store the unanalyzable videos 
    videos_for_manual_review: dict[str, dict[str, bool]] = {}

    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)
        subject_id_number: int = int(re.search("\d+", subject_id).group())

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Skip already created files if we want to skip them 
            tag_task_output_path: str = os.path.join(processing_dir, subject_id, activity_name, "tag_task_start_end.mat")
            if(os.path.exists(tag_task_output_path) and overwrite_existing is False):
                continue 
            
            # Generate the path to the world video 
            world_video_path: str = os.path.join(processing_dir, subject_id, activity_name, "GKA", "W.avi")
            assert os.path.exists(world_video_path), f"Problem with: {world_video_path}"

            # Otherwise, let's make the tag_task struct 
            tag_task_start_end: dict[str, np.ndarray] = {}

            if(verbose is True):
                print(f"Input: {world_video_path}")
                print(f"Output: {tag_task_output_path}")

            # Let's now try to find the event marker in the video 
            event_marker_last_frame: int | None = None
            for keyword in ("place", "front", "camera"):
                event_marker_last_frame = video_io.find_events(world_video_path, 
                                                                           search_text=keyword,
                                                                           min_conf=min_confidence,
                                                                           frame_step=frame_step,
                                                                           verbose=verbose
                                                                      )
                if(event_marker_last_frame is not None):
                    break 

            # If we couldn't find the text, mark it for manual review 
            if(event_marker_last_frame is None):
                warnings.warn(f"Could not find start/end for video: {world_video_path}")

                if(subject_id_number not in videos_for_manual_review):
                    videos_for_manual_review[subject_id_number] = {activity_name: True}
                else:
                    videos_for_manual_review[subject_id_number][activity_name] = True
                
                continue 

            # Otherwise, we can make and save the struct 
            tag_start, tag_end = 1, event_marker_last_frame + tag_end_offset
            tag_task_start_end["tag"] = np.array([tag_start, tag_end])

            world_camera_fps: float = video_io.inspect_video_FPS(world_video_path)
            video_frame_count: int = video_io.inspect_video_frame_count(world_video_path)

            task_start: int = event_marker_last_frame + task_start_offset
            task_end: int = int(np.ceil(task_start + (world_camera_fps * task_length_seconds)))
            assert task_end < video_frame_count, f"Task end: {task_end} > video frame count {video_frame_count} for video @: {world_video_path}"
            
            tag_task_start_end["task"] = np.array([task_start, task_end])

            # Output the .mat file 
            scipy.io.savemat(tag_task_output_path, {"tag_task_start_end": tag_task_start_end})


    return videos_for_manual_review



def verify_virtually_foveated_video_integrity(src_dir: str="/Volumes/FLIC_processing/NEWscriptedIndoorOutdoorVideos2026",
                                              subjects_to_skip: Iterable=set(), 
                                              activities_to_skip: Iterable=set(),
                                              verbose: bool=False, 
                                              target_length_seconds: int= 4 * 60
                                             ) -> None:
    
    # First, let's find all of the subjects in this experiment 
    subject_paths: list[str] = natsorted([os.path.join(src_dir, subject_name) 
                                          for subject_name in os.listdir(src_dir) 
                                          if re.fullmatch(r"FLIC_\d+", subject_name) 
                                          and os.path.isdir(os.path.join(src_dir, subject_name))
                                         ]
                                        ) 
    assert len(subject_paths) > 0, f"No subject directories found in: {src_dir}" 


    # Now, let's iterate over all the subject paths 
    subject_iterator: Iterable = range(len(subject_paths)) if verbose is False else tqdm(range(len(subject_paths)), desc="Processing Subjects", leave=True)
    for subject_num in subject_iterator:
        # Retrieve the subject path and subject name
        subject_path: str = subject_paths[subject_num]
        subject_id: str = os.path.basename(subject_path)
        subject_id_number: int = int(re.search("\d+", subject_id).group())

        # Skip subjects we want ot skip 
        if(subject_id_number in subjects_to_skip):
            continue

        # Iterate over the activites for this subject 
        activites_paths: list[str] = [os.path.join(subject_path, filename) for filename in natsorted(os.listdir(subject_path))
                                      if os.path.isdir(os.path.join(subject_path, filename))
                                     ]
        activities_iterator: Iterable = range(len(activites_paths)) if verbose is False else tqdm(range(len(activites_paths)), desc="Processing Activities", leave=False)
        for activity_num in activities_iterator:
            # Retrieve the activity path and activity name
            activity_path: str = activites_paths[activity_num]
            activity_name: str = os.path.basename(activity_path)

            # Skip acticvities we want to skip 
            if(activity_name in activities_to_skip):
                continue

            # Find all the virtually foveated videos in this directory 
            virtually_foveated_video_names: list = [filename for filename in os.listdir(activity_path)
                                                    if filename.endswith(".avi")
                                                   ]

            # Assert we have the correct number of videos. Should be 3 for all activities (tag_foveated, task_projection, task_foveated)
            assert len(virtually_foveated_video_names) >= (3 if activity_name != "gazeCalibration" else 2), f"Problem with: {activity_path}"

            # Assert the projection and the foveated video have the proper length 
            for video_name in [name for name in virtually_foveated_video_names if "task" in name]:
                video_path: str = os.path.join(activity_path, video_name)

                # Find the FPS of the video and its length 
                video_fps: float = video_io.inspect_video_FPS(video_path)
                video_num_frames: int = video_io.inspect_video_frame_count(video_path)
                video_length: float = video_fps * video_num_frames

                assert video_length >= target_length_seconds, f"Video: {video_path} has length: {video_length}s which is less than target: {target_length_seconds}s"


def main():
    pass 

if(__name__ == "__main__"):
    main()  