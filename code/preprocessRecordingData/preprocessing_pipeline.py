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
def generate_world_videos(src_dir: str="/Volumes/FLIC_raw/scriptedIndoorVideos", 
                          dst_dir: str="/Volumes/FLIC_processing/scriptedIndoorVideos", 
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
def generate_egocentric_mapper_results(src_dir: str="/Volumes/FLIC_raw/scriptedIndoorVideos", 
                                       dst_dir: str="/Volumes/FLIC_processing/scriptedIndoorVideos",
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
def generate_virtually_foveated_videos(src_dir: str="/Volumes/FLIC_raw/scriptedIndoorVideos", 
                                       dst_dir: str="/Volumes/FLIC_processing/scriptedIndoorVideos",
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
            temp_output_dir: str = os.path.join(os.path.expanduser("~/Desktop"), "temp_output_dir")
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

                    shutil.move(temp_output_filepath, output_filepath)

                    # Delete the temporary dir 
                    shutil.rmtree(temp_output_dir)

    # Close the MATLAB engine 
    eng.close() 
    
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
def generate_spds(src_dir: str="/Volumes/FLIC_processing/scriptedIndoorVideos", 
                  dst_dir: str="/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorVideos",
                  overwrite_existing: bool=False,
                  activities_to_skip: Iterable=set(), 
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

            if(verbose is True):
                    print("Input: ")
                    print(f"\t Subject id: {subject_id}")
                    print(f"\t Subject id number: {subject_id_number}")
                    print(f"\t Activity: {activity_name}")

                    print("Output: ")
                    print(f"\t Output dir: {output_dir}")

            # Generate the SPDs and save at targeted output path 
            eng.processSPDs(src_dir, 
                            output_dir, 
                            "subjects", [subject_id_number], 
                            "activities", [activity_name],
                            "verbose", verbose,
                            "overwrite_existing", overwrite_existing,
                            "save_figures", True, 
                            nargout=0
                           )
    # Close the MATLAB engine
    eng.close() 

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
def unpack_neon_recordings(src_dir: str="/Volumes/FLIC_raw/scriptedIndoorVideos",
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


            # Unzip the file
            with zipfile.ZipFile(neon_recording_zip, 'r') as zip_ref:
                zip_ref.extractall(neon_recording_output_dir)

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
def rename_world_recordings(src_dir: str="/Volumes/FLIC_raw/scriptedIndoorVideos",
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
def verify_neon_integrity(src_dir: str="/Volumes/FLIC_raw/scriptedIndoorVideos",
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

"""
def verify_foveated_integrity(src_dir: str="/Volumes/FLIC_raw/scriptedIndoorVideos") -> dict[str, str]:
   
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


    return 
"""

def main():
    pass 

if(__name__ == "__main__"):
    main()  