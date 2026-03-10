import os
import shutil
from natsort import natsorted
import re
from tqdm.auto import tqdm
from typing import Iterable, Literal 
import pathlib
import sys


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

# TODO: Have to edit a bunch of MATLAB code to work nicely with making a dynamic src/dst, so right now it only works with these 
def generate_virtually_foveated_videos(src_dir: str="/Volumes/FLIC_raw/scriptedIndoorVideos", 
                                       dst_dir: str="/Volumes/FLIC_processing/scriptedIndoorVideos",
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

            # Define a temporary output location (to guard against permission issues)
            temp_output_dir: str = os.path.join("~/Desktop", subject_id, activity_name+"temp")
            output_dir: str = os.path.join(dst_dir, subject_id, activity_name)
            os.makedirs(temp_output_dir) # Not okay for this to exist 
            os.makedirs(output_dir, exist_ok=True) # Okay for this to exist 

            # Generate april tag and task for this subjecft/video
            for video_type in ("tag", "task"):    
                if(verbose is True):
                    print("Input: ")
                    print(f"\t Subject id: {subject_id}")
                    print(f"\t Subject id number: {subject_id_number}")
                    print(f"\t Activity: {activity_name}")
                    print(f"\t Video type: {video_type}")

                    print("Output: ")
                    print(f"\t Output dir: {output_dir}")

            
                eng.generateVirtuallyFoveatedVideos([subject_id_number], 
                                                    "output_dir", temp_output_dir, 
                                                    "activity", activity_name, 
                                                    "video_type", video_type, 
                                                    "overwrite_existing", overwrite_existing,
                                                    "verbose", verbose,
                                                    nargout=0
                                                )

                # Move the temp output to the target 
                temp_output_filenames: list[str] = os.listdir(temp_output_dir)
                assert len(temp_output_filenames) == 1, f"Found {len(temp_output_filenames)} temp output files. There should only be 1"
                temp_output_filepath: str = os.path.join(temp_output_dir, temp_output_filenames[0])

                shutil.move(temp_output_filepath, output_dir)

                # Delete the temporary dir 
                shutil.rmtree(temp_output_dir)

    # Close the MATLAB engine 
    eng.close() 
    
    return 

def generate_spds(src_dir: str="/Volumes/FLIC_processing/scriptedIndoorVideos", 
                  dst_dir: str="/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorVideos",
                  overwrite_existing: bool=False,
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


def main():
    pass 

if(__name__ == "__main__"):
    main()  