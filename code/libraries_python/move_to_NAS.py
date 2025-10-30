import argparse
import pathlib
import platform
from tqdm import tqdm
import os
import shutil
import re
import dill
import collections

# Define the list of activities for a given experiment 
experiment_name_set: set[str] = {"scriptedIndoorOutdoor"}
activities_dict: dict[str, str] = { "gazecalibration": "gazeCalibration",
                                    "lunch": "lunch", 
                                    "work": "work", 
                                    "chat": "chat", 
                                    "phone": "phone", 
                                    "walkindoor": "walkIndoor", 
                                    "walkoutdoor": "walkOutdoor", 
                                    "grocery": "grocery", 
                                    "cemetery": "cemetery", 
                                    "walkbiopond": "walkBiopond", 
                                    "sitbiopond": "sitBiopond"
                                  }
modes_dict: dict[str, str] = {"sf": "spatialFrequency", 
                              "tf": "temporalFrequency"
                             }

"""Parse the commandline arguments, returning the source path and the experiment name"""
def parse_args() -> str:
    # Initialize the arg parser 
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description="Move recording files to the NAS")

    # Add arguments 
    parser.add_argument("experiment_name", type="str", choices=experiment_name_set, help="Name of the experiment this recording is for")
    parser.add_argument("recording_path", type=str, help="Path to the recording file")

    # Parse the arguments 
    args: object = parser.parse_args() 

    return args.experiment_name, args.recording_path

"""Given a recording name for the scripted indoor/outdoor experiment
   parse it into the DropBox hierarchical file structure
   E.G. : FLIC_2002_cemetery_sf_session1
          /session_1/cemetary/FLIC_2002/sf

"""
def scripted_indoor_outdoor_recording_name_to_filestructure(recording_name: str) -> str:
    # Let's find the subject ID, activity, mode, and session number 

    # First, we will tokenize the string 
    tokens: list[str] = recording_name.split("_")
    
    # Subject ID is the first two tokens
    project_name: str = tokens[0].upper()
    subject_id: str = tokens[1]
    assert all(char.isnumeric() for char in subject_id), f"SubjectID: {subject_id} has non-numeric characters in it. Is this what you intended?"
    project_subject_ID: str = "_".join([project_name, subject_id])

    # Activity is the next token 
    activity: str = tokens[2].lower()
    assert activity in activities_dict, f"Activity: {activity} unrecognized. Potentially a typo? Valid activites are {activities_dict.keys()}"

    # Mode is the next token 
    mode: str = tokens[3]
    assert mode in modes_dict, f"Mode: {mode} unrecognized. Potentially a typo? Valid modes are: {modes_dict.keys()}"

    # Find the session number 
    session_num: str = "session_" + "".join([char for char in tokens[-1] if char.isnumeric()])

    # Construct a path from this 
    return os.path.join(project_subject_ID, activities_dict[activity], modes_dict[mode])

"""Given a path to a recording directory and an experiment name, 
   move the file to the NAS.
   Returns the path we sent the files to 
"""
def move_to_NAS(experiment_name: str, recording_path: str) -> str:
    # Determine what operating system we are on 
    operating_system: str = platform.system()

    # Parse the input path of a recording file 
    # and the experiment name 
    assert experiment_name in experiment_name_set, f"Experiment name: {experiment_name} unrecognized. Do you have a typo? Valid names are: {experiment_name_set}"
    assert os.path.exists(recording_path), f"Recording path {recording_path} does not exist. Do you have a typo?"
    assert os.path.isdir(recording_path), f"Recording path {recording_path} is not a recording folder."
    
    # Assert that there is a config file going along with the recording 
    assert "config.pkl" in os.listdir(recording_path), f"Recordign path {recording_path} missing config file"

    # Parse the source recording filename into the hierarchical structure 
    # we will use to put it on the NAS
    local_filestructure: str = scripted_indoor_outdoor_recording_name_to_filestructure(os.path.basename(recording_path.rstrip("/")))

    # If we are on MacOS 
    network_drive_path: str = ""
    if(operating_system == "Darwin"):
        # Construct the path to the network drive and ensure it exists 
        network_drive_path = "/Volumes/Aguirre-Brainard Lab NAS Shared/FLIC_raw"
    else:
        raise NotImplementedError()

    # Assert that we have a path to the raw video folder on the network drive 
    assert os.path.exists(network_drive_path), f"Network drive path: {network_drive_path} does not exist. Are you connected on the same Network?"

    # Next, we will append the light logger project and the experiment name 
    network_drive_path = os.path.join(network_drive_path, "lightLogger", experiment_name)
    assert os.path.exists(network_drive_path), f"Network drive path: {network_drive_path} does not exist. Are you connected on the same Network?"

    # Next, we will append the local filestructure 
    video_output_path: str = os.path.join(network_drive_path, local_filestructure)

    # Make the directory structure if it doesn't already exist 
    os.makedirs(video_output_path, exist_ok=True)

    # Let's get the chunk numbers from the files, so that we make sure there are not more than 1 video 
    # in this folder 
    chunk_numbers: list[int] = [ int(re.search(r"chunk_\d+", filename).group().split("_")[1])
                                 for filename in os.listdir(recording_path)
                                 if "config" not in filename 
                                 and "metadata" not in filename
                               ]
    # Let's get the frequency of the numbers 
    chunk_num_frequency_dict: collections.Counter = collections.Counter(chunk_numbers)

    # Assert some chunk started 
    assert 0 in chunk_num_frequency_dict, f"0 not in {recording_path}'s chunk nums. Something has gone very wrong"
    
    # Load in the config file. This will tell us the minimum number of chunk 0s we should have 
    config_dict: dict | None = None
    with open(os.path.join(recording_path, "config.pkl"), "rb") as f:
        config_dict = dill.load(f)

    # We check 2x here since metadata files also have chunk number in 
    assert chunk_num_frequency_dict[0] == len(config_dict["sensors"]), f"Multiple recordings in {recording_path}"
    

    # Next, we will move the files to the NAS 
    filenames: list[str] = [filename for filename in os.listdir(recording_path)]
    for file_num in tqdm(range(len(filenames))):
        # Retrieve the source and dest filepaths 
        src: str = os.path.join(recording_path, filenames[file_num]) 
        dst: str = os.path.join(video_output_path, filenames[file_num])

        # Move the file to the target location 
        shutil.copy(src, dst)

    return video_output_path

"""Verify all the data is present for a given participant"""
def verify_data(path_to_participant_folder: str) -> None:
    # Gather the recording folders for this participant 
    for folder in os.listdir(path_to_participant_folder):
        folder_path: str = os.path.join(path_to_participant_folder, folder)
        spatial_resolution_path: str = os.path.join(folder_path, "spatialFrequency")
        temporal_resolution_path: str = os.path.join(folder_path, "temporalFrequency")

        # Assert that the recordings exist for all folders, except gazecal which only uses temporal
        assert os.path.exists(temporal_resolution_path) and len(os.listdir(temporal_resolution_path)) > 0, f"Problem with {temporal_resolution_path}. Either does not exist or no content"
        if(folder == "gazeCalibration"):
            continue
        assert os.path.exists(spatial_resolution_path) and len(os.listdir(spatial_resolution_path)) > 0, f"Problem with {spatial_resolution_path}. Either does not exist or no content"

    return 

def main():
    # Retrieve the experiment name and recording path from 
    # commandline 
    experiment_name, recording_path = parse_args()

    # Move the single recording to the NAS
    move_to_NAS(experiment_name, recording_path)


if(__name__ == "__main__"):
    main() 