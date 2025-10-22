import argparse
import pathlib
import platform
from tqdm import tqdm
import os
import shutil

# Define the list of activities for a given experiment 
experiment_name_dict: dict = {"scriptedIndoorOutdoor"}
activities_dict: dict = {""}
modes_dict: dict = {"sf": "spatial_frequency", 
                    "tf": "temporal_frequency"
                   }

"""Parse the commandline arguments, returning the source path and the experiment name"""
def parse_args() -> str:
    # Initialize the arg parser 
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description="Move recording files to the NAS")

    # Add arguments 
    parser.add_argument("experiment_name", type="str", choices=experiment_name_dict.keys(), help="Name of the experiment this recording is for")
    parser.add_argument("recording_path", type=str, help="Path to the recording file")

    # Parse the arguments 
    args: object = parser.parse_args() 

    return args.experiment_name, args.recording_path

"""Given a recording name for the scripted indoor/outdoor experiment
   parse it into the DropBox hierarchical file structure
   E.G. : FLIC_2002_cemetery_sf_session1
          /session_1/cemetary/FLIC_2002/sf

"""
def scripted_indoor_outdoor_recording_name_to_filestructure(experiment_name: str, recording_name: str) -> str:
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
    return os.path.join(session_num, activity, project_subject_ID, modes_dict[mode])

"""Given a path to a recording directory and an experiment name, 
   move the file to the NAS.
   Returns the path we sent the files to 
"""
def move_to_NAS(experiment_name: str, recording_path: str) -> str:
    # Determine what operating system we are on 
    operating_system: str = platform.system()

    # Parse the input path of a recording file 
    # and the experiment name 
    assert experiment_name in experiment_name_dict, f"Experiment name: {experiment_name} unrecognized. Do you have a typo? Valid names are: {experiment_name_dict.keys()}"
    assert os.path.exists(recording_path), f"Recording path {recording_path} does not exist. Do you have a typo?"
    assert os.path.isdir(recording_path), f"Recording path {recording_path} is not a recording folder."

    # Parse the source recording filename into the hierarchical structure 
    # we will use to put it on the NAS
    local_filestructure: str = scripted_indoor_outdoor_recording_name_to_filestructure(experiment_name, 
                                                                                       os.path.basename(recording_path.rstrip("/"))
                                                                                      )

    # If we are on MacOS 
    network_dive_path: str = ""
    if(operating_system == "Darwin"):
        # Construct the path to the network drive and ensure it exists 
        network_drive_path: str = "/Volumes/Aguirre-Brainard Lab NAS Shared/FLIC_raw"
    else:
        raise NotImplementedError()

    # Assert that we have a path to the raw video folder on the network drive 
    assert os.path.exists(network_drive_path), f"Network drivep path: {network_drive_path} does not exist. Are you connected on the same Network?"

    # Next, we will append the light logger project and the experiment name 
    network_dive_path = os.path.join(network_dive_path, "lightLogger", experiment_name)
    assert os.path.exists(network_drive_path), f"Network drivep path: {network_drive_path} does not exist. Are you connected on the same Network?"

    # Next, we will append the local filestructure 
    video_output_path: str = os.path.join(network_dive_path, local_filestructure)

    # Make the directory structure if it doesn't already exist 
    #os.makedirs(video_output_path, exist_ok=True)


    # Next, we will move the files to the NAS 
    filenames: list[str] = [filename for filename in os.listdir(recording_path)]
    for file_num in tqdm(range(len(filenames))):
        # Retrieve the source and dest filepaths 
        src: str = os.path.join(recording_path, filenames[file_num]) 
        dst: str = os.path.join(video_output_path, filenames[file_num])

        # Move the file to the target location 
        # shutil.move(src, dst)

    return video_output_path

def main():
    # Retrieve the experiment name and recording path from 
    # commandline 
    experiment_name, recording_path = parse_args()

    # Move the single recording to the NAS
    move_to_NAS(experiment_name, recording_path)


if(__name__ == "__main__"):
    main() 