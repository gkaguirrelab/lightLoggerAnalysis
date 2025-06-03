function generate_playable_videos(recording_path, output_dir, apply_digital_gain, fill_missing_frames, draw_pupil_ROI, debayer_images, num_chunks_to_use) 
% Generate a directory of playable videos from a recording directory 
% of the light logger 
%
% Syntax:
%   generate_playable_videos(recording_path, output_dir, apply_digital_gain, fill_missing_frames, debayer_images, num_chunks_to_use) 
%
% Description:
%   Generates a directory of playable videos, one for each camera sensor
%   on the lightlogger, from a recording directory made with the lightlogger.
%   Optionally applies digital gain, fills in missing frames, and debayers
%   images when generating videos. 
%
% Inputs:
%   recording_path        - String. The path to the directory
%                           containing the chunks of the recording.  
%
%   output_dir            - String. The directory where the new
%                           directory of videos will be output. 
%    
%   apply_digital_gain    - Logical. Whether or not to apply 
%                           digital gain to the frames 
%                           when they are constructed into the video
%
%   fill_missing_frames   - Logical. Whether or not to fill 
%                           in the missing frames in videos 
%                           with black dummy frames 
%
%   draw_pupil_ROI        - Logical. Whether or not to 
%                           highlight the region that is 
%                           used for the pupil AGC 
%
%
%   debayer_images        - Logical. Whether or not to debayer 
%                           the images when constructing them 
%                           into a video
%
%   num_chunks_to_use     - Numeric. The number of chunks to 
%                           use for generating the video from 
%                           each sensor. Default is all. 
%
%
% Outputs:
%
%   NONE                
%
% Examples:
%{
    path_to_experiment = '/Volumes/EXTERNAL1/fixedWalkingPAGC';
    output_dir = "./"; 
    apply_digital_gain = true; 
    fill_missing_frames = true; 
    draw_pupil_ROI = false; 
    debayer_images = false; 
    generate_playable_videos(path_to_experiment, output_dir, apply_digital_gain, fill_missing_frames, draw_pupil_ROI, debayer_images) 
%}

    arguments
        recording_path {mustBeText}; % The path to the folder of the recording 
        output_dir {mustBeText}; % The directory in which a new folder containing videos from each sensor will go 
        apply_digital_gain {mustBeNumericOrLogical} = false; % Whether or not to apply digital gain when generating the videos 
        fill_missing_frames {mustBeNumericOrLogical} = false; % Whether or not to fill in the missing frames in videos with black dummy frames 
        draw_pupil_ROI {mustBeNumericOrLogical} = false; % Whether or not to highlight the region that is used for the pupil AGC 
        debayer_images {mustBeNumericOrLogical} = false; % Whether or not to debayer the images when constructing them into a video 
        num_chunks_to_use {mustBeNumeric} = inf; % The number of chunks to use for generating the video from each sensor. 
    end

    % First, save the current directory the user is working in 
    cwd = pwd(); 

    % Next, find the path to the current file 
    [filedir, ~, ~] = fileparts(mfilename("fullpath"));

    % cd into the directory where this file lives 
    cd(filedir);

    % Import the Python utility function 
    Pi_util = py.importlib.import_module("Pi_util");

    % Call the Python helper function 
    Pi_util.generate_playable_videos(recording_path, output_dir, apply_digital_gain, fill_missing_frames, draw_pupil_ROI, debayer_images, num_chunks_to_use);

end