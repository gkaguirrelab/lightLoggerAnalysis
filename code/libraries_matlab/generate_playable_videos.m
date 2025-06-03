function generate_playable_videos(recording_path, output_dir,...
                                  apply_digital_gain, fill_missing_frames, draw_pupil_ROI, debayer_images,...
                                  password,...
                                  time_ranges,... 
                                  chunk_ranges..., 
                                  Pi_util_path
                                 ) 
% Generate a directory of playable videos from a recording directory 
% of the light logger 
%
% Syntax:
%   generate_playable_videos(recording_path, output_dir,...
%                            apply_digital_gain, fill_missing_frames, draw_pupil_ROI, debayer_images,...
%                            password,...
%                            time_ranges,... 
%                            chunk_ranges..., 
%                            Pi_util_path
%                           ) 
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
%   password               - String. Represents the password 
%                           used to encrypt the data (if encrypted)
%
%   time_ranges           - Struct. Represents the selected 
%                           range of time to parse. Form is 
%                           a struct with fields (WPM)
%                           where each field is a tuple 
%                           of a range [start, end). 
%
%   chunk_ranges          - Struct. Represents the selected 
%                           range of chunks to parse. Form is 
%                           a struct with fields (WPM)
%                           where each field is a tuple 
%                           of a range [start, end). Note:
%                           the use of this and time_ranges is
%                           exclusively OR. 0-indexed. 
%
%   matlab_analysis_libraries_path  - String. Path to the utilized MATLAB 
%                                     helper functions (chunk_dict_to_matlab)
% 
%   Pi_util_path           - String. Path to the Pi_util.py helper file
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
        time_ranges = false; % The timestamps to splice out of a video in the form [start, end) (relative to start of the video)
        chunk_ranges = false; % The chunk numbers to splice out of a video in the form [start, end] (0-indexed)
        password {mustBeText} = "1234"; % The password needed to decrypt encrypted + compressed files (.blosc files)
        matlab_analysis_libraries_path {mustBeText} = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))), "libraries_matlab"); % Path to the utilized MATLAB helper functions (chunk_dict_to_matlab)
        Pi_util_path {mustBeText} = fullfile(fileparts(mfilename("fullpath")), 'Pi_util');  % Path to the Pi_util.py helper file
    end

    % Append path to the Python file importer
    matlab_libraries_path = matlab_analysis_libraries_path;
    addpath(matlab_libraries_path); 

    % Import the Python utility function 
    Pi_util = import_pyfile(Pi_util_path);

    % Call the Python helper function 
    Pi_util.generate_playable_videos(recording_path, output_dir, apply_digital_gain,...
                                     fill_missing_frames, draw_pupil_ROI, debayer_images,..
                                     password,...
                                     time_ranges, chunk_ranges...
                                    );

end