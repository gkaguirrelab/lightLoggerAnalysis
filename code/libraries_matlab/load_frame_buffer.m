function [v, m] = load_frame_buffer(frame_buffer_path, metadata_path, contains_agc_metadata, password)
% Utility function to read in a numpy single frame buffer and optionally 
% its associated metadata into MATLAB as pure MATLAB type. 
%
% Syntax:
%  [v, m] = load_frame_buffer(frame_buffer_path, metadata_path, contains_agc_metadata, password)
%
% Description:
%   Given a path to a numpy array representing a frame buffer for 
%   a given sensor, as well as optionally its associated metadata 
%   matrix path, load them into MATLAB as pure MATLAB type. 
%   Accepts either a pure .npy file as input (or a compressed + encrypted blosc)
%   file. 
%
%
% Inputs:
%   frame_buffer_path    - String. The path to the numpy array
%                          representing a buffer of frames/readings                            
%
%   metadata_buffer_path  - String. The path to the numpy array
%                          representing the matrix of metadata 
%                          for the frames. Optional. 
%   
%   contains_agc_metadata - Bool. Whether or not the metadata 
%                           matrix contains data from an AGC. 
%
%   password             - String. Password used to decrypt 
%                          encrypted data. 
%   
% Outputs:
%
%   v                    - Array. 3D double array of 
%                          the shape (frames, rows, cols)
%
%   m                    - Array. Either single or 
%                          multi-dimensional array.
%                          The first col is timestamp.
%                          the rest are AGC metadata. 
%                          For order, consult specific sensor
%                          _util.py file. 
%                              
%
% Examples:
%{
    path_to_frame_buffer = '/Volumes/T7 Shield/quickOutdoorIndoor/world_chunk_17.npy'; 
    path_to_metadata_mat = '/Volumes/T7 Shield/quickOutdoorIndoor/world_chunk_17_metadata.npy'
    contains_agc_metadata = true; 
    [v, m] = load_frame_buffer(path_to_frame_buffer, path_to_metadata_mat, contains_agc_metadata); 
%}
    
    arguments
        frame_buffer_path {mustBeText}; % Path to a specific frame buffer 
        metadata_path {mustBeText} = ""; % Path to the metadata matrix for this frame buffer 
        contains_agc_metadata {mustBeNumericOrLogical} = false; % Whether or not the metadata buffer contains AGC data. Should usually be TRUE for world
        password {mustBeText} = "1234"; % Password used to decrypt encrpyted files
    end 

    % Let's load in the Python helper file
    addpath(getpref("lightLoggerAnalysis", "light_logger_libraries_matlab"));
    Pi_util = import_pyfile(getpref("lightLoggerAnalysis", Pi_util_path)); 

    % Load the frame buffer
    v_m = cell(Pi_util.load_frame_buffer(frame_buffer_path, metadata_path, contains_agc_metadata, password)); 
    
    % Splice out the frame buffer and the metadata
    [v, m] = v_m{:}; 

    return ; 
end 