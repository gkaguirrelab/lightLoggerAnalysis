%% SCRIPT TO LOAD RAW FRAME BUFFERS AND CONVERT TO CHUNK-LIKE STRUCTURE
%
% Description:
%   This script loads raw frame buffer video data from the world (W) camera,
%   and repackages each buffer into a structure that mimics the format of
%   standard preprocessed 'chunks'. The resulting `chunks` cell array is
%   intended to be compatible with downstream analysis scripts that expect
%   each chunk to contain video data in the form:
%
%   chunks{i}.W.v --> 3D double array (frames × rows × cols)
%
%   PLEASE NOTE!!
%   This script does NOT populate the `.P` or `.M` sensor fields. It also
%   omits the `.t` (timestamp) subfield entirely. This structure serves as 
%   a partial/decoy version of the true `chunks` object and does not capture 
%   the full metadata or sensor diversity present in the raw buffer files.
%
% Usage:
%   path_to_video = 'XXXX'
%   MIN = X     (or define in script)
%   MAX = X     (or define in script)


% sensor_chunks is a Python dict
sensor_chunks = Pi_util.group_sensors_files(path_to_video); 

% Convert to struct with fields: W, P, M
sensor_chunks = struct(sensor_chunks);

% Convert world camera list to MATLAB cell array of 1×2 tuples
sensor_chunks.W = cell(sensor_chunks.W);

% Initialize output cell array
chunks = {};

% !! SELECT INDICIES HERE !!
for i = MIN : MAX
    % Extract metadata and frame buffer paths from tuple
    paths = string(sensor_chunks.W{i});
    metadata_path = paths(1);
    frame_buffer_path = paths(2);

    % Load video data (v) from raw frame buffer
    [v, ~] = load_frame_buffer(frame_buffer_path, metadata_path);

    % Construct quasi-chunk
    chunk = struct();
    chunk.W = struct('v', v);

    % Add to chunks cell array
    chunks{end+1} = chunk;
end
