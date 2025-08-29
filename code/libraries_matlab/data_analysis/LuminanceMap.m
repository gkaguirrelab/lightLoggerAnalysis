function [outputArg1,outputArg2] = LuminanceMap(inputArg1,inputArg2)
%LuminanceMap
%   Using the minispect values to estimate overall illuminance, this code
%   estiamtes the luminance of each pixel

% find paths to data file and load a chunk. This is for testing, pass in
% data later
drop_box_base_dir = getpref('combiExperiments','dropboxBaseDir');
data_path = '/FLIC_data/lightLogger'
load([drop_box_base_dir data_path '/sophia_in_wild_converted_for_sam/chunk_0.mat')

end