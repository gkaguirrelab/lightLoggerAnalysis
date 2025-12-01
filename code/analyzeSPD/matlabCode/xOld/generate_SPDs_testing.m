clear all; 
close all;
subjectIDs = {"FLIC_2001"}
virtually_foveated_video_paths = {"/Users/zacharykelly/Desktop/FLIC_2001_walkIndoor_virtuallyFoveatedVideo.avi"}

for ii = 1:numel(subjectIDs)
    % Retrieve the subjectID and video path 
    subjectID = subjectIDs{ii}; 
    video_path = virtually_foveated_video_paths{ii};

    [whole_video_mean_slopeMap, whole_video_mean_aucMap, whole_video_frequencies] = mapSPDs(video_path, "doPlot", true, "num_frames_to_process", [59121, 59121+1500]); 

    % Save to the desktop these results 
    mapSPDResults.mean_slopeMap = whole_video_mean_slopeMap; 
    mapSPDResults.mean_aucMap = whole_video_mean_aucMap; 

    save_path = sprintf("/Users/zacharykelly/Desktop/%s_virtuallyFoveated_SPD.mat", subjectID);  
    save(save_path, "mapSPDResults"); 

end 
