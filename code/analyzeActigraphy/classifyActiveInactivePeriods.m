function classifyActiveInactivePeriods()
    arguments 

    end 


    %% Load Utility Libraries & Data
    persistent world_util ms_util;
    if isempty(world_util) || isempty(ms_util) || options.force_recalc
        world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path"));
        ms_util = import_pyfile(getpref("lightLoggerAnalysis", "ms_util_path")); 
    end
    
    % Load in the actigraphy data if needed 
    [IMUdata, eyeStateData, blinkData, gazeData] = load_actigraphy_data(IMUdata, eyeStateData, blinkData, gazeData);


end 
