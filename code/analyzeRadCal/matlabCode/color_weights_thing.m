function observed_and_predicted = color_weights_thing(RGB_means)
% Color weights thing.
%
% Syntax:
%   observed_and_predicted = color_weights_thing(RGB_means)
%
% Description:
%   This function color weights thing.
% Inputs:
%   RGB_means                - Input used by the function.
%
% Outputs:
%   observed_and_predicted   - Output produced by the function.
%
% Examples:
%{
    observed_and_predicted = color_weights_thing(RGB_means)
%}

    arguments 
        RGB_means;  % Contrast Gamma Means: [102.6, 154.2, 181.9]
                    % Sphere Means: [123.1, 108.5, 44.9]
    end 

    % Then, we need to copy the way the MS did this code before
    pr670_measurement_results = load("/Users/zacharykelly/Downloads/CloudySkySPD_37degSolarElevation.mat");
    radiance = pr670_measurement_results.radiance; 

    sphereSPDs = radiance; 
    
    assert(all(wls == intersecting_wls)); 

    predictedCounts = sphereSPDs' * intersecting_wls_values; 
    observedCounts = RGB_means;   

    channel_labels = {"R", "G", "B"}; 

    % Define channel colors (RGB)
    channel_colors = [
        1 0 0;   % Red
        0 1 0;   % Green
        0 0 1    % Blue
    ];

    % Normalize to B channel
    predictedRatios = predictedCounts ./ predictedCounts(end);
    observedRatios = observedCounts ./ observedCounts(end);

    observed_and_predicted = [predictedRatios ./ observedRatios];

end 
