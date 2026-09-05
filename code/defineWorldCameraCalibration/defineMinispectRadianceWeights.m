clear

dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'minispectRadianceCalibration',...
    'W2P2M2.mat');
load(dataFileName,'MSLinearityData');

calibration_metadata = MSLinearityData.metadata.ms_linearity;
measurements = MSLinearityData.readings.ms_linearity;

n_ndfs_to_plot = 4;

% Plotting defaults
set(groot, 'DefaultAxesFontSize', 16, 'DefaultTextFontSize', 16, ...
    'DefaultAxesColor', 'w', 'DefaultFigureColor', 'w');

% Path logic
asSensitivityFile = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'ASM7341_spectralSensitivity.mat');

% Hardcoded properties
chip = 'ASM7341';
limits = [-1, 5];
n_detector_channels = size(measurements{1,1}.M.v.AS, 2);

background = calibration_metadata.background;
background_scalars = calibration_metadata.background_scalars;
n_settings_levels = numel(background_scalars);

colorList = ndf_color_list(numel(calibration_metadata.NDFs));
countsAndCalculatedRadiance = cell(1, numel(calibration_metadata.NDFs));

% Process Spectral Sensitivities once 
referenceSourceS = calibration_metadata.cal_files{1}.rawData.S;
referenceWls = SToWls(referenceSourceS);

load(asSensitivityFile, 'T');
minispectS = WlsToS(T.wl);
minispectP_rel = T{:,2:end};
detectorP_rel = nan(referenceSourceS(3), size(minispectP_rel,2));

for jj = 1:size(minispectP_rel,2)
    detectorP_rel(:,jj) = interp1(SToWls(minispectS), minispectP_rel(:,jj), referenceWls);
end

vec = detectorP_rel';
peakSensitivityValue = nan(1, n_detector_channels);
idxVals = nan(1, n_detector_channels);
lambdaMax = nan(1, n_detector_channels);
for ii = 1:n_detector_channels
    [peakSensitivityValue(ii), idxVals(ii)] = max(vec(ii,:));
    lambdaMax(ii) = referenceWls(idxVals(ii));
    detectorP_rel(:,ii) = detectorP_rel(:,ii) ./ peakSensitivityValue(ii);
end

% Iterate over the NDF levels
for nn = 1:numel(calibration_metadata.NDFs)
    cal = calibration_metadata.cal_files{nn};
    sourceP_abs = cal.processedData.P_device;
    
    % Extract and average detector counts
    detectorCounts = extract_detector_counts(nn, measurements);
    detectorCounts = squeeze(mean(detectorCounts, [2, 3]));

    calculatedRadiance = nan(n_settings_levels, n_detector_channels);
    for ss = 1:n_settings_levels
        source_settings = background * background_scalars(ss);
        sphereSPDs = sourceP_abs * source_settings';
        calculatedRadiance(ss,:) = sphereSPDs' * detectorP_rel;
    end
    countsAndCalculatedRadiance{nn} = {detectorCounts, calculatedRadiance};
end 

% Plot linearity directly into a single 3x4 grid figure
logLogFitCoefficients = nan(n_detector_channels, 2);
fitCountLo = 10;
fitCountHigh = 2^16 - 11;

allChannelsFig = figure('Name', sprintf("MS_Linearity_%s_All_Channels", chip));
tl = tiledlayout(allChannelsFig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for ch = 1:n_detector_channels
    ax = nexttile(tl);
    hold(ax, 'on');

    measured_all_NDF_this_chip = nan(n_ndfs_to_plot, n_settings_levels);
    predicted_all_NDF_this_chip = nan(n_ndfs_to_plot, n_settings_levels);

    for nn = 1:n_ndfs_to_plot
        detectorCounts = countsAndCalculatedRadiance{nn}{1};
        calculatedRadiance = countsAndCalculatedRadiance{nn}{2};

        predicted_all_NDF_this_chip(nn,:) = log10(calculatedRadiance(:, ch)).';
        measured_all_NDF_this_chip(nn,:) = log10(detectorCounts(:, ch)).';

        h = scatter(ax, ...
            predicted_all_NDF_this_chip(nn,:), measured_all_NDF_this_chip(nn,:), ...
            'o', 'MarkerFaceColor', colorList(nn,:), ...
            'DisplayName', sprintf('NDF-%g', calibration_metadata.NDFs(nn)));
        h.MarkerFaceAlpha = 0.4;
        h.MarkerEdgeAlpha = 0.4;
    end

    fitX = predicted_all_NDF_this_chip(:);
    fitY = measured_all_NDF_this_chip(:);
    validFitIdx = isfinite(fitX) & isfinite(fitY) & ...
        fitY >= log10(fitCountLo) & fitY <= log10(fitCountHigh);

    fitIntercept = mean(fitY(validFitIdx) - fitX(validFitIdx));
    logLogFitCoefficients(ch,:) = [1, fitIntercept];
    fitXRange = [min(fitX(validFitIdx)), max(fitX(validFitIdx))];
    
    plot(ax, fitXRange, polyval(logLogFitCoefficients(ch,:), fitXRange), ...
        '-k', 'LineWidth', 1.5, 'DisplayName', 'Linear fit');
    plot(ax, [limits(1), limits(2)], [limits(1), limits(2)], ':k', "DisplayName", "IdentityLine");

    title(ax, sprintf("Ch %d | b = %+.3f", ch, logLogFitCoefficients(ch,2)));
    set(ax, 'box', 'off', 'color', 'none');
    
    if ch == 1
        lgd = legend(ax);
        lgd.Layout.Tile = 12;
    end
    hold(ax, 'off');
    xlim([-5 5]);
end

title(tl, sprintf("MS Linearity | %s | All Channels", chip));
xlabel(tl, sprintf('%s predicted radiance [log W/m^2/sr]', chip));
ylabel(tl, sprintf('%s measured counts [log]', chip));

fitObj.coeff = logLogFitCoefficients;
fitObj.lambdaMax = lambdaMax;

% Save the values that relate camera score to average scene radiance
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'minispectRadianceWeights.mat');
readme = ['Created by "defineMinispectRadianceWeights_2.m"\n'...
    'fitObj -- foo.\n'];
save(saveFileName, 'fitObj');

%% LOCAL FUNCTIONS

function colorList = ndf_color_list(num_NDFs)
baseColors = [
    0.6350, 0.0780, 0.1840;  
    0.8500, 0.3250, 0.0980;  
    0.9290, 0.6940, 0.1250;  
    0.4660, 0.6740, 0.1880;  
    0.3010, 0.7450, 0.9330;  
    0, 0.4470, 0.7410;       
    0.4940, 0.1840, 0.5560   
];

if(num_NDFs <= size(baseColors, 1))
    colorList = baseColors(1:num_NDFs, :);
    return;
end
colorList = interp1(linspace(1, num_NDFs, size(baseColors, 1)), baseColors, 1:num_NDFs, 'linear');
end

function counts_mat = extract_detector_counts(NDF_num, measurements)
[~, num_settings_levels, n_measures] = size(measurements);
min_num_readings = inf;
n_channels = size(measurements{NDF_num, 1, 1}.M.v.AS, 2);

for ss = 1:num_settings_levels
    for nn = 1:n_measures
        measurement = measurements{NDF_num, ss, nn};
        measurement_counts = measurement.M.v.AS;

        if isempty(measurement_counts)
            warning("NDF: %d | Settings Level: %d | Measurement: %d has no readings. Generating synthetic datapoint.", NDF_num-1, ss, nn);
            
            good_candidates_idx = setdiff(1:n_measures, nn);
            min_shared_candidate_readings = inf;
            
            for c_idx = good_candidates_idx
                c_counts = measurements{NDF_num, ss, c_idx}.M.v.AS;
                if isempty(c_counts)
                    error("Recording for settings level %d has multiple measurements without readings", ss);
                end
                min_shared_candidate_readings = min(min_shared_candidate_readings, size(c_counts, 1));
            end
            
            synthetic_example = zeros(min_shared_candidate_readings, n_channels);
            for c_idx = good_candidates_idx
                c_counts = measurements{NDF_num, ss, c_idx}.M.v.AS;
                synthetic_example = synthetic_example + c_counts(1:min_shared_candidate_readings, :);
            end
            
            synthetic_example = synthetic_example ./ numel(good_candidates_idx);
            measurement_counts = synthetic_example;
            
            measurement.M.v.AS = measurement_counts;
            measurements{NDF_num, ss, nn} = measurement;
        end
        min_num_readings = min(min_num_readings, size(measurement_counts, 1));
    end
end

counts_mat = nan(num_settings_levels, n_measures, min_num_readings, n_channels);
for ss = 1:num_settings_levels
    for nn = 1:n_measures
        readings = measurements{NDF_num, ss, nn}.M.v.AS;
        counts_mat(ss, nn, 1:min_num_readings, :) = readings(1:min_num_readings, :);
    end
end
end