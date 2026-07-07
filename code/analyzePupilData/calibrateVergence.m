%% calibrateVergence_centerOnly.m
% Compare candidate vergence measures during center fixation at 1 meter.

clear; close all; clc

%% Recording folders
recordings = {
    'Walk Indoor', '/Users/sophiamirabal/Downloads/PupilData_Analysis/walkIndoor_timeseriesData/flic_1034_walkindoor_1-86ecaaf3';
    'Read',        '/Users/sophiamirabal/Downloads/PupilData_Analysis/read_timeseriesData/flic_1034_read_1-a5ad020a';
    'Sit Biopond', '/Users/sophiamirabal/Downloads/PupilData_Analysis/sitBiopond_timeseriesData/flic_1034_sitbiopond_1-57b4e2f7';
    'Chat',        '/Users/sophiamirabal/Downloads/PupilData_Analysis/chat_timeseriesData/flic_1034_chat_1-294a661b';
};

%% Center fixation timestamps, seconds from recording start
centerTimes = {
    'Walk Indoor', 22.691;
    'Read',        14.080;
    'Sit Biopond', 14.712;
    'Chat',        16.540;
};

windowHalfWidth = 0.15; % 300 ms total
smoothWindow = 15;

summary = table();

%% Expected vergence at 1 m, using estimated/measured IED
IED_mm = 63; % replace with participant's measured IED if available
expectedVergence_1m_deg = 2 * atan2d(IED_mm/2, 1000);

%% Plot 1: signed horizontal difference
figure; hold on

for r = 1:size(recordings,1)

    recName = recordings{r,1};
    recDir  = recordings{r,2};

    eye = readtable(fullfile(recDir, '3d_eye_states.csv'), ...
        'VariableNamingRule', 'preserve');

    t = (eye.("timestamp [ns]") - eye.("timestamp [ns]")(1)) * 1e-9;

    centerTime = centerTimes{r,2};
    idxCenter = t >= centerTime - windowHalfWidth & ...
                t <= centerTime + windowHalfWidth;

    tCenter = t(idxCenter);
    tCenterZeroed = tCenter - centerTime;

    %% Optical-axis vectors
    vLeft = [ ...
        eye.("optical axis left x"), ...
        eye.("optical axis left y"), ...
        eye.("optical axis left z")];

    vRight = [ ...
        eye.("optical axis right x"), ...
        eye.("optical axis right y"), ...
        eye.("optical axis right z")];

    %% Horizontal and vertical angles from optical-axis vectors
    hLeft = atan2d(vLeft(:,1), vLeft(:,3));
    hRight = atan2d(vRight(:,1), vRight(:,3));

    vAngLeft = atan2d(vLeft(:,2), vLeft(:,3));
    vAngRight = atan2d(vRight(:,2), vRight(:,3));

    signedHorizDiff = hLeft - hRight;
    absHorizDiff = abs(signedHorizDiff);

    verticalDiff = vAngLeft - vAngRight;
    absVerticalDiff = abs(verticalDiff);

    %% True 3D angular separation between optical axes
    dotLR = sum(vLeft .* vRight, 2);
    normL = vecnorm(vLeft, 2, 2);
    normR = vecnorm(vRight, 2, 2);

    cosTheta = dotLR ./ (normL .* normR);

    % Numerical safety
    cosTheta = max(min(cosTheta, 1), -1);

    vergence3D = acosd(cosTheta);

    %% Center-window data
    signedHorizCenter = signedHorizDiff(idxCenter);
    absHorizCenter = absHorizDiff(idxCenter);
    verticalCenter = verticalDiff(idxCenter);
    absVerticalCenter = absVerticalDiff(idxCenter);
    vergence3DCenter = vergence3D(idxCenter);

    hLeftCenter = hLeft(idxCenter);
    hRightCenter = hRight(idxCenter);

    hLeft_s = smoothdata(hLeftCenter, 'movmedian', smoothWindow, 'omitnan');
    hRight_s = smoothdata(hRightCenter, 'movmedian', smoothWindow, 'omitnan');
    signedHoriz_s = smoothdata(signedHorizCenter, 'movmedian', smoothWindow, 'omitnan');
    vergence3D_s = smoothdata(vergence3DCenter, 'movmedian', smoothWindow, 'omitnan');

    %% Plot signed horizontal difference
    plot(tCenterZeroed, signedHoriz_s, 'LineWidth', 1.5);

    %% Print diagnostic info
    fprintf('\n%s\n', recName);
    fprintf('Median left horizontal angle: %.3f deg\n', median(hLeftCenter, 'omitnan'));
    fprintf('Median right horizontal angle: %.3f deg\n', median(hRightCenter, 'omitnan'));
    fprintf('Median signed horizontal diff: %.3f deg\n', median(signedHorizCenter, 'omitnan'));
    fprintf('Median absolute horizontal diff: %.3f deg\n', median(absHorizCenter, 'omitnan'));
    fprintf('Median 3D optical-axis separation: %.3f deg\n', median(vergence3DCenter, 'omitnan'));

    %% Summary table
    newRow = table( ...
        string(recName), ...
        centerTime, ...
        centerTime - windowHalfWidth, ...
        centerTime + windowHalfWidth, ...
        median(hLeftCenter, 'omitnan'), ...
        median(hRightCenter, 'omitnan'), ...
        median(signedHorizCenter, 'omitnan'), ...
        median(absHorizCenter, 'omitnan'), ...
        median(verticalCenter, 'omitnan'), ...
        median(absVerticalCenter, 'omitnan'), ...
        median(vergence3DCenter, 'omitnan'), ...
        std(signedHorizCenter, 'omitnan'), ...
        std(vergence3DCenter, 'omitnan'), ...
        sum(~isnan(signedHorizCenter)), ...
        'VariableNames', {'Recording', ...
                          'CenterTime_s', ...
                          'StartTime_s', ...
                          'EndTime_s', ...
                          'MedianLeftHorizAngleDeg', ...
                          'MedianRightHorizAngleDeg', ...
                          'MedianSignedHorizDiffDeg', ...
                          'MedianAbsHorizDiffDeg', ...
                          'MedianSignedVerticalDiffDeg', ...
                          'MedianAbsVerticalDiffDeg', ...
                          'Median3DAngleBetweenEyesDeg', ...
                          'SDSignedHorizDiffDeg', ...
                          'SD3DAngleBetweenEyesDeg', ...
                          'NSamples'} ...
    );

    summary = [summary; newRow];

end

yline(0, '--');
yline(expectedVergence_1m_deg, 'k--', '+ expected 1 m');
yline(-expectedVergence_1m_deg, 'k--', '- expected 1 m');

xlabel('Time relative to center fixation (s)');
ylabel('Signed horizontal angle difference, Left - Right (deg)');
title('Signed Horizontal Angle Difference During Center Fixation at 1 m');
legend(recordings(:,1), 'Location', 'best');
grid on

%% Plot 2: 3D optical-axis separation
figure; hold on

for r = 1:size(recordings,1)

    recDir = recordings{r,2};
    eye = readtable(fullfile(recDir, '3d_eye_states.csv'), ...
        'VariableNamingRule', 'preserve');

    t = (eye.("timestamp [ns]") - eye.("timestamp [ns]")(1)) * 1e-9;
    centerTime = centerTimes{r,2};
    idxCenter = t >= centerTime - windowHalfWidth & ...
                t <= centerTime + windowHalfWidth;

    tCenter = t(idxCenter);
    tCenterZeroed = tCenter - centerTime;

    vLeft = [eye.("optical axis left x"), ...
             eye.("optical axis left y"), ...
             eye.("optical axis left z")];

    vRight = [eye.("optical axis right x"), ...
              eye.("optical axis right y"), ...
              eye.("optical axis right z")];

    dotLR = sum(vLeft .* vRight, 2);
    normL = vecnorm(vLeft, 2, 2);
    normR = vecnorm(vRight, 2, 2);
    cosTheta = dotLR ./ (normL .* normR);
    cosTheta = max(min(cosTheta, 1), -1);

    vergence3D = acosd(cosTheta);
    vergence3DCenter = vergence3D(idxCenter);
    vergence3D_s = smoothdata(vergence3DCenter, ...
        'movmedian', smoothWindow, 'omitnan');

    plot(tCenterZeroed, vergence3D_s, 'LineWidth', 1.5);

end

yline(expectedVergence_1m_deg, 'k--', 'Expected 1 m');

xlabel('Time relative to center fixation (s)');
ylabel('3D angle between optical axes (deg)');
title('3D Optical-Axis Separation During Center Fixation at 1 m');
legend(recordings(:,1), 'Location', 'best');
grid on

%% Print and save summary
disp(summary)

outFile = 'vergenceCalibration_centerOnly_fullDiagnostics_summary.csv';
writetable(summary, outFile);

fprintf('\nExpected 1 m vergence for IED %.1f mm: %.3f deg\n', ...
    IED_mm, expectedVergence_1m_deg);

fprintf('Saved summary to: %s\n', outFile);

%% Plot 3: separate left and right horizontal eye angles

figure; hold on

taskColors = lines(size(recordings,1));

for r = 1:size(recordings,1)

    recDir = recordings{r,2};

    eye = readtable(fullfile(recDir, '3d_eye_states.csv'), ...
        'VariableNamingRule', 'preserve');

    t = (eye.("timestamp [ns]") - eye.("timestamp [ns]")(1)) * 1e-9;

    centerTime = centerTimes{r,2};
    idxCenter = t >= centerTime - windowHalfWidth & ...
                t <= centerTime + windowHalfWidth;

    tCenter = t(idxCenter);
    tCenterZeroed = tCenter - centerTime;

    hLeft = atan2d( ...
        eye.("optical axis left x"), ...
        eye.("optical axis left z"));

    hRight = atan2d( ...
        eye.("optical axis right x"), ...
        eye.("optical axis right z"));

    hLeftCenter = hLeft(idxCenter);
    hRightCenter = hRight(idxCenter);

    hLeft_s = smoothdata(hLeftCenter, ...
        'movmedian', smoothWindow, 'omitnan');

    hRight_s = smoothdata(hRightCenter, ...
        'movmedian', smoothWindow, 'omitnan');

    thisColor = taskColors(r,:);

    plot(tCenterZeroed, hLeft_s, '--', ...
        'Color', thisColor, ...
        'LineWidth', 1.4);

    plot(tCenterZeroed, hRight_s, '-', ...
        'Color', thisColor, ...
        'LineWidth', 1.4);

end

yline(0, '--');

xlabel('Time relative to center fixation (s)');
ylabel('Horizontal eye angle (deg)');
title('Separate Left and Right Horizontal Eye Angles During Center Fixation');

legend({'Walk left','Walk right', ...
        'Read left','Read right', ...
        'Biopond left','Biopond right', ...
        'Chat left','Chat right'}, ...
        'Location', 'best');

grid on