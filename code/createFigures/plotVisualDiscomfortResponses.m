%
% This script parses and plots graphs for measurements of visual discomfort 
% taken from both Migraine and Control participants.
%  

filename = "FLIC_VisualDiscomfort_Data.xlsx";
T = readtable(filename, "VariableNamingRule", "preserve");

ids = string(T{:,1});
envNames = T.Properties.VariableNames(2:end);

% data matrix
data = T{:,2:end};

% group participants by ID
isControl = startsWith(ids, "FLIC_0");
isMigraine = startsWith(ids, "FLIC_1");

nControl = sum(isControl);
nMigraine = sum(isMigraine);

figure;
hold on;

nEnv = length(envNames);
offset = 0.18;

for e = 1:nEnv
    controlVals = data(isControl,e);
    migraineVals = data(isMigraine,e);

    controlVals = controlVals(~isnan(controlVals));
    migraineVals = migraineVals(~isnan(migraineVals));

    % x positions for each group
    xControl = (e-offset) * ones(size(controlVals));
    xMigraine = (e+offset) * ones(size(migraineVals));

    % individual participant dots
    plot(xControl, controlVals, "o", ...
        "Color", "b", "MarkerFaceColor", "b", "MarkerSize", 5);

    plot(xMigraine, migraineVals, "o", ...
        "Color", "r", "MarkerFaceColor", "r", "MarkerSize", 5);

    % mean and standard deviation
    meanC = mean(controlVals);
    semC = std(controlVals) / sqrt(length(controlVals));

    meanM = mean(migraineVals);
    semM = std(migraineVals) / sqrt(length(migraineVals));

    % mean +/- SD
    errorbar(e-offset, meanC, semC, semC, ...
        "o", "Color", "b", "MarkerFaceColor", "b", ...
        "MarkerSize", 9, "LineWidth", 2, "CapSize", 10, ...
        "LineStyle", "none");

    errorbar(e+offset, meanM, semM, semM, ...
        "o", "Color", "r", "MarkerFaceColor", "r", ...
        "MarkerSize", 9, "LineWidth", 2, "CapSize", 10, ...
        "LineStyle", "none");
end

% axis formatting
ylim([1 10]);
yticks(1:10);

xlim([0.5 nEnv+0.5]);
xticks(1:nEnv);
xticklabels(envNames);
xtickangle(35);

ax = gca;
ax.XAxis.FontSize = 14;

ylabel("Response value","FontSize",16);
xlabel("Environment","FontSize",16);
title("Visual Discomfort by Environment and Participant Group","FontSize",18);

legend({sprintf('Control (n = %d)',nControl), ...
        sprintf('Migraine (n = %d)',nMigraine)}, ...
        'Location','best', 'FontSize', 14);

hold off;
