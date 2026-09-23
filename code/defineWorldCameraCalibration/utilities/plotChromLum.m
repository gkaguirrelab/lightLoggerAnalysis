function [photopicLuminanceCdM2_Y,chromaticity_xy] = plotChromLum(spd,S)

% Start with the Matlab chromaticity diagram
plotChromaticity;
hold on

% Load the XYZ fundamentals
load('T_xyz1931.mat','T_xyz1931','S_xyz1931');
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
xyYLocus = XYZToxyY(T_xyz);

% Calculate the luminance and the chromaticities
photopicLuminanceCdM2_Y = T_xyz(2,:)*spd;
chromaticity_xy = (T_xyz(1:2,:)*spd/sum(T_xyz*spd));

% Plot the locus of the spectra
plot(chromaticity_xy(1), chromaticity_xy(2), 'o','MarkerEdgeColor','w','MarkerFaceColor', [0.5 0.5 0.5],'LineWidth', 2, 'MarkerSize', 10);

% Add text to give the xy chromaticity values
text(chromaticity_xy(1), chromaticity_xy(2)+0.1,sprintf('[%2.2f, %2.2f]',chromaticity_xy(1:2)))

% Labels
xlabel('x chromaticity');
ylabel('y chromaticity');
title(sprintf('Luminance %2.1f cd/m^2',photopicLuminanceCdM2_Y))
axis square

end