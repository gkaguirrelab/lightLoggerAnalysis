fps = 120;
totalTime = 129.37; %sec
dotTimeFromEnd = [39, 35, 32, 28, 25, 21, 18, 14, 11, 8];
dotTimeFromEnd = 40.833;
dotTimeSec = totalTime - dotTimeFromEnd; % seconds
dotFrame = dotTimeSec.*fps;

pupil_features{1,2}.data{dotFrame}