function AGCSettings = cameraScoreToAGCSettings(cameraScore)

AGCSettings.exposure = 1;
AGCSettings.Again = 1;
AGCSettings.Dgain = 1;

AGCSettings.exposure = cameraScore;

if AGCSettings.exposure > 8333
    AGCSettings.exposure = 8333;
    AGCSettings.Again = cameraScore / 8333;
end

if AGCSettings.Again > 10.666
    AGCSettings.Again = 10.666;
    AGCSettings.Dgain = cameraScore / 8333 / 10.666;
end

end