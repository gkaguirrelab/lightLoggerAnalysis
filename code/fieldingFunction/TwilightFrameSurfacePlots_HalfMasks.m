%% TwilightFrameSurfacePlots_HalfMasks.m

clear; close all; clc

load("selected_twilight_frames.mat")   % frames, times_sec, frames_idx, fps

frames = double(frames);
nFrames = size(frames, 1);

if nFrames ~= 20
    error('Expected 20 frames, but found %d frames.', nFrames)
end

[~, H, W] = size(frames);
frames_nan = frames;

%% Normalize original uncropped frames using the same 0-1 scale
globalMax = max(frames(:), [], 'omitnan');

for k = 1:nFrames
    Ik = squeeze(frames(k,:,:));
    Ik_norm = Ik ./ globalMax;

    figure
    imagesc(Ik_norm)
    axis image
    clim([0 1])

    cb = colorbar;
    cb.Label.String = 'Normalized Intensity';

    title(sprintf('Original uncropped frame %d at %.1f sec', k, times_sec(k)))
end

%% Apply half masks

for k = 1:nFrames
    Ik = squeeze(frames(k,:,:));

    if k <= 5
        Ik(:, round(W/2):end) = NaN;
        maskType = 'right half NaN';

    elseif k <= 10
        Ik(round(H/2):end, :) = NaN;
        maskType = 'bottom half NaN';

    elseif k <= 15
        Ik(:, 1:round(W/2)) = NaN;
        maskType = 'left half NaN';

    else
        Ik(1:round(H/2), :) = NaN;
        maskType = 'top half NaN';
    end

    frames_nan(k,:,:) = Ik;

    figure
    imagesc(Ik ./ globalMax)
    axis image
    clim([0 1])

    cb = colorbar;
    cb.Label.String = 'Normalized Intensity';

    title(sprintf('Masked frame %d at %.1f sec: %s', k, times_sec(k), maskType))
end

%% NaN-aware average and median

avg_img = squeeze(mean(frames_nan, 1, 'omitnan'));
med_img = squeeze(median(frames_nan, 1, 'omitnan'));

avg_norm = avg_img ./ max(avg_img(:), [], 'omitnan');
med_norm = med_img ./ max(med_img(:), [], 'omitnan');

valid_count = squeeze(sum(~isnan(frames_nan), 1));

figure
imagesc(valid_count)
axis image
cb = colorbar;
cb.Label.String = 'Number of Contributing Frames';
title('Number of valid frames contributing at each pixel')

figure
imagesc(avg_norm)
axis image
clim([0 1])
cb = colorbar;
cb.Label.String = 'Normalized Intensity';
title('Normalized NaN-aware average twilight image')

figure
surf(avg_norm, 'EdgeColor', 'none')
view(3)
clim([0 1])
cb = colorbar;
cb.Label.String = 'Normalized Intensity';
zlabel('Normalized Intensity')
title('Normalized average twilight surface')

figure
imagesc(med_norm)
axis image
clim([0 1])
cb = colorbar;
cb.Label.String = 'Normalized Intensity';
title('Normalized NaN-aware median twilight image')

figure
surf(med_norm, 'EdgeColor', 'none')
view(3)
clim([0 1])
cb = colorbar;
cb.Label.String = 'Normalized Intensity';
zlabel('Normalized Intensity')
title('Normalized median twilight surface')

save("twilight_halfmasked_average.mat", ...
    "frames_nan", "avg_img", "med_img", ...
    "avg_norm", "med_norm", "valid_count", ...
    "times_sec", "frames_idx", "fps")
