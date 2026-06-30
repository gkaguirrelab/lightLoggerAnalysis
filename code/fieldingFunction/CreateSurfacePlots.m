%% TwilightFrameSurfacePlots_NoMask.m

clear; close all; clc

load("selected_twilight_frames.mat")

frames = double(frames);

globalCLim = [min(frames(:)) max(frames(:))];

nFrames = size(frames,1);

%% Display all extracted frames

for k = 1:nFrames

    Ik = squeeze(frames(k,:,:));

    figure
    imagesc(Ik)
    axis image
    clim(globalCLim)

    cb = colorbar;
    cb.Label.String = 'Pixel Intensity';

    title(sprintf('Frame %d at %.1f sec', ...
        k, times_sec(k)))

end

%% Average and median images

avg_img = squeeze(mean(frames,1));
med_img = squeeze(median(frames,1));

%% Plot average image

figure
imagesc(avg_img)
axis image
clim(globalCLim)

cb = colorbar;
cb.Label.String = 'Pixel Intensity';

title('Average image')

%% Plot average surface

figure
surf(avg_img,'EdgeColor','none')
view(3)

cb = colorbar;
cb.Label.String = 'Pixel Intensity';

xlabel('X Pixel')
ylabel('Y Pixel')
zlabel('Pixel Intensity')

title('Average image surface')

%% Plot median image

figure
imagesc(med_img)
axis image
clim(globalCLim)

cb = colorbar;
cb.Label.String = 'Pixel Intensity';

title('Median image')

%% Plot median surface

figure
surf(med_img,'EdgeColor','none')
view(3)

cb = colorbar;
cb.Label.String = 'Pixel Intensity';

xlabel('X Pixel')
ylabel('Y Pixel')
zlabel('Pixel Intensity')

title('Median image surface')

%% Save results

save("planetarium_average.mat", ...
    "avg_img", "med_img", ...
    "times_sec", "frames_idx", "fps")