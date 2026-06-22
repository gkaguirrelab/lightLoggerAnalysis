%% TwilightFrameSurfacePlots.m

clear; close all; clc

load("selected_twilight_frames.mat")   % frames, times_sec, frames_idx, fps

frames = double(frames);
nFrames = size(frames, 1);

frames_nan = frames;

for k = 1:nFrames
    Ik = squeeze(frames(k,:,:));

    hFig = figure;
    imagesc(Ik)
    axis image
    colorbar
    title(sprintf('Frame %d at %.1f sec: draw objects to mask', k, times_sec(k)))

    mask_total = false(size(Ik));

    while true
        roi = drawfreehand;

        answer = questdlg('Keep this region?', ...
                          'Mask region', ...
                          'Keep','Redo','Done with image','Keep');

        if strcmp(answer, 'Keep')
            mask_total = mask_total | createMask(roi);

            hold on
            plot(roi.Position(:,1), roi.Position(:,2), 'r', 'LineWidth', 1.5)

        elseif strcmp(answer, 'Redo')
            delete(roi)
            continue

        elseif strcmp(answer, 'Done with image') || isempty(answer)
            delete(roi)
            break
        end

        answer2 = questdlg('Draw another region on this same image?', ...
                           'Continue masking?', ...
                           'Yes','No','Yes');

        if strcmp(answer2, 'No') || isempty(answer2)
            break
        end
    end

    Ik(mask_total) = NaN;
    frames_nan(k,:,:) = Ik;

    close(hFig)
end

avg_img = squeeze(mean(frames_nan, 1, 'omitnan'));
med_img = squeeze(median(frames_nan, 1, 'omitnan'));

figure
imagesc(avg_img)
axis image
colorbar
title('NaN-aware average twilight image')

figure
surf(avg_img, 'EdgeColor', 'none')
view(3)
colorbar
title('NaN-aware average twilight surface')

figure
imagesc(med_img)
axis image
colorbar
title('NaN-aware median twilight image')

figure
surf(med_img, 'EdgeColor', 'none')
view(3)
colorbar
title('NaN-aware median twilight surface')

save("twilight_masked_average.mat", ...
    "frames_nan", "avg_img", "med_img", "times_sec", "frames_idx", "fps")