function plotReconstructionStages(imageStages)

% How many stages did we get?
nStages = length(imageStages);

% Stage labels (we may or may not have the source, depending upon whether
% we are working with actual data or a simulation)
stages = {'raw','linear','flat','equalRGB','radiance','demosaic','impute','source'};

% Plot colors for the histograms
channelColor = {'r','g','b'};

fig = figure;

% Enable custom data tips to show unscaled original values from UserData
dcm = datacursormode(fig);
set(dcm, 'UpdateFcn', @imageDatatipCallback);

tiledlayout(2,nStages,"TileSpacing","tight");
for ss = 1:nStages

    nexttile(ss)

    % Get the image and some basic stats
    I = imageStages{ss};
    isinfI = isinf(I);
    nInf = sum(isinfI(:));
    maxVal = max(I(~isinfI));

    % Show the image uniformly using imshow with a visually scaled copy,
    % while storing the original unscaled image in UserData.
    if ndims(I) == 2

        % Compute display range upper bound (85th percentile of non-inf values)
        validVals = I(~isinf(I) & ~isnan(I));
        if isempty(validVals)
            p85 = 1;
        else
            p85 = prctile(validVals, 85);
        end
        if p85 == 0 || isnan(p85)
            p85 = 1;
        end

        % Normalize display copy for grayscale background
        normI = I / p85;
        normI(normI < 0) = 0;
        normI(normI > 1) = 1;
        normI(isnan(normI)) = 0;

        % Build RGB image for visualization
        R = normI;
        G = normI;
        B = normI;

        isZeros = (I == 0);

        % Color inf points red [1, 0, 0] and zeros blue [0, 0, 1]
        R(isinfI) = 1;   G(isinfI) = 0;   B(isinfI) = 0;
        R(isZeros) = 0;  G(isZeros) = 0;  B(isZeros) = 1;

        RGB = cat(3, R, G, B);

        hIm = imshow(RGB);
    else
        % Handle 3D RGB images with display scaling and solid red inf points
        validVals = I(~isinf(I) & ~isnan(I));
        if isempty(validVals)
            p85 = 1;
        else
            p85 = prctile(validVals, 85);
        end
        if p85 == 0 || isnan(p85)
            p85 = 1;
        end

        normI = I / p85;
        normI(normI < 0) = 0;
        normI(normI > 1) = 1;
        normI(isnan(normI)) = 0;

        R = normI(:,:,1);
        G = normI(:,:,2);
        B = normI(:,:,3);

        isinf3D = any(isinf(I), 3);

        % Color inf points red [1, 0, 0]
        R(isinf3D) = 1;
        G(isinf3D) = 0;
        B(isinf3D) = 0;

        RGB = cat(3, R, G, B);
        hIm = imshow(RGB);
    end

    % Retain original floating-point values in UserData for the custom data tip
    hIm.UserData = I;

    title([sprintf('%d. ',ss),stages{ss}]);

    % Show a histogram (handles both 2D Bayer and 3D RGB images)
    nexttile(ss+nStages)
    edges = (0:255)/255;
    minVals = zeros(1, 3);
    maxVals = zeros(1, 3);

    if ndims(I) == 2
        I_hist = I;
        I_hist(isinfI) = nan;
        [rgbIdx{1}, rgbIdx{2}, rgbIdx{3}] = returnBayerIndices(I_hist, 'BGGR');
        
        for cc = 1:3
            vec = I_hist(rgbIdx{cc});
            minVals(cc) = round(min(vec(:)));
            maxVals(cc) = round(max(vec(:)));

            vec(isnan(vec)) = maxVal;
            vec = vec(:) / maxVal;
            N = histcounts(vec, edges);
            N = N ./ length(rgbIdx{cc});

            plot(edges(1:end-1)*100, N*100, ['.-' channelColor{cc}]);
            hold on
        end
    else
        for cc = 1:3
            vec = I(:,:,cc);
            isinf_cc = isinf(vec);
            vec(isinf_cc) = nan;

            minVals(cc) = round(min(vec(:)));
            maxVals(cc) = round(max(vec(:)));

            vec(isnan(vec)) = maxVal;
            vec = vec(:) / maxVal;
            N = histcounts(vec, edges);
            N = N ./ numel(I(:,:,cc));

            plot(edges(1:end-1)*100, N*100, ['.-' channelColor{cc}]);
            hold on
        end
    end

    % Add min and max text to the upper left corner
    textStr = {sprintf('min = [%d, %d, %d]', minVals(1), minVals(2), minVals(3)), ...
        sprintf('max = [%d, %d, %d]', maxVals(1), maxVals(2), maxVals(3)),...
        sprintf('ceil, floor = [%d, %d]', nInf, sum(I(:) == 0))};
    text(0.25, 0.95, textStr, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);

    if ss == 1
        ylabel('Percentage of pixels');
        xlabel('Percentage max value');
        a = gca();
        a.TickDir = "out";
    else
        axis off
    end
    box off
end

end

% Callback function to override data tip text with original unscaled values
function txt = imageDatatipCallback(~, event_obj)
    pos = event_obj.Position;
    c = round(pos(1));
    r = round(pos(2));
    target = event_obj.Target;
    
    if isgraphics(target, 'image') && ~isempty(target.UserData)
        img = target.UserData;
        sz = size(img);
        if r >= 1 && r <= sz(1) && c >= 1 && c <= sz(2)
            if ndims(img) == 3
                val = squeeze(img(r, c, :));
                valStr = sprintf('[%g, %g, %g]', val(1), val(2), val(3));
            else
                valStr = num2str(img(r, c));
            end
            txt = { ...
                sprintf('X: %g', pos(1)), ...
                sprintf('Y: %g', pos(2)), ...
                sprintf('Value: %s', valStr) ...
            };
            return;
        end
    end
    txt = {sprintf('X: %g', pos(1)), sprintf('Y: %g', pos(2))};
end