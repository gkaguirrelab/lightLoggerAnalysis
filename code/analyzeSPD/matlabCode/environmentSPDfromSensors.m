function [sourceSPD,detectorSPD,modelS,sourceWeights,fVal] = environmentSPDfromSensors( chunks, options )
% SPD from natural and artificial light sources to fit minispect weights
%
% Syntax:
%   [estimatedSPD,enviroS,sourceWeights,fVal] = environmentSPDfromMiniSpect( detectorWeights )
%
% Description:
%   We wish to reconstruct the environmental spectral power distribution
%   from the weights on a set of narrow-band channels from the minispect.
%   To constrain this inverse mapping, we implememnt a forward model of the
%   environmental SPDs that might be created from linear combinations of
%   daylight and artificial light sources. We use a non-linear optimization
%   to find the light source combinations that yields an SPD predicted to
%   produce a set of minispect weights that best matches the observed
%   values.
%
% Inputs:
%   none
%   foo                   - Scalar. Foo foo foo foo foo foo foo foo foo foo
%                           foo foo foo foo foo foo foo foo foo foo foo foo
%                           foo foo foo
%
% Optional key/value pairs:
%   none
%  'bar'                  - Scalar. Bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar
%
% Outputs:
%   none
%   baz                   - Cell. Baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz
%
% Examples:
%{
    load('/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_data/LightLogger_Recordings/quickOutdoorIndoor.mat');
    environmentSPDfromSensors( chunks);
%}

arguments
    chunks (1,:) {iscell}
    options.modelSet (1,:) char {mustBeMember(options.modelSet,{'CombiLED','polynomial','environment'})} = 'environment'
    options.cal (1,:) struct = [];
end

% Concatenate the data from the chunks
asm7341_counts = [];
asm7341_time = [];
imx219_counts = [];
imx219_time = [];

for cc = 1:length(chunks)
    t = chunks{cc}.M.t;
    v = chunks{cc}.M.v.AS;
    asm7341_counts = [asm7341_counts v'];
    asm7341_time = [asm7341_time t];

    t = chunks{cc}.W.t;
    v = chunks{cc}.W.v;
    imx219_counts = [imx219_counts v];
    imx219_time = [imx219_time t];

end


% Load or create the model set
switch options.modelSet
    case 'environment'
        enviroSPDPath = fullfile(tbLocateProjectSilent('combiExperiments'),'data','CIEDaylightComponents.mat');
        load(enviroSPDPath,'CIEDaylightComponents_T');
        modelS = WlsToS(CIEDaylightComponents_T.wls);
        modelP = CIEDaylightComponents_T{:,2:end};
        % Set the bounds. For the daylight components, it is possible to
        % have negative loadings on the 2nd and 3rd component. Otherwise,
        % all components are bound by zero.
        lb = [0 -Inf -Inf];
        ub = [Inf Inf Inf];
        x0 = [1 1 1];
    case 'polynomial'
        degree = 10;
        modelS = [300,1,501];
        shapefcn = polyBasis('chebyshev',degree);
        modelP = zeros(modelS(3),degree+1);
        for ii = 0:modelS(3)-1
            idx = (ii-modelS(3)/2)/(modelS(3)/2);
            modelP(ii+1,:) = [1 shapefcn(idx)];            
        end
        lb = -inf(1,degree+1);
        lb(1) = 0;
        ub = inf(1,degree+1);
        x0 = zeros(1,degree+1);
    case 'CombiLED'
        modelS = options.cal.rawData.S;
        modelP = options.cal.processedData.P_device;
        lb = [0 0 0 0 0 0 0 0];
        ub = [Inf Inf Inf Inf Inf Inf Inf Inf];
        x0 = [1 1 1 1 1 1 1 1];
end

% Define the model wls for plotting later
modelWls = SToWls(modelS);

% Load the SM7341 SPDs and reformat to the model wavelength sampling
asm7341_spdPath = fullfile(tbLocateProjectSilent('combiExperiments'),'data','ASM7341_spectralSensitivity.mat');
load(asm7341_spdPath,'T');
asm7341_S = WlsToS(T.wl);
asm7341_P = T{:,2:end};
asm7341_P = asm7341_P(:,1:8);
asm7341_meanWeight = mean(max(asm7341_P));
asm7341_P = asm7341_P/asm7341_meanWeight;
asm7341_nChannels = size(asm7341_P,2);
asm7341_P_resamp = [];
for ii = 1:asm7341_nChannels
    asm7341_P_resamp(:,ii) = interp1(SToWls(asm7341_S),asm7341_P(:,ii),SToWls(modelS));
end

% Load the IMX219 SPDs and reformat to the model wavelength sampling
imx219_spdPath = fullfile(tbLocateProjectSilent('combiExperiments'),'data','IMX219_spectralSensitivity.mat');
load(imx219_spdPath,'T');
imx219_S = WlsToS(T.wls);
imx219_P = T{:,2:end};
imx219_nChannels = size(imx219_P,2);
imx219_P_resamp = [];
for ii = 1:imx219_nChannels
    imx219_P_resamp(:,ii) = interp1(SToWls(imx219_S),imx219_P(:,ii),SToWls(modelS));
end

% Create a forward model of the sensor weights based upon a combination
% of souces / components. 
mySPD = @(x) modelP*x';
myASM_weights = @(x) mySPD(x)'*asm7341_P_resamp;
myIMX_weights = @(x) mySPD(x)'*imx219_P_resamp;
myContrastUnits = @(v) (v-mean(v))/mean(v);
mySPDwiggle = @(x) max([0,std(myContrastUnits(mySPD(x)))]);
myNonlcon = @(x) negativeWeights(x,mySPD,mySPDwiggle);


% Set the fmincon options
options = optimset('fmincon');
options.Display = 'off';

figure
plot(asm7341_time,asm7341_counts);

% Loop over the measures
figure
for mm = 1:length(asm7341_time)

    % Get these asm7341 sensor weights
    asm7341_weights = asm7341_counts(1:8,mm);

    % Define the objective function
    myObj = @(x) norm(asm7341_weights-myASM_weights(x)');

    % Find the source weights that best fit the observed sensor
    % weights
    [sourceWeights,fVal] = fmincon(myObj,x0,[],[],[],[],lb,ub,[],options);

    % Make this solution the start for the next search (if it is a good
    % solution)
    if mySPDwiggle(sourceWeights)<0.5
        x0 = sourceWeights;
    else
        % Define an x0 by simple linear regression
        x1 = (modelP\(asm7341_weights'*asm7341_P_resamp')')';
        x1 = x0;
        [sourceWeights,fVal] = fmincon(myObj,x1,[],[],[],[],lb,ub,myNonlcon,options);
        if mySPDwiggle(sourceWeights)<0.5
            x0 = sourceWeights;
        end
    end
    
    sourceSPD = mySPD(sourceWeights);
    detectorSPD = (myASM_weights(sourceWeights)*asm7341_P_resamp')';

    % Plot this
    subplot(1,2,1)
    yyaxis left
    sourceSPD(sourceSPD<=0)=1e-4;
    plot(modelWls,asm7341_weights'.*asm7341_P_resamp,'-r');
    yyaxis right
    plot(modelWls,(sourceSPD),'-k','LineWidth',2)
    
%    ylim([-1 3]);
    title(sprintf('mySPDwiggle: %2.2f sd',mySPDwiggle(sourceWeights)))
    subplot(1,2,2)
    plot(log10(asm7341_weights),log10(myASM_weights(sourceWeights)),'*');
    xlim([-1 5]);
    ylim([-1 5]);
    hold on
    refline(1,0);
    hold off
    title(sprintf('%2.2f secs',asm7341_time(mm)))
    if mySPDwiggle(sourceWeights)>0.75
        foo=1;
    end
    pause
end

end


%% LOCAL FUNCTIONS

function [c, ceq] = negativeWeights(x,mySPD,mySPDwiggle)
spd = mySPD(x);
c = sum(spd<=0);
ceq = [];
% if mySPDwiggle(x)>0.5
%     ceq = mySPDwiggle(x);
% else
%     ceq = 0;
% end
end