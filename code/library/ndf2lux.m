function lux = ndf2lux(cal_filepath, options)
    arguments 
        cal_filepath {mustBeText}
        options.force_recalc (1, 1) {mustBeNumericOrLogical} = false; 
    end 

    % Load the XYZ fundamentals
    persistent T_xyz1931 S_xyz1931 
    fundamentals_not_defined = any([isempty(T_xyz1931), isempty(S_xyz1931)]); 

    % Load fundamentals if not loaded in before or forced to refresh 
    if(fundamentals_not_defined || options.force_recalc) 
        load('T_xyz1931.mat','T_xyz1931','S_xyz1931');
    end 

   % Load in the cal file
   cal_data = load(cal_filepath);
   if(~isfield(cal_data, 'cals') || isempty(cal_data.cals))
       error('ndf2lux:MissingCals', ...
             'Calibration file "%s" must contain a non-empty cals variable.', ...
             cal_filepath);
   end
   
   % Retrieve the ending cal
   cal = cal_data.cals{end};
   if(~isfield(cal, 'rawData') || ...
      ~isfield(cal.rawData, 'S') || ...
      ~isfield(cal.rawData, 'gammaCurveMeanMeasurements'))
       error('ndf2lux:InvalidCal', ...
             'Calibration file "%s" must contain rawData.S and rawData.gammaCurveMeanMeasurements.', ...
             cal_filepath);
   end

   S = cal.rawData.S;
   spd = cal.rawData.gammaCurveMeanMeasurements;

   % Calculate the luminance in lux
   T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
   xyYLocus = XYZToxyY(T_xyz);
   lux = T_xyz(2,:)*spd;
   
   return; 


end
