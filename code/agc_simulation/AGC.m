function [gain, exposure] = AGC(s, gain, exposure, speedSetting, fps)
% Description:
% The MATLAB implementation of our custom AGC
%
% Inputs:
%   s                       - Double. Represents the mean intensity 
%                             value of a given frame
%
%   gain                    - Double. Represents the current gain setting
%                           that led to s 
%
%   exposure                - Double. Represents the current exposure setting
%                           that led to s
%
%   speedSetting            - Double. Represents the current speed setting 
%                           of how to adjust the AGC
% Examples:
%{
	s = 123; 
    gain = 1.0; 
    exposure = 52; 
    speedSetting = 0.99;
	[gain, exposure] = AGC(s, gain, exposure, speedSetting); 
%}


        signalTarget = 127;
        gainRange = [1 10.666];
        exposureRange = [37,floor(1e6/fps)];
        signalRange = [0,255];
        
        % Calculate the adjustment
        correction = 1+(signalTarget-s)/signalTarget;

        fprintf('Correction_1 %f\n', correction)

        % Set the speed
        speed = speedSetting;
    
        % Move quickly if we are pegged at the signal range
        if s == signalRange(1) || s == signalRange(2)
            speed = speedSetting^3;
        end
    
        % Move quickly if we are close to the destination
        if abs(correction - 1)<0.25
            speed = speedSetting^2;
        end
    
        % Correct the correction
        correction = 1 + ((1-speed) * (correction-1));

        fprintf('Correction_2 %f\n', correction)

    
        % If correction > 1, it means we need to turn up gain or exposure.
        if correction > 1
            % First choice is to turn up exposure
            if exposure < exposureRange(2)
                exposure = exposure * correction;
                exposure = min([exposure,exposureRange(2)]);
                exposure = max([exposure,exposureRange(1)]);
            else
                gain = gain * correction;
                gain = min([gain,gainRange(2)]);
                gain = max([gain,gainRange(1)]);
            end
        end
    
        % If correction < 1, it means we need to turn down gain or exposure.
        if correction < 1
            % First choice is to turn down gain
            if gain > gainRange(1)
                gain = gain * correction;
                gain = min([gain,gainRange(2)]);
                gain = max([gain,gainRange(1)]);
            else
                exposure = exposure * correction;
                exposure = min([exposure,exposureRange(2)]);
                exposure = max([exposure,exposureRange(1)]);
            end
        end
    
end