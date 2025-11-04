function NDF = str2ndf(string)
% Converts the string representation of an NDF level to its numeric equivalent
%
% Syntax:
%   NDF = str2ndf(string)
%
% Description:
%   Converts a string representation of an NDF level,
%   according to the lab conventions, back to its numeric 
%   representation. 
%
% Inputs:
%   string                - String. The string representation of the NDF 
%                           level.             
%
% Outputs:
%   NDF                   - Double. The numeric representation of the 
%                           NDF level
%
% Examples:
%{
    string = '0x2';
    NDF = str2ndf(string);
%}
        
    NDF = str2double(strrep(string, 'x','.'));
    
end