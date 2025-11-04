function string = ndf2str(ndf)
% Converts a numeric NDF level to its string representation
%
% Syntax:
%   string = ndf2str(ndf)
%
% Description:
%  Converts a numeric NDF level representation to its string representation,
%  abiding by label conventions of using 0x to represent decimal values.
%
% Inputs:
%   ndf                   - Int/Float. Represents the NDF's numeric representation.                
%
% Outputs:
%   string                - String. Represents the NDF's string representation 
%                           with lab convention in mind. 
%
% Examples:
%{
    NDF = 0.2; 
    string = ndf2str(NDF);
%}
    
    string = strrep(num2str(ndf), '.','x');

end