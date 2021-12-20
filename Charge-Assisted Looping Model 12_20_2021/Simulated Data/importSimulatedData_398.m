function simulatedData_398 = importSimulatedData_398()
%%  Summary
%
%   Inputs:
%       None
%
%   Outputs:
%       simulatedData_398: Table containing the start sites and
%       loop circumferences of hand-counted, 2/3 loop flowers of 398 nm
%       long DNA.
%
%   Parent functions:
%       Main
%
%   Child functions:
%       None
%
%%  Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["ImageNumber", "Flower2Loop", "Flower3Loop", "Other"];
opts.VariableTypes = ["double", "double", "double", "categorical"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Other", "EmptyFieldRule", "auto");

% Import the data
simulatedData_398 = readtable('Flower Start Sites 25 Protamines.csv', opts);

end