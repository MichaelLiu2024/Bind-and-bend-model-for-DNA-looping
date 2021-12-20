function experimentalData_217_398_1023 = importExperimentalData_217_398_1023()
%%  Summary
%
%   Inputs:
%       None
%
%   Outputs:
%       experimentalData_217_398_1023: Table containing the start sites and
%       loop circumferences of experimental molecules.
%
%   Parent functions:
%       Main
%
%   Child functions:
%       None
%
%%  Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 10);

%   Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";

%   Specify column names and types
    opts.VariableNames = ["circumferences_217", "circumferences_398", "circumferences_1023", "circumferences_398_2loop", "circumferences_398_3loop", "startSites_217", "startSites_398", "startSites_398_2loop", "startSites_398_3loop", "startSites_1023"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

%   Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

%   Import the data
    experimentalData_217_398_1023 = readtable('experimentalData_217_398_1023.txt', opts);

%   Convert diameters to fractional circumferences
    experimentalData_217_398_1023.circumferences_217 = experimentalData_217_398_1023.circumferences_217 * pi/217;
    experimentalData_217_398_1023.circumferences_398 = experimentalData_217_398_1023.circumferences_398 * pi/398;
    experimentalData_217_398_1023.circumferences_1023 = experimentalData_217_398_1023.circumferences_1023 * pi/1023;
    experimentalData_217_398_1023.circumferences_398_2loop = experimentalData_217_398_1023.circumferences_398 * pi/398;
    experimentalData_217_398_1023.circumferences_398_3loop = experimentalData_217_398_1023.circumferences_398 * pi/398;
end