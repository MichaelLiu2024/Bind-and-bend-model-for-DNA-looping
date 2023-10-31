function experimentalData_105 = importExperimentalData_105()
%%  Summary
%
%   Inputs:
%       None
%
%   Outputs:
%       experimentalData_105: Structure containing the fitted midlines to
%       AFM images of 105nm DNA molecules with various concentrations of
%       protamine.
%
%   Parent functions:
%       Main
%
%   Child functions:
%       None
%
%%  Function

%   Load all experimental data
    allData = load('experimentalData_105.mat').DNAprop;
    
%   Contour lengths in nm
    contourLengths = [allData.DNA_length];
    
%   Protamine concentrations in micromolar
    protConcs = [allData.prot_conc];
    
%   Fitted midlines in pixels
    allMidlines = {allData.midline};
    
%   Number of points in the midlines
    numPts = cellfun(@(x) size(x, 1), allMidlines);

%   Maximum number of points in the midlines
	maxNumPts = max(numPts);
    
%   Nanometers per pixel
    pixelSize = 3.90625;
    
%   Convert midlines to complex planar coordinates
    convertToPlanarCoords = @(arrayCoords) cellfun(@(x) pixelSize*(x(:, 1) + 1i*x(:, 2))', arrayCoords, 'UniformOutput', false);
    
%   Pad with NaNs
    convertToArray = @(cell) reshape(cell2mat(cellfun(@(x) [x(1:end), NaN(1, maxNumPts - length(x))], cell, 'UniformOutput', false)), maxNumPts, [])';

%   Take only 105 nm molecules that have more than 3 points on the midline
    long105 = numPts > 3 & contourLengths == 105;
    
%   Store final data
    experimentalData_105.prot0 = convertToArray(convertToPlanarCoords({allMidlines{long105 & protConcs == 0}}));
    experimentalData_105.prot2 = convertToArray(convertToPlanarCoords({allMidlines{long105 & protConcs == 0.2}}));
    experimentalData_105.prot6 = convertToArray(convertToPlanarCoords({allMidlines{long105 & protConcs == 0.6}}));
    experimentalData_105.prot20 = convertToArray(convertToPlanarCoords({allMidlines{long105 & protConcs == 2}}));
end
