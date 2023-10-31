function planarCoords = centerPlanarCoords(planarCoords)
%%  Summary
%
%   Inputs:
%       planarCoords: Coordinates of the DNA to be centered
%
%   Outputs:
%       planarCoords: Centered and aligned DNA coordinates
%
%   Parent functions:
%       plotDNAProperties
%
%   Child functions:
%       None
%
%%  Function

%   Shift to origin
    planarCoords = planarCoords - planarCoords(1);
    
%   Rotate the DNA such that the first segment lies on the x axis
    planarCoords = planarCoords * exp(1i * -angle(planarCoords(2)));
end
