function arcLength = computeArcLength(planarCoords)
%%  Summary
%
%   Inputs:
%       planarCoords: Piecewise linear curve. The program computes the
%       total arc length of this piecewise segment of DNA
%
%   Outputs:
%       arcLength: The total arc length of planarCoords
%
%   Parent functions:
%       computeStartSiteLoopCirc
%
%   Child functions:
%       None
%
%%  Function

%   Represent the DNA segments as complex vectors
    complexVectors = diff(planarCoords);
    
%   Compute the norm of each complex vector (i.e. the length of each DNA segment)
    segmentLengths = abs(complexVectors);
    
%   The arc length is the sum of all segment lengths
    arcLength = sum(segmentLengths);
end