function crossData = computeCrossoverData(planarCoords)
%%  Summary
%
%   Inputs:
%       planarCoords: Described in parent function.
%
%   Outputs:
%       crossData: A vector of four coordinates, which represent the points
%       just before, at, just after, and at the intersection point, if it
%       exists. This data is needed later to compute the start site, loop
%       circumference, and loop orientation
%
%   Parent functions:
%       performTrial
%
%   Child functions:
%       computeCrossoverPoint
%       nchoose2
%
%%  Function

%   Compute all combinations of segments that can intersect
    segmentMatrix = nchoose2(1:length(planarCoords)-1);

%   Replace each segment in the segmentMatrix with its two complex
%   coordinate endpoints
    coordinateMatrix(:,[2,4]) = planarCoords(segmentMatrix+1);
    coordinateMatrix(:,[1,3]) = planarCoords(segmentMatrix);

%   Calculate the DNA crossover point
    crossData = computeCrossoverPoint(coordinateMatrix);
end