function planarCoords = computePlanarCoords(DNA)
%%  Summary
%
%   Inputs:
%       DNA: Described in parent function.
%       AGENT: Described in parent function.
%
%   Outputs:
%       planarCoords: Described in parent function.
%
%   Parent functions:
%       performTrial
%
%   Child functions:
%       None
%
%%  Function

%   Standard angles and lengths of each DNA segment
    segmentAngles = cumsum(DNA.vertexAngles(1:end-1));
    segmentLengths = diff(DNA.linearCoords);
    
%   Represent the DNA segments with complex vectors
    segmentVectors = segmentLengths .* exp(1i*segmentAngles);
    
%   Represent each DNA endpoint and protamine binding site with a complex coordinate
    planarCoords = cumsum([0, segmentVectors]);
end

