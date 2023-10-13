function orientation = determineOrientation(DNA)
%%  Summary
%
%   Inputs:
%       DNA: Described in parent function.
%
%   Outputs:
%       orientation: Whether the condensing agent bends the DNA clockwise
%       (-1) or counterclockwise (1).
%
%   Parent functions:
%       performTrial
%
%   Child functions:
%       computeCrossoverData
%
%%  Function

%   Calculate the crossover data for the DNA
    crossData = computeCrossoverData(DNA.planarCoords);
    
%   If the DNA is single looped and if the loop is counterclockwise
%   oriented, the orientation is 1. Otherwise, it is -1. This is the most
%   important factor for DNA orientation
    if(size(crossData, 1) == 1)
        if(ccwQ(crossData(1), crossData(2), crossData(3)))
            orientation = 1;
        else
            orientation = -1;
        end
% Else, take the largest bend angle
    else
        orientation = sign(DNA.vertexAngles(find(abs(DNA.vertexAngles)==max(abs(DNA.vertexAngles)))));
    end
end
