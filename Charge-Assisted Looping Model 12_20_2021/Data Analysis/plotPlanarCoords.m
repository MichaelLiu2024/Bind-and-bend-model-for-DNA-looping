function plotPlanarCoords(DNA, crossoverCoords)
%%  Summary
%
%   Inputs:
%       DNA: Described in parent function.
%       crossoverCoords: Either crossoverData or closenessData, described
%       in their respective computation functions.
%
%   Outputs:
%       None (figures)
%
%   Parent functions:
%       performTrial
%
%   Child functions:
%       None
%
%%  Function

%   Create the figure
    figure;
    hold on
    axis equal

%   Plot the DNA contour
    plot(DNA.planarCoords)
    
%   Plot the protamine locations
    plot(DNA.planarCoords(DNA.occupiedVertices), 'ok', 'MarkerSize', 2)
    
%   Plot the crossover points or the closest DNA contacts
    plot(crossoverCoords, 'or')
end