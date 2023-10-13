function [startSite, loopCirc] = computeStartSiteLoopCirc(DNA, intersectionData)
%%  Summary
%
%   Inputs:
%       DNA: Described in parent function.
%       intersectionData: Either crossoverData or closenessData, described
%       in their respective computation functions.
%
%   Outputs:
%       startSite: Described in parent function.
%       loopCirc: Described in parent function.
%
%   Parent functions:
%       performTrial
%
%   Child functions:
%       computeArcLength
%
%%  Function

%   Extract the complex coordinates
    z1 = intersectionData(1);
    zI1 = intersectionData(2);
    z3 = intersectionData(3);
    zI2 = intersectionData(4);

%   Find the indicies of z1 and z3
    indexZ1 = find(DNA.planarCoords == z1);
    indexZ3 = find(DNA.planarCoords == z3);

%   Determine the complex coordinates that make up the two end segments and the DNA loop
    leftEnd = [DNA.planarCoords(1:indexZ1), zI1];
    loop = [zI1, DNA.planarCoords(indexZ1+1:indexZ3), zI2];
    rightEnd = [zI2, DNA.planarCoords(indexZ3+1:end)];
    
%   Calculate the start site and loop circumference
    startSite = min([computeArcLength(leftEnd), computeArcLength(rightEnd)]);
    loopCirc = computeArcLength(loop);
end