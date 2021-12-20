function closenessData = computeClosenessData(DNA)
%%  Summary
%
%   Inputs:
%       DNA: Described in parent function.
%
%   Outputs:
%       closenessData: A vector of two coordinates. They represent the
%       first and second close contacts on the DNA, ordered starting at the
%       origin. The function computes the closest DNA contacts that are
%       both far enough from each other along the DNA contour, and close
%       enough to each other in the plane.
%
%   Parent functions:
%       performTrial
%
%   Child functions:
%       None
%
%%  Function

%   Calculate linear coordinates
    linearCoords = cumsum([0, abs(diff(DNA.planarCoords))]);

%   Compute all combinations of DNA verticies that can be close to each
%   other
    vertexMatrix = nchoose2(1:length(linearCoords));
    
%   Store the complex coordinates of these vertices
    planarMatrix = DNA.planarCoords(vertexMatrix);
    
%   Determine the linear distances along the DNA contour between each pair
%   of vertices
    linearDistances = diff(linearCoords(vertexMatrix),1,2);
    
%   Determine the planar (absolute) distance between each pair of vertices
    planarDistances = abs(diff(planarMatrix,1,2));
    
%   To be considered as a close contact, two points must be at least 3
%   persistence lengths away from each other along the DNA contour, and
%   they must be less than thresholdDistance away from each other in the
%   plane
    validQ = (linearDistances > DNA.kLinear * DNA.persistenceLength) & (planarDistances < DNA.kPlanar * DNA.persistenceLength);

%   Find the index of the minimum distance pair among the valid pairs
    [~,minIndex] = min(planarDistances(validQ));

%   List the valid pairs of DNA vertices
	validComplexCoords = planarMatrix(validQ,:);
    
%   Extract the coordinates of the two closest contacts on the DNA
    closenessData = validComplexCoords(minIndex,:);
end