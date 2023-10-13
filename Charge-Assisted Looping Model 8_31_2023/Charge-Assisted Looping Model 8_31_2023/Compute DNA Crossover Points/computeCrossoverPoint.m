function crossData = computeCrossoverPoint(coordinateMatrix)
%%  Summary
%
%   Inputs:
%       coordinateMatrix: A matrix of all combinations of segments in the
%       DNA. The function extracts the endpoints of the given segments,
%       performs a quick test to determine which segments intersect, and
%       calculates the intersection point, if it exists.
%
%   Outputs:
%       crossData: Described in parent function.
%
%   Parent functions:
%       computeCrossoverData
%
%   Child functions:
%       ccwQ
%
%%  Function

%   Extract the endpoints of the left segment
    z1 = coordinateMatrix(:,1);
    z2 = coordinateMatrix(:,2);
    
%   Extract the endpoints of the right segment
    z3 = coordinateMatrix(:,3);
    z4 = coordinateMatrix(:,4);

%   Quick test to find segments that intersect at non-endpoints
    intersectionIndices = find(ccwQ(z1,z3,z4)~=ccwQ(z2,z3,z4) & ccwQ(z1,z2,z3)~=ccwQ(z1,z2,z4) & z2~=z3);

%   Store the crossover data
    crossData = NaN(length(intersectionIndices),3);

%   Calculate the intersection point and associated data
    for ind = 1:length(intersectionIndices)
%       For ease of notation
        i = intersectionIndices(ind);
        
%       Solve for the intersection point with matrices
        A = [real(z2(i)-z1(i)),-real(z4(i)-z3(i));imag(z2(i)-z1(i)),-imag(z4(i)-z3(i))];
        B = [real(z3(i))-real(z1(i));imag(z3(i))-imag(z1(i))];
        X = linsolve(A,B);

%       Store the needed crossover data
        crossData(ind,:) = [z1(i),z1(i)+(z2(i)-z1(i))*X(1),z3(i)];
    end