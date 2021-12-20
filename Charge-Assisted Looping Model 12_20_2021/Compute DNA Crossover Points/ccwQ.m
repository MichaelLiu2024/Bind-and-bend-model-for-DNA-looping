function ccwQ = ccwQ(z1,z2,z3)
%%  Summary
%
%   Inputs:
%       z1, z2, z3: Three ordered points in the complex plane.
%
%   Outputs:
%       ccwQ: Whether or not the three ordered points z1, z2, and z3 are
%       oriented counterclockwise (1) or clockwise (0)
%
%   Parent functions:
%       computeCrossoverPoint
%
%   Child functions:
%       None
%
%%  Function

%   Two segments
    v2 = z2-z1;
    v3 = z3-z1;
    
%   2D cross product
    ccwQ = real(v2).*imag(v3) - real(v3).*imag(v2) > 0;
end

