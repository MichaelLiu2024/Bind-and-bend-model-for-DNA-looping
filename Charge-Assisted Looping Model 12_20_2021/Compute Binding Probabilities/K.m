function K = K(r, C)
%%  Summary
%
%   Inputs:
%       r: Described in parent function.
%       C: Described in parent function.
%
%   Outputs:
%       K: A vector of dielectric constants.
%
%   Parent functions:
%       uElec
%
%   Child functions:
%       None
%
%%  Function

%   Compute the distance dependent dielectric constants.
    K = C.D - (C.D-1)/2 * ((r/C.H*C.C).^2 + 2*(r/C.H*C.C) + 2) .* exp(-(r/C.H*C.C));
end