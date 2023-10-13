function uElec = uElec(r, effectiveCharges, C)
%%  Summary
%
%   Inputs:
%       r: Vector of distances.
%       effectiveCharges: Numerator of Coulomb's Law.
%       C: Described in parent function.
%
%   Outputs:
%       uElec: A vector of electrostatic energies.
%
%   Parent functions:
%       computeBindingProbabilities
%
%   Child functions:
%       K
%
%%  Function

%   Compute electrostatic energy
    uElec = 1./(4*pi*K(r, C)*C.E0) .* effectiveCharges./(r*1e-9);
end
