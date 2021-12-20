function potential = calculatePotential(allPoints)
%   All point pairs for computing the total energy of the system
    pointPairs = nchoose2(1:size(allPoints, 2));
    
%   Pairwise distances
    r = sqrt(sum(reshape(diff(reshape(allPoints(:, pointPairs), [], 2), 1, 2) .^ 2, 3, [])));
    
%   Total system energy
    potential = sum(1./r);
end