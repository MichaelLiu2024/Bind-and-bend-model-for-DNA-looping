function optimalParameters = optimizeEllipses(n, C, l)
%   Global optimization object
    globalSearch = GlobalSearch;
    
%   Objective function
    potentialFunction = @(mAngles) calculatePotential(generateEllipses(mAngles(1), n, C, l, [0, 0, 0; mAngles(2:4)]));
    
%   Optimization problem
    optimProblem = createOptimProblem('fmincon', 'x0', zeros(1, 4), 'objective', potentialFunction, 'lb', zeros(1, 4), 'ub', [1, repelem(2*pi, 3)]);
    
%   Run solver
    optimalParameters = run(globalSearch, optimProblem);
end

