function ellipsePoints = generateEllipsePoints(m, n, C, l)
%   Complete elliptic integral of the second kind
    [~, E] = ellipke(m);
%   Semi-major axis of DNA
    a = C/4/E;
%   Semi-minor axis of DNA
    b = sqrt((1-m)*a^2);

%   Each point on the DNA
    k = linspace(1, n, n);
    
%   Calculate the central angle theta of each point on the DNA. The initial guess for theta assumes that m = 0
    options = optimoptions('fsolve', 'JacobPattern', speye(n), 'Algorithm', 'trust-region-reflective', 'PrecondBandWidth', 0, 'Display', 'None');
    theta = fsolve(@(phi) ellipticE(phi, m) - C/n/a * k, C/n/a * k, options);
    
%   Convert theta to equidistant points on the ellipse
    ellipsePoints = [(l+a) + a*sin(theta); b*cos(theta); zeros(1, length(theta))];
end