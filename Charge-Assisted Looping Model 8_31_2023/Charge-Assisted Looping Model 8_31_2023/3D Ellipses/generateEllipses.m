function allPoints = generateEllipses(m, n, C, l, angles)
%   The number of ellipses is determined by the number of rows in angles
    numEllipses = size(angles, 1);

%   All points on the ellipses
    allPoints = NaN(3, n*numEllipses);
    
%   The default (unrotated) ellipse lies in the xy plane with major axis
%   coincident with the x axis
    initialEllipsePoints = generateEllipsePoints(m, n, C, l);
    
%   Generate each ellipse
    for i = 1:numEllipses
        allPoints(:, (1:n) + n*(i-1)) = rotateEllipsePoints(initialEllipsePoints, angles(i, :));
    end
    
%   Plot all point charges (debug)
%   figure
%   plot3(allPoints(1, :), allPoints(2, :), allPoints(3, :), 'or')
%   axis equal
end