function ellipsePoints = rotateEllipsePoints(ellipsePoints, angles)
%   Rotation matrices
    rX = [1, 0, 0; 0, cos(angles(1)), -sin(angles(1)); 0, sin(angles(1)), cos(angles(1))];
    rY = [cos(angles(2)), 0, sin(angles(2)); 0, 1, 0; -sin(angles(2)), 0, cos(angles(2))];
    rZ = [cos(angles(3)), -sin(angles(3)), 0; sin(angles(3)), cos(angles(3)), 0; 0, 0, 1];

%   Rotate the ellipsePoints
    ellipsePoints = rZ * rY * rX * ellipsePoints;
end