function [handles] = combo_calc(planarCoords)
% Extract x and y coordinates from planarCoords
x_norm = rmmissing(real(planarCoords));
y_norm = rmmissing(imag(planarCoords));

seg = length(x_norm) - 2;

% ---------------------------------------------------------------
% Check knots spline increments over the fibril
% INTERVAL FUNCTION
for i = 1:(seg+1)
    dist = ( x_norm(i) - x_norm(i+1) ).^2 + ( y_norm(i) - y_norm(i+1) ).^2;
    intervals(i) = sqrt(dist);
end
handles.intervals = intervals;

% MIDPOINT-FLUCT
% --------------------------------------------------------------------
% -------Calculate the deviation of the fibril from secant midpoints

midpointX = 0;
midpointY = 0;
sec_length = 0;
lastpoint = length(x_norm);

for n = 2:lastpoint - 1
    for j = 1:lastpoint - n
        % coordinates of secant midpoint
        midpointX = (x_norm(j) + x_norm(j+n))/2;
        midpointY = (y_norm(j) + y_norm(j+n))/2;
        %secant length associated to previous midpoint
        sec_length = sqrt (( y_norm(j+n) - y_norm(j) ).^2 +...
            (( x_norm(j+n) - x_norm(j) ).^2 ));
        shortest_dist = 1e+10;
        for i = 1:lastpoint
            % dist is the distance between secant midpoint defined above and
            % each knot/point of the fibril spline
            dist = (( midpointY - y_norm(i) ).^2 +...
                (( midpointX - x_norm(i) ).^2));
            dist = sqrt (dist);
            % Select the closest knot of the fibril spline to the secant
            % midpoint and keeps the distance that separate them
            if (dist <= shortest_dist)
                shortest_dist = dist;
                closest_i = i;
            end
        end
        % refine shortest dist
        if closest_i ~= 1 && closest_i <= lastpoint -1
            ci = closest_i;
            Avec = sqrt ( (x_norm(ci) - midpointX).^2 + (y_norm(ci) - midpointY).^2 );
            Bvec = sqrt ( (x_norm(ci+1) - midpointX).^2 + (y_norm(ci+1) - midpointY).^2 );
            Cvec = sqrt ( (x_norm(ci+1) - x_norm(ci)).^2 + (y_norm(ci+1) - y_norm(ci)).^2 );
            cos_theta1 = (Avec.^2 + Cvec.^2 - Bvec.^2) ./ (2 * Avec * Cvec);
            theta1 = acos(cos_theta1);
            shortest_dist1 = Avec * sin(theta1);
            Bvec2 = sqrt ( (x_norm(ci-1) - midpointX).^2 + (y_norm(ci-1) - midpointY).^2 );
            Cvec2 = sqrt ( (x_norm(ci-1) - x_norm(ci)).^2 + (y_norm(ci-1) - y_norm(ci)).^2 );
            cos_theta2 = (Avec.^2 + Cvec2.^2 - Bvec2.^2) ./ (2 * Avec * Cvec2);
            theta2 = acos(cos_theta2);
            shortest_dist2 = Avec * sin(theta2);
            if shortest_dist1 <= shortest_dist2
                shortest_dist = shortest_dist1;
            else
                shortest_dist = shortest_dist2;
            end
        end
        % stores the result of each loop (one loop corresponding to one value
        % of j; ie the first of the two knots that belong to the secant, and
        % the loop being reproduced for each knot over the fibril spline)
        delta(j)= shortest_dist;
        secant(j) = sec_length;
    end
    deltas(:,n-1) = delta;
    secants(:,n-1) = secant;
end
handles.deltas = deltas;
handles.secants = secants;

% pick the 2 matrices and make one in 2 dimensions
A = deltas;
B = secants;
p = length(A);
j = 0;
for n = 0:(p - 1)
    for k = 1:(p - n)
        u = k + (n*p - j);  % jump to next column (-j in terms of rows)
        C(u, 1) = A(k, n + 1 );
        C(u, 2) = B(k, n + 1 );
    end
    j = j + n;
end
deviat_single = C;
handles.deviat_single = deviat_single;

%TANTAN-COREL
%--------------------------------------------------------------------
%-------Calculate the decay of tangent-tangent correlations over the fibril

intervals = handles.intervals;

for j = 1:ceil(lastpoint/2)
    
    for i = 1:lastpoint - j - 1
        % defines coordinates of each tangent vector along the spline
        spX = x_norm(i+1) - x_norm(i);
        spY = y_norm(i+1) - y_norm(i);
        spnextX = x_norm(i+j+1) - x_norm(i+j);
        spnextY = y_norm(i+j+1) - y_norm(i+j);
        % define tangent vectors to the spline at the i^th point
        tan_i = [ spX, spY ];
        tan_k = [ spnextX, spnextY ]; %k is the (i+j)th point of the fibril
        % scalar product of the normalized vectors
        scal_prod = dot(tan_i, tan_k) / ( (norm(tan_i)) * (norm(tan_k)) );
        % Return the tangent-tangent correl = cos of normalized vector for each
        % knot and the corresponding length over the fibril between the 2 knots
        tantan_corel(i) = scal_prod;
        % increment the value of the contour length by adding the value of one
        % more segment at each loop
        contour_length(i) = sum(intervals(1,i:i+j-1));
        
    end
    D(:,j) = tantan_corel;
    E(:,j) = contour_length;
    cosines(j) = mean(tantan_corel);
    contours(j) = mean(contour_length);
end

handles.cosines = cosines;
handles.contours = contours;

% %get exponential decay
% g = fittype('a*exp(-x/2/b)');
% f0 = fit(contours,cosines,g,'Start',[1 0.02]);
% handles.Lp=f0.b;

% pick the 2 matrices and make one in 2 dimensions

p = size(D,2);
j = 0;
for n = 0:(p - 1)
    for k = 1:(p - n)
        u = k + (n*p - j);  
        F(u, 1) = D(k, n + 1 );
        F(u, 2) = E(k, n + 1 );
    end
    j = j + n;
end
corel_single = F;
handles.corel_single = corel_single;

%CONTOUR-ENDEND
%--------------------------------------------------------------------
%-------Calculate the contour length vs. the end-to-end length

for j = 1:lastpoint - 2
    for i = 1:lastpoint - j - 1
        %secant length associated to previous midpoint
        end2end_length(i) = sqrt (( y_norm(i+j) - y_norm(i) ).^2 +...
            (( x_norm(i+j) - x_norm(i) ).^2 ));
        % increment the value of the contour length by adding the value of 
        % one more segment at each loop
        contour_length(i) = sum(intervals(1,i:i+j-1));
    end
    G(:,j) = end2end_length;
    H(:,j) = contour_length;
    end2end(j)=(mean(end2end_length))^2;
    contour(j)=mean(contour_length);
end
handles.end2end = end2end;
handles.contour = contour;

% pick the 2 matrices and make one in 2 dimensions
p = length(G);
j = 0;
for n = 0:(p - 1)
    for k = 1:(p - n)
        u = k + (n*p - j);  
        I(u, 1) = G(k, n + 1 );
        I(u, 2) = H(k, n + 1 );
    end
    j = j + n;
end
% end2cont_single = I;
% handles.end2cont_single = end2cont_single;

% extract the maximum contour length calculated, ie. the fibril length
goodcontour = max(I(:,2));
handles.goodcontour = goodcontour ;

% --------------------------------------------------------------------
% BELOW is the code that fill in the structure with all the elements, pertaining 
% to one fibril, calculated above in this very same function

% IMPORTANT to keep this

% handles.mat_length_struct(i).contourL = handles.goodcontour;
% handles.mat_intervals_struct(i).meanitv_fib = handles.mean_itv;
% handles.mat_sd_struct(i).meansd_fib = handles.sd_itv;
% 
% handles.bsplines_struct(i).bspline = handles.bspline_norm_coord;
% handles.bsplines2_struct(i).bspline = handles.bspline_norm_coord2;
% 
% handles.deviations_struct(i).deviafib = handles.deviat_single;
% handles.correlations_struct(i).corelfib = handles.corel_single;
% handles.wormlike_struct(i).wormfib = handles.end2cont_single;
end
