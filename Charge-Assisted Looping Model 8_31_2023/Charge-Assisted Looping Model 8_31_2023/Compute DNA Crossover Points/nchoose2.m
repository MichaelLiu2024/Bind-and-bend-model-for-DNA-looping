function Y = nchoose2(X)
%%  Summary
%
%   NCHOOSE2 - all combinations of two elements
%       Y = NCHOOSE2(X) returns all combinations of two elements of the array X.
%  
%   All credit due to:
%   (c) Jos van der Geest
%   email: samelinoa@gmail.com
%   http://www.mathworks.uk/matlabcentral/fileexchange/authors/10584
%
%%  Function

%   By creating an (N*(N-1)/2)-by-2 index matrix I the output can be
%   retrieved directly. This index matrix I equals nchoosek(1:N, 2).
%   N is the number of elements if X. We create I step-by-step using
%   left-hand indexing.
                                  % Example for N = 4 ->
    V = numel(X)-1:-1:2 ;         % V : 3 2
    R = cumsum([1 V], 2) ;        % R : 1 4 6
    
%   Step 1 - create I, filling the two columns (c1 and c2)
    I(R,2) = [0 -V] + 1 ;         % -> c1: 0  0  0  0  0  0
                                  %    c2: 1  0  0 -2  0 -1
%   Step 2                              
    I(R,1) = 1 ;                  % -> c1: 1  0  0  1  0  1
                                  %    c2: 1  0  0 -2  0 -1
%   Step 3
    I(:,2) = I(:,2) + 1 ;         % -> c1: 1  0  0  1  0  1
                                  %    c2: 2  1  1 -1  1  0
%   Step 4
    I = cumsum(I, 1) ;            % -> c1: 1  1  1  2  2  3
                                  %    c2: 2  3  4  3  4  4
%   Now we use I to index directly into X to create the output
    Y = X(I) ;                    
end