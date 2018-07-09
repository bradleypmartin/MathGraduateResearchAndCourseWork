function [RBFFDweights example] = ...
    pluspolyweightsdzeroaltpoint(A,RHS,xpositions,ypositions,xeval,yeval,...
    polydegree,normfactor)

% Function PLUSPOLYWEIGHTSDXDY takes the "A" matrix for RBF-FD weight
% computation for dx or dy operators, augments it with up to POLYDEGREE 2D
% polynomials, and evaluates the RBF-FD stencil weights (as described in
% Fornberg and Lehto (2011).

% INPUTS:

% A - the RBF-FD interpolation matrix.
% RHS - the RHS of the RBF-FD weight evaluation (for dx or dy)
%       WITHOUT polynomials added.
% xpositions - (1,n) row vector of x-coordinates for stencil nodes.  Note:
%              the x-coordinate of the evaluation point is in (1,1).
% ypositions - same as above, for y-coordinates.
% polydegree - all support polynomials UP TO AND INCLUDING this degree are
%              added.
% normfactor - if desired, you can add a normalization factor to the x and
%              y variables in the polynomials.  Adding this will change
%              variables from x and y to x*normfactor and y*normfactor (in
%              the polynomials only).
% dxdyswitch - if the entry here is 0 (FALSE), the function helps create a
%              dx operator matrix.  If the entry is 1 (TRUE) the function
%              creates a dy matrix.

% OUTPUTS:

% RBFFDweights - this is a (1,n) row vector of stencil weights
%                corresponding to the n stencil nodes input into the
%                function.

numpolys = (polydegree+1)*(polydegree+2)/2;
n = size(A,1);
augmatrix = zeros(n + numpolys);

augPblock = ones(n,numpolys);
augRHSblock = ones(numpolys,1);

xpositionsmod = normfactor*(xpositions-xeval)';
ypositionsmod = normfactor*(ypositions-yeval)';

counter1 = 0;
counter2 = 1;

for j = 2:numpolys
    
    augPblock(:,j)=xpositionsmod.^(counter2-counter1).*...
        ypositionsmod.^(counter1);

    counter1 = counter1+1;

    if counter1 > counter2
        counter1 = 0;
        counter2 = counter2+1;
    end

end

augRHSblock(1,1) = 1;

counter1 = 0;
counter2 = 1;

for m = 2:numpolys
    
    augRHSblock(m,1)= ...
        0.^(counter2-counter1).*...
        0.^(counter1);             
    
    counter1 = counter1+1;
    
    if counter1 > counter2
        counter1 = 0;
        counter2 = counter2+1;
    end
end

augmatrix(1:n,1:n) = A;
augmatrix(1:n,(n+1):(n+numpolys)) = augPblock;
augmatrix((n+1):(n+numpolys),1:n) = augPblock';

newRHS = [RHS; augRHSblock];

RBFFDweightstemp = (augmatrix\newRHS)';

RBFFDweights = RBFFDweightstemp(1,1:n);

example = [augmatrix newRHS];

%augmatrix

%newRHS