function c = weights(z,x,m)
% Calculates FD weights. The parameters are:
% z location where approximations are to be accurate,
% x vector with x-coordinates for grid points,
% m highest derivative that we want to find weights for
% c array size m+1,lentgh(x) containing (as output) in
% successive rows the weights for derivatives 0,1,...,m.
n = length(x); c = zeros(m+2,n); c(2,1) = 1; x1 = x(ones(1,n),:);
A = x1'-x1; b = cumprod([ones(n,1),A],2); rm = ...
    cumsum(ones(m+2,n-1))-1; d = diag(b); d(1:n-1) = d(1:n-1)./d(2:n);
for i = 2:n
    mn = min(i,m+1);
    c(2:mn+1,i) = d(i-1)*(rm(1:mn,1).*c(1:mn,i-1)-(x(i-1)-z)*...
        c(2:mn+1,i-1));
    c(2:mn+1,1:i-1) = ((x(i)-z)*c(2:mn+1,1:i-1)-rm(1:mn,1:i-1).* ...
        c(1:mn,1:i-1))./(x(i)-x1(1:mn,1:i-1));
end
c(1,:) = [];
end

