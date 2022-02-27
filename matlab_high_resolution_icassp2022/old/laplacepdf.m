function [y] = laplacepdf(x,mu,std);
%
% generates laplacian pdf
% [y] = laplacepdf(x,mu,nu)
%
% input:
% x - positions to evaluate PDF at
% mu - average
% nu - standard deviation
% output:
% y - laplacian pdf
%
% E. Hetland, Univ Michigan, GeoSci468-W10
% ehetland@umich.edu
%

flp=0;
if size(x,1)>1
  x = x';
  flp=1;
end

lam = std/sqrt(2);

neg = find(x-mu<0);
pos = find(x-mu>=0);
y = zeros(size(x));
y(neg) = fliplr(exppdf(fliplr(abs(x(neg)-mu)),lam));
y(pos) = exppdf(x(pos)-mu,lam);

y = y./2;

return