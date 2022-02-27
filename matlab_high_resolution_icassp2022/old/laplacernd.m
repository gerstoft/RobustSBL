function [noise] = laplacernd(mu,std,r,c);
%
% generates laplacian random variables
% [noise] = laplacernd(mu,nu,i,j)
%
% input:
% mu - average
% nu - standard deviation
% i - # of rows
% j - # of columns
%
% output:
% noise - ixj matrix of laplacian random numbers
%
% E. Hetland, Univ Michigan, GeoSci468-W10
% ehetland@umich.edu
%

lam = std/sqrt(2);

if c~=1
  noise = [exprnd(lam,r,floor(c/2)) ...
	   -exprnd(lam,r,ceil(c/2))];
  for k=1:size(noise,1)
    noise(k,:) = noise(k,randperm(size(noise,2)));
  end
else
  noise = [exprnd(lam,floor(r/2),c); ...
	   -exprnd(lam,ceil(r/2),c)];
  noise = noise(randperm(size(noise,1)),:);
end

noise = noise+mu;

return