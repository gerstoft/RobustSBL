function [noise] = laplacelike_rand(mu,sigma,r,c);
% noise = laplacelike_rand(mu,std)
% noise = laplacelike_rand(mu,std,[r,c])
% noise = laplacelike_rand(mu,std,r,c)
%
% generates a matrix of complex-valued Laplace-like random variables
%
% input:
% mu    - mean (complex-valued)
%           mu may be a scalar
%           mu may be a column vector with r rows
%           mu may be a matrix of size (r x c)
% sigma - standard deviation of noise (real-valued positive scalar)
% r     - # of rows
% c     - # of columns
%
% output:
% noise - r x c matrix of complex-valued Laplace-like random numbers
%
% examples:
% noise = laplacelike_rand(0,1,1e6,1); mean(noise), var(noise),
% noise = laplacelike_rand(0,1,[1e6,1]); mean(noise), var(noise),
% noise = laplacelike_rand(zeros(1e6,1),1); mean(noise), var(noise),

% cfm
% Dec 08, 2018

if nargin == 2,
    [r,c] = size(mu);
elseif nargin == 3,
    c = r(2);
    r = r(1);
end

if size(mu) == [1,1],
    noise = mu * ones(r,c);
elseif size(mu) == [r,1],
    noise = mu * ones(1,c);
elseif size(mu) == [r,c],
    noise = mu;
else
    error(sprintf('size(mu)=[%d,%d], r=%d, c=%d',size(mu),r,c));
end

if real(sigma) <= 0
    error('sigma = %f <= 0')
end

A = 2;              % shape parameter of gamma distribution
B = sigma/sqrt(6);  % scale parameter of gamma distribution

R   = gamrnd(A,B,[r,c]); % magnitude is gamma distributed
PHI = 2*pi*rand(r,c);    % phase is uniformly distributed

noise = noise + R .* complex(cos(PHI),sin(PHI));

% END