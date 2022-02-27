% L1LinearRegression: Calculates L-1 multiple linear regression by IRLS
% by Will Dwinnell
% http://matlabdatamining.blogspot.com/2007/10/l-1-linear-regression.html
%
% B = L1LinearRegression(X,Y)
%
% B  = discovered linear coefficients
% X  = independent variables
% Y  = dependent variable
%
% Note 1: An intercept term is assumed (do not append a unit column).
% Note 2: a.k.a. LAD, LAE, LAR, LAV, least absolute, etc. regression
%
% Last modified: Mar-27-2009
%
% Example:
% n = 30;  rand('twister',2323);  randn('state',1865);
% BTrue = [1 -2];
% X = sort(rand(n,1));
% Y = BTrue(1) + X * BTrue(2) + 0.05 * randn(n,1);
% Y(1) = -0.5;
% Y(end) = 0.65;
%
% B = L1LinearRegression(X,Y);
% BS = [ones(n,1) X] \ Y;
%
% figure
% plot(X,Y,'bs',X,B(1) + X * B(2:end),'r-',X,BS(1) + X * BS(2:end),'k-')
% axis square
% grid on
% legend('Training Data','L-1 Regression','Least Squares Regression')

function B = L1LinearRegression(X,Y,Bstart)

% Determine size of predictor data
[n m] = size(X);

if nargin == 2, % modified by cfm
    % Initialize with least-squares fit
    % B       = [ones(n,1) X] \ Y;  % Least squares regression  % modified by cfm
    B       = X \ Y;  % Least squares regression
else
    B = Bstart; % modified by cfm
end
BOld    = B;
BOld(1) = BOld(1) + 1e-5;  % Force divergence

% Repeat until convergence
while (max(abs(B - BOld)) > 1e-6)
    % Move old coefficients
    BOld = B;
    
    % Calculate new observation weights (based on residuals from old coefficients)
    % W = sqrt(1 ./ max(abs((BOld(1) + (X * BOld(2:end))) - Y),1e-6));  % Floor to avoid division by zero
    W = sqrt(1 ./ max(abs((X * BOld) - Y),1e-6));  % Floor to avoid division by zero % modified by cfm
    
    % Calculate new coefficients
    % B = (repmat(W,[1 m+1]) .* [ones(n,1) X]) \ (W .* Y);
    B = (repmat(W,[1 m]) .* X) \ (W .* Y); % modified by cfm
end


% EOF


