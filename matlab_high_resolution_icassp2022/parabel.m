function [Y_opt,X_opt] = parabel(X,Y)
%
%  Find Extremum of interpolating parabel through
%  three points with coordinates X(i),Y(i).
%  Bug/Feature: X is assumed to be equidistant
%
DX = X(2)-X(1);

B = (Y(3)-Y(1))/(2*DX);
A = (Y(3)-2*Y(2)+Y(1))/(2*DX*DX);

%
% Parabel:  Y = Y2 + B*(X-X2) + A*(X-X2)^2/2
%
X_off = -B/(2*A);
X_opt = X(2) + X_off;
Y_opt = Y(2) + B*X_off + A*X_off^2;
return;
