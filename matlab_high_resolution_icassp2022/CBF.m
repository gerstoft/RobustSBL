function [DOA,Am] = CBF(Y) 

[N,L] = size(Y);
M = 18001; % high resolution dictionary, 0.01 deg azimuth resolution
dphi=180/(M-1);
phi_vec = -90:dphi:90;
A = exp(-1j*pi*(0:(N-1))'*sind(phi_vec)); % ULA matrix

RY = (1/L)*Y*(Y'); 
% Other option I tried when N > L is to compute M-estimate of scatter: 
% RY = MVT(Y.','nu',nu,'scaling','gaussian');
% and use it in place of SCM  above but this gave worse results.
CBF = real(sum(conj(A).*(RY*A)));
[~,indx] = max(CBF);
Am = A(:,indx);        % only active replicas
DOA = phi_vec(indx);
end
