close all
clear all

model = 'Laplace-like'
% model = 'Gaussian'

sigma = 1.000; %   0 dB standard deviation of the additive complex noise
sigma = 0.707; %  -3 dB
%sigma = 0.500; %  -6 dB
%sigma = 0.250; % -12 dB
%sigma = 0.125; % -18 dB

N = 20;  % no. of sensors in ULA
d = 0.5; % sensor spacing of ULA measured in wavelengths

% array steering matrix of size N x M
M = 1801;
dtheta=180/(M-1);
phi_vec = [-90:dtheta:90];

for m=1:M
    A(:,m)   = exp(-1j*2*pi*d*[0:N-1]'*sind(phi_vec(m)))/sqrt(N);
end

m_src = (M+1)/2; % broadside source
m_src = (M-1)/4+1; % 45 deg source
a_src = A(:,m_src);
x_src =  exp(1j*pi/4); % source amplitude, magnitude = 1

L = 1e4; % number of snapshots


SNR = 10*log10(norm(a_src*x_src,'fro')^2/sigma^2);

switch model,
    case 'Laplace-like',
        y = laplacelike_rand(a_src*x_src, sigma, N, L);
    case 'Gaussian',
        y = (a_src*ones(1,L)) * x_src + sigma * complex(randn(N,L),randn(N,L))/sqrt(2);
    otherwise,
        error(['unknown model ', model]);
end

% evaluate CBF and phase-only beampatterns for the all noisy snapshots

DOA_CBF   = zeros(1,L);
DOA_phase = zeros(1,L);
sidelobe_level_CBF   = zeros(1,L);
sidelobe_level_phase = zeros(1,L);

for ell=1:L
    
    b_CBF   = abs(A'*y(:,ell)).^2;
    [maxval,index] = max(b_CBF);
    b_CBF          = b_CBF/maxval;
    DOA_CBF(ell)   = phi_vec(index);
    
    b_phase = abs(A'*(y(:,ell)./abs(y(:,ell)))).^2;
    [maxval,index] = max(b_phase);
    b_phase        = b_phase/maxval;
    DOA_phase(ell) = phi_vec(index);
    
    % evaluation of maximum sidelobe levels
    phi_tol = 360/(pi*N); % half beamwidth (from-maximum-to-first-zero) of CBF
    sidelobe_level_CBF(ell)   = maximum_sidelobe_level(10*log10(b_CBF),phi_vec, phi_vec(m_src),phi_tol);
    sidelobe_level_phase(ell) = maximum_sidelobe_level(10*log10(b_phase),phi_vec, phi_vec(m_src),phi_tol);
    
    %fprintf('maximum sidelobe level CBF    = %f dB\n', sidelobe_level_CBF(ell));
    %fprintf('maximum sidelobe level phase  = %f dB\n', sidelobe_level_phase(ell));
    
    % for ell=1000, 2000,...
    % make a pretty plot of both beamformers with indicated sidelobe levels
    if rem(ell,100)==0
        fprintf('%4d ',ell);
        if rem(ell,1000)==0
            fprintf('\n');
            figure(1); clf;
            pretty_plot_both_beamformers;
            drawnow;
        end
        
    end
end

fprintf('\n');

fprintf('mean(DOA_CBF)   = %f degrees\n', mean(DOA_CBF));
fprintf('median(DOA_CBF) = %f degrees\n', median(DOA_CBF));
fprintf('std(DOA_CBF)    = %f degrees\n', std(DOA_CBF));
fprintf('mad(DOA_CBF,1)  = %f degrees\n', mad(DOA_CBF,1));
fprintf('RMSE(DOA_CBF)   = %f degrees\n', sqrt(mean((DOA_CBF-phi_vec(m_src)).^2)));
fprintf('median(sidelobe_level_CBF) = %f dB\n',median(sidelobe_level_CBF));
fprintf('mean(sidelobe_level_CBF)   = %f dB\n\n',mean(sidelobe_level_CBF));

fprintf('mean(DOA_phase)   = %f degrees\n', mean(DOA_phase));
fprintf('median(DOA_phase) = %f degrees\n', median(DOA_phase));
fprintf('std(DOA_phase)    = %f degrees\n', std(DOA_phase));
fprintf('mad(DOA_phase,1)  = %f degrees\n', mad(DOA_phase,1));
fprintf('RMSE(DOA_phase)   = %f degrees\n', sqrt(mean((DOA_phase-phi_vec(m_src)).^2)));
fprintf('median(sidelobe_level_phase) = %f dB\n',median(sidelobe_level_phase));
fprintf('mean(sidelobe_level_phase)   = %f dB\n\n',mean(sidelobe_level_phase));


% testing of phase-only estimator for complex source amplitude, Equation (17)

x_CBF = a_src' * y;

c_phase = 1.0 ./ sum(1.0 ./ abs(y),1);  % Equation (20)
sign_y = y ./ abs(y);                   % Equation (21)
x_phase = c_phase .* (a_src' * sign_y); % Equation (19)

% Experimental stuff, closed form solution to Equation (30)
x_phase2 = complex(zeros(size(x_phase)),zeros(size(x_phase)));
for l=1:L
    aa = real(a_src' * sign_y(:,l));
    bb = imag(a_src' * sign_y(:,l));
    cc = real(a_src .* conj(sign_y(:,l)));
    dd = imag(a_src .* conj(sign_y(:,l)));
    uu = aa ./ ( 1.0 ./ c_phase(l) - sum(cc.^2 ./abs(y(:,l))));
    vv = bb ./ ( 1.0 ./ c_phase(l) - sum(dd.^2 ./abs(y(:,l))));
    x_phase2(l) = complex(uu,vv);
end

abs(mean(x_CBF) - x_src).^2;

abs(mean(x_phase) - x_src).^2;

abs(mean(x_phase2) - x_src).^2;

