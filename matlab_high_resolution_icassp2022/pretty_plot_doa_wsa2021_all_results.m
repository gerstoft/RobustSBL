% pretty_plot_doa_wsa2021_all_results.m

degrees = 180.0/pi;
close all

for icase=[6,7,8]
figure(icase);
switch(icase)
  case 6, datacase='MMV_L=25_Gaussian_noise'
  case 7, datacase='MMV_L=25_Complex-Student2p1_noise'
  case 8, datacase='MMV_L=25_epscontepsilon=5e-2_lambda=10_noise'
        end
M = load(['./figures/doa_wsa2021_SBL4-phase-only_for_' datacase]);
if ~exist('M.algorithm')
    M.algorithm = 'SBL4 phase-only';
end
N = load(['./figures/doa_wsa2021_SBL4_for_' datacase]);
if ~exist('N.algorithm')
    N.algorithm = 'SBL4';
end
% P = load('./figures/doa_wsa2021_SBL4-tM-estim2p1_for_MMV_L=25_Gaussian_noise.mat');
P = load(['./figures/doa_wsa2021_SBL4-tM-estim-Esa2p1_for_' datacase]);
if ~exist('P.algorithm')
    P.algorithm = sprintf('SBL4 t_\\nu{}M-estimator (\\nu = %3.1f)',P.nu_algorithm);
end
Q = load(['./figures/doa_wsa2021_SBL4-M-estim-Huber0p9_for_' datacase]);
if ~exist('P.algorithm')
    Q.algorithm = sprintf('SBL4 Huber M-estimator (q = %3.1f)',Q.nu_algorithm);
end
R = load(['./figures/doa_wsa2021_CBF_for_' datacase]);
if ~exist('Q.algorithm')
    R.algorithm = 'CBF';
end
G = load(['./figures/doa_wsa2021_SBL4-tM-estim-b-free2p1_for_' datacase]);
if ~exist('G.algorithm')
    G.algorithm = sprintf('SBL4 t_\\nu{}M-estimator (\\nu = %3.1f, b=1)',G.nu_algorithm);
end

% Cramer-Rao Lower Bound for single source and uniform linear array, Gaussian noise
%M.ASNR_lin = 10.^([-10 -5 0 M.SNR]/10) * M.antenna_array.N;
M.ASNR_lin = 10.^((M.SNR-10*log10(20))/10) * M.antenna_array.N;
M.CRLB = (1.0 / M.L) * ( 1.0 ./ (M.ASNR_lin) + 1.0 ./ (M.ASNR_lin).^2 )  * ( 6.0 / (M.antenna_array.N^2-1));
M.CRLB = M.CRLB / (pi^2 * (1-cosd(M.DOA_src)^2)); % correction by Esa 2021-10-03
% M.ASNR = 10*log10(M.ASNR_lin) - 10*log10(sqrt(M.antenna_array.N));
M.ASNR = 10*log10(M.ASNR_lin);

if 1
    icase,
    switch(icase)
        case 6, % 'Gaussian',
            SNR_shift = 0.0; % Gaussian noise variance is not same level as Heteroscedastic
        % case 'Laplace-like',
        %     SNR_shift = -10*log10(N*L); % Laplace-like noise variance is not same level as Heteroscedastic
        % case 'Heteroscedastic',
        %     SNR_shift = 0.0;
        case 7, % 'Complex-Student', % method in Ollila & Koivunen, PIMRC 2003
            SNR_shift = 0.0;        
        case 8, % 'epscont', % method in Ollila & Koivunen, PIMRC 2003
            M.epsilon_model, M.lambda_model,
            SNR_shift = 0.0;
            % SNR_shift = 10.0 * log10(1.0 - M.epsilon_model + M.epsilon_model * (M.lambda_model)^2);        
        otherwise,
            error('unknown model')
    end
    
    if SNR_shift ~=0,
        warning(['pretty plot shifts SNR by ',num2str(SNR_shift),' dB'])
    end
    
end

M.SNR = M.SNR + SNR_shift;
N.SNR = N.SNR + SNR_shift;
P.SNR = P.SNR + SNR_shift;

Q.SNR = Q.SNR + SNR_shift;
R.SNR = R.SNR + SNR_shift;
M.ASNR = M.ASNR + SNR_shift;
G.SNR = G.SNR + SNR_shift;
semilogy(M.SNR,M.RMSE_DOA_SBL4,'b-',N.SNR,N.RMSE_DOA_SBL4,'k*-',P.SNR,P.RMSE_DOA_SBL4,'ro-',G.SNR,G.RMSE_DOA_SBL4,'r+-',Q.SNR,Q.RMSE_DOA_SBL4,'mx-',R.SNR,R.RMSE_DOA_SBL4,'-',M.ASNR,sqrt(M.CRLB)*degrees,'--')
legend(M.algorithm,N.algorithm,P.algorithm,G.algorithm,Q.algorithm,R.algorithm,'CRB for Gaussian noise')

% semilogy(M.SNR,M.RMSE_DOA_SBL4,'b-',N.SNR,N.RMSE_DOA_SBL4,'k*-',P.SNR,P.RMSE_DOA_SBL4,'ro-',Q.SNR,Q.RMSE_DOA_SBL4,'mx-',M.ASNR,sqrt(M.CRLB)*degrees,'--')
% legend(M.algorithm,N.algorithm,P.algorithm,Q.algorithm,'CRB for Gaussian noise')

% semilogy(M.SNR,M.RMSE_DOA_SBL4,'b-',N.SNR,N.RMSE_DOA_SBL4,'k*-',P.SNR,P.RMSE_DOA_SBL4,'ro-')
% legend(M.algorithm,N.algorithm,P.algorithm)

if (M.L>1)
    title(sprintf('DOA estimates (%d DOAs, Simulation with %s noise, MMV, L=%d)',M.Number_of_DOAs,M.model,M.L));
else
    title(sprintf('DOA estimates (%d DOAs, Simulation with %s noise, SMV)',M.Number_of_DOAs,M.model));
end
xlabel('ASNR (dB)')
ylabel('RMSE of DOA Estimate (^\circ)')
grid on
% axis([min(M.SNR) max(M.SNR) 1e-2 1e+2])
axis([0 30 0.01 5])

print(['fig_RMSE_of_DOA_for_' datacase],'-dpng');
savefig(['fig_RMSE_of_DOA_for_' datacase]);
end