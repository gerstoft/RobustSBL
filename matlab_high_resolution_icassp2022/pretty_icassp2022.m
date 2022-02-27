% pretty_icassp2022.m

degrees = 180.0/pi;
close all
figure; set(gcf,'position',[50,200,3000,700]);
for icase=[1,2,3]
subplot(1,3,icase,'pos', [0.04+0.32*(icase-1)    0.1100    0.31    0.8150]);
%subplot(1,3,icase);
switch(icase)
    case 1, datatitle='a) Gaussian noise'
    case 2, datatitle='b) MVT noise'
    case 3, datatitle='c) \epsilon contaminated noise'
end

switch(icase)
  case 1, datacase='MMV_L=25_Gaussian_noise'
  case 2, datacase='MMV_L=25_Complex-Student2p1_noise'
  case 3, datacase='MMV_L=25_epscontepsilon=5e-2_lambda=10_noise'
end

M = load(['./figures/doa_wsa2021_SBL4-phase-only_for_' datacase]);
if ~exist('M.algorithm')
    M.algorithm = 'SBL4 phase-only';
end
N = load(['./figures/doa_wsa2021_SBL4_for_' datacase]);
if ~exist('N.algorithm')
    N.algorithm = 'SBL4-G';
end
% P = load('./figures/doa_wsa2021_SBL4-tM-estim2p1_for_MMV_L=25_Gaussian_noise.mat');
P = load(['./figures/doa_wsa2021_SBL4-tM-estim-Esa2p1_for_' datacase]);
if ~exist('P.algorithm')
    P.algorithm = sprintf('SBL4-T (\\nu: %3.1f)',P.nu_algorithm);
    % P.algorithm = 'SBL4-T';
end
Q = load(['./figures/doa_wsa2021_SBL4-M-estim-Huber0p9_for_' datacase]);
if ~exist('P.algorithm')
    Q.algorithm = sprintf('SBL4-H (q: %3.1f)',Q.nu_algorithm);
    % Q.algorithm = 'SBL4-H';
end
R = load(['./figures/doa_wsa2021_CBF_for_' datacase]);
if ~exist('Q.algorithm')
    R.algorithm = 'CBF';
end
G = load(['./figures/doa_wsa2021_SBL4-tM-estim-b-free2p1_for_' datacase]);
if ~exist('P.algorithm')
    G.algorithm = sprintf('SBL4-T B-free (\\nu = %3.1f)',G.nu_algorithm);
end

% Cramer-Rao Lower Bound for single source and uniform linear array, Gaussian noise
M.ASNR_lin = 10.^(M.SNR/10) * M.antenna_array.N;
M.CRLB = (1.0 / M.L) * ( 1.0 ./ (M.ASNR_lin) + 1.0 ./ (M.ASNR_lin).^2 )  * ( 6.0 / (M.antenna_array.N^2-1)); % VanTrees Eq. (8.130)
M.CRLB = M.CRLB / (pi^2 * (1-cosd(M.DOA_src)^2)); % correction by Esa 2021-10-03
% M.ASNR = 10*log10(M.ASNR_lin) - 10*log10(sqrt(M.antenna_array.N));
M.ASNR = 10*log10(M.ASNR_lin);

%semilogy(M.SNR,M.RMSE_DOA_SBL4,'b-',N.SNR,N.RMSE_DOA_SBL4,'k*-',...
%    P.SNR,P.RMSE_DOA_SBL4,'ro-',Q.SNR,Q.RMSE_DOA_SBL4,'mx-',R.SNR,R.RMSE_DOA_SBL4,'-',M.ASNR,sqrt(M.CRLB)*degrees,'--','linewidth',2)

ULA_flag = 1; % chose 0 or 1, if =1 then evaluate CRLB analytically for ULA with half wavelength spacing
Gamma_src = diag(abs(M.x_src).^2); % source covariance matrix
M.CRLB = CRLB_VanTrees8p106(M.A, M.phi_vec, M.m_src, Gamma_src, M.ASNR_lin, M.L,ULA_flag);

semilogy(M.SNR,M.RMSE_DOA_SBL4,'b-',N.SNR,N.RMSE_DOA_SBL4,'k*-',...
    P.SNR,P.RMSE_DOA_SBL4,'ro-',Q.SNR,Q.RMSE_DOA_SBL4,'mx-',R.SNR,R.RMSE_DOA_SBL4,'-',M.ASNR,sqrt(M.CRLB)*degrees,'--', ...
    M.ASNR,sqrt(M.CRLB),'gx-','linewidth',3)


title(datatitle)
set(gca,'fontsize',28)
xlabel('ASNR (dB)')
if icase==1
    ylabel('RMSE of DOA Estimate (^\circ)')
legend(M.algorithm,N.algorithm,P.algorithm,Q.algorithm,R.algorithm,'CRB','CRB Eq.(8.106)')
else
    set(gca,'yticklabel', [])
end 
grid on
axis([min(M.SNR) max(M.SNR) 1e-2 1e+2])

%print(['fig_RMSE_of_DOA_for_' datacase],'-dpng');
%savefig(['fig_RMSE_of_DOA_for_' datacase]);
end
%print('fig_RMSE_DOA_icassp2022','-dpng');
%savefig(['fig_RMSE_DOA_icassp2022']);