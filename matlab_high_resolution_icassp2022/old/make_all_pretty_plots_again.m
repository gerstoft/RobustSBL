
if 1
clear all
% load figures_for_Gaussian_noise/doa_laplace_performance_versus_snr_L=5e4_Gaussian.mat
load figures_for_Gaussian_noise/doa_laplace_performance_versus_snr_with_SBL3_L=5e4_Gaussian.mat
pretty_plot_doa_laplace_performance_versus_snr;
% ! mv -f fig_*.eps fig_*.pdf figures_for_Gaussian_noise/
keyboard
end

if 1
clear all
load figures_for_Laplace-like_noise/doa_laplace_performance_versus_snr_L=5e4_Laplace-like.mat
clear MMV
% load figures_for_Laplace-like_noise/doa_laplace_performance_versus_snr_with_SBL3_L=5e4_Laplace-like.mat
% load figures_for_Laplace-like_noise/doa_laplace_performance_versus_snr_L=5e4_Laplace-like_old.mat % correct 0-24 dB, no SBL3
pretty_plot_doa_laplace_performance_versus_snr;
% ! mv -f fig_*.eps fig_*.pdf figures_for_Laplace-like_noise/
keyboard
end

if 1
clear all
% load figures_for_Heteroscedastic_noise/doa_laplace_performance_versus_snr_L=5e4_Heteroscedastic.mat
load figures_for_Heteroscedastic_noise/doa_laplace_performance_versus_snr_with_SBL3_L=5e4_Heteroscedastic
pretty_plot_doa_laplace_performance_versus_snr;
% ! mv -f fig_*.eps fig_*.pdf figures_for_Heteroscedastic_noise/
end
