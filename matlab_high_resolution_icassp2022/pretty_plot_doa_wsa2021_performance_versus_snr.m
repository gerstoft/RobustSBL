% pretty_plot_doa_wsa2021_performance_versus_snr.m

if 0
    model,
    switch(model)
        case 'Gaussian',
            SNR_shift = -10*log10(N*L); % Gaussian noise variance is not same level as Heteroscedastic
        case 'Laplace-like',
            SNR_shift = -10*log10(N*L); % Laplace-like noise variance is not same level as Heteroscedastic
        case 'Heteroscedastic',
            SNR_shift = 0.0;
        case 'Complex-Student', % method in Ollila & Koivunen, PIMRC 2003
            SNR_shift = 0.0;        
        otherwise,
            error('unknown model')
    end
    
    if SNR_shift ~=0,
        warning(['pretty plot shifts SNR by ',num2str(SNR_shift)])
    end
    
    SNR = SNR + SNR_shift;
end

figure(2);
semilogy(SNR,RMSE_DOA_SBL4,'o-')
legend(sprintf('%s%s',algorithm,nu_algorithm_string));

if (L>1)
    title(sprintf('DOA estimates (%d DOAs, Simulation with %s %s noise, L=%d MMV)',Number_of_DOAs,model,nu_model_string,L));
else
    title(sprintf('DOA estimates (%d DOAs, Simulation with %s %s noise, L=1 SMV)',Number_of_DOAs,model,nu_model_string));
end
xlabel('SNR (dB)')
ylabel('RMSE of DOA Estimate (^\circ)')
grid on
axis([min(SNR) max(SNR) 1e-2 1e+2])
print -dpng fig_RMSE_of_DOA.png
print -dpdf  fig_RMSE_of_DOA.pdf

