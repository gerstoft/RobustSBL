%pretty_plot_all_beamformers

switch(Number_of_DOAs)
    case 1,
        if sidelobe_level_CBF(ell) > sidelobe_level_phase(ell)
            plot(phi_vec , 10*log10(b_CBF), 'b', phi_vec , 10*log10(b_phase), 'r', ...
                phi_vec , 20*log10(b_mle), 'k', ...
                phi_vec(m_src(1)) * [1 1],[-20 0],'k:', ...
                [phi_vec(1),phi_vec(M)], sidelobe_level_CBF(ell)   * [1 1], 'b--', ...
                [phi_vec(1),phi_vec(M)], sidelobe_level_phase(ell) * [1 1], 'r--');
            legend('conventional beamformer','phase only beamformer','MLE for Laplace-like noise',...
                'true DOA',...
                'sidelobe level CBF',...
                'sidelobelevel phase-only',...
                'Location','southeast');
        else
            plot(phi_vec , 10*log10(b_CBF), 'b', phi_vec , 10*log10(b_phase), 'r', ...
                phi_vec , 20*log10(b_mle), 'k', ...
                phi_vec(m_src(1)) * [1 1],[-20 0],'k:', ...
                [phi_vec(1),phi_vec(M)], sidelobe_level_phase(ell) * [1 1], 'r--', ...
                [phi_vec(1),phi_vec(M)], sidelobe_level_CBF(ell)   * [1 1], 'b--');
            legend('conventional beamformer','phase only beamformer','MLE for Laplace-like noise',...
                'true DOA',...
                'sidelobelevel phase-only',...
                'sidelobe level CBF',...
                'Location','southeast');
        end
        hold on
        plot(DOA_CBF(ell),0,'b*',DOA_phase(ell),0,'r*',DOA_mle(ell),0,'k*');
        title(sprintf('beampattern(%d) N=%d, M=%d, d=%g, %s noise, SNR=%3d dB',ell,N,M,d,model,round(SNR(isnr))))
        xlabel('DOA (degrees)')
        ylabel('10*log10(beampattern)')
        axis([-90 90 -20 0])
        
        % decide which beamformers succeed or fail
        % and write a comment into the plot
        if (abs(DOA_CBF(ell)-phi_vec(m_src)) > phi_tol) && (abs(DOA_phase(ell)-phi_vec(m_src)) > phi_tol) && ...
                (abs(DOA_mle(ell)-phi_vec(m_src)) > phi_tol)
            text(-80,-18,'all beamformers fail');
        elseif (abs(DOA_CBF(ell)-phi_vec(m_src)) < phi_tol) && (abs(DOA_phase(ell)-phi_vec(m_src)) < phi_tol) && ...
                (abs(DOA_mle(ell)-phi_vec(m_src)) < phi_tol)
            text(-80,-18,'all beamformers succeed');
        else
            if abs(DOA_CBF(ell)-phi_vec(m_src)) > phi_tol
                text(-80,-17,'CBF fails');
            else
                text(-80,-17,'CBF succeeds');
            end
            
            if abs(DOA_phase(ell)-phi_vec(m_src)) > phi_tol
                text(-80,-18,'phase-only fails');
            else
                text(-80,-18,'phase-only succeeds');
            end
            
            if abs(DOA_mle(ell)-phi_vec(m_src)) > phi_tol
                text(-80,-19,'MLE Laplace fails');
            else
                text(-80,-19,'MLE Laplace succeeds');
            end
        end
        
    case 3,
        % do nothing
        fprintf(1,'DOA_mle(:,%d) = [ %d; %d; %d ]\n', ell, DOA_mle(1:Number_of_DOAs,ell));
        semilogy(b_mle)
        title(sprintf('LAD(%d) N=%d, M=%d, d=%g, %s noise, SNR=%3d dB',ell,N,M,d,model,round(SNR(isnr))))
        xlabel('index of Active Set')
        ylabel('LAD objective function')
    otherwise,
        error(['Number of DOAs = ',num2str(Number_of_DOAs)]);
end % switch

