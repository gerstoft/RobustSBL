
% evaluate SBL_v4 type algorithms and compare with Table 1 and phase-only processing

clear all
for imodel=2
    switch(imodel)
        case 2, model = 'Gaussian',                      nu_model_string = '';
        case 1, model = 'epscont', epsilon_model = 0.05, lambda_model  = 10.0, nu_model=0.05, nu_model_string = 'epsilon=5e-2_lambda=10', noise_enhancement = (1-epsilon_model) + epsilon_model*lambda_model^2,
        case 3, model = 'Complex-Student', nu_model = 2.1, nu_model_string = '2p1';
            % model = 'Heteroscedastic'   % Case III
            % model = 'Laplace-like',                  nu_model_string = '';
    end
    
for ialgorithm = 4 %[0,1,3,4,5]
    close all
    
    switch(ialgorithm)
        case 0, algorithm = 'CBF',                               nu_algorithm_string = '';
        case 1, algorithm = 'SBL4',                              nu_algorithm_string = '';
       % case 2, algorithm = 'SBL4-tM-estim', nu_algorithm = 2.1, nu_algorithm_string = '2p1';
        case 3, algorithm = 'SBL4-phase-only',                   nu_algorithm_string = '';
        case 4, algorithm = 'SBL4-tM-estim-Esa', nu_algorithm = 2.1, nu_algorithm_string = '2p1';
        case 5, algorithm = 'SBL4-M-estim-Huber',  nu_algorithm = 0.9, nu_algorithm_string = '0p9';
        otherwise, ialgorithm, error('unknown ialgorithm');
    end
    
    freq   =  2.0E+03,    % frequency (Hz)
    
    c0     = 343;         % speed of sound (m/s) in dry air at 20�C
    lambda = c0/freq;     % wavelength (m)
    wavenum= 2*pi/lambda; % wave number (rad/m)
    
    antenna_array.type = 'ULA',  % or 'ULA'
    
    switch(antenna_array.type)
        case 'UCA', % definition of uniform circular array geometry
            theta  = 30.0,        % elevation angle of incoming wave [degrees]
            antenna_array.N = 8;  % no. of antenna elements in UCA
            antenna_array.radius = 0.05; % 5cm radius == 10cm diameter
            for n=0:(antenna_array.N-1)
                phi = 2*pi*n / antenna_array.N; % polar angle of n-th sensor coordinates [rad]
                antenna_array.x(n+1) = antenna_array.radius * cos(phi);
                antenna_array.y(n+1) = antenna_array.radius * sin(phi);
                antenna_array.z(n+1) = 0.0;
            end
            antenna_array.d = norm([antenna_array.x(2) - antenna_array.x(1);
                antenna_array.y(2) - antenna_array.y(1);
                antenna_array.z(2) - antenna_array.z(1);
                ]) / lambda, % sensor spacing of UCA measured in wavelengths
            % array steering matrix of size N x M for all azimuth at one single elevation theta (degrees)
            M = 361;
            dphi=360/(M-1);
            phi_vec = [-180:dphi:180]; % sweep over all azimuth angles (degrees)
        case 'ULA', % definition of uniform linear array geometry
            theta  = 0.0,          % irrelevant elevation angle of incoming wave [degrees]
            antenna_array.N = 20;  % no. of sensors in ULA
            antenna_array.d = 0.5; % sensor spacing of ULA measured in wavelengths
            for n=0:(antenna_array.N-1)
                antenna_array.x(n+1) = n * antenna_array.d * lambda;
                antenna_array.y(n+1) = 0.0;
                antenna_array.z(n+1) = 0.0;
            end
            % array steering matrix of size N x M for all azimuth, the elevation is irrelevant.
            % M = 181;   % standard dictionary, 1 deg azimuth resolution
            M = 18001; % high resolution dictionary, 0.01 deg azimuth resolution
            dphi=180/(M-1);
            phi_vec = [-90:dphi:90];
    end
    
    for m=1:M
        kvec = wavenum * [sind(phi_vec(m))*cosd(theta);cosd(phi_vec(m))*cosd(theta);sind(theta)];
        A(:,m)   = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]); % normalization to |a_n|=1 or ||a||_2^2 = N.
    end
    
    
    L = 25,   % number of snapshots used by DOA estimator
    
    if(L>1)
        MMV = 1; % MMV is true, multiple measurement vector is true
    else
        MMV = 0; % MMV is false, single measurement vector
    end
    
    % noise standard deviation sigma
    if 1
        % choose coarse steps for quicker simulation
        sigma_vec = sqrt(antenna_array.N) * [0.125 0.250 0.500 0.600 0.707 1.000 2.000];
        sigma_vec = sqrt(antenna_array.N) * [0.125 1];
    else
        % choose fine steps for high accuracy plots
        sigma_vec = sqrt(antenna_array.N) * [0.0313 0.0442 0.0625 0.0884 0.1250 0.1768 0.2500 0.3536 0.5000 0.7071 1.0000 1.4142 2.0000 2.8284];
    end
    
    
    Number_of_DOAs = 1, % number of sources, use 1 or 3
    
    switch(Number_of_DOAs)
        case 1,
            DOA_src = -45;          % -45 deg source on the grid
            x_src   = exp(1j*pi/4); % source amplitude, magnitude = 1
        case 3,
            DOA_src = [-3; 2; 75];  % Gerstoft2015.pdf see FIGURE 8
            x_src   = (10 .^([12; 22; 20]/20)) * exp(1j*pi/4);
            x_src   = x_src / norm(x_src);
        otherwise,
            error('this Number_of_DOAs is not implemented')
    end
    
    Mcombinations  = nchoosek( M,  Number_of_DOAs); % number of all combinations of hypothetical DOAs
    Active_Set     = nchoosek(1:M, Number_of_DOAs); % the m-th active set is: Active_Set(m,:)
    fprintf(1,'All combinations of DOAs    = %d\n',Mcombinations);
    
    % clean the active set, get rid of all the bad condition numbers
    Icomp=zeros(1,Mcombinations);
    for m=1:Mcombinations
        Asmall=(A(:,Active_Set(m,:)));
        if cond(Asmall) < 1.0E+02, Icomp(m)=1; end
    end
    Active_Set = Active_Set(Icomp==1, :);
    Mcombinations = size(Active_Set, 1);
    fprintf(1,'reduced no. of combinations = %d\n',Mcombinations);
    
    for k=1:Number_of_DOAs
        m_src(k) = find(phi_vec == DOA_src(k));
    end
    
    a_src = A(:,m_src);
    
    SNR = 10*log10(norm(a_src*x_src,'fro')^2 ./ (sigma_vec.^2)),
    
    LL = 1e3, % number of array data vector observations
    
    mean_DOA_SBL4   = zeros(size(sigma_vec));
    median_DOA_SBL4 = zeros(size(sigma_vec));
    std_DOA_SBL4    = zeros(size(sigma_vec));
    mad_DOA_SBL4    = zeros(size(sigma_vec));
    RMSE_DOA_SBL4   = zeros(size(sigma_vec));
    median_sidelobe_level_SBL4 = zeros(size(sigma_vec));
    mean_sidelobe_level_SBL4   = zeros(size(sigma_vec));
    RMSE_x_SBL4     = zeros(size(sigma_vec));
    
    options = SBLSet;
    
    tic;
    
    for isnr = 1:length(sigma_vec)
        sigma = sigma_vec(isnr);
        
        % SNR(isnr) = 10*log10(norm(a_src*x_src,'fro')^2/sigma^2);
        
 %      rng(isnr, 'twister'); % mix the random number seeds
        rng(1, 'twister');    % always the same random number seed
        switch(model)
            case 'Laplace-like',
                if 0 
                   % deterministric source
                   y = laplacelike_rand(a_src*x_src, sigma, antenna_array.N, LL);
                else
                    % stochastic source
                    error('laplacelike_rand for stochastic source not yet implemented')
                end
            case 'Gaussian',
                noise_realization = sigma * complex(randn(antenna_array.N,LL),randn(antenna_array.N,LL))/sqrt(2);
                if 0, 
                    % deterministric source
                    y = (a_src * x_src * ones(1,LL)) + noise_realization;
                else
                    % stochastic source
                    s_src = complex(randn(1,LL),randn(1,LL))/sqrt(2);
                    y = ( a_src * s_src ) + noise_realization;
                end       
            case 'epscont',
                % epsilon = 0.001;
                % lambda  = 3.0;
                noise_realization = epscont(antenna_array.N * LL, epsilon_model, 0.0, lambda_model, sigma);
                noise_realization = reshape(noise_realization, antenna_array.N, LL);
                if 0,
                   y = (a_src * x_src * ones(1,LL)) + noise_realization;
                else
                    s_src = complex(randn(1,LL),randn(1,LL))/sqrt(2);
                    y = ( a_src * s_src ) + noise_realization;
                end
            case 'Complex-Student', % method in Ollila & Koivunen, PIMRC 2003
                noise_realization = sigma * complex(randn(antenna_array.N,LL),randn(antenna_array.N,LL))/sqrt(2);
                if 0, 
                    % deterministric source
                    y = (a_src * x_src * ones(1,LL)) + noise_realization;
                else
                    % stochastic source
                    s_src = complex(randn(1,LL),randn(1,LL))/sqrt(2);
                    y = ( a_src * s_src ) + noise_realization;
                end       
                % nu_model = 5, % choose number of degrees of freedom of chi^2 distribution
                s = ones(antenna_array.N, 1) * chi2rnd(nu_model * ones(1, LL));
                y =  y ./ sqrt(s/nu_model);
            case 'Heteroscedastic',
                noise_realization = sigma * complex(randn(antenna_array.N,LL),randn(antenna_array.N,LL))/sqrt(2);
                std_dev = 10.^(-1.0+2.0*rand(antenna_array.N,LL));
                std_dev = std_dev/sqrt(sum(sum(std_dev.^2))/(antenna_array.N*LL));
                if 0,
                   y = (a_src*ones(1,LL)) * x_src +  std_dev .* noise_realization;
                else
                    % stochastic source
                    s_src = complex(randn(1,LL),randn(1,LL))/sqrt(2);
                    y = ( a_src * s_src ) + std_dev .* noise_realization;
                end
            otherwise,
                error(['unknown model ', model]);
        end
        
        
        % evaluate SBL_v4 and compare with CBF, phase-only, and MLE
        
        DOA_SBL4  = zeros(1,LL);
        sidelobe_level_SBL4 = zeros(1,LL);
        
        MSE_x_SBL4   = 0;
        
        for ell = 1 : LL/L
            gamma = zeros(M,1);
            Y = y(:,ell+[0:(L-1)]);
            switch(algorithm)
                case 'CBF',                [gamma,report]  = SBL_CBF(A,Y,options);
                case 'SBL4',               [gamma,report]  = SBL_v4(A, Y, options);
                case 'SBL4-tM-estim',      [gamma,report]  = SBL_v4_tM(A, Y, nu_algorithm, options);
                case 'SBL4-phase-only',    [gamma,report]  = SBL_v4(A, Y./abs(Y), options);
                case 'SBL4-tM-estim-Esa',  [gamma,report]  = SBL_v4_tM_esa(A, Y, nu_algorithm, Number_of_DOAs);
                case 'SBL4-M-estim-Huber', [gamma,sigc,report] = SBL_v4_M_icassp( A , Y , 'Huber', nu_algorithm, Number_of_DOAs);
                otherwise, error('unknown algorithm');
            end
            
            b_SBL4 = gamma;
            [maxval,index] = max(b_SBL4);
            
%            if 0
%                b_SBL4        = b_SBL4/maxval;
%                DOA_SBL4_ell =
%                find_interpolated_max(phi_vec,b_SBL4,index); % with interpolation in vicinity of maximum
%            else
                 DOA_SBL4_ell = phi_vec(index); % without interpolation
%            end
             DOA_SBL4(1,ell) = DOA_SBL4_ell; 
            
            kvec  = wavenum * [sind(DOA_SBL4_ell)*cosd(theta);cosd(DOA_SBL4_ell)*cosd(theta);sind(theta)];
            a_SBL4 = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]).';  % offgrid SBL4 steering vector estimate
            
            a = conj(a_SBL4);
            % x_SBL4 =  gamma(index) * a_SBL4' * ((gamma(index)*a_SBL4*a_SBL4'+sigma*sigma*eye(N)) \ y(:,ell));
            x_SBL4 =  mean(a_SBL4' * Y/antenna_array.N); % here I use just the CBF, but this is not SBL3 should be corrected later
            
            MSE_x_SBL4 = MSE_x_SBL4 + abs(x_SBL4 - x_src).^2;
            % keyboard;
            % evaluation of maximum sidelobe levels
            phi_tol = 360/(pi*antenna_array.N); % half beamwidth (from-maximum-to-first-zero) of CBF
            sidelobe_level_SBL4(ell)   = maximum_sidelobe_level(20*log10(b_SBL4),phi_vec, phi_vec(m_src),phi_tol);
            
            % for ell=1000, 2000,...
            % make a pretty plot of all beamformers with indicated sidelobe levels
%             if rem(ell,100)==1
%                 fprintf('%4d ',ell);
%                 if rem(ell,1000)==1
%                     fprintf('\n');
%                     figure(1); clf;
%                     % pretty_plot_all_beamformers;
%                     plot(phi_vec , 10*log10(b_SBL4), 'b', phi_vec(m_src), 10*log10(max(b_SBL4))+3.0,'ro',DOA_SBL4_ell,10*log10(max(b_SBL4)),'kx');
%                     axis([min(phi_vec) max(phi_vec) -30 10])
%                     title(sprintf('%s%s gamma vector estimate, with %s %s noise, L=%d, SNR= %d dB',algorithm,nu_algorithm_string,model,nu_model_string,L,round(SNR(isnr))))
%                     xlabel('DOA (^\circ)')
%                     ylabel('gamma magnitude (dB)')
%                     drawnow;
%                 end
%                 
%             end
        end
        
        mean_DOA_SBL4(isnr)   = mean(DOA_SBL4(1:LL/L));
        median_DOA_SBL4(isnr) = median(DOA_SBL4(1:LL/L));
        std_DOA_SBL4(isnr)    = std(DOA_SBL4(1:LL/L));
        mad_DOA_SBL4(isnr)    = mad(DOA_SBL4(1:LL/L),1);
        RMSE_DOA_SBL4(isnr)   = sqrt(mean((DOA_SBL4(1:LL/L)-phi_vec(m_src)).^2));
        median_sidelobe_level_SBL4(isnr) = median(sidelobe_level_SBL4(1:LL/L));
        mean_sidelobe_level_SBL4(isnr)   = mean(sidelobe_level_SBL4(1:LL/L));
        RMSE_x_SBL4(isnr)     = sqrt(MSE_x_SBL4/L);
        
        fprintf('\n%s noise model: SNR = %f dB\n',model,SNR(isnr));
        
        fprintf('RMSE(DOA_SBL4)  = %f degrees\n', RMSE_DOA_SBL4(isnr));
        fprintf('RMSE x_SBL4     = %f \n\n',RMSE_x_SBL4(isnr));
%        save
        pretty_plot_doa_wsa2021_performance_versus_snr;
        elapsed_time = toc;
        elapsed_time_hrs = round(elapsed_time/360)/10;
        percentage_completed = round(100 * isnr/length(sigma_vec));
        remaining_time_hrs = (1.0 - isnr/length(sigma_vec)) * (elapsed_time/3600) / (isnr/length(sigma_vec));
        fprintf('elapsed time = %f hrs, completed = %d percent\n',elapsed_time_hrs, percentage_completed);
        fprintf('estimated remaining time = %f hrs\n', remaining_time_hrs);
   %     fprintf('predicted end of job = %s\n---\n', datestr(now + remaining_time_hrs/24));
    end
    toc
%    pretty_plot_doa_wsa2021_performance_versus_snr;
    
%     if(L>1) % multiple measurement vector
%         eval(sprintf('! mv -f fig_RMSE_of_DOA.pdf fig_RMSE_of_DOA_%s%s_MMV_L=%d_for_%s_%s_noise.pdf',algorithm,nu_algorithm_string,L,model,nu_model_string));
% %        eval(sprintf('! mv -f fig_RMSE_of_DOA.eps fig_RMSE_of_DOA_%s%s_MMV_L=%d_for_%s_%s_noise.eps',algorithm,nu_algorithm_string,L,model,nu_model_string));
%         save(sprintf('doa_wsa2021_%s%s_for_MMV_L=%d_%s%s_noise.mat',algorithm,nu_algorithm_string,L,model,nu_model_string));
%     else % single measurement vector
%         eval(sprintf('! mv -f fig_RMSE_of_DOA.pdf fig_RMSE_of_DOA_%s%s_SMV_L=%d_for_%s_%s_noise.pdf',algorithm,nu_algorithm_string,L,model,nu_model_string));
% %        eval(sprintf('! mv -f fig_RMSE_of_DOA.eps fig_RMSE_of_DOA_%s%s_SMV_L=%d_for_%s_%s_noise.eps',algorithm,nu_algorithm_string,L,model,nu_model_string));
%         save(sprintf('doa_wsa2021_%s%s_for_SMV_L=%d_%s%s_noise.mat',algorithm,nu_algorithm_string,L,model,nu_model_string));
%     end
    
end
end