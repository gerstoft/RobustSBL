close all
clear all

model = 'Laplace-like'
% model = 'Gaussian'
% model = 'Heteroscedastic'   % Case III

N = 20;  % no. of sensors in ULA
d = 0.5; % sensor spacing of ULA measured in wavelengths

MMV = 1,

L = 4,   % number of snapshots used by DOA estimator

LL = 5e4, % number of array data vector observations, use 1e3 for debugging, 5e4 for publication quality plots

Number_of_DOAs = 1, % number of sources, use 1 or 3

c_phase_version = 2, % select 1 (c_phase=1) or 2 (c_phase depends on data magnitude)

% noise standard deviation sigma
if 0
    % choose coarse steps for quicker simulation
    sigma_vec = sqrt(N) * [0.125 0.250 0.500 0.600 0.707 1.000 2.000];
else
    % choose fine steps for high accuracy plots
    sigma_vec = sqrt(N) * [0.0625 0.088 0.125 0.177 0.250 0.354 0.500 0.550 0.600 0.650 0.707 0.750 1.000 1.414 2.000];
end

% array steering matrix of size N x M
M = 181;
dtheta=180/(M-1);
phi_vec = [-90:dtheta:90];

for m=1:M
    A(:,m)   = exp(-1j*2*pi*d*[0:N-1]'*sind(phi_vec(m))); % same normalization to |a_n|=1 as in the LaTeX code
end

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


mean_DOA_CBF   = zeros(size(SNR));
median_DOA_CBF = zeros(size(SNR));
std_DOA_CBF    = zeros(size(SNR));
mad_DOA_CBF    = zeros(size(SNR));
RMSE_DOA_CBF   = zeros(size(SNR));
median_sidelobe_level_CBF = zeros(size(SNR));
mean_sidelobe_level_CBF   = zeros(size(SNR));
RMSE_x_CBF     = zeros(size(SNR));

mean_DOA_phase   = zeros(size(SNR));
median_DOA_phase = zeros(size(SNR));
std_DOA_phase    = zeros(size(SNR));
mad_DOA_phase    = zeros(size(SNR));
RMSE_DOA_phase   = zeros(size(SNR));
median_sidelobe_level_phase = zeros(size(SNR));
mean_sidelobe_level_phase   = zeros(size(SNR));
RMSE_x_phase     = zeros(size(SNR));

mean_DOA_mle   = zeros(size(SNR));
median_DOA_mle = zeros(size(SNR));
std_DOA_mle    = zeros(size(SNR));
mad_DOA_mle    = zeros(size(SNR));
RMSE_DOA_mle   = zeros(size(SNR));
median_sidelobe_level_mle = zeros(size(SNR));
mean_sidelobe_level_mle   = zeros(size(SNR));
RMSE_x_mle     = zeros(size(SNR));

tic;

for isnr = 1:length(sigma_vec)
    sigma = sigma_vec(isnr);
        
    rng(isnr, 'twister'); % random number seed
    switch model,
        case 'Laplace-like',
            y = laplacelike_rand(a_src*x_src, sigma, N, LL);
        case 'Gaussian',
            y = (a_src * x_src * ones(1,LL)) + sigma * complex(randn(N,LL),randn(N,LL))/sqrt(2);
        case 'Heteroscedastic',
            std_dev = 10.^(-1.0+2.0*rand(N,LL));
            std_dev = std_dev/sqrt(sum(sum(std_dev.^2))/(N*LL));
            y = (a_src * x_src * ones(1,LL)) + sigma * (std_dev .* complex(randn(N,LL),randn(N,LL)))/sqrt(2);            
        otherwise,
            error(['unknown model ', model]);
    end
        
    % evaluate CBF and phase-only beampatterns for the all noisy snapshots
    
    DOA_CBF   = zeros(Number_of_DOAs,LL);
    DOA_phase = zeros(Number_of_DOAs,LL);
    DOA_mle   = zeros(Number_of_DOAs,LL);
    sidelobe_level_CBF   = zeros(1,LL);
    sidelobe_level_phase = zeros(1,LL);
    sidelobe_level_mle   = zeros(1,LL);
    
    MSE_x_CBF   = 0;
    MSE_x_phase = 0;
    MSE_x_mle   = 0;
    
    for ell=1:L:LL

        if Number_of_DOAs == 1,
            b_CBF   = sum(abs(A'*y(:,ell + [0:(L-1)])).^2,2); % multiple snapshot CBF
            [maxval,index]   = max(b_CBF);
            b_CBF            = b_CBF/maxval;
            DOA_CBF(1, ell)   = find_interpolated_max(phi_vec,b_CBF,index);
            a_CBF            = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_CBF(ell))); % offgrid CBF steering vector estimate
            x_CBF            = mean(a_CBF' * y(:,ell + [0:(L-1)])/N); % complex source amplitude estimate from CBF
            MSE_x_CBF        = MSE_x_CBF + abs(x_CBF - x_src)^2;
        else
            b_CBF            = NaN(M,L); % no conventional beamformer result for multiple sources
            DOA_CBF(1:Number_of_DOAs, ell) = NaN(Number_of_DOAs,1);
            MSE_x_CBF        = MSE_x_CBF + NaN;
        end
        
        if c_phase_version == 1
            c_phase = ones(1,L); % HERE WE SET c_phase equal to one, special check 2019-03-03 cfm
        else
            c_phase = 1.0 ./ sum(1.0 ./ abs(y(:,ell + [0:(L-1)])),1);   % Equation (21)
        end
        
        sign_y  = y(:,ell + [0:(L-1)]) ./ abs(y(:,ell + [0:(L-1)]));% Equation (22)
        c_phase_sign_y = repmat(c_phase,N,1) .* sign_y;             % Equation (21+22)
        
        if Number_of_DOAs == 1,
            b_phase = sum(abs(A'*c_phase_sign_y).^2,2); % multiple snapshot phase-only beamformer
            [maxval,index]   = max(b_phase);
            b_phase          = b_phase/maxval;
            DOA_phase(1,ell) = find_interpolated_max(phi_vec,b_phase,index);
            a_phase          = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_phase(ell))); % offgrid phase-only steering vector estimate
            x_phase          = mean(a_phase'*c_phase_sign_y);
            MSE_x_phase      = MSE_x_phase + abs(x_phase - x_src).^2;
        else
            b_phase          = NaN(M,L); % no phase-only beamformer result for multiple sources 
            DOA_phase(1:Number_of_DOAs,ell) = NaN;
            MSE_x_phase      = MSE_x_phase + NaN;
        end
        
        b_mle = zeros(Mcombinations,1);
        parfor m=1:Mcombinations
%             if rem(m,1000)==0
%                fprintf(1,'ell=%d, m=%d \n',ell,m);
%             end
            a = A(:,Active_Set(m,:)); % this is a matrix of size (N x Number_of_DOAs)
            x_hat1 = a \ c_phase_sign_y; % Eq. (19), also Eq. (38)
            x_hat2 = x_hat1;
            for ell2=0:(L-1)
               x_hat2(1:Number_of_DOAs, ell2+1) = L1LinearRegression( a, y(:,ell + ell2),x_hat1(:,ell2+1));
            end
            b_mle(m) = sum(sum(abs(y(:,ell + [0:(L-1)]) - a * x_hat2))); % Laplace-like Likelihood function in Equation (27)
        end
        b_mle = 1.0 ./ (b_mle .* b_mle);
        [maxval,index] = max(b_mle);
        %keyboard;
        b_mle        = b_mle/maxval;
        if Number_of_DOAs == 1,
           DOA_mle(1,ell) = find_interpolated_max(phi_vec,b_mle,index);
        else
           DOA_mle(1:Number_of_DOAs, ell) = phi_vec(Active_Set(index,:)).'; % for multiple sources, the interpolation is not implemented
        end
        a_mle        = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_mle(:,ell).')); % offgrid MLE steering vector estimate        
        x_hat1 = a_mle \ c_phase_sign_y; % Eq. (19), also Eq. (38)
        x_hat2 = x_hat1;
        for ell2=0:(L-1)
           x_hat2(1:Number_of_DOAs, ell2+1) = L1LinearRegression( a_mle, y(:,ell + ell2),x_hat1(:,ell2+1));
        end
        x_hat2 = mean(x_hat2);
        MSE_x_mle = MSE_x_mle + abs(x_hat2 - x_src).^2;
% keyboard;
        if Number_of_DOAs == 1,
            % evaluation of maximum sidelobe levels
            phi_tol = 360/(pi*N); % half beamwidth (from-maximum-to-first-zero) of CBF
            sidelobe_level_CBF(ell)   = maximum_sidelobe_level(10*log10(b_CBF),phi_vec, phi_vec(m_src),phi_tol);
            sidelobe_level_phase(ell) = maximum_sidelobe_level(10*log10(b_phase),phi_vec, phi_vec(m_src),phi_tol);
            sidelobe_level_mle(ell)   = maximum_sidelobe_level(20*log10(b_mle),phi_vec, phi_vec(m_src),phi_tol);
        else
            sidelobe_level_CBF(ell)   = NaN; % maximum_sidelobe_level(10*log10(b_CBF),phi_vec, phi_vec(m_src),phi_tol);
            sidelobe_level_phase(ell) = NaN; % maximum_sidelobe_level(10*log10(b_phase),phi_vec, phi_vec(m_src),phi_tol);
            sidelobe_level_mle(ell)   = NaN; % maximum_sidelobe_level(20*log10(b_mle),phi_vec, phi_vec(m_src),phi_tol);            
        end
      
        % for ell=1000, 2000,...
        % make a pretty plot of both beamformers with indicated sidelobe levels
        if rem(ell,100)==1
            fprintf('%4d ',ell);
            if rem(ell,1000)==1
                fprintf('\n');
                figure(1); clf;
                pretty_plot_all_beamformers;
                drawnow;
            end
        end
    end % end of parfor ell=1:L:LL
    
    mean_DOA_CBF(isnr)   = mean(DOA_CBF(1:L:LL));
    median_DOA_CBF(isnr) = median(DOA_CBF(1:L:LL));
    std_DOA_CBF(isnr)    = std(DOA_CBF(1:L:LL));
    mad_DOA_CBF(isnr)    = mad(DOA_CBF(1:L:LL),1);
    RMSE_DOA_CBF(isnr)   = sqrt(mean((DOA_CBF(1:L:LL)-phi_vec(m_src)).^2));
    median_sidelobe_level_CBF(isnr) = median(sidelobe_level_CBF(1:L:LL));
    mean_sidelobe_level_CBF(isnr)   = mean(sidelobe_level_CBF(1:L:LL));
    RMSE_x_CBF(isnr)     = sqrt(MSE_x_CBF/(LL/L));

    mean_DOA_phase(isnr)   = mean(DOA_phase(1:L:LL));
    median_DOA_phase(isnr) = median(DOA_phase(1:L:LL));
    std_DOA_phase(isnr)    = std(DOA_phase(1:L:LL));
    mad_DOA_phase(isnr)    = mad(DOA_phase(1:L:LL),1);
    RMSE_DOA_phase(isnr)   = sqrt(mean((DOA_phase(1:L:LL)-phi_vec(m_src)).^2));
    median_sidelobe_level_phase(isnr) = median(sidelobe_level_phase(1:L:LL));
    mean_sidelobe_level_phase(isnr)   = mean(sidelobe_level_phase(1:L:LL));
    RMSE_x_phase(isnr)     = sqrt(MSE_x_phase/(LL/L));

    mean_DOA_mle(isnr)   = mean(DOA_mle(1:L:LL));
    median_DOA_mle(isnr) = median(DOA_mle(1:L:LL));
    std_DOA_mle(isnr)    = std(DOA_mle(1:L:LL));
    mad_DOA_mle(isnr)    = mad(DOA_mle(1:L:LL),1);
    RMSE_DOA_mle(isnr)   = sqrt(mean((DOA_mle(1:L:LL)-phi_vec(m_src)).^2));
    median_sidelobe_level_mle(isnr) = median(sidelobe_level_mle(1:L:LL));
    mean_sidelobe_level_mle(isnr)   = mean(sidelobe_level_mle(1:L:LL));
    RMSE_x_mle(isnr)     = sqrt(MSE_x_mle/(LL/L));
    
    fprintf('\n%s noise model: SNR = %f dB\n',model,SNR(isnr));
    
    fprintf('mean(DOA_CBF)   = %f degrees\n', mean_DOA_CBF(isnr));
    fprintf('median(DOA_CBF) = %f degrees\n', median_DOA_CBF(isnr));
    fprintf('std(DOA_CBF)    = %f degrees\n', std_DOA_CBF(isnr));
    fprintf('mad(DOA_CBF,1)  = %f degrees\n', mad_DOA_CBF(isnr));
    fprintf('RMSE(DOA_CBF)   = %f degrees\n', RMSE_DOA_CBF(isnr));
    fprintf('median(sidelobe_level_CBF) = %f dB\n',  median_sidelobe_level_CBF(isnr));
    fprintf('mean(sidelobe_level_CBF)   = %f dB\n',mean_sidelobe_level_CBF(isnr));
    fprintf('RMSE x_CBF      = %f \n\n',RMSE_x_CBF(isnr));
    
    fprintf('mean(DOA_phase)   = %f degrees\n', mean_DOA_phase(isnr));
    fprintf('median(DOA_phase) = %f degrees\n', median_DOA_phase(isnr));
    fprintf('std(DOA_phase)    = %f degrees\n', std_DOA_phase(isnr));
    fprintf('mad(DOA_phase,1)  = %f degrees\n', mad_DOA_phase(isnr));
    fprintf('RMSE(DOA_phase)   = %f degrees\n', RMSE_DOA_phase(isnr));
    fprintf('median(sidelobe_level_phase) = %f dB\n',median_sidelobe_level_phase(isnr));
    fprintf('mean(sidelobe_level_phase)   = %f dB\n',mean_sidelobe_level_phase(isnr));
    fprintf('RMSE x_phase      = %f \n\n',RMSE_x_phase(isnr));

    fprintf('mean(DOA_mle)   = %f degrees\n', mean_DOA_mle(isnr));
    fprintf('median(DOA_mle) = %f degrees\n', median_DOA_mle(isnr));
    fprintf('std(DOA_mle)    = %f degrees\n', std_DOA_mle(isnr));
    fprintf('mad(DOA_mle,1)  = %f degrees\n', mad_DOA_mle(isnr));
    fprintf('RMSE(DOA_mle)   = %f degrees\n', RMSE_DOA_mle(isnr));
    fprintf('median(sidelobe_level_mle) = %f dB\n',median_sidelobe_level_mle(isnr));
    fprintf('mean(sidelobe_level_mle)   = %f dB\n',mean_sidelobe_level_mle(isnr));
    fprintf('RMSE x_hat2     = %f \n\n',RMSE_x_mle(isnr));

    save
    pretty_plot_doa_laplace_performance_versus_snr;
    elapsed_time = toc;
    fprintf('elapsed time = %f hrs, completed = %d percent\n---\n',round(elapsed_time/360)/10, round(100 * isnr/length(sigma_vec)));
end

