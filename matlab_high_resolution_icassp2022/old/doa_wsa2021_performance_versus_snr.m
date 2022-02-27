close all
clear all

% model = 'Laplace-like'
model = 'Gaussian'
% model = 'Heteroscedastic'   % Case III

freq   =  2.0E+03,    % frequency (Hz)

c0     = 343;         % speed of sound (m/s) in dry air at 20°C
lambda = c0/freq;     % wavelength (m)
wavenum= 2*pi/lambda; % wave number (rad/m)

antenna_array.type = 'UCA',  % or 'ULA'

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
        M = 181;
        dphi=180/(M-1);
        phi_vec = [-90:dphi:90];
end

for m=1:M
    kvec = wavenum * [sind(phi_vec(m))*cosd(theta);cosd(phi_vec(m))*cosd(theta);sind(theta)];
    A(:,m)   = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]); % normalization to |a_n|=1 or ||a||_2^2 = N.
end


c_phase_version = 2, % select 1 (c_phase=1) or 2 (c_phase depends on data magnitude)

MMV = 0;

L = 1,   % number of snapshots used by DOA estimator


% noise standard deviation sigma
if 0
    % choose coarse steps for quicker simulation
    sigma_vec = sqrt(antenna_array.N) * [0.125 0.250 0.500 0.600 0.707 1.000 2.000];
else
    % choose fine steps for high accuracy plots
    sigma_vec = sqrt(antenna_array.N) * [0.0625 0.088 0.125 0.177 0.250 0.354 0.500 0.550 0.600 0.650 0.707 0.750 1.000 1.414 2.000];
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

LL = 1e4, % number of array data vector observations

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
            y = laplacelike_rand(a_src*x_src, sigma, antenna_array.N, LL);
        case 'Gaussian',
            y = (a_src * x_src * ones(1,LL)) + sigma * complex(randn(antenna_array.N,LL),randn(antenna_array.N,LL))/sqrt(2);
        case 'Heteroscedastic',
            std_dev = 10.^(-1.0+2.0*rand(antenna_array.N,LL));
            std_dev = std_dev/sqrt(sum(sum(std_dev.^2))/(antenna_array.N*LL));
            y = (a_src*ones(1,LL)) * x_src + sigma * (std_dev .* complex(randn(antenna_array.N,LL),randn(antenna_array.N,LL)))/sqrt(2);            
        otherwise,
            error(['unknown model ', model]);
    end
        
    % evaluate CBF and phase-only beampatterns for the all noisy snapshots
    
    DOA_CBF = zeros(Number_of_DOAs,LL);
    DOA_phase = zeros(1,LL);
    DOA_mle   = zeros(1,LL);
    sidelobe_level_CBF   = zeros(1,LL);
    sidelobe_level_phase = zeros(1,LL);
    sidelobe_level_mle   = zeros(1,LL);
    
    MSE_x_CBF   = 0;
    MSE_x_phase = 0;
    MSE_x_mle   = 0;
    
    
    parfor ell=1:L:LL % parfor!!!
        
        if Number_of_DOAs == 1,
            b_CBF   = sum(abs(A'*y(:,ell + [0:(L-1)])).^2,2); % multiple snapshot CBF
            [maxval,index]   = max(b_CBF);
            b_CBF            = b_CBF/maxval;
            DOA_CBF(ell)   =  find_interpolated_max(phi_vec,b_CBF,index);  % CBF offgrid DOA estimate: azimuth in degrees
            
            kvec  = wavenum * [sind(DOA_CBF(ell))*cosd(theta);cosd(DOA_CBF(ell))*cosd(theta);sind(theta)];
            a_CBF = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]).';  % offgrid CBF steering vector estimate
            x_CBF     = mean(a_CBF' * y(:,ell + [0:(L-1)])/antenna_array.N); % complex source amplitude estimate from CBF
            MSE_x_CBF = MSE_x_CBF + abs(x_CBF - x_src)^2;    
        else
            b_CBF            = NaN(M,L); % no conventional beamformer result for multiple sources
            % DOA_CBF(:,ell)   = NaN(Number_of_DOAs,1); % parfor doesn't like this
            DOA_CBF(ell)   = NaN; % this is buggy
            MSE_x_CBF        = MSE_x_CBF + NaN;
        end    
        
        if c_phase_version == 1
            c_phase = ones(1,L); % HERE WE SET c_phase equal to one, special check 2019-03-03 cfm
        else
            c_phase = 1.0 ./ sum(1.0 ./ abs(y(:,ell + [0:(L-1)])),1);   % Equation (21)
        end
        
        sign_y  = y(:,ell + [0:(L-1)]) ./ abs(y(:,ell + [0:(L-1)]));% Equation (22)
        c_phase_sign_y = repmat(c_phase,antenna_array.N,1) .* sign_y;             % Equation (21+22)
        
        b_phase = sum(abs(A'*sign_y).^2,2); % multiple snapshot phase-only beamformer (23) without c_phase(!)
        [maxval,index] = max(b_phase);
        b_phase        = b_phase/maxval;
        DOA_phase(ell) = find_interpolated_max(phi_vec,b_phase,index); % phase-only offgrid DOA estimate: azimuth in degrees
        
        kvec           = wavenum * [sind(DOA_phase(ell))*cosd(theta);cosd(DOA_phase(ell))*cosd(theta);sind(theta)];
        a_phase        = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]).'; % offgrid phase-only steering vector estimate 
        x_phase        = mean(a_phase'*c_phase_sign_y);
        MSE_x_phase    = MSE_x_phase + abs(x_phase - x_src).^2;
        
        b_mle = zeros(size(b_CBF));
        for m=1:M
            a = A(:,m);
            x_hat1 = a' * c_phase_sign_y; % Eq. (19), also Eq. (38)
            for ell2=0:(L-1)
               x_hat2(ell2+1) = L1LinearRegression( a, y(:,ell + ell2),x_hat1(ell2+1));
            end
            b_mle(m) = sum(sum(abs(y(:,ell + [0:(L-1)]) - a * x_hat2))); % Laplace-like Likelihood function in Equation (27)
        end
        b_mle = 1.0 ./ (b_mle .* b_mle);
        [maxval,index] = max(b_mle);
        b_mle        = b_mle/maxval;
        DOA_mle(ell) = find_interpolated_max(phi_vec,b_mle,index); % MLE offgrid DOA estimate: azimuth in degrees
        
        kvec         = wavenum * [sind(DOA_mle(ell))*cosd(theta);cosd(DOA_mle(ell))*cosd(theta);sind(theta)];
        a_mle        = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]).'; % offgrid MLE steering vector estimate  
        
        a = conj(a_mle);
        x_hat1 = a.' * c_phase_sign_y; % Eq. (19), also Eq. (38)

        for ell2=0:(L-1)
           x_hat2(ell2+1) = L1LinearRegression( a_mle, y(:,ell + ell2),x_hat1(ell2+1));
        end
        x_hat2 = mean(x_hat2);
        MSE_x_mle = MSE_x_mle + abs(x_hat2 - x_src).^2;
% keyboard;
        % evaluation of maximum sidelobe levels
        phi_tol = 360/(pi*antenna_array.N); % ULA half beamwidth (from-maximum-to-first-zero) of CBF
        sidelobe_level_CBF(ell)   = maximum_sidelobe_level(10*log10(b_CBF),phi_vec, phi_vec(m_src),phi_tol);
        sidelobe_level_phase(ell) = maximum_sidelobe_level(10*log10(b_phase),phi_vec, phi_vec(m_src),phi_tol);
        sidelobe_level_mle(ell)   = maximum_sidelobe_level(20*log10(b_mle),phi_vec, phi_vec(m_src),phi_tol);
        
      
        % for ell=1000, 2000,...
        % make a pretty plot of both beamformers with indicated sidelobe levels
%         if rem(ell,100)==1
%             fprintf('%4d ',ell);
%             if rem(ell,1000)==1
%                 fprintf('\n');
%                 figure(1); clf;
%                 pretty_plot_all_beamformers;
%                 drawnow;
%             end
%         end
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
    pretty_plot_doa_wsa2021_performance_versus_snr;
    elapsed_time = toc;
    elapsed_time_hrs = round(elapsed_time/360)/10;
    percentage_completed = round(100 * isnr/length(sigma_vec));
    remaining_time_hrs = (1.0 - isnr/length(sigma_vec)) * (elapsed_time/3600) / (isnr/length(sigma_vec));
    fprintf('elapsed time = %f hrs, completed = %d percent\n',elapsed_time_hrs, percentage_completed);
    fprintf('estimated remaining time = %f hrs\n', remaining_time_hrs);
    fprintf('predicted end of job = %s\n---\n', datestr(now + remaining_time_hrs/24));
end

