clear; clc;
close all;

dbstop if error;
% rng(777)

addpath([cd,'/_common'])

%% Noise model
model = 'Gaussian',                      nu_model_string = '';
% model = 'epscont', epsilon_model = 0.05, lambda_model  = 10.0, nu_model_string = 'epsilon=5e-2_lambda=10', 
%                 noise_enhancement = (1-epsilon_model) + epsilon_model*lambda_model^2,
% model = 'Complex-Student', nu_model = 2.1, nu_model_string = '2p1';

%% Environment parameters
freq   =  2.0E+03;    % frequency (Hz)
c0     =  343;        % speed of sound (m/s) in dry air at 20°C
lambda =  c0/freq;    % wavelength (m)
wavenum=  2*pi/lambda;% wave number (rad/m)

%% Array configuration
antenna_array.type = 'ULA';  % or 'UCA'
switch(antenna_array.type)
    case 'ULA' % definition of uniform linear array geometry
        theta  = 0.0;          % irrelevant elevation angle of incoming wave [degrees]
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
    case 'UCA' % definition of uniform circular array geometry
        theta  = 30.0;        % elevation angle of incoming wave [degrees]
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
            ]) / lambda; % sensor spacing of UCA measured in wavelengths
        % array steering matrix of size N x M for all azimuth at one single elevation theta (degrees)
        M = 361;
        dphi=360/(M-1);
        phi_vec = [-180:dphi:180]; % sweep over all azimuth angles (degrees)
end

% Design/steering matrix (Sensing matrix)
sensingMatrix = zeros(antenna_array.N,M);
for m=1:M
    kvec = wavenum * [sind(phi_vec(m))*cosd(theta);cosd(phi_vec(m))*cosd(theta);sind(theta)];
    sensingMatrix(:,m)   = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]); % normalization to |a_n|=1 or ||a||_2^2 = N.
end

%% Number of sensors / grid-points / snapshots
Nsensor     = antenna_array.N;  % number of sensors
Ntheta      = M;                % number of angular-search grid
Nsnapshot   = 25;               % number of snapshots

%% Simulation parameters
% noise standard deviation sigma

% SNRs = 30:-3:-9;
% SNR  = SNRs(6);
SNR = [14.9081800077602]; % dB

% number of sources
Number_of_DOAs = 2  ;
switch(Number_of_DOAs)
    case 1
        DOA_src = -45;          % -45 deg source on the grid
%         x_src   = exp(1j*pi/4); % source amplitude, magnitude = 1

        % YP: Based on 'generate_signal' function 'stochastic source'
        %     To generate Complex Gaussian source amplitude ~ CN(0,1), it
        %     should be one.
        x_src   = ones(Number_of_DOAs,1);
    case 2
        DOA_src = [-45; -25];
% %         x_src   = (10 .^([12; 20]/20)) * exp(1j*pi/4);
%         x_src   = (10 .^([20; 20]/20)) * exp(1j*pi/4); %-------------------------Equal strength
%         x_src   = x_src / norm(x_src);

        % YP: Based on 'generate_signal' function 'stochastic source'
        %     To generate Complex Gaussian source amplitude ~ CN(0,1), it
        %     should be one.
        x_src   = ones(Number_of_DOAs,1);
    case 3
        DOA_src = [-3; 2; 75];  % Gerstoft2015.pdf see FIGURE 8
%         x_src   = (10 .^([12; 22; 20]/20)) * exp(1j*pi/4);
%         x_src   = x_src / norm(x_src);

        % YP: Based on 'generate_signal' function 'stochastic source'
        %     To generate Complex Gaussian source amplitude ~ CN(0,1), it
        %     should be one.
        x_src   = ones(Number_of_DOAs,1);
    otherwise
        error('this Number_of_DOAs is not implemented')
end

% Steering vectors for true sources
for k=1:Number_of_DOAs
    m_src(k) = find(phi_vec == DOA_src(k));
end
a_src = sensingMatrix(:,m_src);

NmonteCarlo = 1;
LSnapshot = Nsnapshot * NmonteCarlo; % Number of array data vector observations "Large"


%% loop over various SNR levels
for isnr = 6 %1:length(SNRs)
    % for isnr = 1:length(SNRs)

    % Noise modeling
    sigma = 1 * norm(a_src*x_src,'fro') / (10^(SNR/20));
%     SNR_gen = 10*log10(norm(a_src*x_src,'fro')^2 ./ (sigma.^2));
%     % check the generated SNR

    switch(model)
        case 'Laplace-like'
            y = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                            sigma,model);
        case 'Gaussian'
            y = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                            sigma,model);
        case 'epscont'
            y = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                            sigma,model,epsilon_model,lambda_model);
        case 'Complex-Student' % method in Ollila & Koivunen, PIMRC 2003
            y = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                            sigma,model,nu_model);
        case 'Heteroscedastic'
            y = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                            sigma,model);
        otherwise
            error(['unknown model ', model]);
    end

% evaluate SBL
    options = SBLSet;
    options.Nsource = Number_of_DOAs;
%     options.gamma_range=10^-3;

    for n_monteCarlo = 1:NmonteCarlo  % parfor loop over snapshot realizations
        Y = y(:,(n_monteCarlo-1)*Nsnapshot+(1:Nsnapshot));

%% SBL OLD
        % SBL old v1
        % Same as SBL_v4 (Not sure the version date) w/o Esa update
        tic;
        [gamma,report] = SBL_pre_Peter(sensingMatrix, Y, options);
        toc;

        % SBL clean-up v1
        % Almost the same as SBL_v4p00 but clean version
        tic;
        [gamma1,report1] = SBL_v4p11(sensingMatrix, Y, options);
        toc;

        % SBL clean-up v4
        % Esa's initialization
        % SBL_v4p11 Gamma update
        % SBL_v4p12 Sigma Initialization
        tic;
        [gamma2,report2] = SBL_v4p14(sensingMatrix, Y, options);
        toc;

%% Gaussian
        % Same as SBL_v4 w/ Esa update
        % Same as SBL_esa.m w/ convergence error 10^(-3)
        tic;
        disp(['SBL Esa version initialized.']);
        [gamma3,report3] = SBL_pre_Esa(sensingMatrix, Y, 'SBL-G', -90:dphi:90, Number_of_DOAs);
        toc;

        % SBL_v4p11 Robust version
        tic;
        loss_param = inf;
        [gamma30,report30] = SBL_v5p01(sensingMatrix, Y, 'SBL-G', loss_param, options);
        toc;

        % SBL_v4p14 Robust version
        tic;
        loss_param = inf;
        [gamma33,report33] = SBL_v5p04(sensingMatrix, Y, 'SBL-G', loss_param, options);
        toc;

        % SBL_v4p14 Robust version [Current version]
        tic;
        loss_param = inf;
        [gamma34,report34] = SBL_v5p05(sensingMatrix, Y, 'SBL-G', loss_param, options);
        toc;        

%% MVT        
        tic;
        disp(['SBL Esa version initialized.']);
        [gamma4,report4] = SBL_pre_Esa(sensingMatrix, Y, 'SBL-T', -90:dphi:90, Number_of_DOAs);
        toc;

        % Robust SBL used in TSP paper
        % Same as SBL_v4_tM_esa.m w/ convergence error 10^(-3)
        tic;
        disp(['SBL TSP SBL4-T initialized.']);
        nu_algorithm = 2.1;
        [gamma5,~,report5]  = SBL_pre_tsp_MVT(sensingMatrix, Y, nu_algorithm, Number_of_DOAs);
        toc;

        % SBL_v4p11 Robust version
        tic;
        loss_param = 2.1;
        [gamma10,report10] = SBL_v5p01(sensingMatrix, Y, 'SBL-T', loss_param, options);
        toc;

        % SBL_v4p14 Robust version
        tic;
        loss_param = 2.1;
        [gamma13,report13] = SBL_v5p04(sensingMatrix, Y, 'SBL-T', loss_param, options);
        toc;

        % SBL_v4p14 Robust version [Current version]
        tic;
        loss_param = 2.1;
        [gamma14,report14] = SBL_v5p05(sensingMatrix, Y, 'SBL-T', loss_param, options);
        toc;

%% Huber
        tic;
        disp(['SBL Esa version initialized.']);
        [gamma6,report6] = SBL_pre_Esa(sensingMatrix, Y, 'SBL-H', -90:dphi:90, Number_of_DOAs);
        toc;

        % Robust SBL used in ICASSP paper
        % Same as SBL_v4_M_icassp.m w/ convergence error 10^(-3)
        tic;
        disp(['SBL TSP SBL4-H initialized.']);
        nu_algorithm = 0.9;
        [gamma7,~,report7] = SBL_pre_icassp_Huber(sensingMatrix , Y , 'Huber', nu_algorithm, Number_of_DOAs);
        toc;

        % SBL_v4p11 Robust version
        tic;
        loss_param = 0.85;
        [gamma20,report20] = SBL_v5p01(sensingMatrix, Y, 'SBL-H', loss_param, options);
        toc;

        % SBL_v4p14 Robust version
        tic;
        loss_param = 0.85;
        [gamma23,report23] = SBL_v5p04(sensingMatrix, Y, 'SBL-H', loss_param, options);
        toc;

        % SBL_v4p14 Robust version [Current version]
        tic;
        loss_param = 0.85;
        [gamma24,report24] = SBL_v5p05(sensingMatrix, Y, 'SBL-H', loss_param, options);
        toc;


%% Results
        disp([' '])
        disp(['OLD models'])

        [~, Ilocs] = findpeaks(abs(gamma),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv4p00 : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma1),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv4p11 : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma2),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv4p14 : ',num2str(sqrt(mean(power(DoA_error,2))))])


        disp([' '])
        disp(['Gaussian models (-G)'])

        [~, Ilocs] = findpeaks(abs(gamma3),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv4_Esa: ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma30),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p01 : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma33),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p04 : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma34),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p05 : ',num2str(sqrt(mean(power(DoA_error,2))))])


        disp([' '])
        disp(['MVT-loss models (-T)'])

        [~, Ilocs] = findpeaks(abs(gamma4),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv4_Esa: ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma5),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBL4-T   : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma10),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p01 : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma13),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p04 : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma14),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p05 : ',num2str(sqrt(mean(power(DoA_error,2))))])


        disp([' '])
        disp(['Huber-loss models (-H)'])

        [~, Ilocs] = findpeaks(abs(gamma6),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv4_Esa: ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma7),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBL4-H   : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma20),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p01 : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma23),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p04 : ',num2str(sqrt(mean(power(DoA_error,2))))])

        [~, Ilocs] = findpeaks(abs(gamma24),'SORTSTR','descend','Npeaks',Number_of_DOAs);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p05 : ',num2str(sqrt(mean(power(DoA_error,2))))])

    end % end of the for-loop

end % end of for isnr=1:length(sigma_vec) loop


%%
% figure,
% plot(phi_vec,gamma1/max(gamma1),'--')
% hold on;
% plot(phi_vec,gamma2/max(gamma2),'--')
% 
% plot(phi_vec,gamma4/max(gamma4))
% plot(phi_vec,gamma5/max(gamma5))
% plot(phi_vec,gamma10/max(gamma10))
% plot(phi_vec,gamma13/max(gamma13))
% 
% plot(phi_vec,gamma6/max(gamma6),'-.')
% plot(phi_vec,gamma7/max(gamma7),'-.')
% plot(phi_vec,gamma20/max(gamma20),'-.')
% plot(phi_vec,gamma23/max(gamma23),'-.')
% hold off;

%%
rmpath([cd,'/_common'])
%% End------------------------------------------------------------------------------------------------------------------------ %%


%% Signal generation
function receivedSignal = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
    sigma,model,model_param1,model_param2)
% function to generate sensor observations
if strcmpi(model,'Laplace-like')
    if 0
        % deterministric source
        receivedSignal = laplacelike_rand(a_src*x_src, sigma, Nsensor, LSnapshot);
    else
        % stochastic source
        error('laplacelike_rand for stochastic source not yet implemented')
    end
elseif (strcmpi(model,'Gaussian') || isempty(model) )
    noise_realization = sigma * complex(randn(Nsensor,LSnapshot),randn(Nsensor,LSnapshot))/sqrt(2);
    if 0
        % deterministric source
        receivedSignal = (a_src * x_src * ones(1,LSnapshot)) + noise_realization;
    else
        % stochastic source
        s_src = x_src .* complex(randn(Number_of_DOAs,LSnapshot),randn(Number_of_DOAs,LSnapshot))/sqrt(2 * Number_of_DOAs);
        receivedSignal = ( a_src * s_src ) + noise_realization;
    end
elseif (strcmpi(model,'epscont') || isempty(model) )
    if nargin < 9
        epsilon_model = 0.05; lambda_model  = 10.0;
    else
        epsilon_model = model_param1; lambda_model  = model_param2;
    end
    noise_realization = epscont(Nsensor,LSnapshot,sigma,epsilon_model, 0.0,lambda_model);
    %old: noise_realization = epscont_old(Nsensor * LSnapshot, epsilon_model, 0.0, lambda_model, sigma);
    %old: noise_realization = reshape(noise_realization, Nsensor, LSnapshot);
    if 0
        % deterministric source
        receivedSignal = (a_src * x_src * ones(1,LSnapshot)) + noise_realization;
    else
        % stochastic source
        s_src = x_src .* complex(randn(Number_of_DOAs,LSnapshot),randn(Number_of_DOAs,LSnapshot))/sqrt(2 * Number_of_DOAs);
        receivedSignal = ( a_src * s_src ) + noise_realization;
    end
elseif (strcmpi(model,'Complex-Student') || isempty(model) ) % method in Ollila & Koivunen, PIMRC 2003
    if nargin < 8
        nu_model = 2.1;
    else
        nu_model = model_param1;
    end
    noise_realization = sigma * complex(randn(Nsensor,LSnapshot),randn(Nsensor,LSnapshot))/sqrt(2);
    if 0
        % deterministric source
        receivedSignal = (a_src * x_src * ones(1,LSnapshot)) + noise_realization;
    else
        % stochastic source
        s_src = x_src .* complex(randn(Number_of_DOAs,LSnapshot),randn(Number_of_DOAs,LSnapshot))/sqrt(2 * Number_of_DOAs);
        receivedSignal = ( a_src * s_src ) + noise_realization;
    end
    % nu_model = 5, % choose number of degrees of freedom of chi^2 distribution
    s = ones(Nsensor, 1) * chi2rnd(nu_model * ones(1, LSnapshot));
    receivedSignal =  receivedSignal ./ sqrt(s/nu_model);
elseif (strcmpi(model,'Heteroscedastic') || isempty(model) )
    noise_realization = sigma * complex(randn(Nsensor,LSnapshot),randn(Nsensor,LSnapshot))/sqrt(2);
    std_dev = 10.^(-1.0+2.0*rand(Nsensor,LSnapshot));
    std_dev = std_dev/sqrt(sum(sum(std_dev.^2))/(Nsensor*LSnapshot));
    if 0
        % deterministric source
        receivedSignal = (a_src * x_src * ones(1,LSnapshot)) + std_dev .* noise_realization;
    else
        % stochastic source
        s_src = x_src .* complex(randn(Number_of_DOAs,LSnapshot),randn(Number_of_DOAs,LSnapshot))/sqrt(2 * Number_of_DOAs);
        receivedSignal = ( a_src * s_src ) + std_dev .* noise_realization;
    end
else
    error(['please specify noise model as a string equal to Laplace-like, ...' ...
        'Gaussian, epscont, Complex-Student or Heteroscedastic\n']);
end
end
