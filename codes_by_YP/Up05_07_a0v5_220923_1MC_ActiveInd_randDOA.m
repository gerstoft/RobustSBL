clear; clc;
close all;

dbstop if error;
% rng(777)

addpath([cd,'/_common'])

% SNRs = 30:-3:-9;
% NmonteCarlo = 100;
% Number_of_DOAs = 2  ;
% function errorDOAcutoff, threshold = 10 [deg.]
% save results: DOAs & errors; struct('theta',theta(Ilocs),'error',DoA_error);
% findpeaks for error calculation, 'Npeaks', Number_of_DOAs "+ 2", minPeakSeparation: 5[deg.] (=5/dphi points)

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
%         M = 361;   % standard dictionary, .5 deg azimuth resolution
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
sensingMatrixD = zeros(antenna_array.N,M); %-- CRB-YP
for m=1:M
    kvec = wavenum * [sind(phi_vec(m))*cosd(theta);cosd(phi_vec(m))*cosd(theta);sind(theta)];
    kvecD = wavenum * [cosd(phi_vec(m))*cosd(theta);-sind(phi_vec(m))*cosd(theta);sind(theta)];
    sensingMatrix(:,m)   = exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]); % normalization to |a_n|=1 or ||a||_2^2 = N.
    sensingMatrixD(:,m)  = (-1j * kvecD.' * [antenna_array.x;antenna_array.y;antenna_array.z])...
        .* exp(-1j * kvec.' * [antenna_array.x;antenna_array.y;antenna_array.z]); % normalization to |a_n|=1 or ||a||_2^2 = N.
end

%% Number of sensors / grid-points / snapshots
Nsensor     = antenna_array.N;  % number of sensors
Ntheta      = M;                % number of angular-search grid
Nsnapshot   = 25;               % number of snapshots

%% Simulation parameters
% noise standard deviation sigma

SNRs = 36:-3:-6;
% SNRs = 30:-3:-9;
% SNRs = [50;45;40];
% SNR  = SNRs(6);
% SNR = [14.9081800077602]; % dB

Number_of_DOAs = 3;

NmonteCarlo = 100;
LSnapshot = Nsnapshot;
% LSnapshot = Nsnapshot * NmonteCarlo; % Number of array data vector observations "Large"


%% loop over various SNR levels

for isnr = 6 %1:length(SNRs)
% for isnr = 1:length(SNRs)

%     rng('default') % YP: We need to be careful where to set 'rng'
    rng(1,'twister')

    SNR  = SNRs(isnr);

% evaluate SBL
    options = SBLSet;
    options.Nsource = Number_of_DOAs;
%     options.gamma_range=10^-3;

    errorDOAseparation = 1; % [deg.]
    errorDOAsepP = floor(errorDOAseparation/dphi) - 1;
    errorDOApeak = Number_of_DOAs + 2;
    errCut = 10; % Maximum RMSE cut-off.

% obtain active indices --------------------------%
    options.activeIndices = 1;
    options.activeIndRepN = 20; % for SBL_v5p07
%-------------------------------------------------%

    for n_monteCarlo = 1
%     for n_monteCarlo = 1:NmonteCarlo  % parfor loop over snapshot realizations

        disp(' ')
        disp(['SNR',num2str(SNR),'#Sim : ',num2str(n_monteCarlo)])

        % number of sources
        x_src   = ones(Number_of_DOAs,1);
        DOA_src = floor(asind(gen_bs(sind(-75), sind(75), Number_of_DOAs, asin(2/Nsensor))));
        DOA_MC(:,n_monteCarlo) = DOA_src;

        % Steering vectors for true sources
        for k=1:Number_of_DOAs
            m_src(k) = find(phi_vec == DOA_src(k));
        end
        a_src = sensingMatrix(:,m_src);

        % Noise modeling
        sigma = 1 * norm(a_src*x_src,'fro') / (10^(SNR/20));
        %     SNR_gen = 10*log10(norm(a_src*x_src,'fro')^2 ./ (sigma.^2));
        %     % check the generated SNR

        switch(model)
            case 'Laplace-like'
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model);
            case 'Gaussian'
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model);
            case 'epscont'
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model,epsilon_model,lambda_model);
            case 'Complex-Student' % method in Ollila & Koivunen, PIMRC 2003
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model,nu_model);
            case 'Heteroscedastic'
                [y,xAmp] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
                    sigma,model);
            otherwise
                error(['unknown model ', model]);
        end

        Y = y;
%         Y = y(:,(n_monteCarlo-1)*Nsnapshot+(1:Nsnapshot));


%% Gaussian
        % SBL_v5p05 Robust version [Current version]
%         tic;
        loss_param = inf;
        [gamma51,report51] = SBL_v5p05(sensingMatrix, Y, 'SBL-G', loss_param, options);
        [gammaInd71,report71] = SBL_v5p07(sensingMatrix, Y, 'SBL-G', loss_param, options, errorDOApeak, errorDOAsepP);
%         toc;
        
% obtain active indices plot ---------------------%
        for itr = 1:2
        if itr==1, actIndvar   = report51; else, actIndvar   = report71; end
        actInd      = reshape(actIndvar.results.activeIndices.',...
            size(actIndvar.results.activeIndices,1)*size(actIndvar.results.activeIndices,2),1);
        actIndItr   = reshape(repmat(1:size(actIndvar.results.activeIndices,1),size(actIndvar.results.activeIndices,2),1),...
            size(actIndvar.results.activeIndices,1)*size(actIndvar.results.activeIndices,2),1);
        if itr==1, figure, set(gcf,'Position',[1,50,560,420]); else, figure, set(gcf,'Position',[1,525,560,420]); end
        scatter(actIndItr,actInd,30,'k*','linewidth',1.0)
        hold on; plot(actIndvar.results.iteration_L1*ones(1,2),[0 length(actIndvar.results.final_iteration.gamma)+1],'r'); hold off;
        set(gca,'fontsize',24,'TickLabelInterpreter','latex')
        xlabel('Iteration','interpreter','latex');
        ylabel('Index No.','interpreter','latex');
        title(['Active indices, SBL ',actIndvar.options.SBL_v])
        axis([0 size(actIndvar.results.activeIndices,1) 0 length(actIndvar.results.final_iteration.gamma)])
        end

%% MVT        
        % SBL_v5p05 Robust version [Current version]
%         tic;
        loss_param = 2.1;
        [gamma52,report52] = SBL_v5p05(sensingMatrix, Y, 'SBL-T', loss_param, options);
        [gammaInd72,report72] = SBL_v5p07(sensingMatrix, Y, 'SBL-T', loss_param, options, errorDOApeak, errorDOAsepP);
%         toc;

% obtain active indices plot ---------------------%
        for itr = 1:2
        if itr==1, actIndvar   = report52; else, actIndvar   = report72; end
        actInd      = reshape(actIndvar.results.activeIndices.',...
            size(actIndvar.results.activeIndices,1)*size(actIndvar.results.activeIndices,2),1);
        actIndItr   = reshape(repmat(1:size(actIndvar.results.activeIndices,1),size(actIndvar.results.activeIndices,2),1),...
            size(actIndvar.results.activeIndices,1)*size(actIndvar.results.activeIndices,2),1);
        if itr==1, figure, set(gcf,'Position',[560,50,560,420]); else, figure, set(gcf,'Position',[560,525,560,420]); end
        scatter(actIndItr,actInd,30,'k*','linewidth',1.0)
        hold on; plot(actIndvar.results.iteration_L1*ones(1,2),[0 length(actIndvar.results.final_iteration.gamma)+1],'r'); hold off;
        set(gca,'fontsize',24,'TickLabelInterpreter','latex')
        xlabel('Iteration','interpreter','latex');
        ylabel('Index No.','interpreter','latex');
        title(['Active indices, SBL ',actIndvar.options.SBL_v])
        axis([0 size(actIndvar.results.activeIndices,1) 0 length(actIndvar.results.final_iteration.gamma)])
        end

%% Huber
        % SBL_v5p05 Robust version [Current version]
%         tic;
        loss_param = 0.9;
        [gamma53,report53] = SBL_v5p05(sensingMatrix, Y, 'SBL-H', loss_param, options);
        [gammaInd73,report73] = SBL_v5p07(sensingMatrix, Y, 'SBL-H', loss_param, options, errorDOApeak, errorDOAsepP);
%         toc;

% obtain active indices plot ---------------------%
        for itr = 1:2
        if itr==1, actIndvar   = report53; else, actIndvar   = report73; end
        actInd      = reshape(actIndvar.results.activeIndices.',...
            size(actIndvar.results.activeIndices,1)*size(actIndvar.results.activeIndices,2),1);
        actIndItr   = reshape(repmat(1:size(actIndvar.results.activeIndices,1),size(actIndvar.results.activeIndices,2),1),...
            size(actIndvar.results.activeIndices,1)*size(actIndvar.results.activeIndices,2),1);
        if itr==1, figure, set(gcf,'Position',[1120,50,560,420]); else, figure, set(gcf,'Position',[1120,525,560,420]); end
        scatter(actIndItr,actInd,30,'k*','linewidth',1.0)
        hold on; plot(actIndvar.results.iteration_L1*ones(1,2),[0 length(actIndvar.results.final_iteration.gamma)+1],'r'); hold off;
        set(gca,'fontsize',24,'TickLabelInterpreter','latex')
        xlabel('Iteration','interpreter','latex');
        ylabel('Index No.','interpreter','latex');
        title(['Active indices, SBL ',actIndvar.options.SBL_v])
        axis([0 size(actIndvar.results.activeIndices,1) 0 length(actIndvar.results.final_iteration.gamma)])
        end


%% Results

        disp([' '])
        disp(['Gaussian models (-G)'])

        [~, Ilocs] = findpeaks(abs(gamma51),'MinPeakDistance',errorDOAsepP,'SORTSTR','descend','Npeaks',errorDOApeak);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,errCut);
        disp(['RMSE SBLv5p05 : ',num2str(sqrt(mean(power(DoA_error,2))))])
%         [~, Ilocs] = findpeaks(abs(gamma71),'MinPeakDistance',errorDOAsepP,'SORTSTR','descend','Npeaks',errorDOApeak);
        DoA_error = errorDOAcutoff(phi_vec(gammaInd71),DOA_src,errCut);
        disp(['RMSE SBLv5p07 : ',num2str(sqrt(mean(power(DoA_error,2))))])


        disp([' '])
        disp(['MVT-loss models (-T)'])

        [~, Ilocs] = findpeaks(abs(gamma52),'MinPeakDistance',errorDOAsepP,'SORTSTR','descend','Npeaks',errorDOApeak);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,errCut);
        disp(['RMSE SBLv5p05 : ',num2str(sqrt(mean(power(DoA_error,2))))])
%         [~, Ilocs] = findpeaks(abs(gamma72),'MinPeakDistance',errorDOAsepP,'SORTSTR','descend','Npeaks',errorDOApeak);
        DoA_error = errorDOAcutoff(phi_vec(gammaInd72),DOA_src,errCut);
        disp(['RMSE SBLv5p07 : ',num2str(sqrt(mean(power(DoA_error,2))))])


        disp([' '])
        disp(['Huber-loss models (-H)'])

        [~, Ilocs] = findpeaks(abs(gamma53),'MinPeakDistance',errorDOAsepP,'SORTSTR','descend','Npeaks',errorDOApeak);
        DoA_error = errorDOAcutoff(phi_vec(Ilocs),DOA_src,10);
        disp(['RMSE SBLv5p05 : ',num2str(sqrt(mean(power(DoA_error,2))))])
%         [~, Ilocs] = findpeaks(abs(gamma73),'MinPeakDistance',errorDOAsepP,'SORTSTR','descend','Npeaks',errorDOApeak);
        DoA_error = errorDOAcutoff(phi_vec(gammaInd73),DOA_src,10);
        disp(['RMSE SBLv5p07 : ',num2str(sqrt(mean(power(DoA_error,2))))])

    end % end of the for-loop
%     saveCharVar = who('outputs*');
%     saveChar = ['save([ model(1), ''mode_'', ''s'', num2str(Number_of_DOAs), ''MC'' , num2str(NmonteCarlo) , ''SNRn'' , num2str(isnr) ], ''SNRs'' , ''NmonteCarlo'' '];
%     for ichar = 1:numel(saveCharVar)
%         saveChar = [saveChar,',''',char(saveCharVar{ichar}),''''];
%     end
%     saveChar = [saveChar,');'];
%     eval(saveChar)

end % end of for isnr=1:length(sigma_vec) loop


%% Figure


%%
rmpath([cd,'/_common'])
%% End------------------------------------------------------------------------------------------------------------------------ %%


%% Signal generation
function [receivedSignal,s_src] = generate_signal(a_src,x_src,Nsensor,LSnapshot,Number_of_DOAs,...
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
