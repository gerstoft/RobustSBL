close all
clear all
clc

%define macro for comparing angles
angDiff = @(a1,a2) (180 - abs(abs(a1 - a2) - 180));


rng(0); % random number seed

% model = 'Laplace-like';
% model = 'Gaussian';
model = 'GaussianOutlier';
% model = 'Heteroscedastic';   % Case III

N       = 20;  % no. of sensors in ULA
d       = 0.5; % sensor spacing of ULA measured in wavelengths

L       = 10;   % number of snapshots used by DOA estimator

% noise standard deviation sigma

sigma   = sqrt(N) * 0.500;


% array steering matrix of size N x M
M       = 181;
dtheta  = 180/(M-1);
phi_vec = [-90:dtheta:90];

for m=1:M
    A(:,m)   = exp(-1j*2*pi*d*[0:N-1]'*sind(phi_vec(m))); % same normalization to |a_n|=1 as in the LaTeX code
end

% m_src = (M+1)/2; % broadside source on the grid
m_src   = (M-1)/4+1; % 45 deg source on the grid
a_src   = A(:,m_src);
x_src   = exp(1j*pi/4); % source amplitude, magnitude = 1

SNR     = 10*log10(norm(a_src*x_src,'fro')^2 ./ (sigma.^2));

% LL = 5e4; % number of array data vector observations
% LL      = 100000; % number of array data vector observations

LL      = 10000; % number of array data vector observations

outlierRadVec   = [0:5:25, 30:10:100];


for oCounter = 1:length(outlierRadVec)
    
    rng(0); % reset random numbers -> trade off bias vs variance
    
    outlierRad = outlierRadVec(oCounter);
        
    disp(['outlier ' num2str(oCounter) ' of ' num2str(length(outlierRadVec))])
    
    switch model
        case 'Laplace-like'
            y       = laplacelike_rand(a_src*x_src, sigma, N, LL);
        case 'Gaussian'
            y       = (a_src*ones(1,LL)) * x_src + sigma * complex(randn(N,LL),randn(N,LL))/sqrt(2);
        case 'GaussianOutlier'
            ySig        = (a_src*ones(1,LL)) * x_src;
            
            noiseVec    = sigma * complex(randn(N,LL),randn(N,LL))/sqrt(2);
            
            %one outlier every Lth observation
            outlierPos  = 1:L:LL;
            cleanPos    = setdiff(1:LL,1:L:LL);
            
            outliersDir     = randsphere(length(outlierPos), N).';  %random unit point-mass
            outliers        = outlierRad*repmat(exp(1j*2*pi*rand(1,size(outliersDir,2))),size(outliersDir,1),1).*outliersDir;  %add complex weight
               
            yO(:,outlierPos)   = ySig(:,outlierPos) +  outliers;
            yO(:,cleanPos)     = ySig(:,cleanPos)   + noiseVec(:,cleanPos);
            
            %reference
            y       = (a_src*ones(1,LL)) * x_src + noiseVec;
            
        case 'Heteroscedastic'
            std_dev = 10.^(-1.0+2.0*rand(N,LL));
            std_dev = std_dev/sqrt(sum(sum(std_dev.^2))/(N*LL));
            y       = (a_src*ones(1,LL)) * x_src + sigma * (std_dev .* complex(randn(N,LL),randn(N,LL)))/sqrt(2);
        otherwise
            error(['unknown model ', model]);
    end
    
    % evaluate CBF and phase-only beampatterns for the all noisy snapshots
    
    DOA_CBF         = zeros(1,LL);
    DOA_phase       = zeros(1,LL);
    DOA_phase_c     = zeros(1,LL);
    DOA_mle         = zeros(1,LL);
    DOA_CBF_O       = zeros(1,LL);
    DOA_phase_O     = zeros(1,LL);
    DOA_phase_c_O   = zeros(1,LL);
    DOA_mle_O       = zeros(1,LL);
    
    
    MSE_x_CBF       = 0;
    MSE_x_phase     = 0;
    MSE_x_phase_c   = 0;
    MSE_x_mle       = 0;
    MSE_x_CBF_O     = 0;
    MSE_x_phase_O   = 0;
    MSE_x_phase_c_O = 0;
    MSE_x_mle_O     = 0;
    
    
    
    for ell=1:L:LL
        
        %% ------- conventional beamforming -------------
        % changed 0:L-1 to 0:L-2 (in not disturbed case) to be consistent with Zoubir
        b_CBF           = sum(abs(A'*y(:,ell + [0:(L-2)])).^2,2); % multiple snapshot CBF
        [maxval,index]  = max(b_CBF);
        b_CBF           = b_CBF/maxval;
        DOA_CBF(ell)    = find_interpolated_max(phi_vec,b_CBF,index);
        a_CBF           = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_CBF(ell))); % offgrid CBF steering vector estimate
        x_CBF           = mean(a_CBF' * y(:,ell + [0:(L-2)])/N); % complex source amplitude estimate from CBF
        MSE_x_CBF       = MSE_x_CBF + abs(x_CBF - x_src)^2;
        
        %-- outliers
        b_CBF           = sum(abs(A'*yO(:,ell + [0:(L-1)])).^2,2); % multiple snapshot CBF
        [maxval,index]  = max(b_CBF);
        b_CBF           = b_CBF/maxval;
        DOA_CBF_O(ell)  = find_interpolated_max(phi_vec,b_CBF,index);
        a_CBF           = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_CBF_O(ell))); % offgrid CBF steering vector estimate
        x_CBF           = mean(a_CBF' * yO(:,ell + [0:(L-1)])/N); % complex source amplitude estimate from CBF
        MSE_x_CBF_O     = MSE_x_CBF_O + abs(x_CBF - x_src)^2;
        
        %% -------- phase only beamforming --------------
        c_phase         = 1./ sum(1./ abs(y(:,ell + [0:(L-2)])),1);   % Equation (21)
        sign_y          = y(:,ell + [0:(L-2)]) ./ abs(y(:,ell + [0:(L-2)]));% Equation (22)
        c_phase_sign_y  = repmat(c_phase,N,1) .* sign_y;             % Equation (21+22)
        
        b_phase         = sum(abs(A'*sign_y).^2,2); % multiple snapshot phase-only beamformer (23) without c_phase(!)
        [maxval,index]  = max(b_phase);
        b_phase         = b_phase/maxval;
        DOA_phase(ell)  = find_interpolated_max(phi_vec,b_phase,index);
        a_phase         = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_phase(ell))); % offgrid phase-only steering vector estimate
        x_phase         = mean(a_phase'*c_phase_sign_y);
        MSE_x_phase     = MSE_x_phase + abs(x_phase - x_src).^2;
        
        %-- outliers
        c_phase         = 1./ sum(1./ abs(yO(:,ell + [0:(L-1)])),1);   % Equation (21)
        sign_y          = yO(:,ell + [0:(L-1)]) ./ abs(yO(:,ell + [0:(L-1)]));% Equation (22)
        c_phase_sign_y  = repmat(c_phase,N,1) .* sign_y;             % Equation (21+22)
        
        b_phase         = sum(abs(A'*sign_y).^2,2); % multiple snapshot phase-only beamformer (23) without c_phase(!)
        [maxval,index]  = max(b_phase);
        b_phase         = b_phase/maxval;
        DOA_phase_O(ell)= find_interpolated_max(phi_vec,b_phase,index);
        a_phase         = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_phase_O(ell))); % offgrid phase-only steering vector estimate
        x_phase         = mean(a_phase'*c_phase_sign_y);
        MSE_x_phase_O   = MSE_x_phase_O + abs(x_phase - x_src).^2;
        
        
        %% -------- phase only beamforming including c--------------
        c_phase         = 1./ sum(1./ abs(y(:,ell + [0:(L-2)])),1);   % Equation (21)
        sign_y          = y(:,ell + [0:(L-2)]) ./ abs(y(:,ell + [0:(L-2)]));% Equation (22)
        c_phase_sign_y  = repmat(c_phase,N,1) .* sign_y;             % Equation (21+22)
        
        b_phase         = sum(abs(A'*c_phase_sign_y).^2,2); % multiple snapshot phase-only beamformer (23) without c_phase(!)
        [maxval,index]  = max(b_phase);
        b_phase         = b_phase/maxval;
        DOA_phase_c(ell)    = find_interpolated_max(phi_vec,b_phase,index);
        a_phase         = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_phase(ell))); % offgrid phase-only steering vector estimate
        x_phase         = mean(a_phase'*c_phase_sign_y);
        MSE_x_phase_c   = MSE_x_phase + abs(x_phase - x_src).^2;
        
        %-- outliers
        c_phase         = 1./ sum(1./ abs(yO(:,ell + [0:(L-1)])),1);   % Equation (21)
        sign_y          = yO(:,ell + [0:(L-1)]) ./ abs(yO(:,ell + [0:(L-1)]));% Equation (22)
        c_phase_sign_y  = repmat(c_phase,N,1) .* sign_y;             % Equation (21+22)
        
        b_phase         = sum(abs(A'*c_phase_sign_y).^2,2); % multiple snapshot phase-only beamformer (23) without c_phase(!)
        [maxval,index]  = max(b_phase);
        b_phase         = b_phase/maxval;
        DOA_phase_c_O(ell)  = find_interpolated_max(phi_vec,b_phase,index);
        a_phase         = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_phase_O(ell))); % offgrid phase-only steering vector estimate
        x_phase         = mean(a_phase'*c_phase_sign_y);
        MSE_x_phase_c_O = MSE_x_phase_O + abs(x_phase - x_src).^2;
        
        
        
        %% --------- LAD beamforming -> MLE for Laplace-like noise ------------
        c_phase         = 1./ sum(1./ abs(y(:,ell + [0:(L-2)])),1);   % Equation (21)
        sign_y          = y(:,ell + [0:(L-2)]) ./ abs(y(:,ell + [0:(L-2)]));% Equation (22)
        c_phase_sign_y  = repmat(c_phase,N,1) .* sign_y;             % Equation (21+22)
        
        b_mle           = zeros(size(b_CBF));
        for m=1:M                                       % for all angles
            a           = A(:,m);
            x_hat1      = a' * c_phase_sign_y; % Eq. (19), also Eq. (38)
            parfor ell2=0:(L-2)
                x_hat2(ell2+1) = L1LinearRegression( a, y(:,ell + ell2),x_hat1(ell2+1));
            end
            b_mle(m)    = sum(sum(abs(y(:,ell + [0:(L-2)]) - a * x_hat2))); % Laplace-like Likelihood function in Equation (27)
        end
        b_mle           = 1./ (b_mle .* b_mle);
        [maxval,index]  = max(b_mle);
        b_mle           = b_mle/maxval;
        DOA_mle(ell)    = find_interpolated_max(phi_vec,b_mle,index);
        a_mle           = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_mle(ell))); % offgrid MLE steering vector estimate
        a               = conj(a_mle);
        x_hat1          = a.' * c_phase_sign_y; % Eq. (19), also Eq. (38)

        parfor ell2=0:(L-2)
            x_hat2(ell2+1) = L1LinearRegression( a_mle, y(:,ell + ell2),x_hat1(ell2+1));
        end
        x_hat2          = mean(x_hat2);
        MSE_x_mle       = MSE_x_mle + abs(x_hat2 - x_src).^2;
        
        %--outliers
        c_phase         = 1./ sum(1./ abs(yO(:,ell + [0:(L-1)])),1);   % Equation (21)
        sign_y          = yO(:,ell + [0:(L-1)]) ./ abs(yO(:,ell + [0:(L-1)]));% Equation (22)
        c_phase_sign_y  = repmat(c_phase,N,1) .* sign_y;             % Equation (21+22)
        
        b_mle           = zeros(size(b_CBF));
        for m=1:M                                       % for all angles
            a           = A(:,m);
            x_hat1      = a' * c_phase_sign_y; % Eq. (19), also Eq. (38)
            parfor ell2=0:(L-1)
                x_hat2(ell2+1) = L1LinearRegression( a, yO(:,ell + ell2),x_hat1(ell2+1));
            end
            b_mle(m)    = sum(sum(abs(yO(:,ell + [0:(L-1)]) - a * x_hat2))); % Laplace-like Likelihood function in Equation (27)
        end
        b_mle           = 1 ./ (b_mle .* b_mle);
        [maxval,index]  = max(b_mle);
        b_mle           = b_mle/maxval;
        DOA_mle_O(ell)  = find_interpolated_max(phi_vec,b_mle,index);
        a_mle           = exp(-1j*2*pi*d*[0:N-1]'*sind(DOA_mle_O(ell))); % offgrid MLE steering vector estimate
        a               = conj(a_mle);
        x_hat1          = a.' * c_phase_sign_y; % Eq. (19), also Eq. (38)

        parfor ell2=0:(L-1)
            x_hat2(ell2+1) = L1LinearRegression( a_mle, yO(:,ell + ell2),x_hat1(ell2+1));
        end
        x_hat2          = mean(x_hat2);
        MSE_x_mle_O     = MSE_x_mle_O + abs(x_hat2 - x_src).^2;
        
        
    end % end of parfor ell=1:L:LL
    
    
    if 0   
        %% CBF
        mean_DOA_CBF            = mean(DOA_CBF(1:L:LL));
        median_DOA_CBF          = median(DOA_CBF(1:L:LL));
        std_DOA_CBF             = std(DOA_CBF(1:L:LL));
        mad_DOA_CBF             = mad(DOA_CBF(1:L:LL),1);
        RMSE_DOA_CBF            = sqrt(mean((DOA_CBF(1:L:LL)-phi_vec(m_src)).^2));
        RMSE_x_CBF              = sqrt(MSE_x_CBF/(LL/L));
        %-outliers
        mean_DOA_CBF_O          = mean(DOA_CBF_O(1:L:LL));
        median_DOA_CBF_O        = median(DOA_CBF_O(1:L:LL));
        std_DOA_CBF_O           = std(DOA_CBF_O(1:L:LL));
        mad_DOA_CBF_O           = mad(DOA_CBF_O(1:L:LL),1);
        RMSE_DOA_CBF_O          = sqrt(mean((DOA_CBF_O(1:L:LL)-phi_vec(m_src)).^2));
        RMSE_x_CBF_O            = sqrt(MSE_x_CBF_O/(LL/L));
        
        %% Phase only
        mean_DOA_phase          = mean(DOA_phase(1:L:LL));
        median_DOA_phase        = median(DOA_phase(1:L:LL));
        std_DOA_phase           = std(DOA_phase(1:L:LL));
        mad_DOA_phase           = mad(DOA_phase(1:L:LL),1);
        RMSE_DOA_phase          = sqrt(mean((DOA_phase(1:L:LL)-phi_vec(m_src)).^2));
        RMSE_x_phase            = sqrt(MSE_x_phase/(LL/L));
        %-outliers
        mean_DOA_phase_O        = mean(DOA_phase_O(1:L:LL));
        median_DOA_phase_O      = median(DOA_phase_O(1:L:LL));
        std_DOA_phase_O         = std(DOA_phase_O(1:L:LL));
        mad_DOA_phase_O         = mad(DOA_phase_O(1:L:LL),1);
        RMSE_DOA_phase_O        = sqrt(mean((DOA_phase_O(1:L:LL)-phi_vec(m_src)).^2));
        RMSE_x_phase_O          = sqrt(MSE_x_phase_O/(LL/L));
        
        %% Phase only including c
        mean_DOA_phase_c        = mean(DOA_phase_c(1:L:LL));
        median_DOA_phase_c      = median(DOA_phase_c(1:L:LL));
        std_DOA_phase_c         = std(DOA_phase_c(1:L:LL));
        mad_DOA_phase_c         = mad(DOA_phase_c(1:L:LL),1);
        RMSE_DOA_phase_c        = sqrt(mean((DOA_phase_c(1:L:LL)-phi_vec(m_src)).^2));
        RMSE_x_phase_c          = sqrt(MSE_x_phase_c/(LL/L));
        %-outliers
        mean_DOA_phase_c_O      = mean(DOA_phase_c_O(1:L:LL));
        median_DOA_phase_c_O    = median(DOA_phase_c_O(1:L:LL));
        std_DOA_phase_c_O       = std(DOA_phase_c_O(1:L:LL));
        mad_DOA_phase_c_O       = mad(DOA_phase_c_O(1:L:LL),1);
        RMSE_DOA_phase_c_O      = sqrt(mean((DOA_phase_c_O(1:L:LL)-phi_vec(m_src)).^2));
        RMSE_x_phase_c_O        = sqrt(MSE_x_phase_c_O/(LL/L));
        
        %% LAD
        mean_DOA_mle            = mean(DOA_mle(1:L:LL));
        median_DOA_mle          = median(DOA_mle(1:L:LL));
        std_DOA_mle             = std(DOA_mle(1:L:LL));
        mad_DOA_mle             = mad(DOA_mle(1:L:LL),1);
        RMSE_DOA_mle            = sqrt(mean((DOA_mle(1:L:LL)-phi_vec(m_src)).^2));
        RMSE_x_mle              = sqrt(MSE_x_mle/(LL/L));
        %-outliers
        mean_DOA_mle_O          = mean(DOA_mle_O(1:L:LL));
        median_DOA_mle_O        = median(DOA_mle_O(1:L:LL));
        std_DOA_mle_O           = std(DOA_mle_O(1:L:LL));
        mad_DOA_mle_O           = mad(DOA_mle_O(1:L:LL),1);
        RMSE_DOA_mle_O          = sqrt(mean((DOA_mle_O(1:L:LL)-phi_vec(m_src)).^2));
        RMSE_x_mle_O            = sqrt(MSE_x_mle_O/(LL/L));
    end
    
    
    %% sensitivity curve = emperical influence function
    
    % angular difference from "First- and Second-Order Characterization of
    % DirectionDispersion and Space Selectivity in the RadioChannel"
   
    
    phasor = false;
    if phasor
        angDiff_CBF         = exp(1j*DOA_CBF_O(1:L:LL)*pi/180)-exp(1j*DOA_CBF(1:L:LL)*pi/180);
        angDiff_phase       = exp(1j*DOA_phase_O(1:L:LL)*pi/180)-exp(1j*DOA_phase(1:L:LL)*pi/180);
        angDiff_phase_c     = exp(1j*DOA_phase_c_O(1:L:LL)*pi/180)-exp(1j*DOA_phase_c(1:L:LL)*pi/180);
        angDiff_mle         = exp(1j*DOA_mle_O(1:L:LL)*pi/180)-exp(1j*DOA_mle(1:L:LL)*pi/180);    
    else
        angDiff_CBF         = angDiff(DOA_CBF_O(1:L:LL),DOA_CBF(1:L:LL));
        angDiff_phase       = angDiff(DOA_phase_O(1:L:LL),DOA_phase(1:L:LL));
        angDiff_phase_c     = angDiff(DOA_phase_c_O(1:L:LL),DOA_phase_c(1:L:LL));
        angDiff_mle         = angDiff(DOA_mle_O(1:L:LL),DOA_mle(1:L:LL));     
    end
    SC_CBF(oCounter)        = L*mean(abs(angDiff_CBF));
    SC_phase(oCounter)      = L*mean(abs(angDiff_phase));
    SC_phase_c(oCounter)    = L*mean(abs(angDiff_phase_c));
    SC_mle(oCounter)        = L*mean(abs(angDiff_mle));
    
    
end

figure(1)
title(['Emperical Influence Function with L=' num2str(L) ' and SNR=' num2str(round(SNR)) 'dB'])
hold on
set(gca,'fontsize',14)
plot(outlierRadVec,SC_CBF,'color',[0.8500    0.3250    0.0980],'linewidth',2)
plot(outlierRadVec,SC_phase,'color', [0    0.4470    0.7410],'linewidth',2)
plot(outlierRadVec,SC_phase_c,'color',[0.4940    0.1840    0.5560],'linewidth',2)
plot(outlierRadVec,SC_mle,'color',[0.9290    0.6940    0.1250],'linewidth',2)
grid on
box on
xlabel('outlier radius','fontsize',16)
ylabel('EIF','fontsize',16)
leg=legend('CBF','phase only','phase only incl. c','LAD');
set(leg,'fontsize',16,'location','northwest')
print('-dpng','EIF_new_longSim')

save saved_workspace_new_phase_long.mat