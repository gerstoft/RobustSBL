%% Esa's simpler version of 'doa_wsa2021_SBL4.m'
% Assumes: 
% - ULA with lambda/2 spacing. 
% - uses parfor to make the computations faster

clearvars;
%Algorithm (robust) parameters
nu_algorithm = 2.1; 

%% Choose Model:
mu=0; % mean of noise
% OPTION 1: Gaussian noise: 
epsilon = 0; lambda = 1;  % Gaussian case
% OPTION 2: epsilon-contaminated noise: 
%epsilon = 0.05; lambda  = 10.0; % Non-Gaussian case
        
%% array steering matrix of size N x M: ULA with half a vawelength spacing
% 1 degree grid spacing 
N = 20;  % no. of sensors in ULA
%M = 181;
M = 18001; % high resolution dictionary, 0.01 deg azimuth resolution
dphi = 180/(M-1);
phi_vec = -90:dphi:90;
A = exp(-1j*pi*(0:(N-1))'*sind(phi_vec)); % ULA matrix
avec = @(theta) exp(-1j*pi*(0:(N-1))'*sind(theta));

%% Model parameters 
L = 25;  % number of snapshots used by DOA estimator
K = 1;   % number of sources
Number_of_DOAs = 1; % number of sources, use 1 or 3
DOA_src = -45;      % -45 deg source on the grid
SNR = 24:-3:-6;     % SNR in dB
sigma_vec = sqrt(10.^(-SNR/10));
x_src   = exp(1j*pi/4); % source amplitude, magnitude = 1    
m_src = find(phi_vec == DOA_src);
a_src = A(:,m_src);

%% Simulation parameters
LL = 350; % number of independent Monte-Carlo trials 

mean_DOA_SBL4   = zeros(length(sigma_vec),5);
RMSE_DOA_SBL4   = zeros(length(sigma_vec),5);
RMSE_x_SBL4     = zeros(length(sigma_vec),3);
    
options = SBLSet; % this is needed for SBL4 (but not for SBL4_v4 only)
%%

% Cramer-Rao Lower Bound for single source and uniform linear array, Gaussian noise
ASNR = 10.^(SNR/10)*N;
CRLB = (1/L)*( 1./ASNR + 1./(ASNR).^2 )  * ( 6/(N^2-1));
degrees = 180.0/pi;
CRLBdeg = (degrees^2)*CRLB/(pi^2*(1-cosd(DOA_src)^2));

%%
    
for isnr = 1:length(sigma_vec)
    
    sigma = sigma_vec(isnr);
    fprintf('%d/%d , SNR= %.3f\n',isnr,length(sigma_vec),SNR(isnr));
    
    %% Loop over snapshots
    DOA_SBL4  = zeros(LL,2);
    AVEd0a = 0; AVEd0 = 0; AVEd1  = 0; AVEd2 = 0; AVEd3 = 0;
    MSEd0a = 0; MSEd0 = 0; MSEd1  = 0; MSEd2 = 0; MSEd3 = 0;
    MSEx1  = 0; MSEx2 = 0; MSEx3 = 0;
        
    %%
    rng('default');
    tic; 
  %  parfor ell=1:LL
    for ell=1:LL
        
        % Generate the L snapshots 
        s = complex(randn(1,L),randn(1,L))/sqrt(2);
        y = (a_src*s) + epscont(N,L,sigma,epsilon,mu,lambda);
        %%
        
        [DOAcbf,tmp] = CBF(y);
        RY = (1/L)*y*(y'); 
        AVEd0a = AVEd0a + DOAcbf;
        MSEd0a = MSEd0a + (DOAcbf - DOA_src)^2;
        
        objfnc = @(theta) (-real(avec(theta)'*RY*avec(theta))); 
        thetamin = max(DOAcbf-2,-90);
        thetamax = min(DOAcbf+2,90);
        DOA0 = fminbnd(objfnc,thetamin,thetamax);
        
        AVEd0 = AVEd0 + DOA0;
        MSEd0 = MSEd0 + (DOA0 - DOA_src)^2;        
        
        %% SBL_v4    
        
        % tic;
        [gamma1,report1] = SBL_v4(A, y, options);   
        % toc;
        b1 = gamma1;
        [maxval1,index1] = max(b1); 
        b1 = b1/maxval1;
      
        DOA1 = phi_vec(index1);
        AVEd1 = AVEd1 + DOA1;
        MSEd1 = MSEd1 + (DOA1 - DOA_src)^2;
        
        a1 = exp(-1i*pi*sind(DOA1)*(0:N-1).');
        x1 =  mean(a1'*y/N); % here I use just the CBF
        MSEx1 = MSEx1 + abs(x1 - x_src).^2;
                      
        %% SBL_v4 (t-loss)
        
        %tic;
        [gamma2,sigc2,report2] = SBL_v4_M_icassp(A,y,'t-loss',nu_algorithm,K);
        %toc;
        b2 = gamma2;
        [maxval2,index2] = max(b2);   
        b2 = b2/maxval2;
        %DOA2 = find_interpolated_max(phi_vec,b2,index2);
        DOA2 = phi_vec(index2);
        AVEd2 = AVEd2 + DOA2;
        MSEd2 = MSEd2 + (DOA2 - DOA_src)^2; 
        
        a2 = exp(-1i*pi*sind(DOA2)*(0:N-1).');
        x2 =  mean(a2'*y/N); % here I use just the CBF

        MSEx2 = MSEx2 + abs(x2 - x_src).^2;
        
        %% SBL_v4 (Huber's-loss)
        
        %tic;
        [gamma3,sigc3,report3] = SBL_v4_M_icassp(A,y,'Huber',0.9,K);
        %toc;
        b3 = gamma3;
        [maxval3,index3] = max(b3);   
        b3 = b3/maxval3;
        %DOA3 = find_interpolated_max(phi_vec,b3,index3);
        DOA3 = phi_vec(index3);
        AVEd3 = AVEd3 + DOA3;
        MSEd3 = MSEd3 + (DOA3 - DOA_src)^2; 
        
        a3 = exp(-1i*pi*sind(DOA3)*(0:N-1).');
        x3 =  mean(a3'*y/N); % here I use just the CBF

        MSEx3 = MSEx3 + abs(x3 - x_src).^2;
        
        %%
        %[index1 index2 index3]
             
    end
    toc;
    mean_DOA_SBL4(isnr,:) = [AVEd0a AVEd0 AVEd1 AVEd2 AVEd3]/LL;
    RMSE_DOA_SBL4(isnr,:)  = sqrt([MSEd0a MSEd0 MSEd1 MSEd2 MSEd3]/LL);
    RMSE_x_SBL4(isnr,:)    = sqrt([MSEx1 MSEx2 MSEx3]/LL);

    %% Print 
    fprintf('\nepsilon = %.2f. SNR = %f dB\n',epsilon,SNR(isnr));
    fprintf('RMSE(DOA_SBL4)  = %.3f %.3f %.3f degrees\n', RMSE_DOA_SBL4(isnr,1),RMSE_DOA_SBL4(isnr,2),RMSE_DOA_SBL4(isnr,3));
    %fprintf('RMSE x_SBL4     = %f \n\n',RMSE_x_SBL4(isnr));
    %%
end
fignro=5;

%% Plot
figure(fignro); clf
semilogy(SNR,RMSE_DOA_SBL4(:,2),'b*-','DisplayName','CBF','LineWidth',2);
hold on;
%semilogy(SNR,RMSE_DOA_SBL4(:,1),'b*-','DisplayName','CBF-grid')
semilogy(SNR,RMSE_DOA_SBL4(:,3),'g*-','DisplayName','SBL4-G','LineWidth',1.5);
semilogy(SNR,RMSE_DOA_SBL4(:,4),'ro-','DisplayName','SBL4-T','LineWidth',1.5);
semilogy(SNR,RMSE_DOA_SBL4(:,5),'k*-','DisplayName','SBL4.H','LineWidth',1.5);
semilogy(SNR,sqrt(CRLBdeg),'-','DisplayName','CRLB','LineWidth',2)
legend('Location','Best','FontSize',18)
title(sprintf('DOA estimates (%d DOAs, \\epsilon=%.2f, \\lambda=%.0f, LL = %d, \\nu=%.2f)',Number_of_DOAs,epsilon,lambda,LL,nu_algorithm),'FontSize',18);
xlabel('SNR (dB)')
ylabel('RMSE of DOA Estimate (dB)')
grid on
xlim([min(SNR) max(SNR)])
set(gca,'FontSize',20,'LineWidth',1.5,'FontName','Helvetica');

%%
figure(fignro+1); clf
hold on;
plot(SNR,mean_DOA_SBL4(:,2),'b*-','DisplayName','CBF','LineWidth',2);
plot(SNR,mean_DOA_SBL4(:,3),'g*-','DisplayName','SBL4-G','LineWidth',1.5);
plot(SNR,mean_DOA_SBL4(:,4),'ro-','DisplayName','SBL4-T','LineWidth',1.5);
plot(SNR,mean_DOA_SBL4(:,5),'k*-','DisplayName','SBL4.H','LineWidth',1.5);
title(sprintf('DOA estimates (%d DOAs, eps=%.2f, LL = %d, nu=%.2f)',Number_of_DOAs,epsilon,LL,nu_algorithm));
xlabel('SNR (dB)')
ylabel('Average DOA Estimate (^\circ)')
grid on
xlim([min(SNR) max(SNR)])

