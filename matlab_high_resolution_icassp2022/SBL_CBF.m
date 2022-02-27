function [ gamma, report] = SBL_CBF( A , Y, options )
%
% function [ gamma , report ] = SBL_CBF( A , Y , Nfreq, Nsource,flag )
% conventional beamformer
%
% Inputs
%
% A - Multiple frequency augmented dictionary < n , m, f>
%     f: number of frequencies
%     n: number of sensors
%     m: number of replicas
%   Note: if f==1, A = < n , m >
%
%
% Y - Multiple snapshot multiple frequency observations < n , L, f>
%     f: number of frequencies
%     n: number of sensors
%     L: number of snapshots
%
% options - see SBLset.m 
%
%
% Outputs
%
% gamma <m , 1> - vector containing source power
%                 1: surfaces found by minimum error norm
%
% report - various report options
%
%--------------------------------------------------------------------------
% Version 1.0:
% Code originally written by cfm
%%
options.SBL_v = 'CBF';

%% slicing
% 
% if options.tic == 1
%     tic
% end

%% Initialize variables
Nfreq     = size(A,3);  % number of frequencies
Nsource = options.Nsource; %number of sources
Nsensor   = size(A,1);% number of sensors
Ntheta    = size(A,2);% number of dictionary entries
Nsnapshot = size(Y,2);% number of snapshots in the data covariance
% noise power initialization
sigc      = ones(Nfreq,1) * options.noisepower.guess;
%sigc      = options.noisepower.guess;
% posterior
%x_post    = zeros( Ntheta, Nsnapshot,Nfreq);
% minimum (global) gamma
gmin_global = realmax;
% reduce broacast
%phi         = options.phi;
% space allocation
errornorm   = zeros(options.convergence.maxiter,1);

% initialize with CBF output (assume single frequency)
%gamma        = Bartlett_processor(squeeze(A), squeeze(Y));
% gamma        = zeros(Ntheta, 1);
% for iF = 1:Nfreq
%     Af = squeeze(A(:,:,iF));
%     Yf = squeeze(Y(:,:,iF));
%     gamma = gamma+sum(abs(Af' * Yf).^2, 2) / Nsnapshot;
% end

gamma        = ones(Ntheta, 1);
gamma_num    = zeros(Ntheta,Nfreq);
gamma_denum  = zeros(Ntheta,Nfreq);

% Sample Covariance Matrix
SCM = zeros(  Nsensor , Nsensor, Nfreq);
for iF = 1 : Nfreq
    SCM(:,:,iF) = squeeze(Y(:,:,iF)) * squeeze(Y(:,:,iF))' / Nsnapshot;
    maxnoise(iF)=real(trace(squeeze(SCM(:,:,iF))))/Nsensor;
    sigc(iF)=maxnoise(iF);
    CBF(:,iF) = real(diag(conj(squeeze(A(:,:,iF)).') * squeeze(SCM(:,:,iF))* squeeze(A(:,:,iF))));   
end
gamma=sum(CBF,2);

%% Report section
% vectors containing errors
report.results.error    = NaN;
% Error when minimum was obtained
report.results.iteration_L1 = NaN;
% General info
report.results.final_iteration.iteration = NaN;
report.results.final_iteration.noisepower = sigc;

% debug output parameters (assuming single frequency)
%report.SigmaYinv = SigmaYinv;
%report.SCM = squeeze(SCM);
% 
% if options.tic == 1
%     report.results.toc = toc;
% else
%     report.results.toc = 0;
% end
% data
report.results.final_iteration.gamma  = gamma  ;
%report.results.final_iteration.x_post = x_post ;
report.options = options;
end