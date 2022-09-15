function [gamma,sigc, Ilocs, phi_vec,report] = SBL_pre_Esa(A,Y,method,phi_vec,K,upar)
%
% function [ gamma , report ] = SBL(A , Y , method, phi_vec, nu, K )
%
% Inputs
%
% A       - Dictionary of size N x M
% Y       - Matrix of N x L (L is the number of snapshots)
% method  - either 'SBL-G', 'SBL-T' or 'SBL-H' (default = 'SBL-T')
% phi_vec - DOA grid  
% K       - number of sources 
% upar    - parameter for u(.) function. Equals degrees of freedom parameter 
%          (default upar = 2.1 for MVT-loss function) and the quantile 0<q<1 
%          for Huber's loss function (default q = 0.85). 
% 
% Outputs
%
% gamma   -  vector containing source power
% sigc    -  noise variance estimate 
% Ilocs   -  index locations on the grid (phi_vec) of found DoA-s 
% phi_vec -  used phi_vec grid 
% report  -  a structure array with some stuff


%% Initialize variables
N = size(A,1);% number of sensors
M = size(A,2);% number of dictionary entries
L = size(Y,2);% number of snapshots in the data covariance

robustSBL = true; % default is to use robust SBL-T method 
if strcmpi(method,'SBL-G')
    ufun = @(t,v) ones(size(t)); % this is actually not used 
    b = 1;
    upar = 1;
    robustSBL = false;
elseif (strcmpi(method,'SBL-T') || isempty(method) )
    ufun = @(t,v) (v+2*N)./(v + 2*t);  % weight function is t-loss 
    if nargin < 6 
       upar = 2.1; % d.o.f v=2.1 is used as the default parameter for t-loss
    end   
    assert( isreal(upar) && upar>=2 && ~isinf(upar),'loss parameter is incorrect');  
    b = tloss_consistency_factor(N,upar,false); % concistency factor    
elseif (strcmpi(method,'SBL-H') || isempty(method) )   
    ufun = @(t,c) ((t<=c) + (c./t).*(t>c)); % weight function u(t)
    if nargin < 6
       q = 0.85;
    else
       q = upar;
    end          
    assert(isreal(q) && isfinite(q) && (0<q) && (q <1),'loss parmeter is incorrect'); 
    upar = chi2inv(q,2*N)/2;  
    b = chi2cdf(2*upar,2*(N+1))+(upar/N)*(1-q);    
else
    error('please specify method as a string equalto SLB-G, SBL-T, or SBL-H\n');
end
const = 1/(b*L);

%% -- Compute initial locations of DoA-s 
anorm = norm(A(:,1)); % assumed that steering vectors all have same length
RY = (1/L)*Y*(Y'); 
SCM = RY;
maxnoise   = real(trace(SCM))/N;
gamma_init = real(sum(conj(A).*(RY*A)))/anorm^4;%  
if K==1 
   [~,indx] = max(gamma_init);
else
    [~, indx] = find_peaks_1D(gamma_init,K);
end
Ilocs  = indx;

%% -- Compute the initial noise variance 
Am = A(:,indx);        % only active replicas
Aplus = pinv(Am);
P_N = eye(N)-Am*Aplus;
sigc = real(trace(P_N*RY))/(N - K);  % noise covariance estimate 

%% -- Compute initial estimates of gamma-s 
gamma = (N/(N-1))*subplus(gamma_init - real(trace(RY)/N)/anorm^2);  % estimates of  gamma 
izeros = gamma==0;
delta_value = (1/10)*min(gamma(~izeros));
gamma(izeros) = delta_value; % assign negative gammas a small positive value 

%%  initial screening ? 
% It could be possible to do some  initial screening (so to use only non-zeros gammas)
% in order to reduce the grid to locations that are a priori sensible. 
% --
% inonzeros = gamma~=0;
% gamma = gamma(inonzeros);
% A =  A(:,inonzeros);
% M = size(A,2);% number of dictionary entries
% phi_vec = phi_vec(inonzeros);
% [~,indx] = max(gamma);
% Ilocs  = indx;
%% Algorithm parameters
status_report = 20; % status report every xx iterations
flag = false;   % print report (false, then do not print)
convergence_error   = 1*10^(-3);
max_iter = 1000; % maximum number of iterations allowed
min_iter = 15;   % solution only accepted at this iteration

errornorm   = zeros(max_iter,1);
deltaval = 1;
gamma_range = 10^-4;

%%
for j1 = 1:max_iter

    gammaOld = gamma;
    IlocsOld = Ilocs;
    SigmaYinv   =  (sigc*eye(N) + bsxfun(@times,gamma,A)*A')\eye(N); 
    
    if robustSBL 
        t = real(sum(conj(Y).*(SigmaYinv*Y))); % norms 
        u = ufun(t,upar);
        RY = const*(Y.*repmat(u,N,1))*Y'; 
    end
    Itheta= find(gamma>max(gamma)*gamma_range);
    
    B =  SigmaYinv*A(:,Itheta); % \Sigma^-1 a_m , m=1,..,M
    gamma_num = subplus(real(sum(conj(B).*(RY*B))));
    %gamma_num = const*(sum(abs((B'*Y).*repmat(sqrt(u),L,1)).^2,2)).'; 
    gamma_denum = subplus(real(sum(conj(A(:,Itheta)).*B)));
    gamma(Itheta) = gamma(Itheta).*(gamma_num./gamma_denum).^(deltaval);
    
    [~, Ilocs] = find_peaks_1D(gamma,K);

    if ~isempty(setdiff(IlocsOld,Ilocs))
        Am = A(:,Ilocs);        % only active replicas
        P_N = eye(N)-Am*pinv(Am);
        sigc = real(trace(P_N*RY))/(N - K);

        sigc        = min(sigc,maxnoise);  %cant be larger than signal+noise.
        sigc        = max(sigc,maxnoise*10^-10); % snr>100 is unlikely larger than signal.
    end
 
    err1 = norm(gamma-gammaOld,1)/norm(gamma,1);
    errornorm(j1) = err1;
    %gamma(Ilocs)

    if j1 > min_iter && (err1 < convergence_error)
        if flag 
            fprintf('Solution converged!\nIteration: %4d. Error: %.7f\n', j1, err1)
        end
        break; % goodbye     
    end
    
    if flag && mod(j1,status_report) == 0
        fprintf('Iteration: %4d. Error: %.7f\n', j1, errornorm(j1))
    end

end

if j1 == max_iter % not convereged
   if flag
       fprintf('Solution not converged. Error: %.6f,', errornorm(j1)); 
   end
end


%% Report section
report.error    = errornorm;
report.iteration = j1;

end

function b=tloss_consistency_factor(p,v,realdata)
%  Function that computes the scaling factor for multivariate t-weight
%  function so that the returned scatter matrix is concistent estimator of
%  the covariance matrix under the assumption that the data is from 
%  Gaussian distribution
%--------------------------------------------------------------------------   

% First try by numerical integration 
if realdata 
     b = tloss_consistency_factor_int(p,v);
else
     b = tloss_consistency_factor_int(2*p,v);        
end

% If integration did not converge, then compute b by MC simulation 
% --
if isnan(b)
   % need to use MC simul to find b
   MCsimul = 100000;
   
   if realdata
        t = chi2rnd(p,1,MCsimul);
        psifun = @(t,v) (p+v)*t./(v+t); % weight function
   else 
        t = (1/2)*chi2rnd(2*p,1,MCsimul);
        psifun = @(t,v) (v+2*p)*t./(v+2*t); % weight function
   end
   b = (1/p)*mean(psifun(t,v));

end

end

function b=tloss_consistency_factor_int(p,v)
% computes the concistency factor b = (1/p) E[|| x ||^2 u_v( ||x||^2)] when
% x ~ N_p(0,I). 

sfun = @(x,p,v)  (x.^(p/2)./(v+ x) .* exp(-x/2));
c = 2^(p/2)*gamma(p/2);
w = warning('off','all');
q = (1/c)*integral(@(x)sfun(x,p,v),0,Inf);
b = ((v+p)/p)*q; % consistency factor  
warning(w)
end

function [pks, locs] = find_peaks_1D(gamma, Nsources)
% This is the code SBLpeaks_1D with different name:
% 
% [pks, locs] = SBLpeaks_1D(gamma, Nsources)
%
% fast alternative for findpeaks in 1D case
%

% output variables
pks = zeros(Nsources,1);
locs = zeros(Nsources,1);

% zero padding on the boundary
gamma_new = zeros(length(gamma)+2,1);
gamma_new(2:end-1) = gamma;

[~, Ilocs]= sort(gamma,'descend');

% current number of peaks found
npeaks = 0;

for ii = 1:length(Ilocs)
    
    % local patch area surrounding the current array entry i.e. (r,c)
    local_patch = gamma_new(Ilocs(ii):Ilocs(ii)+2);
    
    % zero the center
    local_patch(2) = 0;
    
    if sum(sum(gamma(Ilocs(ii)) > local_patch)) == 3
        npeaks = npeaks + 1;
        
        pks(npeaks) = gamma(Ilocs(ii));
        locs(npeaks) = Ilocs(ii);
        
        % if found sufficient peaks, break
        if npeaks == Nsources
            break;
        end
    end
    
end

% if Nsources not found
if npeaks ~= Nsources
    pks(npeaks+1:Nsources) = [];
    locs(npeaks+1:Nsources) = [];
end

end
