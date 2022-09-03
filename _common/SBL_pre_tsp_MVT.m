function [gamma,sigc, report] = SBL_pre_tsp_MVT(A,Y,nu,K)
%
% function [ gamma , report ] = SBL_v4_tM( A , Y , nu, K )
% EXPERIMENTAL FUNCTION FOR WSA2021 Manuscript, see Table 1.
%
% The idea behind SBL is to find a diagonal replica 'covariance' Gamma.
% Minimizing (YY^T / AGA^T + penality) should lead to the correct
% replica selection (up to a bogus scale factor/amplitude).
%
% Inputs
%
% A  - Dictionary of size N x M
% Y  - Matrix of N x L (L is the number of snapshots)
% nu - Degrees of Freedom parameter used in the robust  t-loss function
% K  - number of sources 
% 
% Outputs
%
% gamma  -  vector containing source power
% sigc   -  noise variance estimate 
% report -  a structure array with some stuff


%% Initialize variables
N = size(A,1);% number of sensors
M = size(A,2);% number of dictionary entries
L = size(Y,2);% number of snapshots in the data covariance

ufun = @(t,v) (v+2*N)./(v + 2*t); % weight function is t-loss 
b = tloss_consistency_factor(N,nu,false); % concistency factor      
const = 1/(b*L);

%% Compute initial sigc (noise variance) and gamma 

RY = (1/L)*Y*(Y'); 
% Other option I tried when N > L is to compute M-estimate of scatter: 
% RY = MVT(Y.','nu',nu,'scaling','gaussian');
% and use it in place of SCM  above but this gave worse results.
gamma = real(sum(conj(A).*(RY*A)));%  
[~,indx] = max(gamma);
Am = A(:,indx);        % only active replicas
P_N = eye(N)-Am*Am'/N; % orthogonal projection matrix 
sigc = real(trace(P_N*RY))/(N - K); 
% possible re-update of gamma: ?
% gamma = real(sum(conj(A).*((RY-sigc*eye(N))*A)))/N^2;
% but also this did not work...

%% Algorithm parameters
status_report = 20; % status report every xx iterations
flag = false;   % print report (false, then do not print)
convergence_error   = 10^(-3);
max_iter = 1200; % maximum number of iterations allowed
min_iter = 15;   % solution only accepted at this iteration

errornorm   = zeros(max_iter,1);

%%
for j1 = 1:max_iter

    gammaOld = gamma;

    SigmaYinv   =  (sigc*eye(N) + bsxfun(@times,gamma,A)*A')\eye(N); 
    t = real(sum(conj(Y).*(SigmaYinv*Y))); % norms 
    u = ufun(t,nu);
    RY = const*(Y.*repmat(u,N,1))*Y'; 
    B =  SigmaYinv*A; % \Sigma^-1 a_m , m=1,..,M
    gamma_num = subplus(real(sum(conj(B).*(RY*B))));
    %gamma_num = const*(sum(abs((B'*Y).*repmat(sqrt(u),L,1)).^2,2)).'; 
    gamma_denum = subplus(real(sum(conj(A).*B)));
    gamma = gamma.*(gamma_num./gamma_denum);
    [~, Ilocs] = SBLpeaks_1D(gamma,K);
    Am = A(:,Ilocs);        % only active replicas
    P_N = eye(N)-Am*pinv(Am);
    sigc = real(trace(P_N*RY))/(N - K); 
 
    errornorm(j1) = norm(gamma-gammaOld,1)/norm(gamma,1);

    if j1 > min_iter && (errornorm(j1) < convergence_error)
        if flag 
            fprintf('Solution converged!\nIteration: %4d. Error: %.7f\n', j1, errornorm(j1))
        end
        break; % goodbye     
    elseif j1 == max_iter % not convereged
        if flag 
            fprintf('Solution not converged. Error: %.6f,', errornorm(j1))
        end
    elseif flag && mod(j1,status_report) == 0 
        fprintf('Iteration: %4d. Error: %.7f\n', j1, errornorm(j1))
    end

end


%% Report section
report.error    = errornorm;
report.iteration = j1;

end

function [pks, locs] = SBLpeaks_1D(gamma, Nsources)
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
