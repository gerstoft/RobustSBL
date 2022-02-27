function CRLB = CRLB_VanTrees8p106(A,phi_vec,m_src,Gamma_src,ASNR_lin,L,ULA_flag)
%CRLB = CRLB_VanTrees8p106(A,phi_vec,m_src,gamma_src,ASNR_lin,L,ULA_flag)
%
%Evaluates the trace of the general CRLB expression according to Van Trees Book, Eq.(8.106)
%The sources must be on the grid, i.e. the true steering vectors are some columns in the dictionary matrix A.
%
%A          dictionary matrix (N x M), columns are steering vectors for DOAs in phi_vec. 
%           We need normalization norm(A(:,m),2) == sqrt(N) for all m = 1:M
%
%phi_vec    corresponding vector of DOAs, DOAs need not be equidistant, size must be 1 x M
%
%m_src      indices of true sources, pointers into A and phi_vec, size must be 1 x K
%
%Gamma_src  source covariance matrix, should be Hermitian positive semidefinite, 
%           often just a diagonal matrix, size must be square, K x K
%
%ASNR_lin   vector of Array SNR values for the plot, linear values, not dB
%
%L          number of snapshots used for sample covariance matrix, L should be larger than N
%
%ULA_flag   if ULA_flag == 1 then analytical derivative of plane wave steering vectors is used, 
%                This requires that array geometry is ULA with half wavelength spacing
%                and the phi_vec elements list the DOAs in degrees, not radians.
%           otherwise the derivative of the steering vectors is approximated by difference quotient. 
%                This approximation is valid for any array manifold, but requires a high resolution dictionary
%
% cfm @ UCSD, Feb 21, 2022

[N,M] = size(A); % N = number of sensors, M = dictionary size
K = length(m_src); % K = number of sources
% L = number of snapshots

A_M = A(:, m_src);             % true steering vectors  $A_\mathcal{M}$

if length(phi_vec) ~= M
    error('number of columns in dictionary matrix A does not match length of phi_vec')
end

if size(Gamma_src,1) ~= K || size(Gamma_src,2) ~= K,
    error('problem with size of Gamma_src');
end

sigma2 = N ./ ASNR_lin;


% Van Trees Book, Eq. (8.96) 
% Projection onto noise subspace
I_N = eye(N);
P_perp = I_N - A_M * pinv(A_M);

% Van Trees Book, Eq. (8.97)
if ULA_flag == 1,
  % analytical derivative of ULA steering vectors w.r.t. azimuth
  % this is valid only for ULA with half wavelength spacing
  % phi_vec must be in degrees
  D = complex(zeros(N,K),zeros(N,K));
  for k=1:K
     D(:,k)   = exp(-1j * pi * [0:N-1]' * sind(phi_vec(m_src(k)))); % normalization to |a_n|=1 or ||a||_2^2 = N.
     D(:,k)   = D(:,k) .* ( -1j * pi * [0:N-1]' * cosd(phi_vec(m_src(k))) ) * (pi/180.0);
  end
else
  % central difference to approximate derivative of steering vectors w.r.t. azimuth
  % this is a valid approximation for any antenna array manifold
  m1 = find(m_src == 1); m_src(m1) = 2;   % handle boundary case
  m2 = find(m_src == M); m_src(m2) = M-1; % handle boundary case
  if ~isempty(m1) || ~isempty(m2)
      warning('boundary cases for DOA')
  end
  delta_phi = phi_vec(m_src+1) - phi_vec(m_src-1); % degrees or radians, the unit doesn't matter
  D = (A(:, m_src+1) - A(:, m_src-1)) ./ (ones(N,1) * delta_phi); 
end

% Van Trees Book, Eq. (8.101)
H  = D' * P_perp * D;
HT = H.';

I_K = eye(K);

CRLB = zeros(size(ASNR_lin));

% we now evaluate Van Trees Book, Eq.(8.106) for all elements in the sigma2 vector
for isigma2 = 1:length(sigma2)
    Aux_M = A_M' * A_M * Gamma_src / sigma2(isigma2);
    CRLB_matrix = (sigma2(isigma2)/ (2*L)) * inv(real(Gamma_src * ( inv(I_K + Aux_M)*Aux_M ) .* HT));
    CRLB(isigma2) = real(trace(CRLB_matrix));
end

return
