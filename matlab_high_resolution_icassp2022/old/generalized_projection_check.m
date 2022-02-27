
N=20;
M=3;
A = complex(randn(N,M),randn(N,M));

Sigma = A*diag(randn(M,1).^2)*A'+eye(N);

P1 = A*inv(A'*A)*A';

P2 = A*inv(A'*inv(Sigma)*A)*A'*inv(Sigma);

