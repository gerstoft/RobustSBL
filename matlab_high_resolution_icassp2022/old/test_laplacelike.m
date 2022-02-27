
mu = 1+1j;
sigma = 3;

N =   40;
L = 100000;
y = laplacelike_rand(mu,sigma,N,L);

figure(1)
subplot(2,1,1)
hist(real(y(:)),100);
title('histogram of real part')
subplot(2,1,2)
hist(imag(y(:)),100);
title('histogram of imag part')

mu_hat = mean(y(:));

fprintf('mean     = complex(%g, %g)\n',real(mu_hat),imag(mu_hat));
fprintf('median(real(y)) = %g\n',median(real(y(:))));
fprintf('median(imag(y)) = %g\n\n',median(imag(y(:))));
%
fprintf('variance = %g\n',var(y(:)));
fprintf('E(|y|^2) = %g\n',mean(abs(y(:)).^2));
fprintf('std      = %g\n\n',std(y(:)));
%
fprintf('kurtosis(real(y(:)),1) = %g\n',kurtosis(real(y(:))));
fprintf('kurtosis(imag(y(:)),1) = %g\n\n',kurtosis(imag(y(:))));

figure(2)
x = laplacernd(real(mu),sigma,N,L);
qqplot(x(:),real(y(:)-mu));
title('QQ Plot of real-valued Laplace X vs. Re(Y-E(Y))')

figure(3)
x = abs(x);
qqplot(x(:),abs(y(:)-mu).^(6/5));
title('QQ Plot of real-valued Laplace X vs. |Y-E(Y)|^{1.2}')

figure(4)
plot(real(y(1:1000)),imag(y(1:1000)),'.',real(mu),imag(mu),'ro')
axis equal

% now we check the quantiles of the Laplace-like distribution
a_vec = sigma * [1/8 1/4 1/2 1 2 4 8];
for index = 1:length(a_vec)
    a = a_vec(index);
    count  = sum(abs(y(:)-mu) > a); % how many samples are far away from mu ?
    as = a * sqrt(6)/sigma;
    Ecount = (1+as)*exp(-as) * N * L; % the expected count for Laplace-like
    fprintf('sigma=%f, a=%f, count = %d, E(count)=%d, count/E(count)=%f\n',sigma,a,count,round(Ecount),count/Ecount);
end


