clear;

alpha_vals = linspace( 1, 200, 2000 );
param_vals = linspace( 0.1, 20, 500 );
n = length(alpha_vals);

T_exact = zeros(1,n);
T_am    = zeros(1,n);
T_asymp = zeros(1,n);

for i=1:n
  alpha  = alpha_vals(i);
  k      = 0.2; %param_vals(i);
  beta   = k*alpha;
  lambda = k/(1+k);

  s = 6; %param_vals(i);
  t = 1; %param_vals(i); 

  T_exact(i) = integral( @(l) exp(-t*l) .* transit_density( l, s, alpha, beta, lambda ), 0, s ) ...
             + (1-lambda)*exp(-beta*s) + lambda*exp(-alpha*s)*exp(-t*s);
  T_am(i)    = exp(-lambda*t*s);
  T_asymp(i) = porous_asymp( s, t, alpha, beta, lambda );
end

err_am    = abs( T_exact - T_am    );
err_asymp = abs( T_exact - T_asymp );
clf;
% plot( alpha_vals, T_exact, '-k' );
% hold on
% plot( alpha_vals, T_asymp, '--k' );
% plot( alpha_vals, T_am, '-.k' );
% legend( 'Exact', 'Asymptotic', 'Atomic Mix' )
% xlabel( '\alpha [cm^{-1}]' )
% ylabel( 'Transmission Probability' )
% grid on
semilogy( alpha_vals, err_asymp, '-k' )
hold on
semilogy( alpha_vals, err_am, '--k' )
legend( 'Asymptotic', 'Atomic Mix' )
xlabel( '\alpha [cm^{-1}]' )
ylabel( 'Absolute Error' )
grid on
