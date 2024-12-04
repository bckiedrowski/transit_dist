clear;

alpha_vals = linspace( 1, 200, 2000 );
param_vals = linspace( 0.1, 10, 100 );
n = length(alpha_vals);

T_exact = zeros(1,n);
T_am    = zeros(1,n);
T_asymp = zeros(1,n);

for i=1:n
  alpha  = alpha_vals(i);
  k      = 0.2; %param_vals(i);
  beta   = k*alpha;
  lambda = k/(1+k);

  X = 1;
  t = 6; %param_vals(i); %1;
  
  mu = linspace( 0.01, 1, 100 );
  m  = length(mu);
  psi_exact = zeros(1,m);
  psi_am    = zeros(1,m);
  psi_asymp = zeros(1,m);
  for j=1:m
    s = X/mu(j);
    psi_exact(j) = integral( @(l) exp(-t*l) .* transit_density( l, s, alpha, beta, lambda ), 0, s ) ...
             + (1-lambda)*exp(-beta*s) + lambda*exp(-alpha*s)*exp(-t*s);
    psi_am(j)    = exp(-lambda*t*s);
    psi_asymp(j) = porous_asymp( s, t, alpha, beta, lambda );
  end

  T_exact(i) = trapz( mu, 2*mu.*psi_exact );
  T_am(i)    = trapz( mu, 2*mu.*psi_am    );
  T_asymp(i) = trapz( mu, 2*mu.*psi_asymp );
end

err_am    = abs( T_exact - T_am    );
err_asymp = abs( T_exact - T_asymp );
clf;
% plot( alpha_vals, T_exact, '-k' );
% hold on
% plot( alpha_vals, T_asymp, '--k' );
% plot( alpha_vals, T_am, '-.k' );
% legend( 'Exact', 'Asymptotic', 'Atomic Mix' )
% ylim( [0 1.0] )
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

