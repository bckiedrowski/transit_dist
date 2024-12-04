clear;
clf;

N = 499;

s = 2;
lambda = 2/3;
alpha_vals = linspace( 1, 100, N );
L1_error = zeros(1,N);
for i=1:length(alpha_vals)
   alpha = alpha_vals(i);
   beta  = 5*alpha;
   
   f_exact = @(l) transit_density(    l, s, alpha, beta, lambda );
   f_asymp = @(l) transit_asymptotic( l, s, alpha, beta, lambda );
   
   L1_error(i) = integral( @(l) abs( f_exact(l) - f_asymp(l) ), 0, s );
end
loglog( alpha_vals, L1_error, '-k' )
grid on
xlabel('\alpha [cm^{-1}]')
ylabel('L-1 Norm')