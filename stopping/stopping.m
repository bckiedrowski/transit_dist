clear;

%stopping power
rho_fe = 7.874;
fe_data = load('./fe_stopping.dat');
Edata  = fe_data(:,2);
S = @(E) interp1( Edata, rho_fe * fe_data(:,5), E, 'linear', 'extrap' );

N = 501;
E0 = 100.0;
Evals = logspace(-3,2,N);
range = integral( @(E) 1./S(E), Evals(1), E0 );
l = zeros(1,N);
for i=1:N
  l(i) = integral( @(E) 1./S(E), Evals(i), E0 );  
end
lE = @(E) interp1( Evals, l, E, 'linear', 'extrap' );

k      = 0.25;
lambda = k/(1+k);

%trasmission (beam source)
s_vals     = [ 1 2 4 5 10 20 ] * range;
alpha_vals = linspace( 0.1, 20, 200 ) / range;
T    = zeros(length(alpha_vals),length(s_vals));

for j=1:length(s_vals)
  s = s_vals(j);
  l = linspace(0,s,1001);
  for i=1:length(alpha_vals)
    alpha  = alpha_vals(i);
    beta   = k * alpha;
 
%     [ fA fB ] = transit_density( l, s, alpha, beta );
%     fA = @(x) interp1( l, fA, x, 'linear', 0.0 );
%     fB = @(x) interp1( l, fB, x, 'linear', 0.0 );
%     f  = @(x) lambda * fA(x) + ( 1 - lambda ) * fB(x);
    f = @(x) transit_density( x, s, alpha, beta, lambda );
  
    T(i,j) = integral( @(x) f(x), 0, range ) + (1-lambda)*exp(-beta*s) ...
           + ( s < range ) * lambda*exp(-alpha*s);
  end
end
%plotting
clf;
plot( alpha_vals*range, T(:,1), '-k' );
hold on
plot( alpha_vals*range, T(:,2), '--k' );
plot( alpha_vals*range, T(:,3), '-.k' );
legend( 's = R', 's = 2R', 's = 4R' )
grid on;
xlabel('\alphaR')
ylabel('Transmission Probability')

%trasmission (beam source - exact versus asymptotic)
s_vals     = [ 2 10 ] * range;
alpha_vals = linspace( 0.1, 20, 200 ) / range;
T_exact    = zeros(length(alpha_vals),length(s_vals));
T_asymp    = zeros(length(alpha_vals),length(s_vals));
for j=1:length(s_vals)
  s = s_vals(j);
  l = linspace(0,s,1001);
  for i=1:length(alpha_vals)
    alpha  = alpha_vals(i);
    beta   = k * alpha;

    f = @(x) transit_density( x, s, alpha, beta, lambda );

    T_exact(i,j) = integral( @(x) f(x), 0, range ) + (1-lambda)*exp(-beta*s) ...
           + ( s < range ) * lambda*exp(-alpha*s);
    T_asymp(i,j) = stopping_asymptotic( range, s, alpha, beta, lambda );    
  end
end
clf
semilogy( alpha_vals * range, abs( T_exact(:,1) - T_asymp(:,1) )./T_exact(:,1), '-k' )
hold on
semilogy( alpha_vals * range, abs( T_exact(:,2) - T_asymp(:,2) )./T_exact(:,2), '--k' )
legend( 's = 2R', 's = 10R' )
grid on;
xlabel('\alphaR')
ylabel('Absolute Relative Error of Transmission Probability')

%energy spectrum (beam source)
s = 10*range;
alpha_vals = [ 1 2 5 10 ] / range;
spectrum = zeros(length(alpha_vals),length(Evals));
for i=1:length(alpha_vals)
  alpha = alpha_vals(i) / range;
  beta  = 0.25 * alpha;
  
  f = @(x) transit_density( x, s, alpha, beta, lambda );
  
  fE = @(E) f(lE(E))./S(E);
  %[ integral( @(E) fE(E), 0, E0 ) (1-lambda)*exp(-beta*s) ]
  spectrum(i,:) = fE(Evals)';
end
clf
semilogx( Evals, spectrum(1,:), '-k' )
hold on
semilogx( Evals, spectrum(2,:), '--b' )
semilogx( Evals, spectrum(3,:), '-.k' )
semilogx( Evals, spectrum(4,:), '-r'  )
grid on
xlabel('Transmitted Energy (MeV)')
ylabel('Tramission Spectrum (MeV^{-1})')
legend('\alphaR = 1', '\alphaR = 2', '\alphaR = 5', '\alphaR = 10')

%angular flux  (isotropic flux)
s   = 2*range;
mu  = linspace(0.01,1,100);
psi = zeros(length(mu),length(alpha_vals));

for j=1:length(alpha_vals)
  alpha = alpha_vals(j);
  beta  = 0.25 * alpha;
  for i=1:length(mu)
    f = @(x) transit_density( x, s/mu(i), alpha, beta, lambda );

    psi(i,j) = integral( @(x) f(x), 0, range ) + (1-lambda)*exp(-beta*s/mu(i)) ...
             + ( s/mu(i) < range ) * lambda*exp(-alpha*s/mu(i));
  end
  psi(:,j) = 2*psi(:,j);
  J(j) = trapz( mu, mu'.*psi(:,j) );
end
clf;
hold on
plot( mu, mu'.*psi(:,1), '-k' )
plot( mu, mu'.*psi(:,2), '--b' )
plot( mu, mu'.*psi(:,3), '-.k' )
plot( mu, mu'.*psi(:,4), '-r' )
grid on
xlim( [ 0 1 ] )
xlabel('\mu')
ylabel('Transmitted Angular Flux [cm^{-2}]')
legend('\alphaR = 1', '\alphaR = 2', '\alphaR = 5', '\alphaR = 10', 'Location', 'NorthWest')
%Jout = trapz( mu, mu.*psi )

%partial current (isotropic flux)
alpha_fine = linspace( 0.1, 60, 600 ) / range;
s   = 2*range;
mu  = linspace(0.01,1,100);
psi_exact = zeros(length(mu),length(alpha_fine));
psi_asymp = zeros(length(mu),length(alpha_fine));
J_exact   = zeros(1,length(alpha_fine));
J_asymp   = zeros(1,length(alpha_fine));

for j=1:length(alpha_fine)
  alpha = alpha_fine(j);
  beta  = 0.25 * alpha;
  for i=1:length(mu)
    f = @(x) transit_density( x, s/mu(i), alpha, beta, lambda );

    psi_exact(i,j) = integral( @(x) f(x), 0, range ) + (1-lambda)*exp(-beta*s/mu(i)) ...
             + ( s/mu(i) < range ) * lambda*exp(-alpha*s/mu(i));

    psi_asymp(i,j) = stopping_asymptotic( range, s/mu(i), alpha, beta, lambda );
  end
  psi_exact(:,j) = 2*psi_exact(:,j);
  psi_asymp(:,j) = 2*psi_asymp(:,j);
  J_exact(j) = trapz( mu, mu'.*psi_exact(:,j) );
  J_asymp(j) = trapz( mu, mu'.*psi_asymp(:,j) );
end
clf
plot( alpha_fine * range, J_exact, '-k' )
hold on
plot( alpha_fine * range, J_asymp, '--k' )
grid on
%ylim( [0.6 1.0] )
legend('Exact','Asymptotic')
xlabel('\alphaR')
ylabel('Transmitted Partial Current')
%Jout = trapz( mu, mu.*psi )

% lam = k/(1+k); 
% b  = 0.5*sqrt((1+k)/(k*s));
% ef = erf(sqrt(k*a*s)) + erf(sqrt(a*s));
% ex = ( exp(-k*a*s) - exp(-a*s ) )/sqrt(pi*a);
% CA = 2.0/( ef + b * ex );
% CB = 2.0/( ef - b * ex );
% f
% umin = -sqrt(k*s);
% umax =  sqrt(L) - sqrt(k*(s-L));
% 
% IA = integral( @(u) lam     * CA * ( 1 + b*u ) .* sqrt(a/pi) .* exp(-a*u.^2), umin, umax );
% IB = integral( @(u) (1-lam) * CB * ( 1 - b*u ) .* sqrt(a/pi) .* exp(-a*u.^2), umin, umax );
% 
% RA = 0.5 * CA * ( erf(sqrt(k*a*s)) + erf(sqrt(a)*umax) ...
%     + b * ( exp(-k*a*s) - exp(-a*umax^2) )/sqrt(pi*a) );
% RB = 0.5 * CB * ( erf(sqrt(k*a*s)) + erf(sqrt(a)*umax) ...
%     - b * ( exp(-k*a*s) - exp(-a*umax^2) )/sqrt(pi*a) );
% 
% lam*RA + (1-lam)*RB
% 
% 0.5 * ( lam*CA + (1-lam)*CB ) * ( erf(sqrt(k*a*s)) + erf(sqrt(a)*umax) ) + ...
% 0.5 * ( lam*CA - (1-lam)*CB ) *  b * ( exp(-k*a*s) - exp(-a*umax^2) )/sqrt(pi*a)  


