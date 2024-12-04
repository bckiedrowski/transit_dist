function T = porous_asymp( s, t, alpha, beta, lambda )

  a = alpha;
  k = beta / alpha;

  b    = (2*lambda - 1)*0.5*sqrt((1+k)/(k*s));
  if ( b > 0 ) umin =  max( -sqrt(k*s), -1/b ); else umin = -sqrt(k*s); end;
  if ( b < 0 ) umax =  min(  sqrt(s),   -1/b ); else umax =  sqrt(s);   end; 
  sqa = sqrt(a);
  ef  = erf(sqa*umax) - erf(sqa*umin);
  ex  = ( exp(-a*umax^2) - exp(-a*umin^2) )/sqrt(pi*a);
  C   = 2.0/( ef - b * ex );
  a1  = 2*t*sqrt(k*s*(1+k))/(1+k)^2;
  a2  = t*(k^2 + 2*k*t*s - 1 )/(1+k)^3;
  c1  = b  - a1;
  c2  = a2 - a1*b;
  c3  = a2*b;
    
  Psi = @(u) (c2 + 2*a)/(4*a)*erf(sqrt(a)*u) ...
           - (a*(c1 + c2*u) + c3*( 1 + a*u.^2 ))/(2*sqrt(pi)*a^1.5).*exp(-a*u.^2);
  T = exp(-lambda*t*s) .* ( Psi(umax) - Psi(umin) );
end