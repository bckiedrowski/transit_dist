function T = stopping_asymptotic( range, s, alpha, beta, lambda )
    a    = alpha;
    k    = beta / alpha;
    b    = (2*lambda - 1)*0.5*sqrt((1+k)/(k*s));
    umin = max( -sqrt(k*s), -1/abs(b) );
    umax = min(  sqrt(s),    1/abs(b) );
    
    sqa = sqrt(a);
    ef  = erf(sqa*umax) - erf(sqa*umin);
    ex  = ( exp(-a*umax^2) - exp(-a*umin^2) )/sqrt(pi*a);
    C   = 2.0/( ef - b * ex );
    
    ustar = min( sqrt(range) - sqrt(k*(s-range)), 1/abs(b) );
    
    T = 0.5 * C * ( erf(sqa*ustar) - erf(sqa*umin) - b/sqrt(pi*a) * ( exp(-a*ustar^2) - exp(-a*umin^2) ) );

 %   T = 0.5 * ( lam*CA + (1-lam)*CB ) * ( erf(sqrt(k*a*s)) + erf(sqrt(a)*umax) ) + ...
 %       0.5 * ( lam*CA - (1-lam)*CB ) *  b * ( exp(-k*a*s) - exp(-a*umax^2) )/sqrt(pi*a);
end