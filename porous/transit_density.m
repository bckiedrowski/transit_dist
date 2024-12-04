function f = transit_density( l, s, alpha, beta, lambda )
  n = length(l);
  f = zeros(1,n);

  a = alpha;
  k = beta / alpha;
  for i=1:n
    li = l(i);
    %case for outside the range
    if ( li < 0 || li > s )
      f(i) = 0; continue;
    end
    z  = 2.0*a*sqrt(k*li.*(s-li));
    if ( z <= 700.0 )
      %use exact distribution
      %special cases at endpoints
      if ( li == 0 )
        f(i) = ( lambda*alpha        + ( 1 - lambda )*beta * ( 1 + alpha*s ) )*exp( -beta*s); continue;
      elseif ( li == s )
        f(i) = ( ( 1 - lambda )*beta + lambda*alpha        * ( 1 + beta*s  ) )*exp(-alpha*s); continue;          
      end
      %normal case in interior
      e0 = exp(-a*li - a*k*(s-li));
      c0 = lambda * alpha + ( 1 - lambda ) * beta;
      c1 = sqrt(alpha*beta) * ( (2*lambda - 1)*li + (1-lambda)*s )/sqrt(li*(s-li));
      I0 = besseli(0,z);
      I1 = besseli(1,z);
      
      f(i) = e0 * ( c0*I0 + c1*I1 );      
    else
      %use asymptotic form
      u    = sqrt(li) - sqrt(k*(s-li));
      dudl = ( sqrt(s-li) + sqrt(k*li) )/(2*sqrt(li*(s-li)));

      b    = (2*lambda - 1)*0.5*sqrt((1+k)/(k*s));
      %umin = max( -sqrt(k*s), -1/b );
      %umax = min(  sqrt(s),   -1/b );
      if ( b > 0 ) umin =  max( -sqrt(k*s), -1/b ); else umin = -sqrt(k*s); end;
      if ( b < 0 ) umax =  min(  sqrt(s),   -1/b ); else umax =  sqrt(s);   end;
      if ( u < umin || u > umax )
        f(i) = 0; continue; 
      end
    
      sqa = sqrt(a);
      ef  = erf(sqa*umax) - erf(sqa*umin);
      ex  = ( exp(-a*umax^2) - exp(-a*umin^2) )/sqrt(pi*a);
      C   = 2.0/( ef - b * ex );
    
      f(i) = C * ( 1 + b*u ) *  sqrt(a/pi) .* exp(-a*u^2) .* abs(dudl);
    end
  end
end

% return transit length distribution functions
% fA = function at (l,s) for transit length l starting in zone A
% fB = same but starting in zone B
% function [ fA fB ] = transit_density( l, s, alpha, beta )
%   n = length(l);
%   fA = zeros(1,n);
%   fB = zeros(1,n);
% 
%   a = alpha;
%   k = beta / alpha;
%   for i=1:n
%     li = l(i);
%     %case for outside the range
%     if ( li < 0 || li > s )
%       fA(i) = 0; fB(i) = 0; continue;
%     end
% 
%     z  = 2.0*a*sqrt(k*li.*(s-li));
%     if ( z <= 700.0 )
%       %use exact distribution
%       e0 = exp(-a*li - a*k*(s-li));
%       I0 = besseli(0,z);
%       I1 = besseli(1,z);
% 
%       if ( li < s )
%         fA(i) = e0 .* (   a*I0 + a.*sqrt(k*li./(s-li)).*I1 );
%       else
%         fA(i) = alpha * ( 1 + beta * s ) * exp(-alpha*s);
%       end
%       if ( li > 0 )
%         fB(i) = e0 .* ( k*a*I0 + a.*sqrt(k*(s-li)./li).*I1 );
%       else
%         fB(i) = beta * ( 1 + alpha * s ) * exp(-beta*s);
%       end
%     else
%       %use asymptotic form
%       u    = sqrt(li) - sqrt(k*(s-li));
%       dudl = ( sqrt(s-li) + sqrt(k*li) )/(2*sqrt(li*(s-li)));
%       b  = 0.5*sqrt((1+k)/(k*s));
%       ef = erf(sqrt(k*a*s)) + erf(sqrt(a*s));
%       ex = ( exp(-k*a*s) - exp(-a*s ) )/sqrt(pi*a);
%       CA = 2.0/( ef + b * ex );
%       CB = 2.0/( ef - b * ex );
%       g = sqrt(a/pi) .* exp(-a*u.^2);
% 
%       fA(i) = CA * ( 1 + b*u ) .* g .* dudl;
%       fB(i) = CB * ( 1 - b*u ) .* g .* dudl;
%     end
%   end
% end
