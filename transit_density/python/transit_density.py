def transit_density( l, s, alpha, beta, lmbda ):
  import math
  import scipy.special

  f = []
  k = beta / alpha
  i = 0
  for li in l:
    if li < 0 or li > s:
      f.append( 0.0 )
      continue
    z = 2.0*alpha*math.sqrt(k*li*(s-li))
    if z < 700:
      if li == 0:
        f.append( ( lmbda*alpha        + ( 1 - lmbda )*beta * ( 1 + alpha*s ) )*math.exp( -beta*s) )
        continue
      elif li == s:
        f.append( ( ( 1 - lmbda )*beta + lmbda*alpha        * ( 1 + beta*s  ) )*math.exp(-alpha*s) )
        continue
      e0 = math.exp(-alpha*li - alpha*k*(s-li))
      c0 = lmbda * alpha + ( 1 - lmbda ) * beta
      c1 = math.sqrt(alpha*beta) * ( (2*lmbda - 1)*li + (1-lmbda)*s )/math.sqrt(li*(s-li))
      I0 = scipy.special.i0( z )
      I1 = scipy.special.i1( z )
      
      f.append( e0 * ( c0*I0 + c1*I1 ) )
    else:
      #use asymptotic form
      pi   = math.pi
      u    = math.sqrt(li) - math.sqrt(k*(s-li))
      dudl = ( math.sqrt(s-li) + math.sqrt(k*li) )/(2*math.sqrt(li*(s-li)))
      b    = (2*lmbda - 1)*0.5*math.sqrt((1+k)/(k*s))
      umin = max( -math.sqrt(k*s), -1/b ) if b > 0 else -math.sqrt(k*s)
      umax = min(  math.sqrt(s),   -1/b ) if b < 0 else  math.sqrt(s)
      if ( u < umin or u > umax ):
        f.append( 0.0 )
        continue
      sqa = math.sqrt(alpha);
      ef  = math.erf(sqa*umax) - math.erf(sqa*umin)
      ex  = ( math.exp(-alpha*umax**2) - math.exp(-alpha*umin**2) )/math.sqrt(pi*alpha)
      C   = 2.0/( ef - b * ex )
    
      f.append( C * ( 1 + b*u ) *  math.sqrt(alpha/pi) * math.exp(-alpha*u**2) * abs(dudl) )
  return f
