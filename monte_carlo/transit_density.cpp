#include <cmath>
#include <vector>

//case for evaluation at a single point
double transit_density( const double l, const double s, const double alpha,
                        const double beta, const double lambda ) {
  //case for outside the range
  if ( l < 0 || l > s ) return 0.0;
 
  const auto k = beta / alpha;  
  const auto z = 2.0*alpha*std::sqrt(k*l*(s-l));
  if ( z <= 700.0 ) {
    //use exact distribution
    //special cases at endpoints
    if ( l == 0 ) {
      return ( lambda*alpha        + ( 1 - lambda )*beta * ( 1 + alpha*s ) )*std::exp( -beta*s);
    }
    else if ( l == s ) {
      return ( ( 1 - lambda )*beta + lambda*alpha        * ( 1 + beta*s  ) )*std::exp(-alpha*s);       
    }
    //normal case in interior
    const auto e0 = std::exp(-alpha*l - alpha*k*(s-l));
    const auto c0 = lambda * alpha + ( 1 - lambda ) * beta;
    const auto c1 = std::sqrt(alpha*beta) * ( (2*lambda - 1)*l + (1-lambda)*s )/sqrt(l*(s-l));
    const auto I0 = std::cyl_bessel_i(0,z);
    const auto I1 = std::cyl_bessel_i(1,z);
    
    return e0 * ( c0*I0 + c1*I1 );      
  }
  else {
    //use asymptotic form
    constexpr auto pi = M_PI;
    const auto u    = std::sqrt(l) - std::sqrt(k*(s-l));
    const auto dudl = ( std::sqrt(s-l) + std::sqrt(k*l) )/(2*std::sqrt(l*(s-l)));
    const auto b    = (2*lambda - 1)*0.5*sqrt((1+k)/(k*s));
    //establish non-negative range and return zero if outside range
    const auto umin = b > 0 ? std::max( -std::sqrt(k*s), -1/b ) : -std::sqrt(k*s);
    const auto umax = b < 0 ? std::min(  std::sqrt(s),   -1/b ) :  std::sqrt(s);
    if ( u < umin || u > umax ) return 0.0;
  
    const auto sqa = std::sqrt(alpha);
    const auto ef  = std::erf(sqa*umax) - std::erf(sqa*umin);
    const auto ex  = ( std::exp(-alpha*umax*umax) - std::exp(-alpha*umin*umin) )/std::sqrt(pi*alpha);
    const auto C   = 2.0/( ef - b * ex );
  
    return C * ( 1 + b*u ) * std::sqrt(alpha/pi) * std::exp(-alpha*u*u) * std::fabs(dudl);
  }
}
