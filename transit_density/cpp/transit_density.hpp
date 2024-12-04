#ifndef _TRANSIT_LENGTH_DENSITY_
#define _TRANSIT_LENGTH_DENSITY_
#include <vector>
//case for evaluation at a single point
double transit_density( const double l, const double s, const double alpha,
                        const double beta, const double lambda );
//case for evaluation over a templated container of values of lengths for fixed parameters
template < typename T >
std::vector< double > transit_density( const T& l, const double s, const double alpha,
                                       const double beta, const double lambda )  {
  const auto n = l.size();
  std::vector< double > f( n, 0.0 );
  size_t k = 0;
  for ( auto li : l ) {
    f[k] = transit_density( li, s, alpha, beta, lambda ); ++k;
  }
  return f;
}
#endif
