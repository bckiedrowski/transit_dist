#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <limits>

#include "transit_density.hpp"

static double eps = 1000.0*std::numeric_limits<double>::epsilon();

static std::mt19937_64 rng;
static std::uniform_int_distribution<uint64_t> uint_dist;
double urand() {
  return (double) uint_dist(rng)/UINT64_MAX;
}

int main() {

  constexpr auto nsamples = 1e9;

  constexpr auto L   = 2.0;
  constexpr auto lam = 1.0;
  constexpr auto mu  = 2.0;

  constexpr auto omega = 2.0/3.0;

/*
  auto psi1 = []( const double x ){
    const auto z = 2.0*sqrt(lam*mu*x*(L-x));
    return exp(-lam*x - mu*(L-x))
         * ( lam*std::cyl_bessel_i(0,z) + sqrt(lam*mu*x/(L-x))*bessi1(z) );
  };
  auto psi2 = []( const double x ){
    const auto z = 2.0*sqrt(lam*mu*x*(L-x));
    return exp(-lam*x - mu*(L-x))
         * ( mu*bessi0(z) + sqrt(lam*mu*(L-x)/x)*bessi1(z) );
  };
*/

  constexpr std::array< double, 2 > r = { lam, mu };

  constexpr auto nbins = 100;
  std::array< double, nbins > t1;
  t1.fill( 0.0 );

  auto fA = 0.0;
  auto fB = 0.0;

  for ( auto sample = 0; sample < nsamples ; ++sample ) {
    auto t = 0.0;
    auto x = 0.0;
    auto k = urand() < omega ? 0 : 1;
    auto m = 0;

    while ( x < L ) {
      auto s = -log(urand())/r[k];
      s = x + s < L ? s : L-x;
      x += s + eps;

      if ( k == 0 ) { t += s; }

      k = (k+1) % 2;
      ++m;
    }
    k = floor( t/L * nbins );
    if ( m == 1 ) {
      if ( k == 0 ) fA += 1.0;
      else          fB += 1.0;
    }
    if ( k < nbins && m > 1 ) { t1[k] += 1.0; }
  }

  constexpr auto dx = L/nbins;
  auto x = 0.5*dx;
  for ( auto i = 0 ; i < nbins ; ++i ) {
    const auto calc = t1[i]/nsamples;
    const auto ref = transit_density( x, L, lam, mu, omega );
    std::cout << std::fixed << std::right << std::setw(5) << std::setprecision(2) << x << "  "
      << std::setw(10) << std::setprecision(6) << ref << "   " 
      << calc/dx << " +/- " << std::scientific << std::sqrt( calc * ( 1 - calc ) / nsamples )/dx << "   " 
      << std::scientific << ( calc/dx - ref ) << '\n';
    x += dx;
  }
  fA /= nsamples; fB /= nsamples;
  std::cout << '\n' << ( 1 - omega ) * std::exp(-mu* L) << "   " << fA << " +/- " << std::scientific << std::sqrt( fA*(1-fA)/nsamples ) 
            << '\n' << omega         * std::exp(-lam*L) << "   " << fB << " +/- " << std::scientific << std::sqrt( fB*(1-fB)/nsamples ) << '\n';
  return 0;
}
