#include <cmath>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "hypergeometric_pfq.h"
#include "struve.h"
#include "theo.h"

using namespace std;

/* Constants */

// Magnetic permeability of free space
const double mu_0 = 4.0e-7 * M_PI;


/* Helper functions not exposed */
inline double hyper(double x);
inline double halfspace_int_hyperg(double x1, double x2);
double halfspace_int_struve(double x1, double x2);
double borehole_int_struve(double x1, double x2);


/* Implementation of exposed functions */

gsl_complex z0(const coil_params* coil, const conductor_params* conductor, double frequency) {
  // Set aliases for arguments
  const double& mu_r = conductor->mu_r;
  const double& sigma = conductor->sigma;
  const double& r1 = coil->r1;
  const double& r2 = coil->r2;
  const double& z1 = coil->z1;
  const double& z2 = coil->z2;
  const double& w = coil->windings;

  // Angular frequency
  const double omega = 2.0 * M_PI * frequency;

  // Boundary
  const double b = 40 * r2;

  // Number of summation terms + one
  const int n = 160 + 1;

  // Precomputed quantities
  gsl_complex coefficient = gsl_complex_rect(0.0, omega * 2 * M_PI * mu_0 * gsl_pow_2(w) / gsl_pow_2(r2 - r1) / gsl_pow_2(z2 - z1));
  double q[n];
  gsl_complex p[n];
  double intSquared[n];
  double besJ[n];

  for (unsigned int i = 0; i < n; ++i) {
    q[i] = gsl_sf_bessel_zero_J1(i) / b;
    p[i] = gsl_complex_sqrt(gsl_complex_rect(gsl_pow_2(q[i]), omega * mu_0 * mu_r * sigma));
    intSquared[i] = gsl_pow_2(halfspace_int_hyperg(q[i] * r1, q[i] * r2));
    besJ[i] = intSquared[i] / gsl_pow_2(q[i] * b * gsl_sf_bessel_J0(q[i] * b)) / gsl_pow_5(q[i]);
  }

  double sum = 0;
  for (int i = 1; i < n; ++i) {
    sum += 2 * (q[i] * (z2 - z1) - 1 + exp(q[i] * (z1 - z2))) * besJ[i];
  }

  gsl_complex z0 = gsl_complex_mul_real(coefficient, sum);
  return z0;
}

gsl_complex halfspace_dz(const coil_params* coil, const conductor_params* conductor, double frequency) {
  // Set aliases for arguments
  const double& mu_r = conductor->mu_r;
  const double& sigma = conductor->sigma;
  const double& r1 = coil->r1;
  const double& r2 = coil->r2;
  const double& z1 = coil->z1;
  const double& z2 = coil->z2;
  const double& w = coil->windings;

  // Angular frequency
  const double omega = 2.0 * M_PI * frequency;

  // Boundary
  const double b = 40 * r2;

  // Number of summation terms plus one
  const int n = 160 + 1;

  // Precomputed quantities
  gsl_complex coefficient = gsl_complex_rect(0.0, omega * 2 * M_PI * mu_0 * gsl_pow_2(w) / gsl_pow_2(r2 - r1) / gsl_pow_2(z2 - z1));
  double q[n];
  gsl_complex p[n];
  double intSquared[n];
  double besJ[n];

  for (unsigned int i = 0; i < n; ++i) {
    q[i] = gsl_sf_bessel_zero_J1(i) / b;
    p[i] = gsl_complex_sqrt(gsl_complex_rect(gsl_pow_2(q[i]), omega * mu_0 * mu_r * sigma));

    // Integral
    intSquared[i] = gsl_pow_2(halfspace_int_hyperg(q[i] * r1, q[i] * r2));
    //intSquared[i] = gsl_pow_2(halfspace_int_struve(q[i] * r1, q[i] * r2));

    besJ[i] = intSquared[i] / gsl_pow_2(q[i] * b * gsl_sf_bessel_J0(q[i] * b)) / gsl_pow_5(q[i]);
  }

  gsl_complex sum = gsl_complex_rect(0.0, 0.0);
  for (int i = 1; i < n; ++i) {
    double exponents = gsl_pow_2((exp(-q[i] * z1) - exp(-q[i] * z2)));
    gsl_complex q_minus_p = gsl_complex_sub(gsl_complex_rect(q[i] * mu_r, 0.0), p[i]);
    gsl_complex q_plus_p = gsl_complex_add_real(p[i], q[i] * mu_r);
    gsl_complex qp = gsl_complex_div(q_minus_p, q_plus_p);
    sum = gsl_complex_add(sum, gsl_complex_mul_real(qp, exponents * besJ[i]));
  }

  gsl_complex dz = gsl_complex_mul(coefficient, sum);
  return dz;
}

gsl_complex borehole_dz(const coil_params* coil, const conductor_params* conductor, double frequency) {
/*
  // Set aliases for arguments
  const double& mu_r = conductor->mu_r;
  const double& sigma = conductor->sigma;
  const double& h = conductor->h;
  const double& b = conductor->r;
  const double& r1 = coil->r1;
  const double& r2 = coil->r2;
  const double& z1 = coil->z1;
  const double& z2 = coil->z2;
  const double& w = coil->windings;

  // Angular frequency
  const double omega = 2.0 * M_PI * frequency;

  // Number of summation terms plus one
  const int n = 160 + 1;

  // Precomputed quantities
  gsl_complex coefficient = gsl_complex_rect(0.0, omega * 4 * M_PI * mu_0 * gsl_pow_2(w) / (h * gsl_pow_2(r2 - r1) * gsl_pow_2(z2 - z1)));
  double q[n];
  gsl_complex p[n];
  double intSquared[n];
  double sinSquared[n];
  double besselNom[n];
  double besselDenom[n];

  for (unsigned int i = 0; i < n; ++i) {
    q[i] = gsl_sf_bessel_zero_J1(i) / b;
    p[i] = gsl_complex_sqrt(gsl_complex_rect(gsl_pow_2(q[i]), omega * mu_0 * mu_r * sigma));

    intSquared[i] = gsl_pow_2(borehole_int_struve(q[i] * r1, q[i] * r2));
    sinSquared[i] = gsl_pow_2(sin(q[i] * z1) - sin(q[i] * z2)) / gsl_pow_6(q[i]);

    double K0_q = gsl_sf_bessel_K0(q[i] * b);
    double K0_p = gsl_sf_bessel_K0(p[i] * b);
    double K1_q = gsl_sf_bessel_K1(q[i] * b);
    double K1_p = gsl_sf_bessel_K1(p[i] * b);
    // TODO: need Bessel functions K0, K1 that take complex numbers as arguments
  }
*/

  gsl_complex dummy = gsl_complex_rect(0.0, 0.0);
  return dummy;
}


/* Implementation of helper function */

inline double hyper(double x) {
  return pow(x, 3.0) * hypergeometric_1f2(1.5, 2.5, 2.0, -pow(x, 2.0) / 4.0) / 6.0;
}

// Int(x1, x2) for half-space using hypergeometric function
inline double halfspace_int_hyperg(double x1, double x2) {
  return hyper(x2) - hyper(x1);
}

// Int(x1, x2) for half-space using Struve functions
double halfspace_int_struve(double x1, double x2) {
  double H0_x1 = struve_Hn(0, x1);
  double H0_x2 = struve_Hn(0, x2);
  double H1_x1 = struve_Hn(1, x1);
  double H1_x2 = struve_Hn(1, x2);

  double J0_x1 = gsl_sf_bessel_J0(x1);
  double J0_x2 = gsl_sf_bessel_J0(x2);
  double J1_x1 = gsl_sf_bessel_J1(x1);
  double J1_x2 = gsl_sf_bessel_J1(x2);

  return 0.5 * M_PI * (-x1 * H0_x1 * J1_x1 + x1 * H1_x1 * J0_x1 + x2 * H0_x2 * J1_x2 - x2 * H1_x2 * J0_x2);
}

// Int(x1, x2) for borehole using Struve functions
double borehole_int_struve(double x1, double x2) {
  double L0_x1 = struve_Ln(0, x1);
  double L0_x2 = struve_Ln(0, x2);
  double L1_x1 = struve_Ln(1, x1);
  double L1_x2 = struve_Ln(1, x2);

  double I0_x1 = gsl_sf_bessel_I0(x1);
  double I0_x2 = gsl_sf_bessel_I0(x2);
  double I1_x1 = gsl_sf_bessel_I1(x1);
  double I1_x2 = gsl_sf_bessel_I1(x2);

  return 0.5 * M_PI * (-x1 * L0_x1 * I1_x1 + x1 * L1_x1 * I0_x1 + x2 * L0_x2 * I1_x2 - x2 * L1_x2 * I0_x2);
}
