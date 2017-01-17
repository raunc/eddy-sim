#include <cstdio>
#include <gsl/gsl_complex.h>
#include "lib/theo.h"

using namespace std;

// Constants
const double sigma = 16.45e6;
const double mu_r = 1.0;

double f = 1.0e3;// Hz
double omega = 2.0 * M_PI * f;

// Spatial parameters
const double r1 = 3.015e-3;
const double r2 = 5.46e-3;
const double z1 = 1.32e-3; // Distance between conducting plate and the coil
const double z2 = 4.26e-3;

const double h = 40 * r2; // Boundary from the coil axis
const int Ns = 160;
const double WT = 900.0;// Number of windings

void print_complex(gsl_complex z) {
  printf("%g%+gi", z.dat[0], z.dat[1]);
}

int main(void) {

  // Define coil
  coil_params coil;
  coil.r1 = 3.015e-3;
  coil.r2 = 5.46e-3;
  coil.z1 = 1.32e-3;
  coil.z2 = 4.26e-3;
  coil.windings = 900;

  // Define conducting space
  conductor_params conductor;
  conductor.sigma = 16.45e6;
  conductor.mu_r = 1.0;

  double frequencies[] = {1e3, 10e3, 100e3};

  for (int i = 0; i < 3; ++i) {
    printf("f = %g Hz", frequencies[i]);
    printf("  dz = ");
    print_complex(halfspace_dz(&coil, &conductor, frequencies[i]));
    printf("\n");
  }

  return 0;
}
