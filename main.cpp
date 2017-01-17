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
  coil_params coil_c9;
  coil_c9.r1 = 3.015e-3;
  coil_c9.r2 = 5.46e-3;
  coil_c9.z1 = 1.32e-3;
  coil_c9.z2 = 4.26e-3;
  coil_c9.windings = 900;

  coil_params coil_c20;
  coil_c20.r1 = 1.56e-3;
  coil_c20.r2 = 1.83e-3;
  coil_c20.z1 = 0.68e-3;
  coil_c20.z2 = 4.73e-3;
  coil_c20.windings = 121;

  // Define conducting space
  conductor_params conductor;
  conductor.sigma = 16.45e6;
  conductor.mu_r = 1.0;

  double frequencies[] = {1e3, 10e3, 100e3};

  printf("Coil C9\n");
  for (int i = 0; i < 3; ++i) {
    printf("f = %g Hz", frequencies[i]);
    printf("  dz = ");
    print_complex(halfspace_dz(&coil_c9, &conductor, frequencies[i]));
    printf("\n");
  }

  printf("Coil C20\n");
  for (int i = 0; i < 3; ++i) {
    printf("f = %g Hz", frequencies[i]);
    printf("  dz = ");
    print_complex(halfspace_dz(&coil_c20, &conductor, frequencies[i]));
    printf("\n");
  }

  return 0;
}
