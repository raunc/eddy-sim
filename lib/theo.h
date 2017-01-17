#ifndef EDDY_SIM_THEO_H
#define EDDY_SIM_THEO_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_math.h>

struct coil_params {
  double r1;
  double r2;
  double z1;
  double z2;
  double windings;
};

struct conductor_params {
  double mu_r;
  double sigma;
  double h;
  double r;
};

gsl_complex z0(const coil_params*, const conductor_params*, double frequency);
gsl_complex halfspace_dz(const coil_params*, const conductor_params*, double frequency);
gsl_complex borehole_dz(const coil_params*, const conductor_params*, double frequency);

#endif //EDDY_SIM_THEO_H
