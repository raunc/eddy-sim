#include <cmath>
#include <gsl/gsl_sf_gamma.h>

#include "struve.h"

using namespace std;

double struve_Hn(unsigned int n, double x, double precision) {
  double prev_sum = 0;
  double sum = 0;
  int m = 0;

  do {
    prev_sum = sum;
    sum += pow(-1.0, m) / (gsl_sf_gamma(1.5 + m) * gsl_sf_gamma(1.5 + m + n)) * pow(x / 2.0, 2.0 * m + n + 1);

    ++m;
  } while (abs(sum - prev_sum) > precision);

  return sum;
}

double struve_Ln(unsigned int n, double x, double precision) {
  double prev_sum = 0;
  double sum = 0;
  int m = 0;
  double coefficient = pow(x / 2, n + 1);

  precision /= coefficient;

  do {
    prev_sum = sum;
    sum += pow(x / 2, 2 * m) / (gsl_sf_gamma(1.5 + m) * gsl_sf_gamma(1.5 + m + n));

    ++m;
  } while (abs(sum - prev_sum) > precision);

  return coefficient * sum;
}
