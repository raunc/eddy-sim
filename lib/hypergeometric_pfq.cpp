#include <cmath>

#include <arb.h>
#include <arb_hypgeom.h>
#include <arf.h>

#include "hypergeometric_pfq.h"

double hypergeometric_1f2(double a1, double b1, double b2, double z) {
  // Initialize arb type variables
  arb_t result;
  arb_struct a[1];
  arb_struct b[2];
  arb_t z_arb;
  fmpz_t man;
  fmpz_t exp;

  arb_init(result);
  arb_init(&a[0]);
  arb_init(&b[0]);
  arb_init(&b[1]);
  arb_init(z_arb);
  fmpz_init(man);
  fmpz_init(exp);

  // Set values to arb type variables
  arb_set_d(&a[0], 1.5);
  arb_set_d(&b[0], 2.5);
  arb_set_d(&b[1], 2.0);
  arb_set_d(z_arb, z);

  // Calculate result
  arb_hypgeom_pfq(result, a, 1, b, 2, z_arb, 0, 32);

  // Convert result to double
  arf_get_fmpz_2exp(man, exp, &result->mid);
  double result_double = (double)(*man) * pow(2.0, (double)(*exp));

  // Free memory
  arb_clear(result);
  arb_clear(&a[0]);
  arb_clear(&b[0]);
  arb_clear(&b[1]);
  arb_clear(z_arb);
  fmpz_clear(man);
  fmpz_clear(exp);

  return result_double;
}
