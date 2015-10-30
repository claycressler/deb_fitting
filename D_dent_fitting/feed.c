#include <R.h>
#include <stdio.h>

static double parms[4];

#define Imax parms[0]
#define g parms[1]
#define fh parms[2]
#define L parms[3]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=4;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {
  // state variables
  double F = y[0];

  // balance the equations
  ydot[0] = -Imax * F / (fh + F) * pow(L, g);

}

