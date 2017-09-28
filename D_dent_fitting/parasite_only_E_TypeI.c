#include <R.h>
#include <stdio.h>

static double parms[4];

#define E parms[0]
#define aP parms[1]
#define eP parms[2]
#define P0 parms[3]
#define Perr parms[4]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=5;
  odeparms(&N, parms);
}

/* derivatives */
/* this is for the model where parasites "eat" reserves */
void derivs (int *neq, double *t, double *y, double *ydot) {
  double P = y[0]; // parasites

  // balance the equations
  ydot[0] = eP * aP*E*P;
}

