#include <R.h>
#include <stdio.h>

static double parms[5];

#define E parms[0]
#define aP parms[1]
#define eP parms[2]
#define hP parms[3]
#define P0 parms[4]
#define Perr parms[5]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=6;
  odeparms(&N, parms);
}

/* derivatives */
/* this is for the model where parasites "eat" reserves */
void derivs (int *neq, double *t, double *y, double *ydot) {
  double P0 = y[0]; // immature parasites
  double P1 = y[1]; // mature parasites


  // balance the equations
  ydot[0] = eP * aP*E*P/(hP+P);
}

