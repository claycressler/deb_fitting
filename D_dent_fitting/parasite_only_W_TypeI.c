#include <R.h>
#include <stdio.h>

static double parms[4];

#define aP parms[0]
#define eP parms[1]
#define P0 parms[2]
#define Perr parms[3]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=4;
  odeparms(&N, parms);
}

/* derivatives */
/* this is for the model where parasites "eat" reserves */
void derivs (int *neq, double *t, double *y, double *ydot) {
  double P = y[0]; // parasites

  // what is the weight now? (based on cubic regression on real data)
  double age = *t;
  double W = 0.0018*pow(0.3189 + 0.1266*age - 0.003876*age*age + 0.00004259*age*age*age, 3.0);

  // balance the equations
  ydot[0] = eP * aP*W*P;
}

