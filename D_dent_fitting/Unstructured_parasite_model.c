#include <R.h>
#include <stdio.h>

static double parms[7];

#define theta parms[0]
#define r parms[1]
#define aP parms[2]
#define hP parms[3]
#define b parms[4]
#define P0 parms[5]
#define obs_sd parms[6]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=7;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {
  // balance the equations
  ydot[0] = theta - r*y[0] - aP*y[0]*y[1]/(hP+y[0]);
  ydot[1] = b*aP*y[0]*y[1]/(hP+y[0]);

}

