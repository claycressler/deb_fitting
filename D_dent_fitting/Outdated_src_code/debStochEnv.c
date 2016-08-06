#include <R.h>
#include <stdio.h>

static double parms[16];

#define Imax parms[0]
#define Fh parms[1]
#define g parms[2]
#define rho parms[3]
#define eps parms[4]
#define V parms[5]
#define F0 parms[6]
#define K parms[9]
#define km parms[10]
#define v parms[12]
#define ER parms[11]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=16;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {
  // state variables
  double F = y[0]; // algal concentration
  double E = y[1]; // energy reserves
  double W = y[2]; // structural weight
  double L = pow(W, 1/3); // structural length

  // ingestion
  double ing = Imax * F/(Fh+F) * pow(L,g);
  // mobilization
  double pc = E * (v/L + km) / (1 + K*E/W);

  // balance the equations
  ydot[0] = -ing;
  ydot[1] = rho*eps*V*ing - pc;
  ydot[2] = K*pc - km*W;
  ydot[3] = (1-K)*pc/ER;

}

