/* This code simulates the standard DEB model defined in terms of reserve density W=E/V and structural length L=(V/xi)^(1/q) instead of E and V */
#include <R.h>
#include <stdio.h>

static double parms[17];

#define Imax parms[0]
#define fh parms[1]
#define g parms[2]
#define rho parms[3]
#define eps parms[4]
#define Vol parms[5]
#define F0 parms[6]
#define xi parms[7]
#define q parms[8]
#define K parms[9]
#define km parms[10]
#define ER parms[11]
#define v parms[12]
#define EG parms[13]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=17;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {
  // state variables
  double F = y[0]; // algal concentration
  double W = y[1]; // reserve density
  double L = y[2]; // structural length

  // ingestion
  double ing = Imax * F/(fh+F) * pow(L,g);

  // balance the equations
  ydot[0] = -ing;
  ydot[1] = rho*eps*Vol*ing/(xi*pow(L,q)) - W*v/L;
  ydot[2] = L/q * (K/EG * W * v/L - km) / (1 + K/EG * W);
  ydot[3] = (1-K)/ER * W * xi * pow(L,q) * (v/L + km)/(1 + K/EG * W);

}

