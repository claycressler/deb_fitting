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
  double E = y[1]; // energy reserves
  double V = y[2]; // structural weight
  double L = pow(V/xi, 1/q); // structural length

  // ingestion
  double ing = Imax * F/(fh+F) * pow(L,g);
  // mobilization
  double pc = E * (v/L + km) / (1 + K*E/EG/V);

  // balance the equations
  ydot[0] = -ing;
  ydot[1] = rho*eps*Vol*ing - pc;
  ydot[2] = K*pc/EG - km*V;
  ydot[3] = (1-K)*pc/ER;

}

