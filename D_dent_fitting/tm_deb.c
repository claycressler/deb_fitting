#include <R.h>
#include <stdio.h>

static double parms[12];

#define Imax parms[0]
#define Fh parms[1]
#define g parms[2]
#define rho parms[3]
#define K parms[4]
#define km parms[5]
#define v parms[6]
#define F0 parms[7]
#define E0 parms[8]
#define W0 parms[9]
#define Lerr parms[10]
#define Rerr parms[11]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=12;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {
  // fixed parameters whose values will not ever vary
  double eps = 0.0000000445; // carbon content of algae
  double V = 30; // volume of the container
  double xi = 0.00262; // length-weight regression coefficient
  double q = 2.4; // length-weight regression exponent

  // state variables
  double F = y[0]; // algal concentration
  double E = y[1]; // energy reserves
  double W = y[2]; // structural weight
  double L = pow(W, 1/3); // structural length
  double Lobs = pow((W+E)/xi, 1/q); //observed length

  // ingestion
  double ing = Imax * F/(Fh+F) * pow(Lobs,g);
  // mobilization
  double pc = E * (v/L + km) / (1 + K*E/W);

  // balance the equations
  ydot[0] = -ing;
  ydot[1] = rho*eps*V*ing - pc;
  ydot[2] = K*pc - km*W;
  ydot[3] = (1-K)*pc/(E0+W0);

}
