#include <R.h>
#include <stdio.h>

static double parms[12];

#define rho parms[0]
#define eps parms[1]
#define Imax parms[2]
#define g parms[3]
#define F parms[4]
#define xi parms[5]
#define q parms[6]
#define K parms[7]
#define km parms[8]
#define ER parms[9]
#define v parms[10]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=12;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {
  // state variables
  double E = y[0];
  double W = y[1];
  double L = pow(W, 2/3);
  double Wobs = W+E;
  double Lobs = pow(Wobs/xi, 1/q);

  // mobilization
  double pc = E * (v/L + km) / (1 + K*E/W);

  // balance the equations
  ydot[0] = rho*eps*Imax*pow(Lobs,g)*F - pc;
  ydot[1] = K*pc - km*W;
  ydot[2] = (1-K)*pc/ER;

}

/* when t > t_rep, set C=0 */
/*void event(int *n, double *t, double *y) {
  y[1] = 0;
  } */

