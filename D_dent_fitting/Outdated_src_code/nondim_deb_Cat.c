#include <R.h>
#include <stdio.h>

static double parms[7];
#define g parms[0]
#define q parms[1]
#define alpha parms[2]
#define rho parms[3]
#define phi parms[4]
#define beta parms[5]
#define w_scalar parms[6]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=7;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {
  // state variables
  double f = y[0];
  double w = y[1];
  double l = y[2];

  // balance the equations
  ydot[0] = -f/(1+f)*pow(l,g);
  ydot[1] = rho/w_scalar*phi*f/(1+f)*pow(l,g-q) - alpha*w/l;
  ydot[2] = alpha*l/q*((w/l-1)/(1+w));
  ydot[3] = beta*alpha*pow(l,q)*((w/l+w)/(1+w));

}

