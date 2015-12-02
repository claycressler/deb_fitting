#include <R.h>
#include <stdio.h>

static double parms[6];
#define g parms[0]
#define q parms[1]
#define sigma parms[2]
#define phi parms[3]
#define alpha parms[4]
#define beta parms[5]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=6;
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
  ydot[1] = sigma*phi*f/(1+g)*pow(l,g-q);
  ydot[2] = alpha*l/q*((w/l-q)/(1+w));
  ydot[3] = beta*alpha*pow(l,q)*((w/l+w)/(1+w));

}

