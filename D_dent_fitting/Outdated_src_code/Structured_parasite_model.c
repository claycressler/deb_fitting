#include <R.h>
#include <stdio.h>

static double parms[10];

#define theta parms[0]
#define r parms[1]
#define aC parms[2]
#define hC parms[3]
#define aG parms[4]
#define hG parms[5]
#define b parms[6]
#define m parms[7]
#define t_rep parms[8]
#define obs_sd parms[9]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=10;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {

  // balance the equations
  ydot[0] = theta - r*y[0] - aC*y[0]*y[1]/(hC+y[0]) - aG*y[0]*y[2]/(hG+y[0]);
  ydot[1] = 0;
  ydot[2] = b*aC*y[0]*y[1]/(hC+y[0]) - m*aG*y[0]*y[2]/(hG+y[0]);
  ydot[3] = m*aG*y[0]*y[2]/(hG+y[0]);

}

/* when t > t_rep, set C=0 */
void event(int *n, double *t, double *y) {
  y[1] = 0;
}
