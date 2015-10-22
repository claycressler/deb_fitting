#include <R.h>
#include <stdio.h>

static double parms[16];

#define rho parms[0]
#define Imax parms[1]
#define g parms[2]
#define K parms[3]
#define km parms[4]
#define EG parms[5]
#define ER parms[6]
#define v parms[7]
#define aC parms[8]
#define hC parms[9]
#define aG parms[10]
#define hG parms[11]
#define b parms[12]
#define m parms[13]
#define t_rep parms[14]
#define obs_sd parms[15]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=16;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {
  // state variables
  double E = y[0];
  double L = pow(y[1],2/3);
  double V = y[1];

  // mobilization
  double pc = E * (v/L + km) / (1 + K*E/(EG*V));

  // balance the equations
  ydot[0] = rho*Imax*pow(L,g) - pc;
  ydot[1] = (K*pc - km*EG*V)/EG;
  ydot[2] = (1-K)*pc/ER;

}

/* when t > t_rep, set C=0 */
void event(int *n, double *t, double *y) {
  y[1] = 0;
}

