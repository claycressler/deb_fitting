#include <R.h>
#include <stdio.h>

static double parms[16];

#define K parms[0]
#define km parms[1]
#define eG parms[2]
#define eR parms[3]
#define v parms[4]
#define Rmbar parms[5]
#define pam parms[6]
#define Fh parms[7]
#define eA parms[8]
#define vol parms[9]
#define E0 parms[10]
#define L0 parms[11]
#define Re0 parms[12]
#define R0 parms[13]
#define F0 parms[14]
#define Lsd parms[15]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=16;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {
  double rate[4];

  // state variables: surface area, volume, and reserve density
  double SA;
  double V;
  double W;
  SA = y[1]*y[1];
  V = y[1]*y[1]*y[1];
  W = y[0]/V;

  // fluxes
  double PA;
  double PC;
  PA = y[4]/(Fh+y[4])*SA*pam;
  PC = V*(v/y[1]+km)*W/(1+K/eG*W);

  // compute the transition rates
  rate[0] = PA-PC;
  rate[1] = (K/eG*PC-km*V)/(3*SA);
  rate[2] = (1-K)*PC;
  rate[3] = -PA/eA/vol;

  // balance the equations
  ydot[0] = rate[0];
  ydot[1] = rate[1];
  if (y[2] < Rmbar) {
    ydot[2] = rate[2];
    ydot[3] = 0;
  }
  else {
    ydot[2] = 0;
    ydot[3] = rate[2]/eR;
  }
  ydot[4] = rate[3];

}

/* every two days, food is reset to F0 */
void event(int *n, double *t, double *y) {
  y[4] = F0;
}
