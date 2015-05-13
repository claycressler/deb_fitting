#include <R.h>
#include <stdio.h>

static double parms[15];

#define K parms[0]
#define km parms[1]
#define eG parms[2]
#define eR parms[3]
#define v parms[4]
#define Rmbar parms[5]
#define Imax parms[6]
#define Fh parms[7]
#define eA parms[8]
#define E0 parms[9]
#define L0 parms[10]
#define R0 parms[11]
#define F0 parms[12]
#define Lsd parms[13]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=14;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {

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
  PA = eA*F0;
  PC = V*(v/y[1]+km)*W/(1+K/eG*W);

  // compute the transition rates
  ydot[0] = PA-PC;
  ydot[1] = (K/eG*PC-km*V)/(3*SA);
  ydot[2] = (1-K)*PC;

}

