#include <R.h>
#include <stdio.h>

static double parms[14];

#define alpha   parms[0]
#define kappa   parms[1]
#define beta    parms[2]
#define W_scalar       parms[3]
#define L_scalar       parms[4]
#define Rm_scalar      parms[5]
#define time_scalar    parms[6]
#define f       parms[7]
#define T_M     parms[8]
#define W_0     parms[9]
#define Rm_0    parms[10]
#define L_0     parms[11]
#define L_sd    parms[12]
#define R_sd    parms[13]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=14;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {

  // state variables: energy density, length
  double W = y[0];
  double L = y[3];
  double Lmeas = y[1];

  // balance the equations
  ydot[0] = (f/T_M)*(alpha/(L*L*L))-(W/L);
  ydot[3] = (kappa*W-L)/(3*(1+kappa*W));
  ydot[1] = 0; if(ydot[3]>0 && L>=Lmeas){ydot[1] = ydot[3];}
  ydot[2] = ((1-kappa)*L*L*(1+L)*W)/(1+kappa*W);

}


