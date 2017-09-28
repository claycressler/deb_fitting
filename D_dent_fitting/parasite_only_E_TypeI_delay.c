#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>

static double parms[6];

#define E parms[0]
#define aP parms[1]
#define eP parms[2]
#define tau parms[3]
#define P0 parms[4]
#define Perr parms[5]

/* Interface to dede utility functions in package deSolve */
void lagvalue(double T, int *nr, int N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
  return fun(T, nr, N, ytau);
}

void lagderiv(double T, int *nr, int N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagderiv");
  return fun(T, nr, N, ytau);
}

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=6;
  odeparms(&N, parms);
}

/* derivatives */
/* this is for the model where parasites "eat" reserves */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
  double Pi = y[0]; // immature parasites

  // total number of state variables
  int Nout = 2;
  // nr is an integer that defines the number of the state variable or its derivative whose delay we want
  int nr[2] = {0, 1};
  // array; initialize with default values
  double ytau[2] = {0.0, 0.0};

  // what time are we "lagging from"?
  double T = *t - tau;

  if (*t > tau) {
    lagvalue(T, nr, Nout, ytau);
  }

  // balance the equations
  ydot[0] = eP*aP*E*(Pi - ytau[0]);
  ydot[1] = eP*aP*E*ytau[0];
}

