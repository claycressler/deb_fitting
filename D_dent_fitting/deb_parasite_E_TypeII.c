#include <R.h>
#include <stdio.h>

static double parms[11];

#define rho parms[0]
#define K parms[1]
#define km parms[2]
#define v parms[3]
#define F0 parms[4]
#define Lerr parms[5]
#define aP parms[6]
#define eP parms[7]
#define hP parms[8]
#define P0 parms[9]
#define Perr parms[10]

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=11;
  odeparms(&N, parms);
}

/* derivatives */
/* this is for the model where parasites "eat" reserves */
void derivs (int *neq, double *t, double *y, double *ydot) {
  // fixed parameters whose values will not ever vary
  double eps = 0.00000000816536; // carbon content of algae from Meg's data
  double V = 30; // volume of the container
  double xi = 0.0018; // length-weight regression coefficient from Spencer
  double q = 3; // length-weight regression exponent from Spencer

  // foraging-dependent parameters
  double Fh = 12000; // half-saturation constant
  double Imax = 12162; // ingestion rate
  double g = 1.467; // size-dependence of ingestion
  double a = 3.2e-5; // spore coefficient
  double h = 3.57; // size -dependence of spore dependence

  // state variables
  double F = y[0]; // algal concentration
  double E = y[1]; // energy reserves
  double W = y[2]; // structural weight
  double P = y[3]; // parasites
  double L = pow(W, 1/3); // structural length
  double Lobs = pow(W/xi, 1/q); //observed length


  // ingestion
  double ing = Imax * F/(Fh+F) * pow(Lobs,g) * exp(-a * P/pow(Lobs,h));
  // mobilization
  double pc = E * (v/L + km) / (1 + K*E/W);

  // balance the equations
  ydot[0] = -ing;
  ydot[1] = rho*eps*V*ing - pc - aP*E*P/(hP+E);
  ydot[2] = K*pc - km*W;
  ydot[3] = eP * aP*E/(hP+E)*P;
}

