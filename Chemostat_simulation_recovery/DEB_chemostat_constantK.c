/* pomp model file: DEB_chemostat_constantK */

#include <pomp.h>
#include <R_ext/Rdynload.h>

#define K	(__p[__parindex[0]])
#define km	(__p[__parindex[1]])
#define eG	(__p[__parindex[2]])
#define eR	(__p[__parindex[3]])
#define v	(__p[__parindex[4]])
#define Rmbar	(__p[__parindex[5]])
#define f	(__p[__parindex[6]])
#define E_0	(__p[__parindex[7]])
#define L_0	(__p[__parindex[8]])
#define Re_0	(__p[__parindex[9]])
#define R_0	(__p[__parindex[10]])
#define L_sd	(__p[__parindex[11]])
#define PA_sd	(__p[__parindex[12]])
#define PC_sd	(__p[__parindex[13]])
#define E	(__x[__stateindex[0]])
#define L	(__x[__stateindex[1]])
#define Re	(__x[__stateindex[2]])
#define R	(__x[__stateindex[3]])
#define Lobs	(__y[__obsindex[0]])
#define Robs	(__y[__obsindex[1]])
#define DE	(__f[__stateindex[0]])
#define DL	(__f[__stateindex[1]])
#define DRe	(__f[__stateindex[2]])
#define DR	(__f[__stateindex[3]])
#define TK	(__pt[__parindex[0]])
#define Tkm	(__pt[__parindex[1]])
#define TeG	(__pt[__parindex[2]])
#define TeR	(__pt[__parindex[3]])
#define Tv	(__pt[__parindex[4]])
#define TRmbar	(__pt[__parindex[5]])
#define Tf	(__pt[__parindex[6]])
#define TE_0	(__pt[__parindex[7]])
#define TL_0	(__pt[__parindex[8]])
#define TRe_0	(__pt[__parindex[9]])
#define TR_0	(__pt[__parindex[10]])
#define TL_sd	(__pt[__parindex[11]])
#define TPA_sd	(__pt[__parindex[12]])
#define TPC_sd	(__pt[__parindex[13]])
#define lik	(__lik[0])

void DEB_chemostat_constantK_par_trans (double *__pt, double *__p, int *__parindex)
{

  TK = exp(K)/(1+exp(K));
  Tkm = exp(km);
  TeG = exp(eG);
  TeR = exp(eR);
  Tv = exp(v);
  TRmbar = exp(Rmbar);
  Tf = exp(f);
  TE_0 = exp(E_0);
  TL_0 = exp(L_0);
  TRe_0 = exp(Re_0);
  TR_0 = exp(R_0);
  TL_sd = exp(L_sd);
  TPA_sd = exp(PA_sd);
  TPC_sd = exp(PC_sd);
 
}


void DEB_chemostat_constantK_par_untrans (double *__pt, double *__p, int *__parindex)
{

  TK = log(K/(1-K));
  Tkm = log(km);
  TeG = log(eG);
  TeR = log(eR);
  Tv = log(v);
  TRmbar = log(Rmbar);
  Tf = log(f);
  TE_0 = log(E_0);
  TL_0 = log(L_0);
  TRe_0 = log(Re_0);
  TR_0 = log(R_0);
  TL_sd = log(L_sd);
  TPA_sd = log(PA_sd);
  TPC_sd = log(PC_sd);
 
}


void DEB_chemostat_constantK_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
  Lobs = rnorm(L,L_sd);
  Robs = rpois(R);  
}


void DEB_chemostat_constantK_dmeasure (double *__lik, double *__y, double *__x, double *__p, int give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
  lik = dnorm(Lobs,L,L_sd,give_log) + dpois(Robs,R,give_log);  
}


void DEB_chemostat_constantK_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)
{

  double trans[2]; // transition numbers

  // state variables: surface area, volume, and reserve density
  double SA;
  double V;
  double W;
  SA = L*L;
  V = L*L*L;
  W = E/V;

  // fluxes
  double PA;
  double PC;
  PA = f*SA;
  PC = V*(v/L+km)*W/(1+K/eG*W);

  // stochasticity is in assimilation and mobilization
  trans[0] = rnorm(PA, PA_sd);
  trans[1] = rnorm(PC, PC_sd);

  // balance the equations
  E += trans[0]-trans[1];
  L += (K/eG*trans[1]-km*V)/(3*SA);
  if (Re < Rmbar) {
    Re += (1-K)*trans[1];
    R += 0;
  }
  else {
    Re += 0;
    R += (1-K)*trans[1]/eR;
  }
 
}


void DEB_chemostat_constantK_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{

  double rate[3]; // transition rates

  // state variables: surface area, volume, and reserve density
  double SA;
  double V;
  double W;
  SA = L*L;
  V = L*L*L;
  W = E/V;

  // fluxes
  double PA;
  double PC;
  PA = f*SA;
  PC = V*(v/L+km)*W/(1+K/eG*W);

  // compute the transition rates
  rate[0] = PA-PC;
  rate[1] = (K/eG*PC-km*V)/(3*SA);
  rate[2] = (1-K)*PC;

  // balance the equations
  DE = rate[0];
  DL = rate[1];
  if (Re < Rmbar) {
    DRe = rate[2];
    DR = 0;
  }
  else {
    DRe = 0;
    DR = rate[2]/eR;
  }
 
}

#undef K
#undef km
#undef eG
#undef eR
#undef v
#undef Rmbar
#undef f
#undef E_0
#undef L_0
#undef Re_0
#undef R_0
#undef L_sd
#undef PA_sd
#undef PC_sd
#undef E
#undef L
#undef Re
#undef R
#undef Lobs
#undef Robs
#undef DE
#undef DL
#undef DRe
#undef DR
#undef TK
#undef Tkm
#undef TeG
#undef TeR
#undef Tv
#undef TRmbar
#undef Tf
#undef TE_0
#undef TL_0
#undef TRe_0
#undef TR_0
#undef TL_sd
#undef TPA_sd
#undef TPC_sd
