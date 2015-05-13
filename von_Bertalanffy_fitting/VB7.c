/* pomp model file: VB7 */

#include <pomp.h>
#include <R_ext/Rdynload.h>

#define r	(__p[__parindex[0]])
#define Linf	(__p[__parindex[1]])
#define L_0	(__p[__parindex[2]])
#define L_sd	(__p[__parindex[3]])
#define G_sd	(__p[__parindex[4]])
#define L	(__x[__stateindex[0]])
#define Lobs	(__y[__obsindex[0]])
#define DL	(__f[__stateindex[0]])
#define Tr	(__pt[__parindex[0]])
#define TLinf	(__pt[__parindex[1]])
#define TL_0	(__pt[__parindex[2]])
#define TL_sd	(__pt[__parindex[3]])
#define TG_sd	(__pt[__parindex[4]])
#define lik	(__lik[0])

void VB7_par_trans (double *__pt, double *__p, int *__parindex)
{

  Tr = exp(r);
  TLinf = exp(Linf);
  TL_sd = exp(L_sd);
  TL_0 = exp(L_0);
  TG_sd = exp(G_sd);
 
}


void VB7_par_untrans (double *__pt, double *__p, int *__parindex)
{

  Tr = log(r);
  TLinf = log(Linf);
  TL_sd = log(L_sd);
  TL_0 = log(L_0);
  TG_sd = log(G_sd);
 
}


void VB7_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
  Lobs = rnorm(L,L_sd);  
}


void VB7_dmeasure (double *__lik, double *__y, double *__x, double *__p, int give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
  lik = dnorm(Lobs,L,L_sd,give_log);  
}


void VB7_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)
{
 double rate; //transition rates
  double trans; // transition numbers
  rate = r*(Linf-L);
  trans = rnorm(rate*dt, G_sd*sqrt(dt));
  L += trans;
 
}


void VB7_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
 double rate; //transition rates
  rate = r*(Linf-L);
  DL += rate;
 
}

#undef r
#undef Linf
#undef L_0
#undef L_sd
#undef G_sd
#undef L
#undef Lobs
#undef DL
#undef Tr
#undef TLinf
#undef TL_0
#undef TL_sd
#undef TG_sd
