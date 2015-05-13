## The likelihood of observing length Lobs and births Robs, given that
## L and R were the actual length and births. The length measurement
## error is normally distributed and the number of neonates counted is
## Poisson distributed.
" lik = dnorm(Lobs,L,L_sd,give_log) + dpois(Robs,R,give_log); " -> dmeas

## Simulate a measurement
" Lobs = rnorm(L,L_sd);
  Robs = rpois(R); " -> rmeas

## Encode the deterministic skeletons.
## Model 1: chemostat ingestion with constant kappa
## For this model, I do not need to estimate the maximum surface-area
## specific assimilation rate or the half-saturation constant, because
## the functional response pam*F/(Fh+F) will be constant. So I can
## just esimate that SA-specific functional response as a parameter f.
"
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
" -> skel

## Encode the stochastic model process simulator
"
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
" -> stepfn

## Parameter transformations. On the estimation scale, we the
## parameters to be unconstrained, even though on the natural scale
## they are constrained to be positive (all parameters) and less than
## one (K). Note that the initial conditions for the state variables
## and the standard deviation of the observation error for length
## measurements are also parameters. The transform from the estimation
## scale to the natural scale is given as:
"
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
" -> trans
## The transform from the natural scale to the estimation scale is given as:
"
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
" -> untrans
