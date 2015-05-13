## This R-file implements the versions of the DEB model decided on by
## Cressler and Nelson following a meeting 19 Nov 2012. The model
## differs from the standard DEB model in assuming that there is no
## "maturity maintenance" (k_j). We implement several versions, with
## different assumptions regarding assimilation and starvation.

## DEALING WITH INGESTION: Since the experimental data fed individuals
## a fixed amount of food, we know exactly how much food was ingested
## - all that matters is the rate at which the food is depleted. The
## simplest (but incorrect) assumption is that the assimilation rate
## is a constant F/tau, where F is the food ration and tau is the
## transfer interval. This assumption is incorrect because it assumes
## that the food lasts the entire transfer interval and that
## individuals that differ in size will not differ in assimilation
## rate. A second possible assumption is that food is depleted
## according to a surface-area specific ingestion rate. This is the
## standard DEB assumption. Since ingestion is surface-area specific,
## large individuals will deplete the food faster. If they run out of
## food quickly, then they will enter into starvation conditions and
## the internal dynamics of energy will change qualitatively
## (discussed below). A third possible assumption is that assimilation
## is constant F/(theta*tau), where theta is the fraction of time food
## is present - for (1-theta)*tau, assimilation is 0. We could assume
## further either that theta is a function of size, or that it is
## independent of size.

## DEALING WITH STARVATION: Starvation, in the DEB world, is defined
## by the allocation to growth being insufficient to cover somatic
## maintenance (kappa*p_c < k_m*V). We consider three possible rules
## for starvation.
## 1. Keep kappa constant and burn structure ("shrink") to pay somatic
## maintenance. Reproduction could potentially continue even as the individual was shrinking.
## 2. When kappa*pc < k_m*V, set kappa equal to the value necessary for kappa*p_c to just meet maintenance. Reproduction could potentially continue even after the individual stopped growing.
## 3. Make kappa a function of reserve density (e=E/V) - when e is
## high, hold kappa constant at some minimum value, but increase kappa
## as e gets smaller. The shape and inflection point of this curve
## will be something that we use the data to estimate - for now it
## will be specified as a smooth curve involving the hyperbolic
## tangent function
## kappa(e)=0.5*((Kmin+Kmax)-(Kmax-Kmin)*tanh(w*(e-e_in)), where Kmin
## is the min. value of kappa, Kmax is the max. value of kappa
## (presumably 1), w is a parameter that controls how quickly kappa
## moves from Kmin to Kmax (the larger the faster), and e_in is the
## inflection point. This implies three parameters that need to be
## fitted, which may prove too difficult. For this initial data, I will assume that the inflection point of the curve is at half of maximum reserve density - in DEB theory, the maximum E/V is given by max. surface-area specific assimilation rate divided by energy conductance: p_am/v. I will adopt a similar formalism for now.

## For a first pass of simulation-recovery experiments, I will
## consider only two possibilities for both ingestion and starvation
## (so four models altogether) - ingest1. type II functional response
## with constant food and ingest2. type II functional response with
## food depletion, and starve1. constant kappa and starve2. kappa as a
## function of e. I am not considering ingestion that is
## size-independent because I am not sure how to calculate the maximum
## reserve density when the assimilation rate is not SA-dependent.

deb.ingest1.starve1 <- function(t, y, params) {
  ## State variables
  E <- y[1] ## reserves (mgC)
  L <- y[2] ## structural length (mm)
  Rm <- y[3] ## maturity (mgC)
  Re <- y[4] ## cumulative eggs (#)

  ## Parameters
  K <- params$Kmin ## allocation to growth
  km <- params$km ## somatic maintenance rate (1/day)
  eG <- params$eG ## cost of growth (mgC/mm^3)
  eR <- params$eR ## conversion/cost of egg (mgC/#)
  v <- params$v ## energy conductance (mm/day)
  Rmbar <- params$Rmbar ## energy required for sexual maturity (mgC)
  pam <- params$pam ## maximum surface-area specific assimilation rate (mgC/mm^2/day)
  Fh <- params$Fh ## half-saturation constant of functional response (mgC/ml)
  F0 <- params$F0 ## total food quantity fed at each transfer, although we are pretending that it is a chemostat with that density of food in it (mgC/ml)

  ## Assimilation rate (mgC/day)
  pa <- pam*F0/(Fh+F0)*y[2]^2

  ## Mobilization rate (mgC/day)
  pc <- y[2]^3*(v/y[2]+km)*(y[1]/y[2]^3)/(1+K/eG*y[1]/y[2]^3)

  ## Balance equations
  y1 <- pa-pc
  y2 <- (K*pc-km*eG*y[2]^3)/(3*eG*y[2]^2)
  if (y[3] < Rmbar) {
    y3 <- (1-K)*pc
    y4 <- 0
  }
  else {
    y3 <- 0
    y4 <- (1-K)*pc/eR
  }

  return(c(y1,y2,y3,y4))
}

deb.ingest1.starve2 <- function(t, y, params) {
  ## State variables
  E <- y[1] ## reserves
  L <- y[2] ## structural length
  Rm <- y[3] ## maturity
  Re <- y[4] ## cumulative eggs

  ## Parameters
  Kmin <- params$Kmin ## allocation to growth
  km <- params$km ## somatic maintenance rate
  eG <- params$eG ## cost of growth
  eR <- params$eR ## conversion/cost of egg
  v <- params$v ## energy conductance
  Rmbar <- params$Rmbar ## energy required for sexual maturity
  pam <- params$pam ## maximum surface-area specific assimilation rate
  Fh <- params$Fh ## half-saturation constant of functional response
  F0 <- params$F0 ## total food quantity fed at each transfer, although we are pretending that it is a chemostat with that density of food in it
  w <- params$w ## slope of K(e) plot

  ## Assimilation rate
  pa <- pam*F0/(Fh+F0)*L^2

  ## Calculate reserve density and kappa
  e <- E/L^3
  emax <- pam/v ## maximum reserve density
  ein <- 0.75*emax ## inflection point is at half e_max
  K <- 0.5*((1+Kmin)-(1-Kmin)*tanh(w*(e-ein)))

  ## Mobilization rate
  pc <- L^3*(v/L+km)*eG*E/(eG*L^3+K*E)

  ## Balance equations
  dEdt <- pa-pc
  dLdt <- (K*pc-km*eG*L^3)/(3*eG*L^2)
  if (Rm < Rmbar) {
    dRmdt <- (1-K)*pc
    dRedt <- 0
  }
  else {
    dRmdt <- 0
    dRedt <- (1-K)*pc/eR
  }

  return(c(dEdt,dLdt,dRmdt,dRedt))
}

deb.ingest2.starve1 <- function(t, y, params) {
  ## State variables
  E <- y[1] ## reserves
  L <- y[2] ## structural length
  Rm <- y[3] ## maturity
  Re <- y[4] ## cumulative eggs
  F <- y[5] ## food

  ## Parameters
  K <- params$Kmin ## allocation to growth
  km <- params$km ## somatic maintenance rate
  eG <- params$eG ## cost of growth
  eR <- params$eR ## conversion/cost of egg
  v <- params$v ## energy conductance
  Rmbar <- params$Rmbar ## energy required for sexual maturity
  pam <- params$pam ## maximum surface-area specific assimilation rate
  Fh <- params$Fh ## half-saturation constant of functional response
  ea <- params$ea ## assimilation efficiency
  F0 <- params$F0 ## total food quantity fed at each transfer
  Vol <- params$Vol ## container volume

  ## Assimilation rate
  pa <- pam*F/(Fh+F)*L^2

  ## Mobilization rate
  pc <- L^3*(v/L+km)*eG*E/(eG*L^3+K*E)

  ## Balance equations
  dEdt <- pa-pc
  dLdt <- (K*pc-km*eG*L^3)/(3*eG*L^2)
  if (Rm < Rmbar) {
    dRmdt <- (1-K)*pc
    dRedt <- 0
  }
  else {
    dRmdt <- 0
    dRedt <- (1-K)*pc/eR
  }
  dFdt <- -pa/ea/Vol

  return(c(dEdt,dLdt,dRmdt,dRedt,dFdt))
}

deb.ingest2.starve2 <- function(t, y, params) {
  ## State variables
  E <- y[1] ## reserves
  L <- y[2] ## structural length
  Rm <- y[3] ## maturity
  Re <- y[4] ## cumulative eggs
  F <- y[5] ## food

  ## Parameters
  Kmin <- params$Kmin ## allocation to growth
  km <- params$km ## somatic maintenance rate
  eG <- params$eG ## cost of growth
  eR <- params$eR ## conversion/cost of egg
  v <- params$v ## energy conductance
  Rmbar <- params$Rmbar ## energy required for sexual maturity
  pam <- params$pam ## maximum surface-area specific assimilation rate
  Fh <- params$Fh ## half-saturation constant of functional response
  ea <- params$ea ## assimilation efficiency
  w <- params$w ## slope of K(e) plot
  F0 <- params$F0 ## total food quantity fed at each transfer
  Vol <- params$Vol ## container volume

  ## Assimilation rate
  pa <- pam*F/(Fh+F)*L^2

  ## Calculate reserve density and kappa
  e <- E/L^3
  emax <- pam/v ## maximum reserve density
  ein <- 0.75*emax ## inflection point is at half e_max
  K <- 0.5*((1+Kmin)-(1-Kmin)*tanh(w*(e-ein)))

  ## Mobilization rate
  pc <- L^3*(v/L+km)*eG*E/(eG*L^3+K*E)

  ## Balance equations
  dEdt <- pa-pc
  dLdt <- (K*pc-km*eG*L^3)/(3*eG*L^2)
  if (Rm < Rmbar) {
    dRmdt <- (1-K)*pc
    dRedt <- 0
  }
  else {
    dRmdt <- 0
    dRedt <- (1-K)*pc/eR
  }
  dFdt <- -pa/ea/Vol

  return(c(dEdt,dLdt,dRmdt,dRedt,dFdt))
}
