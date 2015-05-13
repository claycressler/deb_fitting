## This R-function computes the values of the derivatives of the ODE
## system at time t. For simplicity, I focus on reserves, length, and
## food density as the only state variables.
deb.model.batch <- function(t,y,params) {
    ## State variables
    UE <- y[1] ## total reserves
    L <- y[2] ## structural length
    X <- y[3]  ## food density

    ## Parameters
    km <- 0.3314 ## somatic maintenance rate
    v <- 18.1 ## energy conductance
    g <- 10 ## energy investment ratio
    K <- 1585 ## half-saturation constant
    JAm <- 3.8e5 ## max. SA-specific ingestion rate

    ## Calculate the scaled reserve density
    e <- v*UE/L^3

    ## Calculate the assimilation and mobilization fluxes
    SA <- X/(K+X)*L^2
    SC <- L^2*g*e/(g+e)*(1+L*km/v)

    ## Compute the derivatives of the ODE system
    dUEdt <- SA-SC
    dLdt <- 1/3*(v/(g*L^2)*SC-km*L)
    dXdt <- -JAm*SA

    list(c(dUEdt,dLdt,dXdt))
}

## Structural length at birth
L0 <- 0.851

## Reserves at birth, from e=v*UE/L^3=1
UE0 <- 0.851^3/18.1

## Initial food density, assuming 0.05 mgC/daphnid/day, 1.95e-8
## mgC/cell, and a container volume of 80 ml (Coors et al. 2004 low
## food conditions). This is also the amount of food added each day.
X0 <- 0.05/1.95e-8/80

## Vector of initial conditions
y0 <- c(UE0,L0,X0)

## The numerical solver assumes that the state variables are
## continuous and smooth with smooth derivatives. This assumption is
## violated by the batch culture conditions, where a bolus of food is
## added once per day and individuals are transferred to new vials and
## fed a bolus of food every other day. A simple work-around for this
## problem is to solve the ODE for one-day intervals, with
## discontinuous changes in the food state between.
require(deSolve)

## Number of two-day transfers
no.transfers <- 10
## A vector to record the numerical solution results
out <- vector()

for (i in 1:no.transfers) {
    ## Simulate the first day (up to the first feeding)
    ## Times that explicit estimates of state variables are recorded
    times <- seq(0,1,0.1)

    ## Call to the numerical solver
    out1 <- lsoda(y=y0, times=times, func=deb.model.batch)

    ## Record the numerical solution in out - do not record the data
    ## in the final row, as this data serves as the initial conditions
    ## for the next simulation.
    out <- rbind(out, out1[1:(nrow(out1)-1),])

    ## Set the new initial conditions using the final state (the state
    ## variables are recorded in columns 2-4)
    y0 <- out1[nrow(out1),2:4]
    ## Add the bolus of food
    y0[3] <- y0[3]+X0

    ## Simulate the second day (up to the transfer)
    out2 <- lsoda(y=y0, times=times, func=deb.model.batch, parms=NULL)

    ## Record
    out <- rbind(out, out2[1:(nrow(out2)-1),])

    ## Set the new initial conditions
    y0 <- out2[nrow(out2),2:4]
    ## Transfers reset the amount of food to X0
    y0[3] <- X0

    ## Continue until the required number of transfers have been
    ## simulated.
}

## Plot length through time and compare to figure 6 in the appendix.
plot(out[,3], type='l')

## Plot food density through time
plot(out[,4], type='l')

## Food depletion happens so quickly that individuals rapidly enter
## into starvation conditions, depleting their reserves to essentially
## zero and burning structure to meet somatic maintenance
## requirements. Because of this, individuals never experience any
## sustained growth. This differs dramatically from the experimental
## results of Coors et al. 2004 and from the results presented in
## Figure 6 of the appendix. Am I doing something wrong?

## Even at the highest food level, the individual still goes through
## daily periods where structure is burned to pay maintenance.
## Initial food density, assuming 0.2 mgC/daphnid/day.
X0 <- 5/1.95e-8/80
y0 <- c(UE0,L0,X0)
no.transfers <- 10
out.hi <- vector()
for (i in 1:no.transfers) {
    times <- seq(0,1,0.1)
    out1 <- lsoda(y=y0, times=times, func=deb.model.batch)
    out.hi <- rbind(out.hi, out1[1:(nrow(out1)-1),])
    y0 <- out1[nrow(out1),2:4]
    y0[3] <- y0[3]+X0
    out2 <- lsoda(y=y0, times=times, func=deb.model.batch, parms=NULL)
    out.hi <- rbind(out.hi, out2[1:(nrow(out2)-1),])
    y0 <- out2[nrow(out2),2:4]
    y0[3] <- X0
}
plot(out.hi[,3], type='l') ## Length through time
plot(out.hi[,4], type='l') ## Food through time


## However, if I modify the R-function to approximate the continuous
## culture conditions of Sokull-Kluettgen 1998, I can reproduce the
## results shown in Figure 1 of the manuscript.

deb.model.chemostat <- function(t,y,params) {
    ## State variables
    UE <- y[1] ## total reserves
    L <- y[2] ## structural length

    ## Parameters
    X <- params['X'] ## food density
    km <- 0.3314 ## somatic maintenance rate
    v <- 18.1 ## energy conductance
    g <- 10 ## energy investment ratio
    K <- 1585 ## half-saturation constant
    JAm <- 3.8e5 ## max. SA-specific ingestion rate

    ## Calculate the scaled reserve density
    e <- v*UE/L^3

    ## Calculate the assimilation and mobilization fluxes
    SA <- X/(K+X)*L^2
    SC <- L^2*g*e/(g+e)*(1+L*km/v)

    ## Compute the derivatives of the ODE system
    dUEdt <- SA-SC
    dLdt <- 1/3*(v/(g*L^2)*SC-km*L)

    list(c(dUEdt,dLdt))
}

## Initial conditions for reserve and structural length
y0 <- c(0.851^3/18.1, 0.851)

## Low food density of Sokull-Kluettgen 1998
X0 <- 1000
pars <- c(X=X0) ## parameter vector

## Simulate for 42 days
times <- seq(0,42,0.1)

## Call to the solver
out.chemostat <- lsoda(y=y0, times=times, func=deb.model.chemostat, parms=pars)

## Plot length as a function of age
plot(out.chemostat[,c(1,3)], type='l')

## Repeat the analysis at a higher food density
X0 <- 25000
pars <- c(X=X0) ## parameter vector
times <- seq(0,42,0.1)
out.chemostat <- lsoda(y=y0, times=times, func=deb.model.chemostat, parms=pars)
plot(out.chemostat[,c(1,3)], type='l')

## Both of the plots of the continuous culture experimental conditions
## look very similar to Figure 1 of the manuscript. To me, this
## suggests that the discrepancy between my results and your results
## is arising from the feeding submodel, and possibly with the value
## of the parameter {J_Am}.


Dear Dr. Martin -

    I wanted to ask a question about your forthcoming American Naturalist paper, for which I was privileged to serve as a reviewer. As someone attempting to apply dynamic energy budget theory, I was very happy to see your study, and I think it is an important paper that will serve as a catalyst for further a lot of future theory development.

I am also working on the problem of fitting DEB models to growth and reproduction trajectories. I was encouraged to see that the standard DEB model did a good job of fitting the growth and reproduction trajectories from Coors et al. 2004. The data I have were also collected under batch culture conditions. In particular, it was very good to see (in your Appendix Figure 6) that starvation (in particular, burning structure to pay maintenance) did not seem to be a problem.

However, I have run into problems trying to crudely reproduce the results you present in Figure 6 in the Appendix. When I simulate the DEB equations under batch culture conditions at low food, I find that the individual rapidly depletes the food, causing it to have to burn structure to pay maintenance every day, and therefore never has any sustained growth in length. It appears to me that the value of {J_Am} is too high, but I assume I have made a mistake somewhere in my code. Would you be willing to help me understand what I am doing wrong? I have attached the R-code I am using to simulate the DEB model under batch transfer conditions and have tried to annotate it fully.

Thank you in advance, and congratulations on getting this important study published.

Best,
Clay
