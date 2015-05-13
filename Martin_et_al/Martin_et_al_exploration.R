deb.model.chemostat <- function(t,y,params) {
    ## State variables
    UE <- y[1]
    L <- y[2]
    UH <- y[3]
    UR <- y[4]

    ## Parameters
    X <- as.numeric(params['X'])
    kappa <- 0.678
    KR <- 0.95
    km <- 0.3314
    kj <- 0.1921
    UHb <- 0.111
    UHp <- 2.547
    v <- 18.1
    g <- 10
    K <- 1585

    e <- v*UE/L^3
    SA <- X/(K+X)*L^2
    SC <- L^2*g*e/(g+e)*(1+L*km/v)

    dUEdt <- SA-SC
    dLdt <- 1/3*(v/(g*L^2)*SC-km*L)
    if (UH < UHp) {
        dUHdt <- (1-kappa)*SC-kj*UH
        dURdt <- 0
    }
    else {
        dUHdt <- 0
        dURdt <- (1-kappa)*SC-kj*UHp
    }

    list(c(dUEdt,dLdt,dUHdt,dURdt))
}

## This R-function computes the values of the derivatives of the ODE
## system at time t.
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

