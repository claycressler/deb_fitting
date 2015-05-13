## First, draw DEB parameters from a multivariate normal distribution
## with known means and variance-covariance matrix. the MASS package
## has a function mvrnorm for doing this. mvrnorm(n, mu, Sigma), where
## n is the number of samples required
## mu is a vector giving the means of the variables
## Sigma is a positive-definite symmetric matrix specifying the
## covariance matrix of the variables

## Note: the estimated parameters have to be unconstrained, which means that the entries into the covariance matrix have to be unconstrained - does this compromise my ability to construct a positive-definite covariance matrix? Is it guaranteed that a covariance matrix is positive-definite on the natural scale?

## Easy to deal with - basically, first thing inside the Metropolis-Hastings algorithm to do is take the unconstrained variables and put them on the natural scale, draw the deb parameters, simulate, compute the likelihood, then take the next Metropolis-Hastings step. As far as Metropolis-Hastings algorithm is concerned, the parameters are unconstrained, because the constraint is imposed inside the likelihood evaluation. I want to write C-code for all aspects of this - parameter transformation and untransformation, numerical simulation, and likelihood calculation.

## Estimated parameters of the DEB model:
## alpha = nondimensional energy assimilation
## kappa = fractional energy allocation to growth
## beta  = nondimensional energy requirement for maturity
## a     = energy density scalar
## b     = length scalar
## c     = maturity energy scalar
## d     = time scalar
## W0    = energy density initial condition
## L0    = length initial condition

## Fixed parameters of the DEB model:
## R0    = maturity energy initial condition = 0
## F0    = food at each transfer
## T0    = transfer interval

## Draw the parameters from a multivariate normal distribution

## Dimensional parameters
k_M=0.3 #Structural volume per-capita maintenance rate (/day)
e_G=0.0016 #Cost of growth (?)
e_R=0.08 #Cost and conversion to produce an egg (?)
e_A=0.75 #Conversion efficiency of ingested to assimilated food (?)
v=1.9 #Energy conductance (mm/day)
R_M.bar=2 #Maturation threshold (units of eggs)
K=0.75 #Allocation of mobilization to structure (constant allocation variant)  (?)

## Initial conditions
E0=0.02 #Energy at birth (?)
L0obs=0.65 #Observed length at birth (mm)
R0=0 #Reproduction energy at birth
L0=0.65 #Length at birth (mm)

## Feeding parameters
f=0.05 #Known food provided at transfer (mgC)
T_M=2 #Known transfer interval (days)

## non-dimensional parameters
W0.nondim=(E0/(L0^3))/e_G
L0obs.nondim=L0*k_M/v
R0.nondim=0
L0.nondim=L0obs.nondim
alpha=(k_M^2)*e_A/((v^3)*e_G)
kappa=K
beta=R_M.bar*e_R*(k_M^3)/(e_G*(v^3))
a=e_G*k_M
b=v
c=e_G*(v^3)/(e_R*(k_M^2))
d=k_M

## C-code for running the model:
dyn.load('deb_redimensionalized.so')
y0 <- c(W0.nondim,L0obs.nondim,R0.nondim,L0.nondim)
pars <- c(alpha,kappa,beta,f,T_M)
out <- ode(y0,seq(0,75,0.1),func='derivs',dllname='deb_redimensionalized',parms=pars,initfunc='initmod')

## Get the parameters and estimated initial conditions from a
## multivariate normal distribution with means and variance-covariance
## matrix as specified.
## means
M <- c(alpha,kappa,beta,a,b,c,d,W0.nondim,L0.nondim)

## variance-covariance matrix - note that this must be
## positive-definite.

## To be a variance-covariance matrix, the diagonal
## elements (the variances) must be positive, and off-diagonal
## elements come in pairs: the covariance between variables i and j is
## the same as the covariance between variables j and i. This matters
## for the Metropolis-Hastings algorithm, because it is the entries in
## the variance-covariance matrix that are trying to be estimated.

## There is an R-function for generating a random covariance
## matrix. It works in the following way: it first generates a random
## correlation matrix R via the method proposed in Joe (2006), then
## randomly generates variances (sig_1^2, sig_2^2, ...) from an
## interval specified by the argument rangeVar. The covariance matrix
## is then constructed as
## diag(sig_1,...,sig_n)*R*diag(sig_1,...,sig_n).

require(clusterGeneration)
## generate a random correlation matrix based on random partial correlations
R <- rcorrmatrix(length(M))
## generate a vector of variances
sig <- runif(length(M), min=0, max=10)
## generate a positive-definite variance-covariance matrix
V <- diag(sig) %*% R %*% diag(sig)
## Q: why is the diagonal of V not equal to sig? The entries on the diagonal of a covariance matrix are supposed to be the variances.

## sample from the multivariate normal distribution implied
require(MASS)
p <- mvrnorm(n=1, mu=M, Sigma=V)

## So what is actually being estimated are the parameters of the
## correlation matrix and the variances, which are then used to
## construct a positive-definite variance-covariance matrix.

## However, a new problem arises: given a random correlation matrix
## and random variances, I cannot guarantee positive values for the
## parameters. We have, basically, two sets of constraints:
## 1. I want the estimated parameters (the means, variances, and correlation matrix) to be unconstrained on the estimation scale.
## 2. The DEB parameters must be positive after they are drawn from the multivariate normal distribution implied by the means, variances, and correlation matrix. This positivity requirement implies constraints on the variance-covariance matrix that are hard to define.

## I could have two layers of tranformation (that is, transform the means, variances, and correlation matrix once inside the likelihood evaluation, then transform the parameters after drawing them from the multivariate normal), but I fear this would render as meaningless the means, variances, and correlation matrix that were converged upon. On the other hand, is there any constraint I could set that would *guarantee* that the draws from the multivariate normal were positive? I kind of doubt it.

