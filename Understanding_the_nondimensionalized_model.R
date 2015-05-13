require(deSolve)
##PARAMETERS
k_M=0.3
#Structural volume per-capita maintenance rate (/day)
e_G=0.0016
#Cost of growth (?)
e_R=0.08
#Cost and conversion to produce an egg (?)
e_A=0.75
#Conversion efficiency of ingested to assimilated food (?)
v=1.9
#Energy conductance (mm/day)
R_M.bar=2
#Maturation threshold (units of eggs)
Eo=0.02
#Energy at birth (?)
Lo=0.65
#Length at birth (mm)
K=0.75
#Allocation of mobilization to structure (constant allocation variant)  (?)

T_M=2
#Known transfer interval (days)
f=0.05
#Known food provided at transfer (mgC)

Wo.nondim=(Eo/(Lo^3))/e_G
Lo.nondim=Lo*k_M/v
alpha=(k_M^2)*e_A/((v^3)*e_G)
kappa=K
beta=R_M.bar*e_R*(k_M^3)/(e_G*(v^3))
## scalars for re-dimensionalizing (not sure these are correct)
a=e_G*k_M
b=v
c=e_G*(v^3)/(e_R*(k_M^2))
d=k_M

dim.pars <- c(k_M=k_M,e_G=e_G,e_R=e_R,e_A=e_A,v=v,R_M.bar=R_M.bar,K=K,f=f,T_M=T_M)
nondim.pars <- c(alpha=alpha,kappa=kappa,beta=beta,f=f,T_M=T_M)
dim.y0 <- c(Eo/Lo^3,Lo,0)
nondim.y0 <- c(Wo.nondim,Lo.nondim,0)

LH1 = function(t, y, p) {
    k_M = p[1]
    e_G = p[2]
    e_R = p[3]
    e_A = p[4]
    v = p[5]
    R_M.bar = p[6]
    K = p[7]
    f = p[8]
    T_M = p[9]
    W <- y[1]
    L <- y[2]
    R_M <- y[3]

    P_C <- L^3*(v/L+k_M)*W/(1+K/e_G*W)

    dWdt <- (f/T_M)*(e_A/(L^3))-v*W/L
    dLdt <- (K/e_G*P_C-k_M*L^3)/(3*L^2)
    dR_Mdt <- (1-K)*P_C/e_R
    list(c(dWdt,dLdt,dR_Mdt))
}

LH2=function(t,y,p){ #p=c(alpha,kappa,beta,f,T_M)
    alpha=p[1];
    kappa=p[2];
    beta=p[3];
    f=p[4];
    T_M=p[5]
    W=y[1];
    L=y[2];
    R_M=y[3]

    ##//** ODEs **\\
    dWdt=(f/T_M)*(alpha/(L^3))-(W/L)
    dLdt=(kappa*W-L)/(3*(1+kappa*W))
    dR_Mdt=((1-kappa)*L*L*(1+L)*W)/(1+kappa*W)

    list(c(dWdt,dLdt,dR_Mdt))
}

LH3=function(t,y,p){ #p=c(alpha,kappa,beta,f,T_M)
    alpha=p[1];
    kappa=p[2];
    beta=p[3];
    f=p[4];
    T_M=p[5]
    W=y[1];
    L=y[2];
    R_M=y[3]

    ##//** ODEs **\\
    dWdt=(f/T_M)*(alpha/(L^3))-(W/L)
    dLdt=1/3*((1+L)*kappa*W/(1+kappa*W)-L)
    dR_Mdt=((1-kappa)*L*L*(1+L)*W)/(1+kappa*W)

    list(c(dWdt,dLdt,dR_Mdt))
}

out1 <- lsoda(y=dim.y0, times=seq(0,75,0.1), func=LH1, parms=dim.pars)
out2 <- lsoda(y=nondim.y0, times=seq(0,75,0.1), func=LH2, parms=nondim.pars)
## Bill's call to the integrator is slightly different.
out2.Bill <- lsoda(y=c(dim.y0[1]*d/a,dim.y0[2]*d/b,0), times=seq(0,75,0.1)*d, func=LH2, p=nondim.pars)

## The initial conditions for out2 and out2.Bill are identical
c(dim.y0[1]*d/a,dim.y0[2]*d/b,0)
nondim.y0

## However, out2 and out2.Bill do not agree on the final values of the
## nondimensionalized reproduction energy state variable, because out2
## has been integrated for too long. That is what I figured would
## happen - the reason I was investigating it was because I would like
## to be able to use pomp's facilities for fitting, but I can't
## because inside pomp, the length of time to integrate the
## differential equations (for trajectory matching, say) is determined
## by the length of the time series data. On the nondimensional scale,
## the length of time to integrate the ODEs that is equivalent to the
## length of the time series data is different.

## That still leaves the problem of the difference in the state
## equations between Bill's model and my own nondimensionalization. I
## should be able to redimensionalize the output and get either
## out2.Bill or out3 to lie directly atop the output in out1.
out3 <- lsoda(y=nondim.y0, times=seq(0,75,0.1)*d, func=LH3, parms=nondim.pars)

## You should be able to redimensionalize length by multiplying by
## b/d, which is equal to v/km, the reciprocal of the
## nondimensionalization constant for the L equation.
tail(out1[,3])
tail(out2[,3]*b/d)
tail(out2.Bill[,3]*b/d)
tail(out3[,3]*b/d)

## You should be able to nondimensionalize reproduction by multiplying by c/d.
tail(out1[,4])
tail(out2[,4]*c/d) ## this is much larger because of the too-long integration time
tail(out2.Bill[,4]*c/d)
tail(out3[,4]*c/d)
