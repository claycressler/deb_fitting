load('results-2C-AA-stage-three.rda')
p.2C <- array(NA, dim=c(length(results3),16))
for (i in 1:length(results3)) {
  if (length(results3[[i]])>0)
    p.2C[i,] <- as.numeric(results3[[i]][order(results3[[i]][,'loglik']),][1,1:16])
}
colnames(p.2C) <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar','f','T_M','W.0','Lmeas.0','Rm.0','L.0','L.sd','R.sd','loglik')
save(p.2C, file='Clay_2C-AA_best_fit_parameters.rda')

## Compare the two to see how much the likelihoods have improved, and
## if the best fitting parameter sets are still *about* the same.

## What I want is to plot the distributions of the best-fitting
## parameter values for each animal. To make comparison with Bill's
## results easier, I will convert my scalar parameters to the same as
## his scalars (e.g., his 'a' is equal to my W.scalar*time.scalar
## because he parameterized a as W/T), and I will put back on the
## estimation scale.
a <- log(p.2C[,'W.scalar']*p.2C[,'time.scalar'])
b <- log(p.2C[,'L.scalar']*p.2C[,'time.scalar'])
c <- log(p.2C[,'Rm.scalar']*p.2C[,'time.scalar'])
d <- log(p.2C[,'time.scalar'])
alpha <- log(p.2C[,'alpha'])
kappa <- log(p.2C[,'kappa']/(1-p.2C[,'kappa']))
beta <- log(p.2C[,'beta'])
Wo.dim <- log(p.2C[,'W.0']*p.2C[,'W.scalar'])
Lo.dim <- log(p.2C[,'L.0']*p.2C[,'L.scalar'])
p.t.2C <- as.data.frame(array(NA, dim=c(length(a),9)))
p.t.2C[,1] <- Wo.dim
p.t.2C[,2] <- Lo.dim
p.t.2C[,3] <- alpha
p.t.2C[,4] <- kappa
p.t.2C[,5] <- beta
p.t.2C[,6] <- a
p.t.2C[,7] <- b
p.t.2C[,8] <- c
p.t.2C[,9] <- d
colnames(p.t.2C) <- c('Wo.dim','Lo.dim','alpha','kappa','beta','a','b','c','d')
save(p.t.2C, file='Clay_2C-AA_best_fit_parameters_transformed_scale.rda')

pairs(p.t.2C, col=factor(Lo.dim < -5))
dev.copy2pdf(file='Clay_2C-AA_p_correlations.pdf')

par(mfcol=c(3,3))
hist(Wo.dim, breaks=20)
hist(Lo.dim, breaks=20)
hist(alpha, breaks=20)
hist(kappa, breaks=20)
hist(beta, breaks=20)
hist(a, breaks=20)
hist(b, breaks=20)
hist(c, breaks=20)
hist(d, breaks=20)
dev.copy2pdf(file='Clay_2C-AA_p_histograms.pdf')

load('results-4B-AA-stage-three.rda')
p.4B <- array(NA, dim=c(length(results3),16))
for (i in 1:length(results3)) {
  if (is.data.frame(results3[[i]]))
    p.4B[i,] <- as.numeric(results3[[i]][order(results3[[i]]$loglik),][1,1:16])
}
colnames(p.4B) <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar','f','T_M','W.0','Lmeas.0','Rm.0','L.0','L.sd','R.sd','loglik')
save(p.4B, file='Clay_4B-AA_best_fit_parameters.rda')

## What I want is to plot the distributions of the best-fitting
## parameter values for each animal. To make comparison with Bill's
## results easier, I will convert my scalar parameters to the same as
## his scalars (e.g., his 'a' is equal to my W.scalar*time.scalar
## because he parameterized a as W/T), and I will put back on the
## estimation scale.
a <- log(p.4B[,'W.scalar']*p.4B[,'time.scalar'])
b <- log(p.4B[,'L.scalar']*p.4B[,'time.scalar'])
c <- log(p.4B[,'Rm.scalar']*p.4B[,'time.scalar'])
d <- log(p.4B[,'time.scalar'])
alpha <- log(p.4B[,'alpha'])
kappa <- log(p.4B[,'kappa']/(1-p.4B[,'kappa']))
beta <- log(p.4B[,'beta'])
Wo.dim <- log(p.4B[,'W.0']*p.4B[,'W.scalar'])
Lo.dim <- log(p.4B[,'L.0']*p.4B[,'L.scalar'])
p.t.4B <- as.data.frame(array(NA, dim=c(length(a),9)))
p.t.4B[,1] <- Wo.dim
p.t.4B[,2] <- Lo.dim
p.t.4B[,3] <- alpha
p.t.4B[,4] <- kappa
p.t.4B[,5] <- beta
p.t.4B[,6] <- a
p.t.4B[,7] <- b
p.t.4B[,8] <- c
p.t.4B[,9] <- d
colnames(p.t.4B) <- c('Wo.dim','Lo.dim','alpha','kappa','beta','a','b','c','d')
save(p.t.4B, file='Clay_4B-AA_best_fit_parameters_transformed_scale.rda')

pairs(p.t.4B)
dev.copy2pdf(file='Clay_4B-AA_p_correlations.pdf')

par(mfcol=c(3,3))
hist(Wo.dim)
hist(Lo.dim)
hist(alpha)
hist(kappa)
hist(beta)
hist(a)
hist(b)
hist(c)
hist(d)
dev.copy2pdf(file='Clay_4B-AA_p_histograms.pdf')

