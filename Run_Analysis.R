
rm(list = ls())

### Loading required packages
library(mvtnorm)
library(mBvs)
library(corrplot)

## Example Simulated Dataset 
load(file = "ExampleData.RData")

### Similar Design as the PHACS Data 
n = 254
q = 44
p0 = 1
p1 = 1
p_adj = 1
p <- p1+p_adj

Y <- Data[, 1:q]
nonz <- rep(NA, q)
for(j in 1:q) nonz[j] <- sum(Y[,j] != 0)

set.seed(976908)
#####################
## Hyperparameters ##
rho0     <- q + 3 + 1
Psi0    <- diag(3, q)

mu_alpha0     <- 0
mu_alpha    <- rep(0, q)

mu_beta0    <- 0
mu_beta      <- rep(0, q)

a_alpha0    <- 0.7
b_alpha0     <- 0.7

a_alpha        <- rep(0.7, p)
b_alpha     <- rep(0.7, p)

a_beta0        <- 0.7
b_beta0     <- 0.7

a_beta        <- rep(0.7, p)
b_beta         <- rep(0.7, p)

v_beta = rep(10, q)
omega_beta = rep(0.5, p1)
v_alpha = rep(10, q)
omega_alpha = rep(0.5, p0)

##
hyperParams.gen <- list(rho0=rho0, Psi0=Psi0, mu_alpha0=mu_alpha0, mu_alpha=mu_alpha,
mu_beta0=mu_beta0, mu_beta=mu_beta, a_alpha0=a_alpha0, b_alpha0=b_alpha0,
a_alpha=a_alpha, b_alpha=b_alpha, a_beta0=a_beta0, b_beta0=b_beta0,
a_beta=a_beta, b_beta=b_beta, v_beta=v_beta, omega_beta=omega_beta,
v_alpha=v_alpha, omega_alpha=omega_alpha)

###################
## MCMC SETTINGS ##

## Setting for the overall run
##
numReps    <- 1000000
thin       <- 100
burninPerc <- 0.5

## Settings for storage
##
storeV      <-    FALSE
storeW      <-    FALSE

## Tuning parameters for specific updates
##
##  - Generalized model
beta0.prop.var    <- 0.001
alpha.prop.var    <- 0.5
beta.prop.var    <- 0.001
V.prop.var    <- 0.5

##
mcmc.gen <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
storage=list(storeV=storeV, storeW=storeW),
tuning=list(beta0.prop.var=beta0.prop.var, alpha.prop.var=alpha.prop.var,
beta.prop.var=beta.prop.var, V.prop.var=V.prop.var))

#####################
## Starting Values ##

## Generalized model
##
B <- matrix(0.1, p, q, byrow = T)
A <- matrix(0.1, p, q, byrow = T)

V <- matrix(rnorm(n*q, 0, 0.1), n, q)
W <- matrix(rnorm(n*q, 0, 0.1), n, q)

beta0 <- log(as.vector(apply(Y, 2, mean)))
alpha0 <- log(nonz/n / ((n-nonz)/n))

Sigma_V    <- matrix(0, q, q)
diag(Sigma_V) <- 1

R        <- matrix(0, q, q)
diag(R) <- 1

sigSq_alpha0 <- 1
sigSq_alpha <- rep(1, p)
sigSq_beta0 <- 1
sigSq_beta <- rep(1, p)

startValues.gen <- vector("list", 1)
startValues.gen[[1]] <- list(B=B, A=A, V=V, W=W, beta0=beta0, alpha0=alpha0, R=R,
sigSq_alpha0=sigSq_alpha0,
sigSq_alpha=sigSq_alpha, sigSq_beta0=sigSq_beta0, sigSq_beta=sigSq_beta, Sigma_V=Sigma_V)

###################################
## Fitting the generalized model ##
###################################
fit <- mzipBvs(Y, lin.pred, Data, model="generalized", offset=Data$offset, hyperParams.gen, startValues.gen, mcmc.gen)

save(fit, file = "fit.RData")

#### Printing estimates of coefficients
sum.fit <- summary(fit)
round(sum.fit$AlphaDelta$hiv, 2)
round(sum.fit$BetaGamma$hiv, 2)

fit <- fit$chain1

#### For plots
postscript("sim_ZI.eps", width = 5, height = 5, horizontal = F)
hist(apply(fit$alpha0.p, 2, median), breaks = 20, xlim = c(0.2, 0.7), main="(b) Zero-inflation", xlab=expression(hat(alpha)[0][j]))
dev.off()

postscript("sim_OD.eps", width = 5, height = 5, horizontal = F)
hist(diag(apply(fit$Sigma_V.p, c(1,2), median)), xlim = c(1, 3.0), breaks = 20, main="(a) Overdispersion", xlab=expression(hat(Sigma)[V][j][j]))
dev.off()


corVal <- cor(log(Y+1))
R_V <- array(NA, c(q, q, 5000))
for(i in 1:5000)
{
	R_V[,,i] <- cov2cor(fit$Sigma_V.p[,,i])
}
R_V.est <- apply(R_V, c(1,2), median)
R.est <- apply(fit$R.p, c(1,2), median)

ord <- corrMatOrder(corVal, order = "hclust")

postscript("sim_cor_lnY.eps", width = 5, height = 5, horizontal = F)
corrplot(corVal[ord,ord], order = "original", title = expression(paste("(c) Observed")), mar = c(1, 0, 1.5, 0), cl.pos= "b", tl.pos= "n")
dev.off()

postscript("sim_R.eps", width =5, height = 5, horizontal = F)
corrplot(R.est[ord,ord], order = "original", title = expression(paste("(d) Binary part, ", widehat(R))), mar = c(1, 0, 1.5, 0), cl.pos= "b", tl.pos= "n")
dev.off()

postscript("Sim_R_V.eps", width = 5, height = 5, horizontal = F)
corrplot(R_V.est[ord,ord], order = "original", title = expression(paste("(e) Count part, ", widehat(R)[V])), mar = c(1, 0, 1.5, 0), cl.pos= "b", tl.pos= "n")
dev.off()


postscript("true_R.eps", width = 5, height = 5, horizontal = F)
corrplot(true$R[ord,ord], order = "original", title = expression(paste("(f) True R")), mar = c(1, 0, 1.5, 0), cl.pos= "b", tl.pos= "n")
dev.off()

